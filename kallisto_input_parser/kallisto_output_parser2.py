import os
import pandas as pd
import argparse
import numpy as np


def check_target_ids(target_ids, curr_target_ids, file_path):
    """Check if target_ids match the reference target_ids."""
    if len(curr_target_ids) != len(target_ids):
        print(f"WARNING: {file_path} has {len(curr_target_ids)} target_ids, "
              f"expected {len(target_ids)}. Skipping.")
        return False
    if not np.array_equal(target_ids, curr_target_ids):
        print(f"WARNING: {file_path} has different target_id order/values. Skipping.")
        return False
    return True


def main(output_path, meta_path, kal_out_paths):
    meta_path = os.path.abspath(meta_path)
    output_path = os.path.abspath(output_path)

    try:
        meta_df = pd.read_csv(meta_path, usecols=['target_id'])
    except Exception as e:
        print(f"Error reading metadata file: {e}")
        return

    target_ids = meta_df['target_id'].values
    print(f"Metadata loaded: {len(target_ids)} target_ids")

    # collect all abundance files first
    abundance_files = []
    for kal_out_path in kal_out_paths:
        kal_out_path = os.path.abspath(kal_out_path)
        for root, dirs, files in os.walk(kal_out_path):
            if "abundance.tsv" in files:
                sample_id = os.path.basename(root)
                abundance_files.append((sample_id, os.path.join(root, "abundance.tsv")))

    print(f"Found {len(abundance_files)} samples")

    os.makedirs(output_path, exist_ok=True)

    # process in chunks, save each chunk as parquet
    chunk_size = 100  # adjust based on available RAM
    chunks = [abundance_files[i:i+chunk_size] for i in range(0, len(abundance_files), chunk_size)]
    chunk_files = []
    skipped_samples = []

    for chunk_idx, chunk in enumerate(chunks):
        print(f"Processing chunk {chunk_idx+1}/{len(chunks)}")

        series_list = []

        for sample_id, file_path in chunk:
            try:
                curr_df = pd.read_csv(
                    file_path,
                    sep='\t',
                    usecols=['target_id', 'est_counts'],
                    dtype={'target_id': str, 'est_counts': 'float32'}
                )

                # check target_id order and values match reference
                if not check_target_ids(target_ids, curr_df['target_id'].values, file_path):
                    skipped_samples.append(sample_id)
                    continue

                # extract est_counts as series with target_id as index
                s = curr_df.set_index('target_id')['est_counts'].rename(sample_id)
                series_list.append(s)

            except Exception as e:
                print(f"Error reading file {file_path}: {e}")
                skipped_samples.append(sample_id)
                continue

        if not series_list:
            print(f"No valid samples in chunk {chunk_idx+1}, skipping.")
            continue

        # single concat per chunk
        chunk_df = pd.concat(series_list, axis=1)

        # save chunk to parquet
        chunk_file = os.path.join(output_path, f'chunk_{chunk_idx}.parquet')
        chunk_df.to_parquet(chunk_file)
        chunk_files.append(chunk_file)

        del series_list, chunk_df

    if not chunk_files:
        print("No valid data found. Exiting.")
        return

    # reassemble all chunks
    print("Reassembling chunks...")
    result_df = pd.concat(
        [pd.read_parquet(f) for f in chunk_files],
        axis=1
    )

    output_file = os.path.join(output_path, 'kallisto_raw_counts_merged.parquet')
    result_df.to_parquet(output_file)
    print(f"All done! Output saved to {output_file}")
    print(f"Final matrix shape: {result_df.shape[0]} features x {result_df.shape[1]} samples")

    # cleanup chunk files
    for f in chunk_files:
        os.remove(f)
    print("Chunk files removed.")

    # report skipped samples
    if skipped_samples:
        print(f"\nSkipped {len(skipped_samples)} samples:")
        for s in skipped_samples:
            print(f"  - {s}")
    else:
        print("All samples processed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Merge raw counts from kallisto output in a single csv file. Optimized for large inputs.',
        usage='kallisto_input_parser /path/to/kal/out /path/to/meta.csv /path/to/result'
    )
    parser.add_argument('output_path', type=str, help='Path where the merged output CSV will be saved')
    parser.add_argument('meta_path', type=str, help='Path to the metadata file (CSV format)')
    parser.add_argument('kal_out_paths', type=str, nargs='+', help='Paths to the kallisto output directories')

    args = parser.parse_args()
    main(args.output_path, args.meta_path, args.kal_out_paths)