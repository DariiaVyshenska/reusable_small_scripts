import pandas as pd
import json
import argparse
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_keys_csv(keys_file_path):
    try:
        keys_table = pd.read_csv(keys_file_path, usecols=['target_id', 'peptide_seq'])
        if keys_table.empty:
            logging.error('The CSV file with keys data is empty or lacks required columns.')
            return None
        return keys_table
    except Exception as e:
        logging.error(f"Could not read library keys file: {e}")
        return None

def dump_kmer_db(result_file_path, keys_table, kmer_len):
    with open(result_file_path, 'w') as file:
        for row in keys_table.itertuples():
            peptide_id = row[1]
            peptide_seq = row[2]
            curr_hash = {}
            kmers = set(peptide_seq[idx:idx+kmer_len] for idx in range(len(peptide_seq) - kmer_len + 1))
            json_obj = json.dumps({peptide_id: list(kmers)})
            file.write(json_obj + '\n')

def main(kmer_len, keys_file_path, output_path):
    output_path = os.path.abspath(output_path)
    os.makedirs(output_path, exist_ok=True)
    result_file_path = os.path.abspath(f'{output_path}/peptide_kmer_db.jsonl')
    keys_file_path = os.path.abspath(keys_file_path)
    
    if os.path.exists(result_file_path):
        logging.error(f'The file {result_file_path} already exists. Aborting to prevent data loss.')
        return
    
    keys_table = read_keys_csv(keys_file_path)
    
    if keys_table is None:
        return
    
    keys_table.sort_values(by = 'target_id', inplace=True)
    dump_kmer_db(result_file_path, keys_table, kmer_len)
    logging.info('All done!')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script creates a simple Phip-seq library database (.jsonl) of k-mers per peptide of a given size.',
                                   usage='python phipseq-key-kmer-parser.py -k 9 /path/to/keys/table /path/to/dump/output')
    parser.add_argument('-k', '--kmer_len', type=int, default=9, help='(optional) Desired k-mer length. Default is 9.')
    parser.add_argument('keys_file_path', type=str, help='Path to the Phip-seq keys file (CSV format). Must contain "target_id" and "peptide_seq" columns.')
    parser.add_argument('output_path', type=str, help='Path where the .jsonl database will be saved.')
    
    args = parser.parse_args()
    main(args.kmer_len, args.keys_file_path, args.output_path)