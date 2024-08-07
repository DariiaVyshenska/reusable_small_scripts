import argparse
import logging
import os
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def reorder_columns(summary_df):
  determ_col = [
    'SAMPLE_ID', 'NONSYNONYMOUS', 'SYNONYMOUS',  
    'COMPLEX', 'STOPGAIN', 'STOPLOSS', 'TOTAL_SNPS',
    ]
  unmatched_columns = [col for col in summary_df.columns if col not in determ_col]
  new_col_order = [col for col in determ_col if col in summary_df.columns] + unmatched_columns
  
  return summary_df[new_col_order]

def summarize_snps(input_csv_path, output_csv_path):
  snp_data = pd.read_csv(input_csv_path, usecols=['Sample', 'Syn'])
  
  summary_df = snp_data.groupby(['Sample', 'Syn']).size().unstack(fill_value=0).reset_index()
  summary_df.rename_axis(None, axis=1, inplace=True)
  summary_df.rename(columns={"Sample": "sample_id", 
                             "nonsynonymous SNV": "nonsynonymous", 
                             "synonymous SNV": "synonymous"}, 
                    inplace=True)
  summary_df['total_snps'] = summary_df.iloc[:, 1:].sum(axis=1)
  
  summary_df.columns = summary_df.columns.str.upper()
  summary_df['SAMPLE_ID'] = summary_df['SAMPLE_ID'].str.replace('.fastq.gz', '', regex=False)
  ordered_df = reorder_columns(summary_df)
  ordered_df.to_csv(output_csv_path, index = False) 

def main(input_csv_path, output_path):
  try:
    os.makedirs(output_path, exist_ok=True)
    output_csv_path = os.path.join(output_path, 'snp_summary.csv')
    summarize_snps(input_csv_path, output_csv_path)
  except Exception as e:
    logging.error(f"Error parsing visualization csv file: {e}")
  
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Summarizes RAVA visualization.csv file: types of SNPs and counts per sample.',
                                    usage='python main.py <input_csv_path> <output_path>')
  parser.add_argument('input_csv_path', type=str, help='Path to the RAVA generated visualization.csv file.')
  parser.add_argument('output_path', type=str, help='Path where the output summary CSV will be saved')
    
  args = parser.parse_args()
  main(args.input_csv_path, args.output_path)