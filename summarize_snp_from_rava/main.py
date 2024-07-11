import logging
import os
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
  summary_df.to_csv(output_csv_path, index = False)
  

def main(input_csv_path, output_path):
  try:
    output_csv_path = os.path.join(output_path, 'snp_summary.csv')
    summarize_snps(input_csv_path, output_csv_path)
  except Exception as e:
    logging.error(f"Error parsing visualization csv file: {e}")
  
if __name__ == '__main__':
  # arguments - input csv and output folder
  input_csv_path = './test_data/test_visualization.csv'
  output_path = '.'
  main(input_csv_path, output_path)