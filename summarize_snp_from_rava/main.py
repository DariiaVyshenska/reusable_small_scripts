#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd
from utils import format_snp_table
from io_utils import dump_haf_info, dump_haf_summary, dump_all_snp_info
from parsers import get_haf_info

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main(af_threshold, input_csv_path, output_path):
  if not os.path.exists(output_path):
    os.makedirs(output_path)

  try:
    snp_data = pd.read_csv(input_csv_path)
  except FileNotFoundError:
    print(f"Error: The file {input_csv_path} was not found.")
    return
  except Exception as e:
    logging.error(f"Error parsing visualization csv file: {e}")
    return
  
  f_snp_data = format_snp_table(snp_data)
  dump_all_snp_info(f_snp_data, output_path)

  haf_snp_info = get_haf_info(f_snp_data, af_threshold)
  dump_haf_info(haf_snp_info, af_threshold, output_path)
  dump_haf_summary(haf_snp_info, af_threshold, output_path)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Summarize SNPs identified by RAVA above certain AF threshold.',
                                    usage='./main.py [<af_threshold>] <input_csv_path> <output_path>')
  parser.add_argument('af_threshold', type=int, nargs='?', default=25, help='Allelic frequency threshold (default: 25)')
  parser.add_argument('input_csv_path', type=str, help='Path to the RAVA generated visualization.csv file')
  parser.add_argument('output_path', type=str, help='Path where the output summary CSV will be saved')

  args = parser.parse_args()
  main(args.af_threshold, args.input_csv_path, args.output_path)
  logging.info('High alellic frequency SNPs info is successfully extracted!')