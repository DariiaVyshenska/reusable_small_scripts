#!/usr/bin/env python3

import argparse
import os
import glob
from io_utils import extract_data

def main(json_files_path, output_path):
  json_files_paths = glob.glob(os.path.join(json_files_path, "*.json"))
  output_file_path = f"{output_path}/fastp_qc_report.csv"
  extract_data(json_files_paths, output_file_path)
  
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extracts selected statistics from fastp json report files.',
                                  usage='Usage: ./main.py <json_files_path> <output_dir>')
  parser.add_argument('json_files_path', type=str, help='Path to the directory with fastp json report files')
  parser.add_argument('output_path', type=str, help='Path where the output CSV will be saved')
  
  args = parser.parse_args()
  main(args.json_files_path, args.output_path)
  print('Parsing fastp json reports is done!')