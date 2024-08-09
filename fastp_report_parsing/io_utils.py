import os
import csv
from parsers import parse_fastp_stats

def extract_data(json_files, output_file_path):
  header = [
    'SAMPLE_ID',
    'TOTAL_READS', 'DUPLIC_RATE', 
    'INSERT_SIZE_PEAK', 'INSERT_SIZE_UNKNOWN', 
    'GC_CONT_BEFORE_FILT', 'GC_CONT_AFTER_FILT', 
    'PASSED_FILT_READS', 'PASSED_FILT_RATIO', 
    'LOW_QUAL_READS', 'LOW_QUAL_RATIO', 
    'TOO_MANY_N_READS', 'TOO_MANY_N_RATIO', 
    'TOO_SHORT_READS', 'TOO_SHORT_RATIO', 
  ]

  with open(output_file_path, 'w', newline='') as output_file:
    writer = csv.writer(output_file)
    writer.writerow(header)

    for json_file in json_files:
      basename = os.path.basename(json_file)
      sample_id = [basename.replace('_report.json', '')]
      fastp_stats = parse_fastp_stats(json_file)
      current_row = [sample_id] + fastp_stats
      writer.writerow(current_row)
