import argparse
from parsers import get_gene_locations, extract_sample_id, process_vcf_records
from utils import write_csv
from typing import List, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main(vcf_file_path: str, gb_file_path: str, output_path: str) -> None:
  try:
    all_cds_regions = get_gene_locations(gb_file_path)
  except Exception as e:
    logging.error(f"Error getting gene locations: {e}")
    return
  
  sample_id = extract_sample_id(vcf_file_path)
  try:
    data = process_vcf_records(sample_id, vcf_file_path, all_cds_regions)
  except Exception as e:
    logging.error(f"Error processing VCF records: {e}")
    return
  
  try:
    write_csv(sample_id, output_path, data)
  except Exception as e:
    logging.error(f"Error writing output CSV: {e}")
    return
  
  logging.info("All done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts indels from sorted vcf that: AF >= .1 and depth >= 30.',
                                     usage='Usage: python main.py <vcf_path> <gb_path> <output_dir>')
    parser.add_argument('vcf_file_path', type=str, help='Path to the sorted vcf file')
    parser.add_argument('gb_file_path', type=str, help='Path to the GeneBank file')
    parser.add_argument('output_path', type=str, help='Path where the output CSV will be saved')
    
    args = parser.parse_args()
    main(args.vcf_file_path, args.gb_file_path, args.output_path)