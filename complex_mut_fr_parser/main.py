#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:47:33 2024

@author: Dariia Vyshenska
"""
import argparse
from parsers import extract_codon_frequencies
from io_utils import codon_stats_to_csv

def main(bam_file_path, codon_start_pos, output_file_path):
  codon_frequencies = extract_codon_frequencies(bam_file_path, codon_start_pos)
  codon_stats_to_csv(codon_frequencies, output_file_path, bam_file_path, codon_start_pos)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extracts frequencies of complex mutations from INDEXED .bam file.',
                                  usage='Usage: ./main.py <indexed_bam_path> <ref_codon_start_pos> <output_dir>')
  parser.add_argument('indexed_bam_path', type=str, help='Path to the indexed bam file')
  parser.add_argument('ref_codon_start_pos', type=int, help='Codon start position in reference file')
  parser.add_argument('output_path', type=str, help='Path where the output CSV will be saved')
  
  args = parser.parse_args()
  main(args.indexed_bam_path, args.ref_codon_start_pos, args.output_path)
