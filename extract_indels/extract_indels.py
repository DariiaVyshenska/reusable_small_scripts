import os
import vcfpy
import argparse
import csv
from Bio import SeqIO

# notes on docs:
# 1. the CDS is defined by the first position of the indel. if indel spans
# more than one CDS only the first CDS, where the start of the indel is
# will be recorded.



def get_gene_locations(gb_file_path):
  all_cds = []
  with open(gb_file_path, 'r') as file:
    for record in SeqIO.parse(file, 'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                gene_name = feature.qualifiers.get('gene', ['N/A'])[0]
                product = feature.qualifiers.get('product', ['N/A'])[0]
                all_cds.append((start, end, gene_name, product))
  return(all_cds)

def get_cds_info(pos, cds_list):
  for cds in cds_list:
    if cds[0] <= pos <= cds[1]:
      return cds[2], cds[3]
  return None, None


def main(vcf_file_path, gb_file_path, output_path):
  os.makedirs(output_path, exist_ok=True)
  full_indel_csv_filename = os.path.join(output_path, 'indels_af10_dp30_vcf_extraction.csv')

  all_cds_regions = get_gene_locations(gb_file_path)

  sample_id_filename = os.path.basename(vcf_file_path)
  sample_id = os.path.splitext(sample_id_filename)[0]
  sample_id = os.path.splitext(sample_id)[0] # this is due to how RAVA names vcf-s

  vcf_reader = vcfpy.Reader.from_path(vcf_file_path)
  
  
  with open(full_indel_csv_filename, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['SAMPLE_ID', 'POSITION', 'REFERENCE_SEQ', 
                         'ALTERNATIVE_SEQ', 'SEQ_DEPTH', 
                         'ALT_SEQ_DEPTH', 'FREQUENCY', 'CHANGE_TYPE',
                         'GENE_NAME', 'PRODUCT'])
    for record in vcf_reader:
      if record.ALT and len(record.REF) != len(record.ALT[0].value):
        sample_call = record.call_for_sample['Sample1']

        pos = record.POS
        ref = record.REF
        alt = record.ALT[0].value
        dp = sample_call.data.get('DP')
        ad = sample_call.data.get('AD')
        if dp:
          freq = ad / dp * 100
        else:
          freq = 0
        type = ''
        if len(record.REF) > len(record.ALT[0].value):
          mut_type = 'deletion'
        else:
          mut_type = 'insertion'

        if dp >= 30 and freq >= 10:
          gene_name, product = get_cds_info(pos, all_cds_regions)
          curr_row = [sample_id, pos, ref, alt, dp, 
                     ad, freq, mut_type, gene_name,
                     product]
          csv_writer.writerow(curr_row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extracts indels from sorted vcf that: AF >= .1 and depth >= 30.',
                                   usage='extract_indels ./path/to/sorted/vcf ./path/to/gb/file ./output/dir')
  parser.add_argument('vcf_file_path', type=str, help='Path to the sorted vcf file')
  parser.add_argument('gb_file_path', type=str, help='Path to the GeneBank file')
  parser.add_argument('output_path', type=str, help='Path where the output CSV will be saved')
  
  args = parser.parse_args()
  main(args.vcf_file_path, args.gb_file_path, args.output_path)
  print("All done!")
