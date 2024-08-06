import pysam
import pandas as pd
import os
from exceptions import GenomeRefError

def freq_dict_to_csv(freq_dict, output_file_path, bam_file_path, codon_start_pos):
  df = pd.DataFrame(list(freq_dict.items()), columns=['CODON', 'FREQUENCY'])
  df = df.sort_values(by='FREQUENCY', ascending=False)
  print(df)
  
  basename = os.path.basename(bam_file_path)
  smpl_id = basename.replace('.fastq.gz.bam', '')
  output_path = f"{output_file_path}/{smpl_id}_{codon_start_pos}_complex_freqs.csv"
  df.to_csv(output_path, index=False)
  
def open_and_validate_bam(bam_file):
  bam = pysam.AlignmentFile(bam_file, "rb")
  if len(bam.references) != 1:
    raise GenomeRefError('GenomeRefError: input .bam must contain one and only one reference (chromosome)')
  return bam