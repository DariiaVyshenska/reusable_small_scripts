import pysam
from collections import defaultdict
from utils import calc_freq
from io_utils import open_and_validate_bam

def process_pileup(bam, codon_start_pos):
  codon_end_pos = codon_start_pos + 2
  read_ids = set()
  codon_counts = defaultdict(int)
  ref_name = bam.references[0]

  for pileup_column in bam.pileup(reference=ref_name, start=codon_start_pos-1, stop=codon_end_pos): # do I need stop position?  what does it do?
    for pileup_read in pileup_column.pileups:
      if pileup_read.is_del or pileup_read.is_refskip:
        continue
      
      read_pos = pileup_read.query_position
      if read_pos is None:
        continue
      
      read_id = pileup_read.alignment.query_name
      if read_id in read_ids:
        continue
      
      read = pileup_read.alignment.query_sequence
      if read_pos + 2 >= len(read):
        continue
      
      read_ids.add(read_id)
      codon = read[read_pos:read_pos+3]
      codon_counts[codon] += 1
    break
  return codon_counts, len(read_ids)

def extract_codon_frequencies(bam_file, codon_start_pos):
  try:
    bam = open_and_validate_bam(bam_file)
    codon_counts, total_reads = process_pileup(bam, codon_start_pos)
    print(calc_freq(codon_counts, total_reads))
    return calc_freq(codon_counts, total_reads)
  except Exception as e:
    print(e)
  finally:
    bam.close()

