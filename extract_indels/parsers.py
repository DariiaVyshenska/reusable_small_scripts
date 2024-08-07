from Bio import SeqIO
import os
import vcfpy
import logging
from typing import List, Tuple, Optional, Any

STRND_BIAS_THRESH = 0.9
REF_SEQ_LEN_EXTRACT = 15

def get_full_ref_seq(gb_file_path: str) -> str:
  """
  Extracts the full reference nucleotide sequence from a GenBank file.

  Parameters:
      gb_file_path (str): Path to the GenBank file.
    
  Returns:
      str: The full reference nucleotide sequence.
  """
  with open(gb_file_path, 'r') as file:
      for record in SeqIO.parse(file, 'genbank'):
          return str(record.seq)

def get_gene_locations(gb_file_path: str) -> List[Tuple[int, int, str, str]]:
  """
  Extracts gene locations from a GeneBank file.
  
  Parameters:
    gb_file_path (str): Path to the GeneBank file.
      
  Returns:
    List[Tuple[int, int, str, str]]: List of gene locations with start, end, gene name, and product.
  """
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
  return all_cds

def extract_sample_id(vcf_path: str) -> str:
  """
  Extracts the sample ID from a VCF file path.
    
  Parameters:
    vcf_path (str): Path to the VCF file.
      
  Returns:
    str: Extracted sample ID.
  """
  sample_id_filename = os.path.basename(vcf_path)
  sample_id = os.path.splitext(sample_id_filename)[0]
  sample_id = os.path.splitext(sample_id)[0]
  sample_id = os.path.splitext(sample_id)[0]
  return sample_id

def parse_vcf(vcf_file_path: str) -> vcfpy.Reader:
  """
  Parses a VCF file and returns a VCF reader object.
  
  Parameters:
    vcf_file_path (str): Path to the VCF file.
      
  Returns:
    vcfpy.Reader: VCF reader object.
  """
  return vcfpy.Reader.from_path(vcf_file_path)

def get_cds_info(pos: int, cds_list: List[Tuple[int, int, str, str]]) -> Tuple[bool, Tuple[Any]]:
  """
  Gets CDS information for a given position.
    
  Parameters:
    pos (int): Position to check.
    cds_list (List[Tuple[int, int, str, str]]): List of CDS regions.
        
  Returns:
    Tuple[bool, Tuple[Optional[str]]]: Whether the position is within a CDS and the gene name and product (first CDS, if present in multiple ones).
  """
  for cds in cds_list:
    if cds[0] <= pos <= cds[1]:
      return True, (cds[2])
  return False, (None)

def get_del_placement(pos: int, mut_seq: str, full_ref_seq: str) -> str:
  """
  Provides the string representation of an deletion placement in reference.
  Deletion seq is marked by '*' and surrounded by 15 nt reference seq from both sides.

  Parameters:
      pos (int): Position of the mutation in the reference sequence.
      mut_seq (str): Sequence of the mutation from vcf output (column 'REF' value).
      full_ref_seq (str): The full reference nucleotide sequence.

  Returns:
      str: The reference sequence with the deletion marked by asterisks.
  """
  
  result = []
  idx = pos - 1
  mut_len = len(mut_seq)
  result.append(full_ref_seq[pos - REF_SEQ_LEN_EXTRACT:pos])
  result.append(full_ref_seq[idx+1:idx + mut_len])
  result.append(full_ref_seq[idx + mut_len:idx + mut_len + REF_SEQ_LEN_EXTRACT])
  return '*'.join(result)

def get_insert_placement(pos: int, mut_seq: str, full_ref_seq: str) -> str:
  """
  Provides the string representation of an insertion placement in reference.
  Insertion seq is marked by '*' and surrounded by 15 nt reference seq from both sides.

  Parameters:
      pos (int): Position of the mutation in the reference sequence.
      mut_seq (str): Sequence of the mutation (column 'ALT' value).
      full_ref_seq (str): The full reference nucleotide sequence.

  Returns:
      str: The reference sequence with the insertion marked by asterisks.
  """

  result = []
  result.append(full_ref_seq[pos - REF_SEQ_LEN_EXTRACT:pos])
  result.append(mut_seq[1:])
  result.append(full_ref_seq[pos:pos + REF_SEQ_LEN_EXTRACT])
  return '*'.join(result)

def parse_record(
  sample_id: str, 
  record: vcfpy.Record, 
  all_cds_regions: List[Tuple[int, int, str, str]],
  full_ref_seq: str
  ) -> Tuple[bool, Optional[List[Any]]]:
  """
  Parses a VCF record and extracts relevant information.
    
  Parameters:
    sample_id (str): Sample ID.
    record (vcfpy.Record): VCF record.
    all_cds_regions (List[Tuple[int, int, str, str]]): List of gene locations.
    full_ref_seq (str): full reference sequence string.
        
    Returns:
      Tuple[bool, Optional[List[Any]]]: Validity of the record and extracted information.
    """
  if record.ALT and len(record.REF) != len(record.ALT[0].value):
    sample_call = record.call_for_sample['Sample1']
    dp = sample_call.data.get('DP')
    ad = sample_call.data.get('AD')
    freq = (round(ad / dp * 100, 2)) if dp else 0

    if dp >= 30 and freq >= 10:
      pos = record.POS
      ref = record.REF
      alt = record.ALT[0].value
      
      adf = sample_call.data.get('ADF')
      adr = sample_call.data.get('ADR')
      adf_ratio = round(adf / (adf + adr), 2)
      adr_ratio = round(adr / (adf + adr), 2)
      strn_bias_status = (adf_ratio < STRND_BIAS_THRESH and adr_ratio < STRND_BIAS_THRESH)
      mut_type = 'deletion' if len(record.REF) > len(record.ALT[0].value) else 'insertion'
      valid_len = ((len(ref) - 1) % 3 == 0) if mut_type == 'deletion' else (len(alt) - 1) % 3 == 0
      

      placement =(
        get_del_placement(pos, ref, full_ref_seq) 
        if mut_type == 'deletion' 
        else get_insert_placement(pos, alt, full_ref_seq))
      
      _, (product) = get_cds_info(pos, all_cds_regions)
      return True, [sample_id, pos, ref, alt, dp, ad, freq, 
                    adf_ratio, adr_ratio, strn_bias_status, 
                    valid_len, mut_type, product, placement]
  return False, None

def process_vcf_records(
  sample_id: str, 
  vcf_file_path: str, 
  all_cds_regions: List[Tuple[int, int, str, str]],
  full_ref_seq: str
  ) -> List[List[Any]]:
  """
  Processes VCF records and extracts relevant information.
    
  Parameters:
    sample_id (str): Sample ID.
    vcf_file_path (str): Path to the VCF file.
    all_cds_regions (List[Tuple[int, int, str, str]]): List of CDS regions.
        
  Returns:
    List[List[Any]]: Extracted information from VCF records.
  """
  data = []

  vcf_reader = parse_vcf(vcf_file_path)
  for record in vcf_reader:
    valid_record, indel_info = parse_record(sample_id, record, all_cds_regions, full_ref_seq)
    if valid_record:
      data.append(indel_info)

  return data