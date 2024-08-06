from Bio import SeqIO
import os
import vcfpy
import logging
from typing import List, Tuple, Optional, Any

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

def get_cds_info(pos: int, cds_list: List[Tuple[int, int, str, str]]) -> Tuple[bool, Tuple[Any, Any]]:
  """
  Gets CDS information for a given position.
    
  Parameters:
    pos (int): Position to check.
    cds_list (List[Tuple[int, int, str, str]]): List of CDS regions.
        
  Returns:
    Tuple[bool, Tuple[Optional[str], Optional[str]]]: Whether the position is within a CDS and the gene name and product (first CDS, if present in multiple ones).
  """
  for cds in cds_list:
    if cds[0] <= pos <= cds[1]:
      return True, (cds[2], cds[3])
  return False, (None, None)

def parse_record(
  sample_id: str, 
  record: vcfpy.Record, 
  all_cds_regions: List[Tuple[int, int, str, str]]
  ) -> Tuple[bool, Optional[List[Any]]]:
  """
  Parses a VCF record and extracts relevant information.
    
  Parameters:
    sample_id (str): Sample ID.
    record (vcfpy.Record): VCF record.
    all_cds_regions (List[Tuple[int, int, str, str]]): List of gene locations.
        
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
      mut_type = 'deletion' if len(record.REF) > len(record.ALT[0].value) else 'insertion'
      _, (gene_name, product) = get_cds_info(pos, all_cds_regions)
      return True, [sample_id, pos, ref, alt, dp, ad, freq, mut_type, gene_name, product]
  return False, None

def process_vcf_records(
  sample_id: str, 
  vcf_file_path: str, 
  all_cds_regions: List[Tuple[int, int, str, str]]
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
    valid_record, indel_info = parse_record(sample_id, record, all_cds_regions)
    if valid_record:
      data.append(indel_info)

  return data