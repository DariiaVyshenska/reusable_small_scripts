from Bio import SeqIO
import os
import vcfpy
from typing import List, Tuple, Optional, Any, Iterable, Callable

REF_SEQ_LEN_EXTRACT = 15

RecordParser = Callable[
    [str, str, vcfpy.Record, List[Tuple[int, int, str, str]]],
    Optional[List[Any]]
]

def format_arr(arr: Iterable[Any]) -> str:
  return ';'.join(str(el) for el in arr)

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
          product = feature.qualifiers.get('product', ['N/A'])[0]
          all_cds.append((start, end, product))
  return all_cds

def extract_file_name(vcf_path: str) -> str:
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
  reader = vcfpy.Reader.from_path(vcf_file_path)
  if reader is None:
    raise RuntimeError(f'Failed to parse VCF file: {vcf_file_path}')
  return reader

def get_cds_info(
  pos: int, 
  cds_list: List[Tuple[int, int, str, str]]
  ) -> Tuple[bool, Optional[str]]:
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
      return True, cds[2]
  return False, None

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


def parse_legacy_record(
  vcf_sample_id: str,
  sample_id: str, 
  record: vcfpy.Record, 
  all_cds_regions: List[Tuple[int, int, str, str]],
  full_ref_seq: str
  ) -> Optional[List[Any]]:
  """
  Parses a VCF record and extracts relevant information.
    
  Parameters:
    vcf_sample_id (str): Sample ID as it appears in VCF.
    sample_id (str): Sample ID from the file name.
    record (vcfpy.Record): VCF record.
    all_cds_regions (List[Tuple[int, int, str, str]]): List of gene locations.
    full_ref_seq (str): full reference sequence string.
        
    Returns:
      Optional[List[Any]]: Valid record or None.
    """
  if record.ALT and len(record.REF) != len(record.ALT[0].value):
    sample_call = record.call_for_sample[vcf_sample_id]
    dp = sample_call.data.get('DP')
    ad = sample_call.data.get('AD')
    freq = (round(ad / dp * 100, 2)) if dp else 0

    if dp >= 30 and freq >= 10:
      pos = record.POS
      ref = record.REF
      alt = record.ALT[0].value
      
      mut_type = 'deletion' if len(record.REF) > len(record.ALT[0].value) else 'insertion'
      valid_len = ((len(ref) - 1) % 3 == 0) if mut_type == 'deletion' else (len(alt) - 1) % 3 == 0

      placement =(
        get_del_placement(pos, ref, full_ref_seq) 
        if mut_type == 'deletion' 
        else get_insert_placement(pos, alt, full_ref_seq))
      
      _, product = get_cds_info(pos, all_cds_regions)
      return [sample_id, pos, ref, alt, dp, ad, freq, 
              valid_len, mut_type, product, placement]
  return None

def parse_mutect2_record(
  vcf_sample_id: str,
  sample_id: str, 
  record: vcfpy.Record, 
  all_cds_regions: List[Tuple[int, int, str, str]],
  full_ref_seq: str
  ) -> Optional[List[Any]]:
  """
  Parses a VCF record and extracts relevant information.
    
  Parameters:
    vcf_sample_id (str): Sample ID as it appears in VCF.
    sample_id (str): Sample ID from the file name.
    record (vcfpy.Record): VCF record.
    all_cds_regions (List[Tuple[int, int, str, str]]): List of gene locations.
        
    Returns:
      Optional[List[Any]]: Valid record or None.
  """
  types_to_detect = ('DEL', 'INS', 'INDEL')
  alts = record.ALT or []
  if not any(alt.type in types_to_detect for alt in alts):
    return None
  
  if len(alts) > 1:
    raise ValueError(
      f"{sample_id} ({vcf_sample_id}) has multi-allelic INDEL at "
      f"{record.CHROM}:{record.POS} REF={record.REF} ALTs={[a.value for a in alts]}"
    )
  
  sample_call = record.call_for_sample[vcf_sample_id]
  dp = sample_call.data.get('DP')
  if dp is None or dp < 30:
    return None
  
  af = sample_call.data.get('AF')[0]
  if (af is None) or (af < 0.1):
    return None
  
  ad = sample_call.data.get('AD')[1]
  pos = record.POS
  ref = record.REF
  alt = alts[0].value
  mut_type = alts[0].type
  
  validity = 'NA'
  if mut_type == 'DEL':
    validity = (len(ref) - 1) % 3 == 0
  elif mut_type == 'INS':
    validity = (len(alts[0].value) - 1) % 3 == 0

  placement = 'NA'
  if mut_type == 'DEL':
    placement = get_del_placement(pos, ref, full_ref_seq)
  elif mut_type == 'INS':
    placement =get_insert_placement(pos, alt, full_ref_seq)

  _, product = get_cds_info(pos, all_cds_regions)

  return [sample_id, pos, ref, alt, dp, 
                ad, af,
                validity, mut_type, product, placement]


def detect_record_type(vcf_header: vcfpy.Header) -> str:
  return (
    'mutect2' if any(
    ((line.key == 'source') and (line.value == 'Mutect2')) 
    for line in vcf_header.lines
    ) else 'legacy'
  )

def get_record_parser(record_type: str) -> RecordParser:
  return parse_legacy_record if record_type == 'legacy' else parse_mutect2_record

def get_sample_name_from_vcf(vcf_reader: vcfpy.Reader) -> str:
  sample_names = vcf_reader.header.samples.names
  if len(sample_names) != 1:
    raise ValueError(f"Expected exactly 1 sample, found {len(sample_names)}")
  return sample_names[0]
  
def process_vcf_records(
  record_type: str,
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
    full_ref_seq (str): full reference sequence.
        
  Returns:
    List[List[Any]]: Extracted information from VCF records.
  """
  vcf_reader = parse_vcf(vcf_file_path)
  if record_type == 'auto':
    record_type = detect_record_type(vcf_reader.header)
  record_parser = get_record_parser(record_type)
  vcf_sample_id = get_sample_name_from_vcf(vcf_reader)

  data = []
  for record in vcf_reader:
    indel_info = record_parser(vcf_sample_id, sample_id, record, all_cds_regions, full_ref_seq)
    if indel_info is not None:
      data.append(indel_info)

  return data