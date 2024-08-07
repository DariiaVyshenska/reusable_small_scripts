import csv
import os
from typing import List

def write_csv(sample_id: str, output_path: str, data: List) -> None:
  """
  Writes the extracted indel information to a CSV file.

  Parameters:
    sample_id (str): The sample ID to include in the CSV file name.
    output_path (str): The directory where the CSV file will be saved.
    data (List[List[Any]]): The data to be written to the CSV file. Each sublist represents a row.

  Returns:
    None
  """
  os.makedirs(output_path, exist_ok=True)
  file_name = f"{sample_id}_indels_af10_dp30.csv"
  full_indel_csv_filename = os.path.join(output_path, file_name)
  
  with open(full_indel_csv_filename, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['SAMPLE_ID', 'POSITION', 'REFERENCE_SEQ', 'ALTERNATIVE_SEQ', 
                         'SEQ_DEPTH', 'ALT_SEQ_DEPTH', 'FREQUENCY', 'ADF_RATIO', 
                         'ADR_RATIO', 'STRND_BIAS_PASS', 'VALID_CODON_LEN',  
                         'CHANGE_TYPE', 'PRODUCT', 'PLACEMENT'])
    for row in data:
      csv_writer.writerow(row)