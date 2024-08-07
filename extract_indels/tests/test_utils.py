import unittest
import os
import csv
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import write_csv

class TestUtils(unittest.TestCase):
  def setUp(self):
    self.sample_id = 'test_sample'
    self.output_path = 'test_output'
    self.data = [
      [
        "test_sample", "12345", "A", "T", 
        "100", "20", "0.2", "0.1", "0.9", 
        "FALSE", "TRUE", "insertion", "PRODUCT1", 
        "placement1"],
      [
        "test_sample", "67890", "G", "C", 
        "150", "30", "0.3", "0.5", "0.5", 
        "TRUE", "FALSE", "deletion", "PRODUCT2", 
        "placement2"]
    ]
    self.expected_file = os.path.join(self.output_path, f"{self.sample_id}_indels_af10_dp30.csv")
    
  def tearDown(self):
    if os.path.exists(self.expected_file):
      os.remove(self.expected_file)
    if os.path.exists(self.output_path):
      os.rmdir(self.output_path)
  
  def test_write_csv(self):
    write_csv(self.sample_id, self.output_path, self.data)
    
    self.assertTrue(os.path.exists(self.expected_file))
    
    with open(self.expected_file, 'r') as csv_file:
      csv_reader = csv.reader(csv_file)
      header = next(csv_reader)
      self.assertEqual(header, ['SAMPLE_ID', 'POSITION', 'REFERENCE_SEQ', 'ALTERNATIVE_SEQ', 
                                'SEQ_DEPTH', 'ALT_SEQ_DEPTH', 'FREQUENCY', 'ADF_RATIO', 
                                'ADR_RATIO', 'STRND_BIAS_PASS', 'VALID_CODON_LEN', 
                                'CHANGE_TYPE', 'PRODUCT', 'PLACEMENT'])
      
      for idx, row in enumerate(csv_reader):
        self.assertEqual(row, self.data[idx])
       
            

if __name__ == '__main__':
  unittest.main()