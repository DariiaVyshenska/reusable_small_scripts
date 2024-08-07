import unittest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from parsers import (
  extract_sample_id, 
  get_cds_info, 
  get_del_placement,
  get_insert_placement)

class TestParsers(unittest.TestCase):
  def test_extract_sample_id(self):
    self.assertEqual(extract_sample_id('./my/path/test_full.fastq.vcf'), 'test_full')
    self.assertEqual(extract_sample_id('/path/to/test_full.fastq.vcf'), 'test_full')

  def test_get_cds_info(self):
    self.assertEqual(get_cds_info(3, [(1, 5, 'PRODUCT1')]), (True, ('PRODUCT1')))
    self.assertEqual(get_cds_info(5, [(1, 5, 'PRODUCT1')]), (True, ('PRODUCT1')))
    self.assertEqual(get_cds_info(1, [(1, 5, 'PRODUCT1')]), (True, ('PRODUCT1')))
    self.assertEqual(get_cds_info(8, [(1, 5, 'PRODUCT1')]), (False, (None)))
    self.assertEqual(get_cds_info(11, [(1, 5, 'PRODUCT1'), (6, 10, 'PRODUCT2')]), (False, (None)))
    self.assertEqual(get_cds_info(7, []), (False, (None)))
    
  def test_get_del_placement(self):
    self.assertEqual(
      get_del_placement(
        15, 
        'GCTATACATGTCTCTGGGACCAATGGTA', 
        'TTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATABBB'),
      'TTACTTGGTTCCATG*CTATACATGTCTCTGGGACCAATGGTA*CTAAGAGGTTTGATA')
  
  def test_get_insert_placement(self):
    self.assertEqual(
      get_insert_placement(
        37,
        'ABBB',
        'TTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGACT'
      ),
      'TAATGATGGTGTTTA*BBB*TTTTGCTTCCACTGA'
    )

if __name__ == '__main__':
  unittest.main()
