import unittest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from parsers import extract_sample_id, get_cds_info

class TestParsers(unittest.TestCase):
  def test_extract_sample_id(self):
    self.assertEqual(extract_sample_id('./my/path/test_full.fastq.vcf'), 'test_full')
    self.assertEqual(extract_sample_id('/path/to/test_full.fastq.vcf'), 'test_full')

  def test_get_cds_info(self):
    self.assertEqual(get_cds_info(3, [(1, 5, 'GENE1', 'PRODUCT1')]), (True, ('GENE1', 'PRODUCT1')))
    self.assertEqual(get_cds_info(5, [(1, 5, 'GENE1', 'PRODUCT1')]), (True, ('GENE1', 'PRODUCT1')))
    self.assertEqual(get_cds_info(1, [(1, 5, 'GENE1', 'PRODUCT1')]), (True, ('GENE1', 'PRODUCT1')))
    self.assertEqual(get_cds_info(8, [(1, 5, 'GENE1', 'PRODUCT1')]), (False, (None, None)))
    self.assertEqual(get_cds_info(11, [(1, 5, 'GENE1', 'PRODUCT1'), (6, 10, 'GENE2', 'PRODUCT2')]), (False, (None, None)))
    self.assertEqual(get_cds_info(7, []), (False, (None, None)))

if __name__ == '__main__':
  unittest.main()