import unittest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from parsers import extract_codon_frequencies

class TestParsesr(unittest.TestCase):
  def test_extract_codon_frequencies(self):
    curr_output = extract_codon_frequencies('./test_data/coreF44.fastq.gz.bam', 6747)
    self.assertEqual(curr_output['GCT'], 0.4972)
    self.assertEqual(curr_output['ACT'], 0.3437)
    self.assertEqual(curr_output['GTT'], 0.1075)