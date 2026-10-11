""" Unit test cases for ProteinQuantCSV """
from unittest import TestCase
from msptools.ProteinQuantCSV.ProteinQuantCSV import ProteinQuantCSV

class TestProteinQuantCSV(TestCase):
    """ Test cases for ProteinQuantCSV """
    def test_parse_header_map_protein(self):
        """ Ensure the header map is read correctly from protein csv. """
        with open('test/files/test_protein_quant_csv_reheader_csv_2.csv', 'rt') as handle:
            peptide_quant = ProteinQuantCSV.read(handle, 'protein')
        header = [
            'abundance_1',
            'abundance_2'
        ]
        for record in peptide_quant.records:
            self.assertEqual(record.header, header)

    def test_parse_header_map_peptides(self):
        """ Ensure the header map is read correctly from peptide csv. """
        with open('test/files/test_protein_quant_csv_reheader_csv_3.csv', 'rt') as handle:
            peptide_quant = ProteinQuantCSV.read(handle, 'peptide')
        header = [
            'abundance_1',
            'abundance_2'
        ]
        for record in peptide_quant.records:
            self.assertEqual(record.header, header)
