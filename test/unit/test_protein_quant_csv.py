""" Unit test cases for ProteinQuantCSV """
from unittest import TestCase
from msptools.ProteinQuantCSV.ProteinQuantCSV import ProteinQuantCSV
from msptools.ProteinQuantCSV.ProteinQuantCSVMetadata import ProteinQuantCSVMetadata

class TestProteinQuantCSV(TestCase):
    """ Test cases for ProteinQuantCSV """
    def test_protein_metadata_origin_prefix(self):
        """ Preserve the protein abundance prefix when rewriting metadata. """
        with open('test/files/test_protein_quant_csv_reheader_csv_1.csv', 'rt') as handle:
            metadata = ProteinQuantCSVMetadata.parse(handle)
        self.assertEqual(metadata.to_lines()[0],
            '# Protein abundances computed from file CPTAC_LUAD_01_merged.consensusXML'
        )

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
