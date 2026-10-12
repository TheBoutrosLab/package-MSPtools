""" Unit test cases for ProteinQuantCSV """
from unittest import TestCase
from io import StringIO
import numpy as np
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

    def test_openms_abundance_headers_define_plexes(self):
        """ New OpenMS CSVs encode the plex and channel in each column name. """
        source = StringIO(
            '# Peptide abundances computed from file \'merged.consensusXML\'\n'
            '# Parameters (relevant only): best_charge_and_fraction=true\n'
            '"peptide"\t"protein"\t"n_proteins"\t"charge"\t'
            '"abundance|fraction_1|ch1"\t"abundance|fraction_1|ch2"\t'
            '"abundance|fraction_2|ch1"\t"abundance|fraction_2|ch2"\t"fraction"\n'
            '"PEPTIDE"\t"PROTEIN"\t1\t2\t0\t0\t5\t0\t1\n'
        )
        quant = ProteinQuantCSV.read(source, 'peptide')
        self.assertEqual(quant.metadata.get_plex_2_column(-1),
            {'fraction_1': [0, 1], 'fraction_2': [2, 3]})
        quant.remove_zero_plex()
        self.assertTrue(np.isnan(quant.records[0].abundance[0]))
        self.assertTrue(np.isnan(quant.records[0].abundance[1]))
        self.assertEqual(quant.records[0].abundance[2:], [5.0, 0.0])

    def test_openms_protein_headers_without_metadata(self):
        source = StringIO(
            '# Protein abundances computed from file \'merged.consensusXML\'\n'
            '"protein"\t"n_proteins"\t"protein_score"\t"n_peptides"\t'
            '"abundance|sample|ch1"\t"abundance|sample|ch2"\n'
            '"P1"\t1\t0.9\t2\t0\t0\n'
        )
        quant = ProteinQuantCSV.read(source, 'protein')
        self.assertEqual(quant.metadata.get_plex_2_column(-1), {'sample': [0, 1]})
        quant.remove_zero_plex()
        self.assertTrue(all(np.isnan(x) for x in quant.records[0].abundance))
