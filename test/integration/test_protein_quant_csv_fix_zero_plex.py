""" test ProteinQuantCSV fixZeroPlex """
import argparse
from test.integration import IntegrationTestCase
from msptools import cli


class TestProteinQuantCSVFixZeroPlex(IntegrationTestCase):
    """ Test ProteinQuantCSV fixZeroPlex """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.csv')
        args.consensusXML = None
        args.protein_graoup = None
        args.type = 'protein'
        args.separator = '\t'
        args.quoting = 'double'
        args.greedy_group_resolution = False
        args.nocheck_intensity = False
        return args

    def test_protein_quant_csv_fix_zero_plex_protein_case_1(self):
        """ Test ProteinQuantCSV fixZeroPlex """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_fix_zero_plex_protein_input_1.csv'
        args.type = 'protein'
        cli.ProteinQuantCSV.fix_zero_plex.main(args)

    def test_protein_quant_csv_fix_zero_plex_peptide_case_1(self):
        """ Test ProteinQuantCSV fixZeroPlex """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_fix_zero_plex_peptide_input_1.csv'
        args.type = 'peptide'
        cli.ProteinQuantCSV.fix_zero_plex.main(args)
