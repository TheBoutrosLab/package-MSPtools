""" Test ProteinQuantCSVReheader """
import argparse
from os import PathLike
from typing import List
from test.integration import IntegrationTestCase
from msptools import cli


class TestProteinQuantCSVReheader(IntegrationTestCase):
    """ Test ProteinQuantCSVReheader """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.csv')
        args.consensusXML = None
        args.separator = '\t'
        args.quoting = 'double'
        args.header = None
        return args

    def assert_header_matches(self, file:PathLike, header:List[str]):
        """ Assert the header of a protein quant CSV file matches. """
        new_headers = None
        with open(file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                new_headers = [x.strip('"') for x in fields[4:]]
                if new_headers[-1] == 'fraction':
                    new_headers = new_headers[:-1]
                break
        self.assertTrue(new_headers is not None)
        self.assertTrue(new_headers == header)

    def test_protein_quant_csv_reheader_case_1(self):
        """ Test ProteinQuantCSVReheader case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_reheader_csv_1.csv'
        args.consensusXML = \
            'test/files/test_protein_quant_csv_reheader_consensus_xml_1.consensusXML'
        args.experiment_type = 'TMT'
        cli.ProteinQuantCSV.reheader.main(args)
        headers_expected = [
            'reference',
            *[f"CPT000{x}" for x in range(1,10)]
        ]
        self.assert_header_matches(args.output_path, headers_expected)

    def test_protein_quant_csv_reheader_case_2(self):
        """ Test ProteinQuantCSVReheader case 2 with LFQ """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_reheader_csv_2.csv'
        args.consensusXML = None
        args.experiment_type = 'LFQ'
        cli.ProteinQuantCSV.reheader.main(args)
        headers_expected = ['Liver-2', 'Liver-1']
        self.assert_header_matches(args.output_path, headers_expected)

    def test_protein_quant_csv_reheader_case_3(self):
        """ Test ProteinQuantCSVReheader case 3 with LFQ peptide quant """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_reheader_csv_3.csv'
        args.consensusXML = None
        args.experiment_type = 'LFQ'
        cli.ProteinQuantCSV.reheader.main(args)
        headers_expected = ['Liver-2', 'Liver-1']
        self.assert_header_matches(args.output_path, headers_expected)

    def test_protein_quant_csv_reheader_case_4(self):
        """ Test case with only one column """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_reheader_csv_4.csv'
        args.consensusXML = None
        args.experiment_type = 'LFQ'
        args.header = ['CPCG0001']
        cli.ProteinQuantCSV.reheader.main(args)
        headers_expected = args.header
        self.assert_header_matches(args.output_path, headers_expected)

    def test_protein_quant_csv_reheader_fail_header_no_match(self):
        """ Test that the reheader will fail because the headers provided dont match """
        args = self.create_base_args()
        args.input_path = 'test/files/test_protein_quant_csv_reheader_csv_4.csv'
        args.consensusXML = None
        args.experiment_type = 'LFQ'

        args.header = ['CPCG0001', 'CPCG0002']
        with self.assertRaises(ValueError):
            cli.ProteinQuantCSV.reheader.main(args)
