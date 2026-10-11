""" Test ConsensusXML setColumnHeaders """
import argparse
from test.integration import IntegrationTestCase
from msptools import cli


class TestConsensusXMLMADNormalizer(IntegrationTestCase):
    """ Test the msptools.cli.ConsensusXML.mad_normalizer module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.consensusXML')
        args.min_int = 1e-08
        args.virtual_reference = False
        args.virtual_reference_method = 'median'
        return args

    def test_mad_normalize_case1(self):
        """ Test ConsensusXML MADNormalizer case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        cli.ConsensusXML.mad_normalizer.main(args)

    def test_mad_normalize_case2_virtual_reference(self):
        """ Test ConsensusXML MADNormalizer case 2 with virtual reference """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        args.virtual_reference = True
        cli.ConsensusXML.mad_normalizer.main(args)
