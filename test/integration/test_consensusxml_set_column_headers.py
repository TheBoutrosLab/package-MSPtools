""" Test ConsensusXML setColumnHeaders """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import ConsensusMap, ConsensusXMLFile
from msptools import cli


class TestConsensusXMLSetColumnHeaders(IntegrationTestCase):
    """ Test the msptools::ConsensusXML::setDatabasePath module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.consensusXML')
        return args

    def test_set_column_headers_case1(self):
        """ Test IdXML filterProteinGroups case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/test_consensusxml_set_column_headers_1.consensusXML'
        args.column_headers = ['SBP-0001']
        cli.ConsensusXML.set_column_headers.main(args)

        cmap = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, cmap)
        headers = [x.filename for x in cmap.getColumnHeaders().values()]
        self.assertEqual(headers, args.column_headers)
