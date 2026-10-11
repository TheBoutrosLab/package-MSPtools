""" Test IdXML FilterProteinGroups """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import IdXMLFile
from msptools import cli


class TestIdXMLSetDatabasePath(IntegrationTestCase):
    """ Test the msptools::IdXML::setDatabasePath module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.idXML')
        args.database = None
        return args

    def test_set_database_path_case1(self):
        """ Test IdXML filterProteinGroups case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/test_idxml_set_database_path.idXML'
        args.database = 'test/files/test_idxml_set_database_path_database.fa'
        cli.IdXML.set_database_path.main(args)

        peptides = []
        proteins = []
        IdXMLFile().load(args.output_path, proteins, peptides)
        self.assertEqual(proteins[0].getSearchParameters().db, args.database)
