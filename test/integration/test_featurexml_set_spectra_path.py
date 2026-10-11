""" Test FeatureXML setSpectraPath """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import FeatureMap, FeatureXMLFile
from msptools import cli


class TestIdXMLSetDatabasePath(IntegrationTestCase):
    """ Test the msptools::IdXML::setDatabasePath module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.featureXML')
        return args

    def test_set_database_path_case1(self):
        """ Test IdXML filterProteinGroups case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/test_featurexml_set_spectra_path.featureXML'
        args.spectra_path = 'SBP-0001'
        cli.FeatureXML.set_spectra_path.main(args)

        fmap = FeatureMap()
        FeatureXMLFile().load(args.output_path, fmap)
        ms_runs = []
        fmap.getPrimaryMSRunPath(ms_runs)
        ms_runs = {x.decode('utf-8') for x in ms_runs}
        self.assertEqual(ms_runs, {args.spectra_path})
