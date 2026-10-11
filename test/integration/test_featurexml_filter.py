""" Test FeatureXML filter """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import FeatureXMLFile, FeatureMap
from msptools import cli


class TestFeataureXMLFilter(IntegrationTestCase):
    """ Test the msptools::FeatureXML::filter module """
    def test_filter(self):
        """ Test filter """
        args = argparse.Namespace()
        args.input_path = 'test/files/test_featurexml_filter_1.featureXML'
        args.output_path = str(self.work_dir/'output.featureXML')
        cli.FeatureXML.featurexml_filter.main(args)

        filtered_fmap = FeatureMap()
        FeatureXMLFile().load(args.output_path, filtered_fmap)

        self.assertTrue(not any(x.getPeptideIdentifications() for x in filtered_fmap))
        self.assertEqual(filtered_fmap.size(), 5)
