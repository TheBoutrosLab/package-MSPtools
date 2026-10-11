""" Test MzML filter """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import MSExperiment, MzMLFile, IdXMLFile
from msptools import cli


class TestMzMLFilter(IntegrationTestCase):
    """ Test the msptools::MzML::filter module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.consensus_file = None
        args.id_file = None
        args.rt_tolerance = 0.1
        args.mz_tolerance = 0.1
        args.mz_tolerance_unit = 'ppm'
        args.consensus_check_id = 'false'
        args.use_rt_mz = False
        args.output_path = str(self.work_dir/'output.mzML')
        return args

    def test_set_column_headers_case1(self):
        """ Test MzML filter case 1. This mzML has 10 spectra, 5 MS1 and 5 MS2.
        3 MS2 exist in the idXML and should be removed. """
        args = self.create_base_args()
        args.input_path = 'test/files/test_mzml_filter_1.mzML'
        args.id_file = 'test/files/test_mzml_filter_1.idXML'
        cli.MzML.filter.main(args)

        filtered_mzml = MSExperiment()
        MzMLFile().load(args.output_path, filtered_mzml)
        pep_ids = []
        prot_ids = []
        IdXMLFile().load(args.id_file, prot_ids, pep_ids)

        mzml_spectrum = {x.getNativeID() for x in filtered_mzml}
        idxml_spectrum = {x.getMetaValue('spectrum_reference') for x in pep_ids}

        self.assertTrue(not any(x in idxml_spectrum for x in mzml_spectrum))
