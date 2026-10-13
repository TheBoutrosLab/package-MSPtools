""" Test ConsensusXML MADNormalizer """
import argparse
import numpy as np
from pyopenms import ConsensusMap, ConsensusXMLFile
from test.integration import IntegrationTestCase
from msptools import cli
from msptools.algorithm.MADNormalizer import MADNormalizer


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

    def test_missing_physical_reference_stays_missing(self):
        """ Reference presence is determined by its intensity, not a ratio value. """
        ratios = {'plex': np.array([[0., 1., 2.], [0., 2., 3.], [0., 3., 4.]])}
        normalized = MADNormalizer.normalize_mad(
            ratios, {'plex': 0}, {'plex': [True, False, True]}
        )['plex']
        self.assertEqual(normalized[0, 0], 0)
        self.assertTrue(np.isnan(normalized[1, 0]))
        self.assertEqual(normalized[2, 0], 0)

    def test_mad_normalize_case1(self):
        """ Test ConsensusXML MADNormalizer case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        cli.ConsensusXML.mad_normalizer.main(args)

        original = ConsensusMap()
        normalized = ConsensusMap()
        ConsensusXMLFile().load(args.input_path, original)
        ConsensusXMLFile().load(args.output_path, normalized)
        raw = np.array([[f.getIntensity() for f in cf.getFeatureList()]
                        for cf in original])
        result = np.array([[f.getIntensity() for f in cf.getFeatureList()]
                           for cf in normalized])

        self.assertTrue(np.all(result[:, 0] == 0))
        self.assertTrue(np.all(np.isfinite(result[:, 1:])))
        # The global center must use the sample channels, not the reference.
        sample_medians = np.median(np.log2(raw[:, 1:] / raw[:, [0]]), axis=0)
        self.assertTrue(np.allclose(np.median(result[:, 1:], axis=0),
                                    np.median(sample_medians), atol=1e-5))

    def test_mad_normalize_case2_virtual_reference(self):
        """ Test ConsensusXML MADNormalizer case 2 with virtual reference """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        args.virtual_reference = True
        cli.ConsensusXML.mad_normalizer.main(args)
        normalized = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, normalized)
        self.assertTrue(all(np.isfinite(f.getIntensity()) for cf in normalized
                            for f in cf.getFeatureList()))
