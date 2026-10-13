""" Test ConsensusXML MADNormalizer """
import argparse
import numpy as np
from pyopenms import BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, ConsensusXMLFile
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

    def test_absent_physical_reference_handle_stays_absent(self):
        """ An absent reference handle leaves measured sample ratios undefined. """
        cmap = ConsensusMap()
        headers = {}
        for index, description in [(0, 'reference'), (1, 'sample')]:
            header = ColumnHeader()
            header.filename = 'plex'
            header.setMetaValue('channel_description', description)
            headers[index] = header
        cmap.setColumnHeaders(headers)

        for sample_intensity in [100, 200, 400]:
            cf = ConsensusFeature()
            cf.setMetaValue('precursor_intensity', 1000)
            for index, intensity in [(0, 100), (1, sample_intensity)]:
                feature = BaseFeature()
                feature.setIntensity(intensity)
                cf.insert(index, feature)
            cmap.push_back(cf)

        missing_cf = ConsensusFeature()
        missing_cf.setMetaValue('precursor_intensity', 1000)
        sample = BaseFeature()
        sample.setIntensity(300)
        missing_cf.insert(1, sample)
        cmap.push_back(missing_cf)

        normalized = MADNormalizer(cmap).normalize()
        for cf in list(normalized)[:3]:
            intensities = {f.getMapIndex(): f.getIntensity()
                           for f in cf.getFeatureList()}
            self.assertEqual(intensities[0], 0)
            self.assertTrue(np.isfinite(intensities[1]))

        missing_output = list(normalized)[-1]
        features = {f.getMapIndex(): f.getIntensity()
                    for f in missing_output.getFeatureList()}
        self.assertEqual(set(features), {1})
        self.assertTrue(np.isnan(features[1]))
        self.assertEqual(missing_output.getMetaValue('reference_intensity'), 0)

    def test_sparse_channels_keep_reference_and_sample_positions(self):
        """ Absent handles must not shift ratios into another channel. """
        cmap = ConsensusMap()
        headers = {}
        for index in range(4):
            header = ColumnHeader()
            header.filename = 'plex'
            header.setMetaValue('channel_description',
                                'reference' if index == 2 else f'sample{index}')
            headers[index] = header
        cmap.setColumnHeaders(headers)

        rows = [
            {0: 100, 2: 200, 3: 500},
            {0: 160, 2: 200, 3: 400},
            {0: 130, 2: 200, 3: 600},
            {0: 190, 2: 200},
        ]
        for row in rows:
            cf = ConsensusFeature()
            cf.setMetaValue('precursor_intensity', 1000)
            for index, intensity in row.items():
                feature = BaseFeature()
                feature.setIntensity(intensity)
                cf.insert(index, feature)
            cmap.push_back(cf)

        normalized = MADNormalizer(cmap).normalize()
        sample_values = []
        for original, cf in zip(rows, normalized):
            intensities = {f.getMapIndex(): f.getIntensity()
                           for f in cf.getFeatureList()}
            self.assertEqual(set(intensities), set(original))
            self.assertEqual(intensities[2], 0)
            self.assertTrue(np.isfinite(intensities[0]))
            if 3 in original:
                self.assertTrue(np.isfinite(intensities[3]))
                sample_values.append(intensities[3])
        self.assertTrue(any(value != 0 for value in sample_values))

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
