""" Test ConsensusXML MADNormalizer """
import argparse
from contextlib import redirect_stdout
from io import StringIO
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
        args.quasi_intensity = False
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

        for ref_intensity in [0, 1e-9, np.nan]:
            cf = ConsensusFeature()
            cf.setMetaValue('precursor_intensity', 1000)
            for index, intensity in [(0, ref_intensity), (1, 300)]:
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

        normalizer = MADNormalizer(cmap)
        normalized = normalizer.normalize()
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

        summary = normalizer.summary
        self.assertEqual(summary['ConsensusFeatures processed'], 7)
        self.assertEqual(summary['ConsensusFeatures normalized'], 3)
        self.assertEqual(summary['Valid physical references'], 3)
        self.assertEqual(summary['Missing reference handles'], 1)
        self.assertEqual(summary['Zero reference intensities'], 1)
        self.assertEqual(summary['Other references at or below min_int'], 1)
        self.assertEqual(summary['Non-finite reference intensities'], 1)

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
        output = StringIO()
        with redirect_stdout(output):
            cli.ConsensusXML.mad_normalizer.main(args)
        self.assertIn('MAD normalization tally', output.getvalue())
        self.assertRegex(output.getvalue(), r'ConsensusFeatures processed\s+5')
        self.assertRegex(output.getvalue(), r'ConsensusFeatures normalized\s+5')

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

    def test_quasi_intensity_matches_reference_scaled_ratios(self):
        """ The optional output uses each feature's own reference estimate. """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        cli.ConsensusXML.mad_normalizer.main(args)
        ratios = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, ratios)

        args.quasi_intensity = True
        cli.ConsensusXML.mad_normalizer.main(args)
        abundances = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, abundances)

        for ratio_cf, abundance_cf in zip(ratios, abundances):
            ref_int = ratio_cf.getMetaValue('reference_intensity')
            self.assertGreater(ref_int, 0)
            self.assertAlmostEqual(
                abundance_cf.getMetaValue('reference_intensity'), ref_int
            )
            ratio_handles = {f.getMapIndex(): f.getIntensity()
                             for f in ratio_cf.getFeatureList()}
            for handle in abundance_cf.getFeatureList():
                expected = ref_int * np.exp2(ratio_handles[handle.getMapIndex()])
                self.assertAlmostEqual(handle.getIntensity(), expected,
                                       delta=expected * 1e-5)

    def test_quasi_intensity_requires_valid_reference(self):
        """ Only measured zero samples with valid references become zero. """
        cmap = ConsensusMap()
        headers = {}
        for index, description in [(0, 'reference'), (1, 'sample')]:
            header = ColumnHeader()
            header.filename = 'plex'
            header.setMetaValue('channel_description', description)
            headers[index] = header
        cmap.setColumnHeaders(headers)
        for ref, sample in [
            (100, 100), (100, 200), (100, 400),
            (100, 0), (100, 1e-9), (0, 300)
        ]:
            cf = ConsensusFeature()
            cf.setMetaValue('precursor_intensity', 1000)
            for index, intensity in [(0, ref), (1, sample)]:
                feature = BaseFeature()
                feature.setIntensity(intensity)
                cf.insert(index, feature)
            cmap.push_back(cf)
        missing_ref = ConsensusFeature()
        missing_ref.setMetaValue('precursor_intensity', 1000)
        sample = BaseFeature()
        sample.setIntensity(300)
        missing_ref.insert(1, sample)
        cmap.push_back(missing_ref)

        ratios = MADNormalizer(cmap).normalize()
        for cf in list(ratios)[3:5]:
            sample_ratio = next(f.getIntensity() for f in cf.getFeatureList()
                                if f.getMapIndex() == 1)
            self.assertTrue(np.isnan(sample_ratio))

        result = MADNormalizer(cmap, quasi_intensity=True).normalize()
        self.assertTrue(all(np.isfinite(f.getIntensity()) for cf in list(result)[:3]
                            for f in cf.getFeatureList()))
        zero_sample = {f.getMapIndex(): f.getIntensity()
                       for f in list(result)[3].getFeatureList()}
        self.assertGreater(zero_sample[0], 0)
        self.assertEqual(zero_sample[1], 0)
        low_sample = {f.getMapIndex(): f.getIntensity()
                      for f in list(result)[4].getFeatureList()}
        self.assertTrue(np.isnan(low_sample[1]))
        for cf in list(result)[5:]:
            self.assertTrue(all(np.isnan(f.getIntensity())
                                for f in cf.getFeatureList()))

    def test_quasi_intensity_with_virtual_reference(self):
        """ Virtual references can also produce quasi intensities. """
        args = self.create_base_args()
        args.input_path = 'test/files/MADNormalizer_1_input.consensusXML'
        args.virtual_reference = True
        args.quasi_intensity = True
        cli.ConsensusXML.mad_normalizer.main(args)
        output = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, output)
        self.assertTrue(all(np.isfinite(f.getIntensity()) and f.getIntensity() > 0
                            for cf in output for f in cf.getFeatureList()))
