""" Test ConsensusXML MADNormalizer """
import argparse
from contextlib import redirect_stdout
from io import StringIO
import numpy as np
from pyopenms import (
    AASequence, BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap,
    ConsensusXMLFile, PeptideHit, PeptideIdentification, ProteinIdentification
)
from test.integration import IntegrationTestCase
from msptools import cli
from msptools.algorithm.MADNormalizer import MADNormalizer


class TestConsensusXMLMADNormalizer(IntegrationTestCase):
    """ Test the msptools.cli.ConsensusXML.mad_normalizer module """
    @staticmethod
    def make_merged_map(layout='linked'):
        """ Two plexes and fractions with known bridge reference medians. """
        cmap = ConsensusMap()
        headers = {}
        for index, plex, description in [
            (0, 'A', 'reference'), (1, 'A', 'sample'),
            (2, 'B', 'reference'), (3, 'B', 'sample')
        ]:
            header = ColumnHeader()
            header.filename = plex
            header.label = str(index)
            header.setMetaValue('channel_description', description)
            headers[index] = header
        cmap.setColumnHeaders(headers)
        protein_run = ProteinIdentification()
        protein_run.setIdentifier('synthetic_run')
        protein_run.setSearchEngine('synthetic')
        cmap.setProteinIdentifications([protein_run])
        rows = [
            ({0: 100, 1: 150, 2: 200, 3: 300}, 'F1'),
            ({0: 120, 1: 240, 2: 240, 3: 360}, 'F1'),
            ({0: 80, 1: 40, 2: 160, 3: 80}, 'F2'),
            ({0: 90, 1: 135}, 'F2'),
            ({2: 180, 3: 270}, 'F2'),
        ]
        handle_id = 1000
        for feature_id, (values, fraction) in enumerate(rows, start=100):
            plex_rows = [values] if layout == 'linked' else [
                {index: intensity for index, intensity in values.items()
                 if index // 2 == plex_index}
                for plex_index in (0, 1)
            ]
            for plex_values in plex_rows:
                if not plex_values:
                    continue
                cf = ConsensusFeature()
                cf.setUniqueId(feature_id if layout == 'linked' else feature_id * 10
                               + min(plex_values) // 2)
                cf.setMetaValue('tmt_fraction_id', fraction)
                cf.setMetaValue('custom_metadata', 'preserve_me')
                if layout == 'scan':
                    cf.setMetaValue('precursor_intensity', 1000)
                else:
                    for plex in ('A', 'B'):
                        area = sum(intensity for index, intensity in plex_values.items()
                                   if index // 2 == ('A', 'B').index(plex))
                        if area:
                            cf.setMetaValue(f'tmt_ms1_area::{plex}', area)
                hit = PeptideHit()
                hit.setSequence(AASequence.fromString('PEPTIDEK'))
                identification = PeptideIdentification()
                identification.setIdentifier('synthetic_run')
                identification.setHits([hit])
                cf.setPeptideIdentifications([identification])
                for index, intensity in plex_values.items():
                    feature = BaseFeature()
                    feature.setUniqueId(handle_id)
                    feature.setIntensity(intensity)
                    cf.insert(index, feature)
                    handle_id += 1
                cmap.push_back(cf)
        return cmap

    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.consensusXML')
        args.min_int = 1e-08
        args.virtual_reference = False
        args.virtual_reference_method = 'median'
        args.quasi_intensity = False
        args.layout = 'scan'
        args.bridge_scale = False
        args.allow_virtual_bridge = False
        args.scaling_report_path = None
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

    def test_merged_scan_map_uses_plex_specific_ratios(self):
        """ Scan features across plexes share MAD statistics without mixing handles. """
        cmap = self.make_merged_map('scan')
        normalizer = MADNormalizer(cmap)
        output = normalizer.normalize()
        self.assertEqual(output.size(), cmap.size())
        sample_ratios = {'A': [], 'B': []}
        for original, normalized in zip(cmap, output):
            self.assertEqual(original.getUniqueId(), normalized.getUniqueId())
            self.assertEqual(original.getMetaValue('custom_metadata'),
                             normalized.getMetaValue('custom_metadata'))
            self.assertEqual(len(normalized.getPeptideIdentifications()), 1)
            raw = {h.getMapIndex(): h for h in original.getFeatureList()}
            values = {h.getMapIndex(): h for h in normalized.getFeatureList()}
            self.assertEqual(set(raw), set(values))
            for index in raw:
                self.assertEqual(raw[index].getUniqueId(), values[index].getUniqueId())
            ref = next(index for index in values if index in (0, 2))
            self.assertEqual(values[ref].getIntensity(), 0)
            sample = ref + 1
            sample_ratios['A' if ref == 0 else 'B'].append(values[sample].getIntensity())
        self.assertEqual(normalizer.summary['ConsensusFeatures processed'], cmap.size())
        self.assertAlmostEqual(np.median(sample_ratios['A']),
                               np.median(sample_ratios['B']), places=5)

    def test_linked_bridge_scales_areas_and_preserves_ids(self):
        """ Per-plex reference medians set symmetric factors for every feature. """
        cmap = self.make_merged_map()
        ratios = MADNormalizer(cmap, layout='linked').normalize()
        normalizer = MADNormalizer(cmap, layout='linked', quasi_intensity=True,
                                   bridge_scale=True)
        output = normalizer.normalize()
        self.assertEqual([row['plex_scale_factor'] for row in normalizer.scaling_report],
                         [1.5, 0.75])
        self.assertEqual([row['reference_feature_count'] for row in normalizer.scaling_report],
                         [4, 4])
        for source, ratio_cf, result in zip(cmap, ratios, output):
            self.assertEqual(source.getUniqueId(), result.getUniqueId())
            self.assertEqual(source.getMetaValue('tmt_fraction_id'),
                             result.getMetaValue('tmt_fraction_id'))
            for plex in ('A', 'B'):
                key = f'tmt_ms1_area::{plex}'
                if source.metaValueExists(key):
                    self.assertEqual(source.getMetaValue(key), result.getMetaValue(key))
            self.assertEqual(len(result.getPeptideIdentifications()), 1)
            raw = {h.getMapIndex(): h for h in source.getFeatureList()}
            ratio_values = {h.getMapIndex(): h.getIntensity()
                            for h in ratio_cf.getFeatureList()}
            for handle in result.getFeatureList():
                index = handle.getMapIndex()
                plex = 'A' if index < 2 else 'B'
                factor = 1.5 if plex == 'A' else 0.75
                ref_index = 0 if plex == 'A' else 2
                expected = raw[ref_index].getIntensity() * np.exp2(ratio_values[index]) * factor
                self.assertAlmostEqual(handle.getIntensity(), expected,
                                       delta=expected * 1e-5)
                self.assertEqual(handle.getUniqueId(), raw[index].getUniqueId())
        path = self.work_dir/'linked.consensusXML'
        ConsensusXMLFile().store(str(path), output)
        reloaded = ConsensusMap()
        ConsensusXMLFile().load(str(path), reloaded)
        self.assertEqual(
            [h.getUniqueId() for cf in reloaded for h in cf.getFeatureList()],
            [h.getUniqueId() for cf in cmap for h in cf.getFeatureList()]
        )
        self.assertEqual([cf.getUniqueId() for cf in reloaded],
                         [cf.getUniqueId() for cf in cmap])
        self.assertTrue(all(len(cf.getPeptideIdentifications()) == 1
                            and cf.getMetaValue('custom_metadata') == 'preserve_me'
                            for cf in reloaded))

    def test_bridge_restrictions_and_missing_plex_reference(self):
        cmap = self.make_merged_map()
        with self.assertRaisesRegex(ValueError, 'requires linked layout'):
            MADNormalizer(cmap, quasi_intensity=True, bridge_scale=True).normalize()
        with self.assertRaisesRegex(ValueError, 'requires linked layout'):
            MADNormalizer(cmap, layout='linked', bridge_scale=True).normalize()
        with self.assertRaisesRegex(ValueError, 'requires explicit opt-in'):
            MADNormalizer(cmap, layout='linked', virtual_reference=True,
                          quasi_intensity=True, bridge_scale=True).normalize()
        cmap_only_unmatched = ConsensusMap()
        cmap_only_unmatched.setColumnHeaders(cmap.getColumnHeaders())
        for cf in list(cmap)[3:]:
            cmap_only_unmatched.push_back(cf)
        normalizer = MADNormalizer(cmap_only_unmatched, layout='linked',
                                   quasi_intensity=True, bridge_scale=True)
        normalizer.normalize()
        self.assertEqual([row['reference_feature_count'] for row in normalizer.scaling_report],
                         [1, 1])

        headers = cmap_only_unmatched.getColumnHeaders()
        for index, description in [(4, 'reference'), (5, 'sample')]:
            header = ColumnHeader()
            header.filename = 'C'
            header.label = str(index)
            header.setMetaValue('channel_description', description)
            headers[index] = header
        cmap_only_unmatched.setColumnHeaders(headers)
        with self.assertRaisesRegex(ValueError, 'no valid reference area for plex C'):
            MADNormalizer(cmap_only_unmatched, layout='linked',
                          quasi_intensity=True, bridge_scale=True).normalize()

    def test_bridge_uses_partial_overlaps_across_three_plexes(self):
        """ A-B and B-C overlaps suffice without an A-B-C feature. """
        cmap = ConsensusMap()
        headers = {}
        for index, plex in enumerate(('A', 'A', 'B', 'B', 'C', 'C')):
            header = ColumnHeader()
            header.filename = plex
            header.label = str(index)
            header.setMetaValue('channel_description',
                                'reference' if index % 2 == 0 else 'sample')
            headers[index] = header
        cmap.setColumnHeaders(headers)
        rows = [
            {0: 100, 1: 150, 2: 200, 3: 300},
            {0: 120, 1: 240, 2: 240, 3: 360},
            {2: 180, 3: 270, 4: 360, 5: 540},
            {2: 160, 3: 80, 4: 320, 5: 160},
        ]
        for row in rows:
            cf = ConsensusFeature()
            for plex, first_index in [('A', 0), ('B', 2), ('C', 4)]:
                if first_index in row:
                    cf.setMetaValue(f'tmt_ms1_area::{plex}',
                                    row[first_index] + row[first_index + 1])
            for index, intensity in row.items():
                handle = BaseFeature()
                handle.setIntensity(intensity)
                cf.insert(index, handle)
            cmap.push_back(cf)
        normalizer = MADNormalizer(cmap, layout='linked', quasi_intensity=True,
                                   bridge_scale=True)
        output = normalizer.normalize()
        self.assertEqual(output.size(), len(rows))
        self.assertEqual([item['reference_feature_count']
                          for item in normalizer.scaling_report], [2, 4, 2])
        self.assertEqual([item['median_reference_area']
                          for item in normalizer.scaling_report], [110, 190, 340])
        self.assertTrue(np.allclose(
            [item['plex_scale_factor'] for item in normalizer.scaling_report],
            [190 / 110, 1, 190 / 340]
        ))

    def test_virtual_bridge_requires_and_respects_opt_in(self):
        """ Virtual bridge uses per-plex medians only when explicitly allowed. """
        cmap = self.make_merged_map()
        normalizer = MADNormalizer(cmap, layout='linked', virtual_reference=True,
                                   quasi_intensity=True, bridge_scale=True,
                                   allow_virtual_bridge=True,
                                   scaling_report_path=str(self.work_dir/'scaling.tsv'))
        output = normalizer.normalize()
        self.assertEqual(output.size(), cmap.size())
        self.assertEqual([row['plex_scale_factor'] for row in normalizer.scaling_report],
                         [1.5, 0.75])
        self.assertTrue((self.work_dir/'scaling.tsv').read_text().startswith(
            'plex_id\treference_feature_count\tmedian_reference_area'))

    def test_linked_missing_and_invalid_references_are_per_plex(self):
        """ An invalid A reference leaves the valid B group usable. """
        cmap = self.make_merged_map()
        examples = [({1: 300, 2: 200, 3: 300}, 'missing'),
                    ({0: 0, 1: 300, 2: 200, 3: 300}, 'zero'),
                    ({0: np.nan, 1: 300, 2: 200, 3: 300}, 'nonfinite')]
        for index, (values, _) in enumerate(examples, start=200):
            cf = ConsensusFeature()
            cf.setUniqueId(index)
            cf.setMetaValue('tmt_ms1_area::A',
                            sum(value for map_index, value in values.items()
                                if map_index < 2 and np.isfinite(value)))
            cf.setMetaValue('tmt_ms1_area::B',
                            sum(value for map_index, value in values.items()
                                if map_index >= 2 and np.isfinite(value)))
            for map_index, intensity in values.items():
                handle = BaseFeature()
                handle.setUniqueId(index * 10 + map_index)
                handle.setIntensity(intensity)
                cf.insert(map_index, handle)
            cmap.push_back(cf)
        normalizer = MADNormalizer(cmap, layout='linked')
        output = normalizer.normalize()
        for cf in list(output)[-3:]:
            values = {h.getMapIndex(): h.getIntensity() for h in cf.getFeatureList()}
            self.assertTrue(np.isnan(values[1]))
            self.assertEqual(values[2], 0)
            self.assertTrue(np.isfinite(values[3]))
        self.assertEqual(normalizer.summary['Missing reference handles'], 1)
        self.assertEqual(normalizer.summary['Zero reference intensities'], 1)
        self.assertEqual(normalizer.summary['Non-finite reference intensities'], 1)

    def test_zero_mad_uses_median_centering(self):
        """ A constant sample channel remains finite and is reported. """
        cmap = ConsensusMap()
        headers = {}
        for index, description in [(0, 'reference'), (1, 'sample')]:
            header = ColumnHeader()
            header.filename = 'A'
            header.label = str(index)
            header.setMetaValue('channel_description', description)
            headers[index] = header
        cmap.setColumnHeaders(headers)
        for index in range(3):
            cf = ConsensusFeature()
            cf.setMetaValue('precursor_intensity', 1000)
            for map_index, intensity in [(0, 100), (1, 200)]:
                handle = BaseFeature()
                handle.setIntensity(intensity)
                cf.insert(map_index, handle)
            cmap.push_back(cf)
        normalizer = MADNormalizer(cmap)
        output = normalizer.normalize()
        self.assertEqual(normalizer.fallback_channels, [('A', 1)])
        self.assertTrue(all(np.isclose(next(h.getIntensity() for h in cf.getFeatureList()
                                           if h.getMapIndex() == 1), 1)
                            for cf in output))

    def test_merged_map_rejects_ambiguous_reference_and_unknown_index(self):
        cmap = self.make_merged_map()
        headers = cmap.getColumnHeaders()
        headers[1].setMetaValue('channel_description', 'reference')
        cmap.setColumnHeaders(headers)
        with self.assertRaisesRegex(ValueError, 'plex A; found 2'):
            MADNormalizer(cmap, layout='linked').normalize()

        cmap = self.make_merged_map()
        invalid = ConsensusFeature()
        handle = BaseFeature()
        handle.setIntensity(10)
        invalid.insert(99, handle)
        cmap.push_back(invalid)
        with self.assertRaisesRegex(ValueError, 'Unknown column map index 99'):
            MADNormalizer(cmap, layout='linked').normalize()

    def test_linked_cli_writes_scaling_report_and_preserves_parent_ids(self):
        args = self.create_base_args()
        args.input_path = str(self.work_dir/'linked-input.consensusXML')
        args.layout = 'linked'
        args.quasi_intensity = True
        args.bridge_scale = True
        args.scaling_report_path = str(self.work_dir/'scaling.tsv')
        cmap = self.make_merged_map()
        ConsensusXMLFile().store(args.input_path, cmap)
        cli.ConsensusXML.mad_normalizer.main(args)

        output = ConsensusMap()
        ConsensusXMLFile().load(args.output_path, output)
        self.assertEqual([cf.getUniqueId() for cf in output],
                         [cf.getUniqueId() for cf in cmap])
        self.assertEqual(len(output.getColumnHeaders()), 4)
        self.assertEqual(len((self.work_dir/'scaling.tsv').read_text().splitlines()), 3)

    def test_linked_layout_requires_allocated_area_metadata(self):
        """ A scan map must not silently be interpreted as allocated areas. """
        cmap = self.make_merged_map('scan')
        with self.assertRaisesRegex(ValueError, 'lacks tmt_ms1_area::A'):
            MADNormalizer(cmap, layout='linked', quasi_intensity=True).normalize()

        linked = self.make_merged_map()
        first = linked[0]
        first.setMetaValue('tmt_ms1_area::A', np.nan)
        linked[0] = first
        with self.assertRaisesRegex(ValueError, 'invalid tmt_ms1_area::A'):
            MADNormalizer(linked, layout='linked').normalize()

        first.setMetaValue('tmt_ms1_area::A', 5)
        linked[0] = first
        expected = MADNormalizer(self.make_merged_map(), layout='linked').normalize()
        output = MADNormalizer(linked, layout='linked').normalize()
        self.assertEqual(output[0].getMetaValue('tmt_ms1_area::A'), 5)
        self.assertEqual([h.getIntensity() for h in output[0].getFeatureList()],
                         [h.getIntensity() for h in expected[0].getFeatureList()])

    def test_scaling_report_cannot_replace_consensusxml(self):
        """ Reject both path collisions before opening either output. """
        args = self.create_base_args()
        args.input_path = str(self.work_dir/'linked-input.consensusXML')
        args.layout = 'linked'
        args.quasi_intensity = True
        args.bridge_scale = True
        ConsensusXMLFile().store(args.input_path, self.make_merged_map())
        original = (self.work_dir/'linked-input.consensusXML').read_bytes()

        for report_path in (args.input_path, args.output_path):
            args.scaling_report_path = report_path
            with self.subTest(report_path=report_path):
                with self.assertRaisesRegex(ValueError,
                        'Scaling report path must differ'):
                    cli.ConsensusXML.mad_normalizer.main(args)
                self.assertEqual((self.work_dir/'linked-input.consensusXML').read_bytes(),
                                 original)
                self.assertFalse((self.work_dir/'output.consensusXML').exists())
