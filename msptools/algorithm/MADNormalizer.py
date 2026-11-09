""" MAD Normalizer """
from __future__ import annotations
from typing import List, Dict, Tuple
import statistics
import csv
import numpy as np
from pyopenms import ConsensusMap, ConsensusFeature, FeatureHandle, BaseFeature


VIRTUAL_REFERENCE_AVEAGE_FUNS = {
    'mean': statistics.mean,
    'median': statistics.median
}

class MADNormalizer():
    """ MADNormalizer performs cross-plex normalization using a TMT-integrator
    like algorithm.

    For i in 1..n samples, j in 1..p features:

    1. R_ij = log2(int_ij/int_kj)

    int_ij is the intensity of feature j in sample i and int_kj is the intensity
    of feature j in the reference channel of plex k.

    2. M_i = median(R_ij for j in 1..p)

    M_i is the sample median ratio.

    3. M_0 = meidan(M_i for i in 1..n)

    M_0 is the global meidan ratio.

    4. RC_ij = R_ij - M_i

    RC_ij is the centered ratio of feature j in sample i.

    5. MAD_i = median(abs(RC_ij) for j in 1..p)

    Sample MAD.

    6. MAD_0 = median(MAD_i for i in 1..n)

    Global MAD.

    7. RN_ij = RC_ij / MAD_i x MAD_0 + M_0

    RN_ij is the normalized ratio of feature j in sample i.

    8. REF_kj = MS1Int_kj x MS2Int_k_ref / MS2Int_k_sum

    MS1Int_kj is the MS1 (precursor) intensity of feature j in plex k. MS2Int_k_ref
    is the MS2 intensity of the reference channel in plex k for feature j, and
    MS2_k_sum is the sum of MS2 intensity of all channels in plex k for feature j.

    The normalized ratios (RN_ij) are the default normalization result. For
    scan maps, REF_kj is allocated from precursor intensity as above. For
    linked maps, reporter handles already contain allocated areas, so REF_kj
    is the reference handle area. Each output feature stores its per-plex
    reference areas as metadata. Optionally, output handles contain the quasi
    intensity REF_kj * 2**RN_ij instead of the log2 ratio.
    """
    def __init__(self, cmap:ConsensusMap, virtual_reference:bool=False,
            virtual_reference_method:str='median', min_int:float=1e-08,
            column_2_ref:Dict[int,int]=None, column_2_plex:Dict[int, str]=None,
            plex_2_column:Dict[str,List[int]]=None,
            quasi_intensity:bool=False, layout:str='scan', bridge_scale:bool=False,
            allow_virtual_bridge:bool=False, scaling_report_path:str=None):
        """ constructor """
        self.cmap = cmap
        # column index -> reference channel
        self.column_2_ref = column_2_ref or {}
        # column index -> plex
        self.column_2_plex = column_2_plex or {}
        # plex -> column indices
        self.plex_2_columns = plex_2_column or {}
        self.isobaric_type = str
        self.virtual_reference = virtual_reference
        self.virtual_reference_method = virtual_reference_method
        self.min_int = min_int
        self.quasi_intensity = quasi_intensity
        self.layout = layout
        self.bridge_scale = bridge_scale
        self.allow_virtual_bridge = allow_virtual_bridge
        self.scaling_report_path = scaling_report_path
        self.summary = {}
        self.scaling_report = []
        self.fallback_channels = []

    def get_isobaric_type(self) -> str:
        """ Get isobaric experiment type """
        isobaric_types = set()
        for dp in self.cmap.getDataProcessing():
            if dp.getSoftware().getName() == 'IsobaricAnalyzer':
                isobaric_types.add(dp.getMetaValue('parameter: type'))
        if len(isobaric_types) > 1:
            raise ValueError(
                f"Inconsistent isobaric type found: {isobaric_types}"
            )
        return list(isobaric_types)[0]

    def set_channel_plex(self):
        """ Set reference channels """
        # ConsensusMap column index -> plex
        column_2_plex:Dict[int, str] = {}
        # plex -> matching reference columns
        plex_2_refs:Dict[str, List[int]] = {}
        # plex -> ConsensusMap column indices
        plex_2_columns:Dict[str, List[int]] = {}
        for i, column_header in self.cmap.getColumnHeaders().items():
            plex = column_header.filename
            column_2_plex[i] = plex

            if plex not in plex_2_columns:
                plex_2_columns[plex] = []
            plex_2_columns[plex].append(i)

            if not self.virtual_reference \
                    and column_header.getMetaValue('channel_description') == 'reference':
                plex_2_refs.setdefault(plex, []).append(i)

        column_2_ref = {}
        if not self.virtual_reference:
            plex_2_ref:Dict[str, int] = {}
            for plex in plex_2_columns:
                matches = plex_2_refs.get(plex, [])
                if len(matches) != 1:
                    raise ValueError(
                        f"Expected one physical reference for plex {plex}; "
                        f"found {len(matches)}"
                    )
                plex_2_ref[plex] = matches[0]
            for i,plex in column_2_plex.items():
                column_2_ref[i] = plex_2_ref[plex]

        self.column_2_ref = column_2_ref
        self.column_2_plex = column_2_plex
        self.plex_2_columns = plex_2_columns

    def get_median_ratios(self, ratios, override:bool=False):
        """ Get median ratio """
        median_ratios = [
            x.getMetaValue('median_ratio') for x in self.cmap.getColumnHeaders().values()
        ]
        if any(x is None for x in median_ratios) or override:
            return np.median(ratios, axis = 0)
        return np.array(median_ratios, dtype=np.float64)

    def group_handles(self, cf:ConsensusFeature) -> Dict[str, List[FeatureHandle]]:
        """ Partition one feature's handles by their column-header plex. """
        groups:Dict[str, List[FeatureHandle]] = {}
        for handle in cf.getFeatureList():
            index = handle.getMapIndex()
            if index not in self.column_2_plex:
                raise ValueError(f'Unknown column map index {index} in consensus feature '
                                 f'{cf.getUniqueId()}')
            groups.setdefault(self.column_2_plex[index], []).append(handle)
        if self.layout == 'scan' and len(groups) > 1:
            raise ValueError(f'Scan-layout consensus feature {cf.getUniqueId()} '
                             f'contains multiple plexes: {list(groups)}')
        return groups

    def get_plex_and_ref_intensity(self, cf:ConsensusFeature) -> Tuple[str, float]:
        """ Return the reference reporter for a single-plex consensus feature. """
        groups = self.group_handles(cf)
        if len(groups) != 1:
            raise ValueError('Use per-plex groups for a multi-plex consensus feature')
        plex, handles = next(iter(groups.items()))
        return plex, self.get_reference_intensity(plex, handles)

    def get_reference_intensity(self, plex:str, handles:List[FeatureHandle]) -> float:
        """ Return a physical or virtual reference within one plex. """
        if self.virtual_reference:
            values = [h.getIntensity() for h in handles
                      if np.isfinite(h.getIntensity()) and h.getIntensity() > self.min_int]
            if len(values) <= 1:
                return 0
            return VIRTUAL_REFERENCE_AVEAGE_FUNS[self.virtual_reference_method](values)
        ref_index = self.column_2_ref[self.plex_2_columns[plex][0]]
        ref_feature = next((h for h in handles if h.getMapIndex() == ref_index), None)
        if ref_feature is None:
            return 0
        intensity = ref_feature.getIntensity()
        return 0 if intensity < self.min_int else intensity

    def collect_ratios(self):
        """ Collect log2 ratios and reference areas for each feature/plex group. """
        ratios:Dict[str,List[List[float]]] = {}
        ref_intensities:Dict[Tuple[int,str],float] = {}
        index_map:Dict[Tuple[int,str],int] = {}
        reference_validity:Dict[str,List[bool]] = {}
        column_positions = {
            plex: {column: pos for pos, column in enumerate(columns)}
            for plex, columns in self.plex_2_columns.items()
        }
        cf: ConsensusFeature
        for feature_index, cf in enumerate(self.cmap):
            self.summary['ConsensusFeatures processed'] += 1
            for plex, handles in self.group_handles(cf).items():
                if self.layout == 'linked':
                    area_key = f'tmt_ms1_area::{plex}'
                    if not cf.metaValueExists(area_key):
                        raise ValueError(
                            f'Linked-layout consensus feature {cf.getUniqueId()} '
                            f'lacks {area_key}; use scan layout for reporter intensities'
                        )
                    area = cf.getMetaValue(area_key)
                    if not np.isfinite(area) or area < 0:
                        raise ValueError(
                            f'Linked-layout consensus feature {cf.getUniqueId()} '
                            f'has invalid {area_key}: {area}'
                        )
                ref_int = self.get_reference_intensity(plex, handles)
                if not self.virtual_reference:
                    ref_index = self.column_2_ref[self.plex_2_columns[plex][0]]
                    ref_handle = next((h for h in handles if h.getMapIndex() == ref_index), None)
                    if ref_handle is None:
                        state = 'Missing reference handles'
                    else:
                        raw_ref = ref_handle.getIntensity()
                        if not np.isfinite(raw_ref):
                            state = 'Non-finite reference intensities'
                        elif raw_ref == 0:
                            state = 'Zero reference intensities'
                        elif raw_ref <= self.min_int:
                            state = 'Other references at or below min_int'
                        else:
                            state = 'Valid physical references'
                    self.summary[state] += 1
                    reference_validity.setdefault(plex, []).append(
                        bool(np.isfinite(ref_int) and ref_int > self.min_int)
                    )
                intensities = [h.getIntensity() if np.isfinite(h.getIntensity())
                               and h.getIntensity() > self.min_int else 0 for h in handles]
                if self.layout == 'linked':
                    ref_area = ref_int
                else:
                    total = sum(intensities)
                    precursor = cf.getMetaValue('precursor_intensity')
                    ref_area = 0 if total == 0 else precursor * ref_int / total
                ref_intensities[(feature_index, plex)] = ref_area
                cf_ratios = [np.nan] * len(self.plex_2_columns[plex])
                for handle, intensity in zip(handles, intensities):
                    if intensity > 0 and np.isfinite(ref_int) and ref_int > self.min_int:
                        ratio = np.log2(intensity / ref_int)
                        cf_ratios[column_positions[plex][handle.getMapIndex()]] = ratio
                ratios.setdefault(plex, []).append(cf_ratios)
                index_map[(feature_index, plex)] = len(ratios[plex]) - 1
        ratios_np = {k:np.array(v, dtype=np.float64) for k,v in ratios.items()}
        return ratios_np, ref_intensities, index_map, reference_validity

    @staticmethod
    def normalize_mad(ratios, reference_columns:dict[str, int]=None,
            reference_validity:Dict[str,List[bool]]=None, fallback_channels=None):
        """ Perform MAD normalization, excluding physical reference channels. """
        reference_columns = reference_columns or {}
        reference_validity = reference_validity or {}
        masked_ratios = {k:np.ma.masked_invalid(v.copy()) for k,v in ratios.items()}
        for plex, column in reference_columns.items():
            masked_ratios[plex][:, column] = np.ma.masked

        m_i = {k:np.ma.median(v, axis=0) for k,v in masked_ratios.items()}
        m_0 = np.ma.median(np.ma.masked_invalid(np.ma.concatenate(list(m_i.values()))))
        m_0 = 0.0 if np.ma.is_masked(m_0) else float(m_0)
        centered = {k:v-m_i[k] for k,v in masked_ratios.items()}
        mad_i = {k:np.ma.median(np.ma.masked_invalid(np.absolute(v)), axis = 0)
                 for k,v in centered.items()}
        mad_0 = np.ma.median(np.ma.masked_invalid(np.ma.concatenate(list(mad_i.values()))))
        mad_0 = 0.0 if np.ma.is_masked(mad_0) else float(mad_0)
        normalized = {}
        for plex, values in centered.items():
            output = np.full(values.shape, np.nan)
            for column in range(values.shape[1]):
                if column in (reference_columns.get(plex),):
                    continue
                mad = mad_i[plex][column]
                if np.ma.is_masked(mad) or not np.isfinite(mad) or mad <= 0 or mad_0 <= 0:
                    result = values[:, column] + m_0
                    if fallback_channels is not None:
                        fallback_channels.append((plex, column))
                else:
                    result = values[:, column] / mad * mad_0 + m_0
                output[:, column] = result.filled(np.nan)
            normalized[plex] = output
        for plex, column in reference_columns.items():
            normalized[plex][:, column] = np.where(reference_validity[plex], 0.0, np.nan)
        return normalized

    def calculate_bridge_factors(self, ref_intensities):
        """ Align each plex's median of all valid reference areas. """
        plexes = list(self.plex_2_columns)
        medians = {}
        counts = {}
        for plex in plexes:
            values = [area for (_, area_plex), area in ref_intensities.items()
                      if area_plex == plex and np.isfinite(area) and area > 0]
            counts[plex] = len(values)
            if self.bridge_scale and not values:
                raise ValueError(f'Cannot bridge-scale: no valid reference area for plex {plex}')
            medians[plex] = float(np.median(values)) if values else np.nan
        valid_medians = [value for value in medians.values() if np.isfinite(value)]
        target = float(np.median(valid_medians)) if valid_medians else np.nan
        factors = {plex: target / medians[plex] if self.bridge_scale else 1.0
                   for plex in plexes}
        self.scaling_report = [
            {'plex_id': plex, 'reference_feature_count': counts[plex],
             'median_reference_area': medians[plex], 'target_reference_area': target,
             'plex_scale_factor': factors[plex]}
            for plex in plexes
        ]
        if self.scaling_report_path:
            with open(self.scaling_report_path, 'w', newline='', encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(self.scaling_report[0]),
                                        delimiter='\t')
                writer.writeheader()
                writer.writerows(self.scaling_report)
        return factors

    def create_consensus_map(self, ratios, ref_intensities, index_map, factors):
        """ Copy the map while changing only handle intensities and output metadata. """
        cmap = ConsensusMap(self.cmap)
        cmap.clear(False)
        column_positions = {
            plex: {column: pos for pos, column in enumerate(columns)}
            for plex, columns in self.plex_2_columns.items()
        }
        for i, cf in enumerate(self.cmap):
            new_cf = ConsensusFeature(cf)
            new_cf.clear()
            new_cf.setUniqueId(cf.getUniqueId())
            new_cf.setRT(cf.getRT())
            new_cf.setMZ(cf.getMZ())
            new_cf.setIntensity(cf.getIntensity())
            new_cf.setCharge(cf.getCharge())
            new_cf.setWidth(cf.getWidth())
            new_cf.setQuality(cf.getQuality())
            new_cf.setRatios(cf.getRatios())
            new_cf.setPeptideIdentifications(cf.getPeptideIdentifications())
            keys = []
            cf.getKeys(keys)
            for key in keys:
                new_cf.setMetaValue(key, cf.getMetaValue(key))
            has_finite_sample = False
            for feature in cf.getFeatureList():
                new_feature = BaseFeature()
                new_feature.setUniqueId(feature.getUniqueId())
                new_feature.setWidth(feature.getWidth())
                new_feature.setCharge(feature.getCharge())
                new_feature.setMZ(feature.getMZ())
                new_feature.setRT(feature.getRT())

                map_index = feature.getMapIndex()
                plex = self.column_2_plex[map_index]
                ref_int = ref_intensities[i, plex]
                cf_ratios = ratios[plex][index_map[i, plex]]
                ratio = cf_ratios[column_positions[plex][map_index]]
                if self.quasi_intensity:
                    intensity = np.nan
                    if np.isfinite(ref_int) and ref_int > 0:
                        if np.isfinite(ratio):
                            with np.errstate(over='ignore', under='ignore', invalid='ignore'):
                                intensity = ref_int * np.exp2(ratio) * factors[plex]
                            if not np.isfinite(intensity):
                                intensity = np.nan
                        elif feature.getIntensity() == 0:
                            intensity = 0
                else:
                    intensity = ratio
                if np.isfinite(intensity) and (self.virtual_reference or
                        map_index != self.column_2_ref[map_index]):
                    has_finite_sample = True
                new_feature.setIntensity(intensity)
                new_cf.insert(map_index, new_feature)
            if has_finite_sample:
                self.summary['ConsensusFeatures normalized'] += 1
            groups = self.group_handles(cf)
            for plex in groups:
                new_cf.setMetaValue(f'reference_intensity::{plex}', ref_intensities[i, plex])
                new_cf.setMetaValue(f'normalized_reference_intensity::{plex}',
                                    ref_intensities[i, plex] * factors[plex])
                new_cf.setMetaValue(f'plex_scale_factor::{plex}', factors[plex])
            if len(groups) == 1:
                plex = next(iter(groups))
                new_cf.setMetaValue('reference_intensity', ref_intensities[i, plex])
            cmap.push_back(new_cf)
        cmap.ensureUniqueId()
        return cmap

    def normalize(self):
        """ Conduct MAD normalization. """
        if self.layout not in ('scan', 'linked'):
            raise ValueError(f'Unknown MADNormalizer layout: {self.layout}')
        if self.bridge_scale and (self.layout != 'linked' or not self.quasi_intensity):
            raise ValueError('Bridge scaling requires linked layout and quasi-intensity output')
        if self.bridge_scale and self.virtual_reference and not self.allow_virtual_bridge:
            raise ValueError('Virtual-reference bridge scaling requires explicit opt-in')
        self.set_channel_plex()
        self.summary = {
            'ConsensusMaps processed': 1,
            'ConsensusMaps normalized': 1,
            'ConsensusFeatures processed': 0,
            'ConsensusFeatures normalized': 0,
        }
        if not self.virtual_reference:
            self.summary.update({
                'Valid physical references': 0,
                'Missing reference handles': 0,
                'Zero reference intensities': 0,
                'Other references at or below min_int': 0,
                'Non-finite reference intensities': 0,
            })

        # Sep 1: Get log2 ratios for each channel of each feature.
        ratios, ref_intensities, index_map, reference_validity = self.collect_ratios()
        if not ratios:
            raise ValueError('No reporter handles found in the input ConsensusMap')

        # Step 2: MAD normalization
        reference_columns:dict[str, int] = {}
        if not self.virtual_reference:
            for plex, columns in self.plex_2_columns.items():
                if plex not in ratios:
                    continue
                reference_columns[plex] = columns.index(self.column_2_ref[columns[0]])
        fallback_channels = []
        ratios = self.normalize_mad(ratios, reference_columns, reference_validity,
                                    fallback_channels)
        self.fallback_channels = [
            (plex, self.plex_2_columns[plex][position])
            for plex, position in fallback_channels
        ]
        self.summary['Channels with zero or undefined MAD'] = len(fallback_channels)

        factors = self.calculate_bridge_factors(ref_intensities)

        # Step 3: Put together a normalized ConsensusMap
        cmap = self.create_consensus_map(ratios, ref_intensities, index_map, factors)

        return cmap
