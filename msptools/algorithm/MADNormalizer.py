""" MAD Normalizer """
from __future__ import annotations
from typing import List, Dict, Set, Tuple
import statistics
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

    The normalized ratios (RN_ij) are the final normalization result. The
    reference intensities (REF_kj) are saved to the output consensusXML file
    as a meta value of each consensus feature, and can be utilized in protein
    or peptide quantitation to convert ratio values to the absolute abundance scale.
    """
    def __init__(self, cmap:ConsensusMap, virtual_reference:bool=False,
            virtual_reference_method:str='median', min_int:float=1e-08,
            column_2_ref:Dict[int,int]=None, column_2_plex:Dict[int, str]=None,
            plex_2_column:Dict[str,List[int]]=None):
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
        # plex -> reference channel
        plex_2_ref:Dict[str, int] = {}
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
                plex_2_ref[plex] = i

        column_2_ref = {}
        if not self.virtual_reference:
            for i,plex in column_2_plex.items():
                try:
                    ref = plex_2_ref[plex]
                except KeyError as e:
                    raise ValueError(f"Cannot find reference channel for plex {plex}") from e
                column_2_ref[i] = ref

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

    def get_plex_and_ref_intensity(self, cf:ConsensusFeature) -> Tuple[int, float]:
        """ Get the plex ID and reference intensity of a consensus feature """
        virtual_fun = VIRTUAL_REFERENCE_AVEAGE_FUNS[self.virtual_reference_method]
        ref_indices = set()
        plexes:Set[str] = set()
        if self.virtual_reference:
            for feature in cf.getFeatureList():
                index = feature.getMapIndex()
                plexes.add(self.column_2_plex[index])
            if len(plexes) > 1:
                raise ValueError(f"Non-unique plex found: {plexes}")
            plex = list(plexes)[0]
            ints = [f.getIntensity() for f in cf.getFeatureList()]
            ints = [x for x in ints if x > self.min_int]
            if len(ints) <= 1:
                ref_int = 0
            else:
                ref_int = virtual_fun(ints)
            return plex, ref_int

        features:Dict[int, FeatureHandle] = {}
        for feature in cf.getFeatureList():
            index = feature.getMapIndex()
            features[index] = feature
            ref_indices.add(self.column_2_ref[index])
            plexes.add(self.column_2_plex[index])
        if len(ref_indices) > 1:
            raise ValueError(f"Non-unique reference found: {ref_indices}")
        if len(plexes) > 1:
            raise ValueError(f"Non-unique plex found: {plexes}")
        plex = list(plexes)[0]
        ref_index = list(ref_indices)[0]
        ref_feature = features.get(ref_index)
        if ref_feature is None:
            return plex, 0
        ref_int = ref_feature.getIntensity()
        if ref_int < self.min_int:
            ref_int = 0
        return plex, ref_int

    def collect_ratios(self):
        """ Collect log2 ratios for each channel in each feature. """
        ratios:Dict[List[List[float]]] = {}
        ref_intensities:List[float] = []
        index_map:List[Tuple[str,int]] = []
        reference_validity:Dict[str,List[bool]] = {}
        column_positions = {
            plex: {column: pos for pos, column in enumerate(columns)}
            for plex, columns in self.plex_2_columns.items()
        }
        cf: ConsensusFeature
        for cf in self.cmap:
            plex, ms2_intensity_ref = self.get_plex_and_ref_intensity(cf)
            if not self.virtual_reference:
                reference_validity.setdefault(plex, []).append(
                    bool(np.isfinite(ms2_intensity_ref) and ms2_intensity_ref > self.min_int)
                )
            features = cf.getFeatureList()
            intensities = []
            for f in features:
                intensity = f.getIntensity()
                if intensity > self.min_int:
                    intensities.append(intensity)
                else:
                    intensities.append(0)
            ms2_intensity_sum = sum(intensities)
            ms1_intensity = cf.getMetaValue('precursor_intensity')
            if ms2_intensity_sum == 0:
                ref_intensities.append(0)
            else:
                ref_intensities.append(ms1_intensity * ms2_intensity_ref / ms2_intensity_sum)
            if ms2_intensity_ref == 0:
                feature_ratios = [0 for _ in features]
            else:
                feature_ratios = [i / ms2_intensity_ref for i in intensities]
            cf_ratios = [np.nan] * len(self.plex_2_columns[plex])
            for feature, ratio in zip(features, feature_ratios):
                position = column_positions[plex][feature.getMapIndex()]
                cf_ratios[position] = np.log2(ratio) if ratio != 0 else np.nan

            if plex not in ratios:
                ratios[plex] = []
            ratios[plex].append(cf_ratios)
            index_map.append((plex, len(ratios[plex]) - 1))
        ratios_np = {k:np.array(v, dtype=np.float64) for k,v in ratios.items()}
        return ratios_np, ref_intensities, index_map, reference_validity

    @staticmethod
    def normalize_mad(ratios, reference_columns:dict[str, int]=None,
            reference_validity:Dict[str,List[bool]]=None):
        """ Perform MAD normalization, excluding physical reference channels. """
        reference_columns = reference_columns or {}
        reference_validity = reference_validity or {}
        masked_ratios = {k:np.ma.masked_invalid(v.copy()) for k,v in ratios.items()}
        for plex, column in reference_columns.items():
            masked_ratios[plex][:, column] = np.ma.masked

        m_i = {k:np.ma.median(v, axis=0) for k,v in masked_ratios.items()}
        m_0 = np.ma.median(np.ma.masked_invalid(np.ma.concatenate(m_i.values())))
        centered = {k:v-m_i[k] for k,v in masked_ratios.items()}
        mad_i = {k:np.ma.median(np.ma.masked_invalid(np.absolute(v)), axis = 0)
                 for k,v in centered.items()}
        mad_0 = np.ma.median(np.ma.masked_invalid(np.ma.concatenate(mad_i.values())))
        normalized = {k:((v / mad_i[k]) * mad_0 + m_0).filled(np.nan)
                      for k,v in centered.items()}
        for plex, column in reference_columns.items():
            normalized[plex][:, column] = np.where(reference_validity[plex], 0.0, np.nan)
        return normalized

    def create_consensus_map(self, ratios, ref_intensities, index_map):
        """ Create consensus map with the normalized ratios. """
        cmap = ConsensusMap(self.cmap)
        cmap.clear(False)
        column_positions = {
            plex: {column: pos for pos, column in enumerate(columns)}
            for plex, columns in self.plex_2_columns.items()
        }
        for i, cf in enumerate(self.cmap):
            plex, j = index_map[i]
            cf_ratios = ratios[plex][j]
            ref_int = ref_intensities[i]
            new_cf = ConsensusFeature(cf)
            new_cf.clear()
            for feature in cf.getFeatureList():
                new_feature = BaseFeature()
                new_feature.setWidth(feature.getWidth())
                new_feature.setCharge(feature.getCharge())
                new_feature.setMZ(feature.getMZ())
                new_feature.setRT(feature.getRT())

                map_index = feature.getMapIndex()
                ratio = cf_ratios[column_positions[plex][map_index]]
                new_feature.setIntensity(ratio)
                new_cf.insert(map_index, new_feature)
            new_cf.setMetaValue('reference_intensity', ref_int)
            cmap.push_back(new_cf)
        cmap.ensureUniqueId()
        return cmap

    def normalize(self):
        """ Conduct MAD normalization. """
        self.set_channel_plex()

        # Sep 1: Get log2 ratios for each channel of each feature.
        ratios, ref_intensities, index_map, reference_validity = self.collect_ratios()

        # Step 2: MAD normalization
        reference_columns:dict[str, int] = {}
        if not self.virtual_reference:
            for plex, columns in self.plex_2_columns.items():
                reference_columns[plex] = columns.index(self.column_2_ref[columns[0]])
        ratios = self.normalize_mad(ratios, reference_columns, reference_validity)

        # Step 3: Put together a normalized ConsensusMap
        cmap = self.create_consensus_map(ratios, ref_intensities, index_map)

        return cmap
