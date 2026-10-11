""" MzML filter """
import bisect
from typing import Dict, List, Set
from pyopenms import PeakMap, ConsensusMap, Peak2D, PeptideIdentification, MSSpectrum

from msptools import logger


# pylint: disable=R0912
class MzMLFilterTally():
    """ A tally of spectrum processed by MzMLFilter. """
    # pylint: disable=R0903
    def __init__(self):
        """ Constructor """
        self.n_total = 0
        self.n_total_ms1 = 0
        self.n_total_ms2 = 0
        self.n_total_ms3 = 0
        self.n_invalid_spectra = 0
        self.n_match_id_ms2 = 0
        self.n_match_id_ms3 = 0
        self.n_match_rt_and_mz_ms2 = 0
        self.n_match_rt_and_mz_ms3 = 0
        self.n_psm = 0
        self.n_stored = 0

    def log(self):
        """ log """
        logger(f"Total spectra: {self.n_total}")
        logger(f"Total MS1 spectra: {self.n_total_ms1}")
        logger(f"Total MS2 spectra: {self.n_total_ms2}")
        logger(f"Total MS3 spectra: {self.n_total_ms3}")
        logger(f"Spectra with invalid MS level: {self.n_invalid_spectra}")
        logger(f"PSMs loaded from denylist: {self.n_psm}")
        logger(f"MS2 spectra filtered by ID: {self.n_match_id_ms2}")
        logger(f"MS2 spectra filtered by RT & MZ: {self.n_match_rt_and_mz_ms2}")
        logger(f"MS3 spectra filtered by ID: {self.n_match_id_ms3}")
        logger(f"MS3 spectra filtered by RT & MZ: {self.n_match_rt_and_mz_ms3}")
        logger(f"Spectra stored: {self.n_stored}")

class MzMLFilter():
    """ MzMLFilter
    Filter mzML file with a given consensusXML. MS2 spectra that the precursor
    ions match with the consensusXML file are removed.
    """
    def __init__(self):
        """ Constructor """
        self.input_path:str = None
        self.output_path:str = None
        self.consensus_file:str = None
        self.id_file:str = None
        self.rt_tolerance:float = None
        self.mz_tolerance:float = None
        self.mz_tolerance_unit:bool = None
        self.consensus_check_id:bool = None
        self.use_rt_mz:bool = None
        self.exp = PeakMap()
        self.feature_pos:List[Peak2D] = []
        self.feature_ids:Set[str] = set()
        self.spectrum_map:Dict[str,MSSpectrum] = {}
        self.tally = MzMLFilterTally()

    def count_total_spectra(self):
        """ Tally total spectra. """
        for spectrum in self.exp:
            self.tally.n_total += 1
            if spectrum.getMSLevel() == 1:
                self.tally.n_total_ms1 += 1
            elif spectrum.getMSLevel() == 2:
                self.tally.n_total_ms2 += 1
            elif spectrum.getMSLevel() == 3:
                self.tally.n_total_ms3 += 1
            else:
                self.tally.n_invalid_spectra += 1

    def add_feature_id(self, val:str):
        """ Add feature id """
        before_size = len(self.feature_ids)
        self.feature_ids.add(val)
        after_size = len(self.feature_ids)
        if after_size == before_size:
            raise ValueError(f"Duplicated spectrum ID: {val}")

    def make_spectrum_dict(self) -> None:
        """ Make spectrum dict, a mapping from spectrum ID to the actual
        MSSpectrum object. This is to quickly retrieve spectrum by ID. """
        for spectrum in self.exp:
            native_id = spectrum.getNativeID()
            len_before = len(self.spectrum_map)
            self.spectrum_map[native_id] = spectrum
            len_after = len(self.spectrum_map)
            if len_before == len_after:
                raise ValueError(f"Duplicated spectrum ID: {native_id}")

    def get_feature_pos_id(self, pep_ids:List[PeptideIdentification]):
        """ Compute the feature list to be excluded from the mzML file from
        idXML. """
        self.feature_pos = []
        for pep_id in pep_ids:
            if pep_id.metaValueExists('spectrum_reference'):
                spectrum_reference = pep_id.getMetaValue('spectrum_reference')
                self.add_feature_id(spectrum_reference)
            p = Peak2D()
            p.setRT(pep_id.getRT())
            p.setMZ(pep_id.getMZ())
            self.feature_pos.append(p)
        self.feature_pos.sort(key=lambda x:x.getRT())

    def get_feature_pos_consensus(self, cmap:ConsensusMap):
        """ Compute the feature list to be excluded from the mzML file from
        consensusXML. """
        self.feature_pos = []
        for cf in cmap:
            if self.consensus_check_id:
                if not cf.getPeptideIdentifications():
                    continue
                pep_id = cf.getPeptideIdentifications()[0]
                if pep_id.metaValueExists('spectrum_reference'):
                    spectrum_reference = pep_id.getMetaValue('spectrum_reference')
                    self.add_feature_id(spectrum_reference)
                p = Peak2D()
                p.setRT(pep_id.getRT())
                p.setMZ(pep_id.getMZ())
            else:
                if cf.metaValueExists('spectrum_reference'):
                    spectrum_reference = cf.getMetaValue('spectrum_reference')
                    self.add_feature_id(spectrum_reference)
                p = Peak2D()
                p.setRT(cf.getRT())
                p.setMZ(cf.getMZ())

            self.feature_pos.append(p)
        self.feature_pos.sort(key=lambda x:x.getRT())

    def is_matched_by_native_id(self, spectrum:MSSpectrum) -> bool:
        """ Whether the given spectrum can be matched by spectrum reference ID
        (native ID). """
        if spectrum.getMSLevel() == 1:
            raise ValueError('MS1 should not be called.')

        if spectrum.getMSLevel() == 2:
            native_id = spectrum.getNativeID()
            if native_id in self.feature_ids:
                return True

        else: # MS3
            if self.consensus_file is not None and not self.consensus_check_id:
                # In this case self.feature_ids should contain MS3 spectrum IDs.
                return spectrum.getNativeID() in self.feature_ids

            # Peptide matching is done at MS2 level.
            precursor = None
            for precursor in spectrum.getPrecursors():
                if precursor.metaValueExists('ms level') and \
                        int(precursor.getMetaValue('ms level')) == 2:
                    break
            else:
                native_id = spectrum.getNativeID()
                raise ValueError(f'Can not find MS2 precursor for spectrum {native_id}')

            if precursor is None:
                native_id = spectrum.getNativeID()
                raise ValueError(f'Can not find MS2 precursor for spectrum {native_id}')

            if precursor.metaValueExists('spectrum_ref'):
                native_id = precursor.getMetaValue('spectrum_ref')
                if native_id in self.feature_ids:
                    return True
        return False

    def get_precursor_spectrum(self, spectrum:MSSpectrum,
            precursor_ms_level:int=1) -> MSSpectrum:
        """ Get precursor spectrum """
        native_id = spectrum.getNativeID()
        if not spectrum.getPrecursors()[0].metaValueExists('spectrum_ref'):
            raise ValueError(f"Precursor pectrum reference not found for specturm {native_id}")
        precursor_ref = spectrum.getPrecursors()[0].getMetaValue('spectrum_ref')
        precursor = self.spectrum_map.get(precursor_ref)
        if not precursor:
            raise ValueError(f"Failed to find precursor spectrum for {native_id}")
        if precursor.getMSLevel() > precursor_ms_level:
            return self.get_precursor_spectrum(precursor)
        return precursor

    def increment_tally(self, spectrum:MSSpectrum, id_matched:bool):
        """ Increment tally """
        if spectrum.getMSLevel() == 2:
            if id_matched:
                self.tally.n_match_id_ms2 += 1
            else:
                self.tally.n_match_rt_and_mz_ms2 += 1
        elif spectrum.getMSLevel() == 3:
            if id_matched:
                self.tally.n_match_id_ms3 += 1
            else:
                self.tally.n_match_rt_and_mz_ms3 += 1
        else:
            native_id = spectrum.getNativeID()
            raise ValueError(f"MS1 spectrum should not be filtered: {native_id}")


    def filter_by_denylist(self) -> PeakMap:
        """ Filter MS2 based on the denylist. """
        feature_pos_rt = [x.getRT() for x in self.feature_pos]
        list_idx:Set[int] = set()
        for i,spectrum in enumerate(self.exp):
            if spectrum.getMSLevel() == 1:
                continue

            if not spectrum.getPrecursors():
                continue

            if self.is_matched_by_native_id(spectrum):
                self.increment_tally(spectrum, True)
                list_idx.add(i)

            if not self.use_rt_mz:
                continue

            # Can't get MZ from MSSpectrum. But Precursor does have the .getMZ method..
            if spectrum.getMSLevel() == 2:
                ms2_precursor = spectrum
            else:
                ms2_precursor = self.get_precursor_spectrum(spectrum, 2)
            ms1_precursor_ion = ms2_precursor.getPrecursors()[0]
            ms1_precursor_spectrum = self.get_precursor_spectrum(ms2_precursor, 1)
            pc_mz = ms1_precursor_ion.getMZ()
            pc_rt = ms1_precursor_spectrum.getRT()
            if self.mz_tolerance_unit:
                mz_tol_da = pc_mz * self.mz_tolerance * 10 ** (-6)
            else:
                mz_tol_da = self.mz_tolerance
            p_low = bisect.bisect_left(feature_pos_rt, pc_rt - self.rt_tolerance)
            p_high = bisect.bisect_right(feature_pos_rt, pc_rt + self.rt_tolerance)
            for f_it in range(p_low, p_high):
                if pc_mz - mz_tol_da < self.feature_pos[f_it].getMZ() < pc_mz + mz_tol_da:
                    self.increment_tally(spectrum, False)
                    list_idx.add(i)

        exp2 = PeakMap(self.exp)
        exp2.clear(False)
        for i,spectrum in enumerate(self.exp):
            if i not in list_idx:
                exp2.addSpectrum(spectrum)

        return exp2
