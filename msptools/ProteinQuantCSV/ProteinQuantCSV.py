""" IO """
from __future__ import annotations
from typing import TYPE_CHECKING
from abc import ABC, abstractmethod
import csv
import numpy as np
from pyopenms import ConsensusMap, PeptideProteinResolution
from .ProteinQuantCSVMetadata import ProteinQuantCSVMetadata


if TYPE_CHECKING:
    from typing import IO, Iterable, Dict, List, Set, Tuple
    from pyopenms import (
        ConsensusFeature, PeptideIdentification, PeptideHit, AASequence,
        ProteinIdentification
    )

class QuantRecord(ABC):
    """ Record for ProteinQuantCSV """
    def __init__(self, abundance:List[float], header:List[str]):
        """ Constructor """
        self.abundance = abundance
        self.header = header

    @abstractmethod
    def left_attributes(self) -> List[str]:
        """ Fields at the left of the abundance matrix. """

    @abstractmethod
    def right_attributes(self) -> List[str]:
        """ Fields at the right of the abundance matrix. """

    def to_list(self):
        """ to fields """
        return [
            *self.left_attributes(),
            *self.abundance,
            *self.right_attributes()
        ]

class ProteinQuantRecord(QuantRecord):
    """ ProteinQuantRecord """
    def __init__(self, abundance:List[float], header:List[int], protein:str,
            n_proteins:int, protein_score:float, n_peptides:int):
        """ Constructor """
        super().__init__(abundance=abundance, header=header)
        self.protein = protein
        self.n_proteins = n_proteins
        self.protein_score = protein_score
        self.n_peptides = n_peptides

    def left_attributes(self):
        """ Fields at the left of the abundance matrix. """
        return [
            self.protein,
            self.n_proteins,
            self.protein_score,
            self.n_peptides
        ]

    def right_attributes(self):
        """ Fields at the right of the abundance matrix. """
        return []

class PeptideQuantRecord(QuantRecord):
    """ PeptideQuantRecord """
    def __init__(self, abundance:List[float], header:List[int], peptide:str,
            protein:str, n_proteins:int, charge:int, fraction:str):
        """ Constructor """
        super().__init__(abundance=abundance, header=header)
        self.peptide = peptide
        self.protein = protein
        self.n_proteins = n_proteins
        self.charge = charge
        self.fraction=fraction

    def left_attributes(self):
        """ Fields at the left of the abundance matrix. """
        return [
            self.peptide,
            self.protein,
            self.n_proteins,
            self.charge
        ]

    def right_attributes(self):
        """ Fields at the right of the abundance matrix. """
        return [
            self.fraction
        ]

class ProteinQuantCSV:
    """ ProteinQuantCSV """
    def __init__(self, records:List[QuantRecord]=None,
            metadata:ProteinQuantCSVMetadata=None):
        """ constructor """
        self.records = records or []
        self.metadata = metadata

    def remove_zero_plex(self, summary=None):
        """ For features that all channels of a plex are zero, set them to NA """
        plex2column = self.metadata.get_plex_2_column(-1)
        summarization = bool(summary == {})
        if summarization:
            for k in plex2column:
                summary[k] = 0
        for record in self.records:
            for plex, columns in plex2column.items():
                if all(record.abundance[i] == 0 for i in columns):
                    for i in columns:
                        record.abundance[i] = np.nan
                    if summarization:
                        summary[plex] += 1

    def remove_unquantified_peptides(self, cmap:ConsensusMap,
            summary:List[Tuple[str, int]]=None, check_intensity:bool=True):
        """ Set the intensity of unquantified peptides to nan """
        if summary == []:
            summarization = True
            _summary = [0 for _ in self.metadata.header_map]
        else:
            summarization = False

        # get channel to plex map and plex to channel map
        channel2plex = self.metadata.get_column_2_plex(-1)
        plex2channel = self.metadata.get_plex_2_column(-1)

        # create map from peptide to quantified channel
        # map: peptide -> charge -> quantified channels
        peptide_map:Dict[str, Dict[int, Set[int]]] = {}
        cf:ConsensusFeature
        for cf in cmap:
            quantified_channels = [f.getMapIndex() for f in cf.getFeatureList()]

            pept_ids:List[PeptideIdentification] = cf.getPeptideIdentifications()
            if not pept_ids:
                continue

            hits:List[PeptideHit] = pept_ids[0].getHits()
            if not hits:
                continue

            peptide:AASequence = hits[0].getSequence()
            seq = peptide.toString()
            charge = hits[0].getCharge()

            cur = peptide_map.setdefault(seq, {})
            cur = cur.setdefault(charge, set())
            cur.update(quantified_channels)

        # Now fix the intensities
        record:PeptideQuantRecord
        for record in self.records:
            quantified_channels = peptide_map.get(record.peptide, {}).get(record.charge)
            for i, intensity in enumerate(record.abundance):
                plex_channels = plex2channel[channel2plex[i]]
                if all(np.isnan(record.abundance[x]) for x in plex_channels):
                    continue
                if quantified_channels is None:
                    raise ValueError(
                        'No peptide/charge mapping found in consensusXML file:'
                        f' sequence: {record.peptide}, charge: {record.charge}'
                    )
                if i not in quantified_channels:
                    if check_intensity and intensity != 0:
                        raise ValueError(
                            "Peptide has non-zero intensity value in ProteinQuantCSV"
                            " but marked as unquantified in consensusXML file."
                            f" sequence: {record.peptide}, charge: {record.charge}"
                        )
                    record.abundance[i] = np.nan
                    if summarization:
                        _summary[i] += 1

        if summarization and self.records:
            for i, sample in enumerate(self.records[0].header):
                summary.append((sample, _summary[i]))

    def removed_unquantified_proteins(self, cmap:ConsensusMap,
            protein_groups:ProteinIdentification=None, greedy_group_resolution:bool=False,
            summary:List[Tuple[str, int]]=None, check_intensity:bool=True):
        """ Set the intensity of unquantified proteins to nan """
        if summary == []:
            summarization = True
            _summary = [0 for _ in self.metadata.header_map]
        else:
            summarization = False
        # get accession to leader
        accession_to_leader:Dict[str,str] = {}
        if not protein_groups:
            protein_groups = cmap.getProteinIdentifications()[0]

        if greedy_group_resolution:
            ppr = PeptideProteinResolution(False)
            ppr.buildGraph(protein_groups, [])
            ppr.resolveGraph(protein_groups, [])

        for pg in protein_groups.getIndistinguishableProteins():
            accessions = [x.decode('utf-8') for x in pg.accessions]
            for acc in accessions:
                accession_to_leader[acc] = accessions[0]

        # get channel to plex map and plex to channel map
        channel2plex = self.metadata.get_column_2_plex(-1)
        plex2channel = self.metadata.get_plex_2_column(-1)

        # create map from protein to quantified channel
        # map: protein  -> quantified channels
        protein_map:Dict[str, Set[int]] = {}
        cf:ConsensusFeature
        for cf in cmap:
            quantified_channels = [f.getMapIndex() for f in cf.getFeatureList()]

            pept_ids:List[PeptideIdentification] = cf.getPeptideIdentifications()
            if not pept_ids:
                continue

            hits:List[PeptideHit] = pept_ids[0].getHits()
            if not hits:
                continue

            accessions = [x.decode('utf-8') for x in hits[0].extractProteinAccessionsSet()]
            accession = self.get_accession(accessions, accession_to_leader)

            if accession:
                cur = protein_map.setdefault(accession, set())
                cur.update(quantified_channels)

        record:ProteinQuantRecord
        for record in self.records:
            proteins = record.protein.split('/')

            quantified_channels = set()
            for protein in proteins:
                channels = protein_map.get(protein)
                if channels:
                    quantified_channels.update(channels)
            for i, intensity in enumerate(record.abundance):
                plex_channels = plex2channel[channel2plex[i]]
                if all(np.isnan(record.abundance[x]) for x in plex_channels):
                    continue
                if i not in quantified_channels:
                    if check_intensity and intensity != 0 and not np.isnan(intensity):
                        raise ValueError(
                            "Protein has non-zero intensity value in ProteinQuantCSV"
                            " but marked as unquantified in consensusXML file."
                            f" protein: {record.protein}, column: {i}."
                        )
                    record.abundance[i] = np.nan
                    if summarization:
                        _summary[i] += 1

        if summarization and self.records:
            for i, sample in enumerate(self.records[0].header):
                summary.append((sample, _summary[i]))

    @staticmethod
    def get_accession(pep_accessions:Set[str], accession_to_leader:Dict[str,str]):
        """ Get accession """
        if accession_to_leader:
            leaders:Set[str] = set()
            for acc in pep_accessions:
                leader = accession_to_leader.get(acc)
                if leader:
                    leaders.add(leader)
            if not leaders:
                return None
            if len(leaders) == 1:
                return list(leaders)[0]

        elif len(pep_accessions) == 1:
            return list(pep_accessions)[0]

        return None

    def write(self, handle:IO, fieldnames:List[str], delimiter:str, quoting:str):
        """ Write ProteinQuantCSV """
        for line in self.metadata.to_lines():
            handle.write(line + '\n')

        writer = csv.writer(
            handle,
            delimiter=delimiter,
            quoting=quoting
        )
        writer.writerow(fieldnames)

        for record in self.records:
            writer.writerow(record.to_list())

    @staticmethod
    def get_fieldnames(handle:IO, delim:str):
        """ get fieldnames """
        handle.seek(0)
        for line in handle:
            if line.startswith('#'):
                continue
            fieldnames = [x.strip('"').strip("'") for x in line.rstrip().split(delim)]
            break
        handle.seek(0)
        return fieldnames

    @staticmethod
    def get_quoting(quoting:str) -> str:
        """ get quoting from OpenMS stype to python's CSV style """
        if quoting == 'none':
            return csv.QUOTE_NONE
        if quoting == 'double':
            return csv.QUOTE_NONNUMERIC
        if quoting == 'escape':
            return csv.QUOTE_MINIMAL
        raise ValueError(f"Unrecognized value for `quoting`: '{quoting}'")

    @staticmethod
    def parse(handle:IO, format:str) -> Iterable[ProteinQuantRecord]:   # pylint: disable=redefined-builtin
        """ Parse ProteinQuantCSV file """
        format = format.lower()
        supported_formats = ['protein', 'peptide']
        if format not in supported_formats:
            raise ValueError(
                f"Unsupported format: '{format}'. Supported formats are: {supported_formats}"
            )
        header = None
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            fields = [f.strip('"').strip("'") for f in fields]

            start_i = 4
            if format == 'peptide':
                end_i = len(fields) - 1
            else:
                end_i = len(fields)

            if not header:
                header = fields[start_i:end_i]
                continue

            abundance = [float(x) for x in fields[start_i:end_i]]

            if format == 'peptide':
                fraction = fields[-1]
                yield PeptideQuantRecord(
                    abundance=abundance,
                    header=header,
                    peptide=fields[0],
                    protein=fields[1],
                    n_proteins=int(fields[2]),
                    charge=int(fields[3]),
                    fraction=fraction
                )
            else:
                yield ProteinQuantRecord(
                    abundance=abundance,
                    header=header,
                    protein=fields[0],
                    n_proteins=int(fields[1]),
                    protein_score=float(fields[2]),
                    n_peptides=int(fields[3])
                )

    @classmethod
    def read(cls, handle:IO, format:str):   # pylint: disable=redefined-builtin
        """ Read ProteinQuantCSV """
        metadata = ProteinQuantCSVMetadata.parse(handle)
        if not metadata.header_map:
            fieldnames = cls.get_fieldnames(handle, '\t')
            abundance_headers = fieldnames[4:-1] if format.lower() == 'peptide' else fieldnames[4:]
            if abundance_headers and all(name.startswith('abundance|') and
                    '|' in name[len('abundance|'):] for name in abundance_headers):
                metadata.header_map = {
                    i: tuple(name[len('abundance|'):].rsplit('|', 1))
                    for i, name in enumerate(abundance_headers, start=1)
                }
        records = list(cls.parse(handle, format=format))
        return ProteinQuantCSV(records=records, metadata=metadata)
