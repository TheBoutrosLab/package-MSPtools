""" ProteinQuantCSV annotator """
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import re
import csv
from typing import Dict, List, Set, Tuple, IO
from Bio import SeqIO
from pyopenms import (
    IdXMLFile, PeptideIdentification, ProteinIdentification, PeptideHit
)
from moPepGen import constant
from moPepGen.index import IndexDir
from moPepGen.seqvar import GVFMetadata
from moPepGen.aa import VariantPeptideIdentifier as vpi
from moPepGen.aa.VariantPeptideLabel import LabelSourceMapping, VariantSourceSet


ALT_TRANS_TYPES = [constant.SEC_TERMINATION_TYPE, *constant.CODON_REASSIGNMENTS_TYPES]

class PeptideEntry():
    """ Peptide entry """
    def __init__(self, peptide:str, protein:str, unmod_seq:str):
        """ constructor """
        self.peptide = peptide
        self.protein = protein
        self.unmod_seq = unmod_seq

class QuantifiedPeptide(PeptideEntry):
    """ Quantified Peptide, representing a record from the ProteinQuantifer.csv
    table. """
    def __init__(self, peptide:str, protein:str, unmod_seq:str, n_proteins:int,
            charge:int, abundances:Dict[str,float], fraction:int):
        """ """
        super().__init__(
            peptide=peptide,
            protein=protein,
            unmod_seq=unmod_seq
        )
        self.n_proteins = n_proteins
        self.charge = charge
        self.abundances = abundances
        self.fraction = fraction

class IDPeptide(PeptideEntry):
    """ Identified peptide (from idXML) """
    def __init__(self, peptide:str, protein:str, unmod_seq:str, psm_id:str, charge:int, mz:float,
            score:float, score_type:str):
        """ constructor """
        super().__init__(
            peptide=peptide,
            protein=protein,
            unmod_seq=unmod_seq
        )
        self.psm_id = psm_id
        self.charge = charge
        self.mz = mz
        self.score = score
        self.score_type = score_type

@dataclass
class FastaHeaderAnnotation():
    """ Fasta header entry level annotation """
    annotation:list[str] = field(default_factory=list)
    raw_sources:VariantSourceSet = field(default_factory=VariantSourceSet)
    grouped_sources:VariantSourceSet = field(default_factory=VariantSourceSet)

class PeptideAnnotation():
    """ Peptide annotation information """
    def __init__(self, peptide:str, unmod_seq:str, sample_name:str,
            fasta_header_annotation:Dict[str, FastaHeaderAnnotation],
            raw_sources:Set[str]=None,
            grouped_sources:VariantSourceSet=None):
        """ Constructor """
        self.peptide = peptide
        self.unmod_seq = unmod_seq
        self.sample_name = sample_name
        self.fasta_header_annotation = fasta_header_annotation or {}
        self.raw_sources = raw_sources or VariantSourceSet()
        self.grouped_sources = grouped_sources or VariantSourceSet()

    @staticmethod
    def table_headers() -> str:
        """ Table headers """
        headers = [
            'peptide',
            'unmod_seq',
            'sample',
            'fasta_header_entry',
            'fasta_header_annotation',
            'raw_source_entry',
            'grouped_source_entry',
            'raw_source_peptide',
            'grouped_source_peptide'
        ]
        return '\t'.join(headers)

    def to_lines(self) -> List[str]:
        """ Convert to lines to be written to the output tsv file. """
        raw_sources = str(self.raw_sources)
        grouped_sources = str(self.grouped_sources)
        lines:List[str] = []
        for var_label, header_anno in self.fasta_header_annotation.items():
            fields = [
                self.peptide,
                self.unmod_seq,
                self.sample_name,
                var_label,
                '|'.join(header_anno.annotation),
                str(header_anno.raw_sources),
                str(header_anno.grouped_sources),
                raw_sources,
                grouped_sources
            ]
            lines.append('\t'.join(fields))
        return lines

    @staticmethod
    def parse(handle:IO[str]) -> Dict[str, List[PeptideAnnotation]]:
        """ parser """
        # mapping from peptide -> sample -> annotation
        records:dict[str, dict[str, PeptideAnnotation]] = {}
        reader = csv.DictReader(handle, delimiter='\t')
        for rec in reader:
            peptide = rec['peptide']
            sample = rec['sample']
            val = records.setdefault(peptide, {})
            raw_sources_peptide = rec['raw_source_peptide'].split('-')
            grouped_sources_peptide = rec['grouped_source_peptide'].split('-')
            raw_sources_entry = rec['raw_source_entry'].split('-')
            grouped_sources_entry = rec['grouped_source_entry'].split('-')

            val:PeptideAnnotation = val.setdefault(
                sample,
                PeptideAnnotation(
                    peptide = peptide,
                    unmod_seq = rec['unmod_seq'],
                    sample_name = sample,
                    raw_sources=set(raw_sources_peptide),
                    grouped_sources=set(grouped_sources_peptide),
                    fasta_header_annotation={}
                )
            )
            header = FastaHeaderAnnotation(
                annotation=rec['fasta_header_annotation'].split('|'),
                raw_sources=set(raw_sources_entry),
                grouped_sources=set(grouped_sources_entry)
            )
            val.fasta_header_annotation[rec['fasta_header_entry']] = header
        flatten = {}
        for k,v in records.items():
            flatten[k] = list(v.values())
        return flatten

@dataclass
class TranscriptAnnotation():
    """ Transcript annotation """
    transcript_id:str = None
    gene_id:str = None
    gene_symbol:str = None
    gene_biotype:str = None
    transcript_biotype:str = None
    is_protein_coding:bool = False

@dataclass
class ORFAnnoataion():
    """ ORF annotation """
    orf_id:str = None
    start:int = None
    end:int = None

class Tally():
    """ Annotation tally """
    def __init__(self):
        """ constructor """
        self.n_total_peptides = 0
        self.n_peptides_annotated = 0

    def log(self):
        """ Log tally result """
        print(
            f"  Peptides annotated in at least one sample: {self.n_peptides_annotated}",
            flush=True
        )
        print(
            f"  Peptides not annotated: {self.n_total_peptides - self.n_peptides_annotated}",
            flush=True
        )

    @staticmethod
    def warn_unannotated(peptide:str):
        """ Throw a warning for unannotated peptide. """
        print(f"Unannotated: {peptide}")

class VariantPeptideAnnotator():
    """ Quantified peptide annotator """
    def __init__(self, peptide_paths:Dict[str:Dict[str,Path]]=None,
            merged_dbs:Dict[str:Dict[str,Path]]=None,
            individual_dbs:Dict[str:Dict[str,Path]]=None,
            gvfs:Dict[str,List[Path]]=None,
            peptides:List[PeptideEntry]=None,
            encode_headers:Dict[str,List[str]]=None,
            gvf_info:Dict[str, Dict[str,Dict[str,Dict[str,str]]]]=None,
            variant_labels:Dict[str,Set[str]]=None,
            peptide_annotations:Dict[str,Dict[str,PeptideAnnotation]]=None,
            transcript_annotation:Dict[str,TranscriptAnnotation]=None,
            orf_annotation:Dict[str,Dict[str,ORFAnnoataion]]=None,
            coding_txs:Set[str]=None, label_free:bool=False, tally:Tally=None,
            group_map:Dict[str,str]=None, source_order:Dict[str,int]=None,
            label_map:Dict[str,LabelSourceMapping]=None,
            database_source_map:Dict[str,List[str]]=None,
            static_source_order:bool=False
            ):
        """ Constructor """
        # mapping: plex_id -> database_identifier -> quant.csv
        self.peptide_paths = peptide_paths or {}
        # mapping: plex_id -> database_identifier -> merged.fasta
        self.merged_dbs = merged_dbs or {}
        # mapping: plex_id -> sample_name -> individual.fasta
        self.individual_dbs = individual_dbs or {}
        # mapping: sample_name -> gvf files.
        self.gvfs = gvfs or {}

        self.peptides = peptides or []
        # mapping: encode ID -> variant labels
        self.encode_headers = encode_headers or {}
        # mapping: sample index -> transcript -> variant label -> attr key -> attr value
        self.gvf_info = gvf_info or {}
        # mapping: sample index -> variant labels
        self.variant_labels = variant_labels or {}
        # mapping: peptide seq -> sample name -> PeptideAnnotation
        self.peptide_annotations = peptide_annotations or {}

        self.transcript_annotation = transcript_annotation or {}
        self.coding_txs = coding_txs or set()
        # mapping: tx_id -> orf_id -> orf annotations
        self.orf_annotation = orf_annotation or {}

        self.label_free = label_free

        self.tally = tally

        self.group_map = group_map or {}
        self.source_order = source_order or {}
        self.label_map = label_map or {}
        self.database_source_map = database_source_map or {}
        self.static_source_order = static_source_order

    def clear_tally(self):
        """ Clear tally """
        self.tally = Tally()

    def append_order_internal_sources(self):
        """ Add internal sources that are not present in any GTFs, including
        noncoding, sec termination, and codon reassignment. """
        sources = [
            constant.SOURCE_NOVEL_ORF,
            constant.SOURCE_SEC_TERMINATION,
            constant.SOURCE_CODON_REASSIGNMENT
        ]
        for source in sources:
            if source in self.source_order:
                continue
            if source not in self.source_order:
                self.append_order(source)

    def append_order(self, source:str):
        """ Add a group to the end of the order. """
        if source in self.source_order:
            raise ValueError(f'The source {source} already has an order.')
        self.add_source_order(source)

    def set_source_set_level(self):
        """ Set source set level """
        VariantSourceSet.set_levels(self.source_order)

    def load_transcript_annotation(self, path:Path):
        """ Load transcript annotation from the GTF file. """
        index_dir = IndexDir(path)
        coding_txs = index_dir.load_coding_tx()
        self.coding_txs = coding_txs
        anno = index_dir.load_annotation()
        for tx_id in anno.transcripts:
            tx_model = anno.transcripts[tx_id]
            tx_anno = TranscriptAnnotation(
                transcript_id=tx_id, gene_id=tx_model.gene_id,
                gene_symbol=tx_model.gene_name, gene_biotype=tx_model.gene_type,
                is_protein_coding=bool(tx_id in coding_txs)
            )
            self.transcript_annotation[tx_id] = tx_anno

    def add_source_order(self, source:str):
        """ Add a source to source order """
        if source in self.source_order:
            return
        if not self.source_order:
            self.source_order[source] = 0
        else:
            self.source_order[source] = max(self.source_order.values()) + 1

    def load_gvf_annotation(self, sample_id:str):
        """ Load the GVF info. """
        data:Dict[str,Dict[str,str]] = self.gvf_info.setdefault(sample_id, {})
        label_map = self.label_map.setdefault(sample_id, LabelSourceMapping())
        if not self.static_source_order:
            self.source_order = {}
        for gvf_file in self.gvfs[sample_id]:
            with open(gvf_file, 'rt', encoding='UTF-8') as handle:
                metadata = GVFMetadata.parse(handle)
                source = metadata.source
                # if source in self.group_map:
                #     source = self.group_map[source]

                if not self.static_source_order:
                    self.add_source_order(source)

                line:str
                for line in handle:
                    if line.startswith('#'):
                        continue
                    fields = line.rstrip().split('\t')
                    gene_id = fields[0]
                    var_id = fields[2]
                    info = fields[-1]

                    info_map = {}
                    for info_field in info.split(';'):
                        key, val = info_field.split('=')
                        info_map[key] = val

                    tx_id = info_map['TRANSCRIPT_ID']

                    if tx_id not in data:
                        data[tx_id] = {}

                    data[tx_id][var_id] = info_map

                    label_map.add_record(gene_id, var_id, source)

    def load_orf_annotation(self, path:Path):
        """ Load the orf annotation from the noncoding orf sequence FASTA file. """
        for seq in SeqIO.parse(path, 'fasta'):
            fields = seq.description.split('|')
            if len(fields) == 3:
                tx_id, orf_id, orf = fields
            else:
                tx_id, _, orf_id, orf = fields
            start, end = tuple(int(x) for x in orf.split('-'))
            orf_anno = ORFAnnoataion(orf_id=orf_id, start=start, end=end)
            if tx_id in self.orf_annotation:
                self.orf_annotation[tx_id][orf_id] = orf_anno
            else:
                self.orf_annotation[tx_id] = {orf_id: orf_anno}

    def load_individual_fasta_headers(self, plex_id:str, sample_id:str):
        """ Load the FASTA headers from each individual FASTA file. """
        fasta_file = self.individual_dbs[plex_id][sample_id]
        sample_variant_labels = self.variant_labels.setdefault(sample_id, set())
        with open(fasta_file, 'rt', encoding='UTF-8') as handle:
            line:str
            for line in handle:
                if not line.startswith('>'):
                    continue
                for label in line.lstrip('>').rstrip().split(' '):
                    sample_variant_labels.add(label)

    def load_merged_fasta_headers(self, plex_id:str, db_id:str):
        """ Load the original FASTA header before encoding from fasta dict.
        Should be called only for label free. """
        sample_variant_labels = self.variant_labels.setdefault(plex_id, set())
        fasta_file = Path(self.merged_dbs[plex_id][db_id])
        fasta_dict = self.find_fasta_dict(fasta_file)
        with open(fasta_dict, 'rt', encoding='UTF-8') as handle:
            for line in handle:
                _, val = line.rstrip().split('\t')
                for label in val.split(' '):
                    sample_variant_labels.add(label)

    def load_quantified_peptides(self, path:Path):
        """ Load peptide quantity matrix """
        self.peptides = []
        with open(path, 'rt', encoding='UTF-8') as handle:
            line:str
            while True:
                line = handle.readline()
                if not line or not line.startswith('#'):
                    break
            if not line:
                raise ValueError('Input peptide intensity table is empty')
            sample_names = [x.strip('"') for x in line.rstrip().split('\t')[4:-1]]
            for line in handle:
                fields = [x.strip('"') for x in line.rstrip().split('\t')]
                peptide = fields[0]
                unmod_seq = peptide.strip('.')
                unmod_seq = re.subn(r'\(\S+?\)', '', unmod_seq)[0]

                protein = fields[1]
                # X!Tandem's treats [D]|[P] as a valid cleavage site and cannot be turned off,
                # resulting hits to multiple peptides. Ignoring anything that contains the '/'
                # symbol. See #79
                if '/' in protein:
                    continue

                n_proteins = fields[2]
                charge = int(fields[3])
                fraction = int(fields[-1])
                abundances:Dict[str,float] = {}
                for i, sample_name in enumerate(sample_names):
                    abundances[sample_name] = float(fields[4 + i])
                quantified_peptide = QuantifiedPeptide(
                    peptide=peptide,
                    protein=protein,
                    unmod_seq=unmod_seq,
                    n_proteins=n_proteins,
                    charge=charge,
                    abundances=abundances,
                    fraction=fraction
                )
                self.peptides.append(quantified_peptide)

    def load_detected_peptides(self, path:Path):
        """ Load detected peptides from idXML file """
        self.peptides = []
        prot_ids:List[ProteinIdentification] = []
        pept_ids:List[PeptideIdentification] = []
        IdXMLFile().load(str(path), prot_ids, pept_ids)
        for psm in pept_ids:
            if not psm.getHits():
                continue
            hit:PeptideHit = psm.getHits()[0]
            accessions = [x.decode('utf-8') for x in hit.extractProteinAccessionsSet()]
            # Skipping PSMs to a peptide that have multiple protein accessessions.
            # For example X!Tandem treats [D]|[P] as a valid cleavage site. See #79
            if len(accessions) > 1:
                continue
            protein = accessions[0]
            peptide = hit.getSequence().toString()
            unmod_seq = hit.getSequence().toUnmodifiedString()
            charge = hit.getCharge()
            mz = psm.getMZ()
            psm_id = psm.getMetaValue('spectrum_reference')
            score = hit.getScore()
            score_type = psm.getScoreType()
            id_peptide = IDPeptide(
                peptide=peptide,
                protein=protein,
                unmod_seq=unmod_seq,
                charge=charge,
                mz=mz,
                score=score,
                score_type=score_type,
                psm_id=psm_id
            )
            self.peptides.append(id_peptide)

    @staticmethod
    def find_fasta_dict(path:Path) -> Path:
        """ Find the corresponding *.dict file of the given FASTA file. """
        parent_dir = path.parent
        name = path.name
        fasta_dict = parent_dir/(name + '.dict')
        if not fasta_dict.exists():
            raise ValueError(f"Index file could not be found for {path}")
        return fasta_dict

    def load_peptide_path_csv(self, path:Path):
        """ Load the CSV file for ProteinQuantifier peptide CSV. """
        self.peptide_paths = {}
        with open(path, 'rt') as handle:
            for line in handle:
                plex_id, db_id, csv_file = line.rstrip().split(',')
                if plex_id not in self.peptide_paths:
                    self.peptide_paths[plex_id] = {db_id:csv_file}
                else:
                    self.peptide_paths[plex_id][db_id] = csv_file

    def load_merged_db_csv(self, path:Path):
        """  Load the merged database CSV file. """
        self.merged_dbs = {}
        with open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('sample_name'):
                    continue
                plex_id, _, db_id, db_file = line.rstrip().split(',')
                if plex_id not in self.merged_dbs:
                    self.merged_dbs[plex_id] = {db_id:db_file}
                else:
                    self.merged_dbs[plex_id][db_id] = db_file

    def load_individual_db_csv(self, path:Path):
        """ Load the individual database CSV file. """
        self.individual_dbs = {}
        with open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('sample_id'):
                    continue
                sample_id, plex_id, _, db_path = line.rstrip().split(',')
                if plex_id not in self.individual_dbs:
                    self.individual_dbs[plex_id] = {}
                self.individual_dbs[plex_id][sample_id] = db_path

    def load_gvf_csv(self, path:Path):
        """ Load gvf CSV file. """
        self.gvfs = {}
        with open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('sample_name'):
                    continue
                sample_id, gvf_file = line.rstrip().split(',')
                if sample_id not in self.gvfs:
                    self.gvfs[sample_id] = [gvf_file]
                else:
                    self.gvfs[sample_id].append(gvf_file)

    def load_annotation_data(self, plex_id:str, db_id:str):
        """ Load all necessary data of the plex for quantitation.
        1. Load the merged fasta dict as a dict
        2. Load the individual fasta headers as a set
        3. Load the variant annotation from gvf
        """
        ## merged db
        merged_fasta_path = Path(self.merged_dbs[plex_id][db_id])
        merged_fasta_dict = self.find_fasta_dict(merged_fasta_path)
        self.encode_headers = {}
        with open(merged_fasta_dict, 'rt', encoding='UTF-8') as handle:
            for line in handle:
                key, val = line.rstrip().split('\t')
                self.encode_headers[key] = val.split(' ')

        ## individual DB
        self.gvf_info = {}
        self.variant_labels = {}
        if not self.label_free:
            for sample_id in self.individual_dbs[plex_id]:
                self.load_individual_fasta_headers(plex_id, sample_id)
                if sample_id not in self.gvfs:
                    print(f"No GVF found for sample {sample_id}", flush=True)
                else:
                    self.load_gvf_annotation(sample_id)
        else:
            self.load_merged_fasta_headers(plex_id, db_id)
            if plex_id not in self.gvfs:
                print(f"No GVF found for plex {plex_id}", flush=True)
            else:
                self.load_gvf_annotation(plex_id)

        if not self.static_source_order:
            self.append_order_internal_sources()
            self.set_source_set_level()

    @staticmethod
    def dict_to_line_safely(x:Dict[str,str]) -> str:
        """ Convert dict to line and quote values if they contain any special
        characters. """
        special_characters = (',', ';', '|', ' ', '\t')
        fields = []
        for k,v in x.items():
            if isinstance(v, str) and any(x in v for x in special_characters):
                if '"' in v:
                    v = f"'{v}'"
                else:
                    v = f'"{v}"'
            fields.append(f"{k}={v}")
        return ';'.join(fields)

    def get_gvf_info(self, label:vpi.VariantPeptideIdentifier,
            sample_index:str) -> List[str]:
        """ Get GVF annotation. """
        gvf_info:List[Dict[str,str]] = []
        if isinstance(label, vpi.FusionVariantPeptideIdentifier):
            tx_id = label.first_tx_id
            tx_anno = self.transcript_annotation[tx_id]
            gene_id = tx_anno.gene_id
            gene_symbol = tx_anno.gene_symbol
            fusion_id = label.fusion_id
            label_anno = {
                'VARIANT_ID': fusion_id,
                **self.gvf_info[sample_index][tx_id][fusion_id],
                'ALT_TRANS': False
            }
            gvf_info.append(label_anno)
            for x in label.first_variants:
                grouped_source, raw_source = self.get_variant_source(sample_index, gene_id, x)
                label_anno = {
                    'VARIANT_ID': x,
                    **self.gvf_info[sample_index][label.first_tx_id][x],
                    'GENE_ID': gene_id,
                    'GENE_SYMBOL': gene_symbol,
                    'ALT_TRANS': False,
                    'GROUPED_SOURCE': grouped_source,
                    'RAW_SOURCE': raw_source
                }
                gvf_info.append(label_anno)
            for x in label.second_variants:
                second_tx_anno = self.transcript_annotation[label.second_tx_id]
                second_gene_id = second_tx_anno.gene_id
                second_gene_symbol = second_tx_anno.gene_symbol
                grouped_source, raw_source = self.get_variant_source(
                    sample_index, second_gene_id, x
                )
                label_anno = {
                    'VARIANT_ID': x,
                    **self.gvf_info[sample_index][label.second_tx_id][x],
                    'GENE_ID': second_gene_id,
                    'GENE_SYMBOL': second_gene_symbol,
                    'ALT_TRANS': False,
                    'GROUPED_SOURCE': grouped_source,
                    'RAW_SOURCE': raw_source
                }
                gvf_info.append(label_anno)
            for x in label.peptide_variants:
                grouped_source, raw_source = self.get_variant_source(sample_index, gene_id, x)
                label_anno = {
                    'VARIANT_ID': x,
                    'TRANSCRIPT_ID': tx_id,
                    'GENE_ID': gene_id,
                    'GENE_SYMBOL': gene_symbol,
                    'ALT_TRANS': True,
                    'GROUPED_SOURCE': grouped_source,
                    'RAW_SOURCE': raw_source
                }
                gvf_info.append(label_anno)
        else:
            if isinstance(label, vpi.CircRNAVariantPeptideIdentifier):
                # CIRC-ENST00000614167.2-E2-E1
                circ_rna_id = label.circ_rna_id
                tx_id = circ_rna_id.split('-')[1]
                tx_anno = self.transcript_annotation[tx_id]
                gene_id = tx_anno.gene_id
                gene_symbol = tx_anno.gene_symbol
                grouped_source, raw_source = self.get_variant_source(
                    sample_index, gene_id, circ_rna_id
                )
                label_anno = {
                    'VARIANT_ID': circ_rna_id,
                    'TRANSCRIPT_ID': tx_id,
                    'GENE_ID': gene_id,
                    'GENE_SYMBOL': gene_symbol,
                    **self.gvf_info[sample_index][tx_id][circ_rna_id],
                    'ALT_TRANS': False,
                    'GROUPED_SOURCE': grouped_source,
                    'RAW_SOURCE': raw_source
                }
                gvf_info.append(label_anno)
            else:
                tx_id = label.transcript_id
                tx_anno = self.transcript_annotation[tx_id]
                gene_id = tx_anno.gene_id
                gene_symbol = tx_anno.gene_symbol
            for variant_id in label.variant_ids:
                label_anno = {
                    'VARIANT_ID': variant_id
                }
                if variant_id.split('-')[0] in ALT_TRANS_TYPES:
                    label_anno.update({
                        'ALT_TRANS': True
                    })
                else:
                    label_anno.update({
                        **self.gvf_info[sample_index][tx_id][variant_id],
                        'ALT_TRANS': False
                    })
                grouped_source, raw_source = self.get_variant_source(
                    sample_index, gene_id, variant_id
                )
                label_anno.update({
                    'TRANSCRIPT_ID': tx_id,
                    'GENE_ID': gene_id,
                    'GENE_SYMBOL': gene_symbol,
                    'GROUPED_SOURCE': grouped_source,
                    'RAW_SOURCE': raw_source
                })

                gvf_info.append(label_anno)

        return [self.dict_to_line_safely(x) for x in gvf_info]

    def get_noncoding_annotation(self, label:vpi.NovelORFPeptideIdentifier
            ) -> List[str]:
        """ Get noncoding annotation """
        tx_id = label.transcript_id
        tx_anno = self.transcript_annotation[tx_id]
        orf_anno = self.orf_annotation[tx_id][label.orf_id]
        raw_source = constant.SOURCE_NOVEL_ORF
        try:
            grouped_source = self.group_map[raw_source]
        except KeyError:
            grouped_source = raw_source
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'GENE_ID': tx_anno.gene_id,
            'GENE_SYMBOL': tx_anno.gene_symbol,
            'GENE_BIOTYPE': tx_anno.gene_biotype,
            'ORF_ID': orf_anno.orf_id,
            'ORF_START': orf_anno.start,
            'ORF_END': orf_anno.end,
            'GROUPED_SOURCE': grouped_source,
            'RAW_SOURCE': raw_source
        }
        return [self.dict_to_line_safely(attrs)]

    def get_variant_source(self, sample_id:str, gene_id:str, label:str) -> Tuple[str,str]:
        """ Get variant source """
        if label.split('-')[0] == constant.SEC_TERMINATION_TYPE:
            raw_source = constant.SOURCE_SEC_TERMINATION
        elif label.split('-')[0] in constant.CODON_REASSIGNMENTS_TYPES:
            raw_source = constant.SOURCE_CODON_REASSIGNMENT
        else:
            raw_source = self.label_map[sample_id].get_source(gene_id, label)

        try:
            grouped_source = self.group_map[raw_source]
        except KeyError:
            grouped_source = raw_source
        return grouped_source, raw_source

    def get_peptide_sources(self, var_id:vpi.VariantPeptideIdentifier,
            sample_id:str) -> Tuple[VariantSourceSet, VariantSourceSet]:
        """ Get peptide sources """
        grouped_sources = VariantSourceSet()
        raw_sources = VariantSourceSet()
        if isinstance(var_id, vpi.NovelORFPeptideIdentifier):
            tx_id = var_id.transcript_id
            gene_ids = [var_id.gene_id]
            var_ids = {gene_ids[0]: var_id.codon_reassigns}

        elif isinstance(var_id, vpi.CircRNAVariantPeptideIdentifier):
            circ_rna_id = var_id.circ_rna_id
            tx_id = circ_rna_id.split('-', 2)[1]
            gene_ids = [self.transcript_annotation[tx_id].gene_id]
            var_ids = {gene_ids[0]: [circ_rna_id, *var_id.variant_ids]}

        elif isinstance(var_id, vpi.FusionVariantPeptideIdentifier):
            first_tx_id = var_id.first_tx_id
            first_gene_id = self.transcript_annotation[first_tx_id].gene_id
            second_tx_id = var_id.second_tx_id
            second_gene_id = self.transcript_annotation[second_tx_id].gene_id
            gene_ids = [first_gene_id, second_gene_id]
            var_ids = {
                first_gene_id: var_id.first_variants + [var_id.fusion_id],
                second_gene_id: var_id.second_variants
            }
            tx_id = first_tx_id

        elif isinstance(var_id, vpi.BaseVariantPeptideIdentifier):
            gene_id = self.transcript_annotation[var_id.transcript_id].gene_id
            gene_ids = [gene_id]
            var_ids = {gene_id: var_id.variant_ids}
            tx_id = var_id.transcript_id

        if self.transcript_annotation[tx_id].is_protein_coding is False:
            raw_sources.add(constant.SOURCE_NOVEL_ORF)

        for gene_id, _ids in var_ids.items():
            for _id in _ids:
                if _id.split('-')[0] == constant.SEC_TERMINATION_TYPE:
                    raw_sources.add(constant.SOURCE_SEC_TERMINATION)
                elif _id.split('-')[0] in constant.CODON_REASSIGNMENTS_TYPES:
                    raw_sources.add(constant.SOURCE_CODON_REASSIGNMENT)
                else:
                    source = self.label_map[sample_id].get_source(gene_id, _id)
                    raw_sources.add(source)

        for source in raw_sources:
            if source in self.group_map:
                grouped_sources.add(self.group_map[source])
            else:
                grouped_sources.add(source)
        return grouped_sources, raw_sources

    def annotate(self):
        """ Annotate quantified peptides. """
        gene_id_placeholder = 'XXXXX'
        self.peptide_annotations = {}
        self.tally.n_total_peptides = len(self.peptides)
        for peptide in self.peptides:
            seq = peptide.peptide
            fasta_header = self.encode_headers[peptide.protein]
            peptide_labels:List[vpi.VariantPeptideIdentifier] = []
            for group in fasta_header:
                fields = group.split('|')
                if fields[-2].startswith('ORF') and len(fields) == 3:
                    fields.insert(1, gene_id_placeholder)
                    group = '|'.join(fields)
                label = vpi.parse_variant_peptide_id(group, self.coding_txs)
                if isinstance(label[0], vpi.NovelORFPeptideIdentifier) \
                        and label[0].gene_id == gene_id_placeholder:
                    label[0].gene_id = None
                peptide_labels.append(label[0])

            peptide_annotated = False
            for sample_name in self.variant_labels:
                sample_variant_labels = self.variant_labels[sample_name]
                for label in peptide_labels:
                    # FIXME: NovelORF with variants also need noncoding annotation
                    if isinstance(label, vpi.NovelORFPeptideIdentifier):
                        var_anno = self.get_noncoding_annotation(label)
                    else:
                        is_alt_trans_only = \
                            isinstance(label, vpi.BaseVariantPeptideIdentifier) \
                            and all(v.split('-')[0] in ALT_TRANS_TYPES for v in label.variant_ids)
                        if not is_alt_trans_only \
                                and str(label) not in sample_variant_labels:
                            continue
                        try:
                            var_anno = self.get_gvf_info(label, sample_name)
                        except KeyError as e:
                            print(f"Failed to load GVF info for {label}")
                            raise KeyError from e

                    peptide_annotated = True

                    if seq not in self.peptide_annotations:
                        self.peptide_annotations[seq] = {}

                    grouped_sources, raw_sources = self.get_peptide_sources(label, sample_name)
                    header_anno = FastaHeaderAnnotation(
                            annotation=var_anno,
                            raw_sources=raw_sources,
                            grouped_sources=grouped_sources
                        )

                    if sample_name not in self.peptide_annotations[seq]:
                        anno = PeptideAnnotation(
                            peptide=seq,
                            unmod_seq=peptide.unmod_seq,
                            sample_name=sample_name,
                            fasta_header_annotation={str(label): header_anno},
                            raw_sources=raw_sources,
                            grouped_sources=grouped_sources
                        )
                        self.peptide_annotations[seq][sample_name] = anno
                    else:
                        anno:PeptideAnnotation = self.peptide_annotations[seq][sample_name]
                        label_str = str(label)
                        if label_str in anno.fasta_header_annotation:
                            continue
                        anno.fasta_header_annotation[label_str] = header_anno
                        if grouped_sources < anno.grouped_sources:
                            anno.grouped_sources = grouped_sources
                            anno.raw_sources = raw_sources
                        elif grouped_sources == anno.grouped_sources \
                                and raw_sources < anno.raw_sources:
                            anno.raw_sources = raw_sources

            if peptide_annotated:
                self.tally.n_peptides_annotated += 1
            else:
                self.tally.warn_unannotated(peptide.peptide)

    def write(self, path:Path):
        """ Write out to file. """
        with open(path, 'wt', encoding='UTF-8') as handle:
            handle.write(PeptideAnnotation.table_headers() + '\n')
            for sample_anno in self.peptide_annotations.values():
                for anno in sample_anno.values():
                    for line in anno.to_lines():
                        handle.write(line + '\n')
