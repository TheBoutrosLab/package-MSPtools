"""Build a combined target and decoy protein database."""
from collections import Counter
from importlib.resources import files
from pathlib import Path
import os
import tempfile

from Bio import SeqIO
from pyopenms import AASequence, DecoyGenerator, ProteaseDB


class DecoyDatabase:
    """Build a combined FASTA using whole-protein or peptide reversal."""

    def __init__(
            self, input_path, output_path, include_crap=False, decoy_prefix='DECOY_',
            contam_prefix='CONTAM_', enzyme=None
    ):
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)
        self.include_crap = include_crap
        self.decoy_prefix = decoy_prefix
        self.contam_prefix = contam_prefix
        self.enzyme = enzyme
        self.generator = DecoyGenerator()
        self.accessions = set()
        self.counts = {}

    def validate(self):
        """Check paths and settings before opening an output file."""
        if (self.input_path.resolve() == self.output_path.resolve()
                or (self.output_path.exists()
                    and os.path.samefile(self.input_path, self.output_path))):
            raise ValueError('Input and output paths must differ')
        for label, prefix in (('decoy', self.decoy_prefix),
                              ('contaminant', self.contam_prefix)):
            if not prefix or any(c.isspace() or c == '>' for c in prefix):
                raise ValueError(f'The {label} prefix must be nonempty without whitespace or >')
        if (self.decoy_prefix.startswith(self.contam_prefix)
                or self.contam_prefix.startswith(self.decoy_prefix)):
            raise ValueError('Decoy and contaminant prefixes must not overlap')
        if self.enzyme is not None and not ProteaseDB().hasEnzyme(self.enzyme):
            raise ValueError(f'Unknown enzyme: {self.enzyme}')

    def reverse_sequence(self, sequence, accession):
        """Reverse the sequence while preserving its optional terminal stop."""
        # OpenMS converts * to X; preserve an optional terminal stop explicitly.
        protein = sequence[:-1] if sequence.endswith('*') else sequence
        if not protein or any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO' for c in protein):
            raise ValueError(f'Invalid protein sequence for {accession}')
        if self.enzyme is None:
            # Preserve an initial methionine, as in Philosopher's default mode.
            if protein.startswith('M'):
                tail = protein[1:]
                decoy = 'M' + (self.generator.reverseProtein(
                    AASequence.fromString(tail)).toString() if tail else '')
            else:
                decoy = self.generator.reverseProtein(AASequence.fromString(protein)).toString()
        else:
            decoy = self.generator.reversePeptides(
                AASequence.fromString(protein), self.enzyme).toString()
        if Counter(decoy) != Counter(protein):
            raise ValueError(f'Decoy changed amino acid composition for {accession}')
        return decoy + ('*' if sequence.endswith('*') else '')

    @staticmethod
    def write_record(output, accession, description, sequence):
        """Write a FASTA entry with wrapped sequence lines."""
        output.write(f'>{accession}{description}\n')
        for start in range(0, len(sequence), 60):
            output.write(sequence[start:start + 60] + '\n')

    def write_records(self, source, output, contaminant=False):
        """Validate and write each target followed by its decoy."""
        for record in SeqIO.parse(source, 'fasta'):
            if not record.id:
                raise ValueError('Every FASTA record must have an accession')
            if record.id.startswith(self.decoy_prefix):
                raise ValueError(f'Input already contains a decoy accession: {record.id}')
            accession = (self.contam_prefix if contaminant else '') + record.id
            decoy_accession = self.decoy_prefix + accession
            for identifier in (accession, decoy_accession):
                if identifier in self.accessions:
                    raise ValueError(f'Duplicate output accession: {identifier}')
                self.accessions.add(identifier)
            sequence = str(record.seq).upper()
            decoy = self.reverse_sequence(sequence, record.id)
            description = record.description[len(record.id):]
            for identifier, value in ((accession, sequence), (decoy_accession, decoy)):
                self.write_record(output, identifier, description, value)
            self.counts['Contaminants' if contaminant else 'Targets'] += 1
            self.counts['Decoys'] += 1
            self.counts['Decoys identical to their target'] += int(sequence == decoy)

    def build(self):
        """Write the complete database atomically and return its record tally."""
        self.validate()
        self.accessions = set()
        self.counts = {'Targets': 0, 'Contaminants': 0, 'Decoys': 0,
                       'Decoys identical to their target': 0}
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                    dir=self.output_path.parent, prefix=f'.{self.output_path.name}.',
                    delete=False) as output:
                temporary = Path(output.name)
                with self.input_path.open(encoding='utf-8') as source:
                    # SeqIO ignores leading non-FASTA text; reject that explicitly.
                    first = next((line for line in source if line.strip()), '')
                    if not first.startswith('>'):
                        raise ValueError('Input must be a nonempty protein FASTA file')
                    source.seek(0)
                    self.write_records(source, output)
                if self.include_crap:
                    with files('msptools.data').joinpath('crap.fasta').open(
                            encoding='utf-8') as source:
                        self.write_records(source, output, contaminant=True)
            os.replace(temporary, self.output_path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return self.counts.copy()
