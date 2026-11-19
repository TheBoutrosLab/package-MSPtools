"""Build a combined target and decoy protein database."""
from collections import Counter
from importlib.resources import files
from pathlib import Path
import os
import tempfile

from Bio import SeqIO
from pyopenms import AASequence, DecoyGenerator, ProteaseDB


def create_decoy_database(input_path, output_path, include_crap=False,
                          decoy_prefix='DECOY_', contam_prefix='CONTAM_',
                          enzyme='Trypsin'):
    """Write targets and enzyme-specific reversed decoys, returning counts."""
    input_path, output_path = Path(input_path), Path(output_path)
    if (input_path.resolve() == output_path.resolve()
            or (output_path.exists() and os.path.samefile(input_path, output_path))):
        raise ValueError('Input and output paths must differ')
    for label, prefix in (('decoy', decoy_prefix), ('contaminant', contam_prefix)):
        if not prefix or any(c.isspace() or c == '>' for c in prefix):
            raise ValueError(f'The {label} prefix must be nonempty without whitespace or >')
    if decoy_prefix.startswith(contam_prefix) or contam_prefix.startswith(decoy_prefix):
        raise ValueError('Decoy and contaminant prefixes must not overlap')
    if not ProteaseDB().hasEnzyme(enzyme):
        raise ValueError(f'Unknown enzyme: {enzyme}')

    generator = DecoyGenerator()
    accessions = set()
    counts = {'Targets': 0, 'Contaminants': 0, 'Decoys': 0,
              'Decoys identical to their target': 0}

    def write_records(source, output, contaminant=False):
        for record in SeqIO.parse(source, 'fasta'):
            if not record.id:
                raise ValueError('Every FASTA record must have an accession')
            if record.id.startswith(decoy_prefix):
                raise ValueError(f'Input already contains a decoy accession: {record.id}')
            accession = (contam_prefix if contaminant else '') + record.id
            decoy_accession = decoy_prefix + accession
            for identifier in (accession, decoy_accession):
                if identifier in accessions:
                    raise ValueError(f'Duplicate output accession: {identifier}')
                accessions.add(identifier)
            sequence = str(record.seq).upper()
            # OpenMS converts * to X; preserve an optional terminal stop explicitly.
            protein = sequence[:-1] if sequence.endswith('*') else sequence
            if not protein or any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO' for c in protein):
                raise ValueError(f'Invalid protein sequence for {record.id}')
            decoy = generator.reversePeptides(AASequence.fromString(protein), enzyme).toString()
            if Counter(decoy) != Counter(protein):
                raise ValueError(f'Decoy changed amino acid composition for {record.id}')
            if sequence.endswith('*'):
                decoy += '*'
            description = record.description[len(record.id):]
            for identifier, value in ((accession, sequence), (decoy_accession, decoy)):
                output.write(f'>{identifier}{description}\n')
                for start in range(0, len(value), 60):
                    output.write(value[start:start + 60] + '\n')
            counts['Contaminants' if contaminant else 'Targets'] += 1
            counts['Decoys'] += 1
            counts['Decoys identical to their target'] += int(sequence == decoy)

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                dir=output_path.parent, prefix=f'.{output_path.name}.', delete=False) as output:
            temporary = Path(output.name)
            with input_path.open(encoding='utf-8') as source:
                # SeqIO ignores leading non-FASTA text; reject that explicitly.
                first = next((line for line in source if line.strip()), '')
                if not first.startswith('>'):
                    raise ValueError('Input must be a nonempty protein FASTA file')
                source.seek(0)
                write_records(source, output)
            if include_crap:
                with files('msptools.data').joinpath('crap.fasta').open(
                        encoding='utf-8') as source:
                    write_records(source, output, contaminant=True)
        os.replace(temporary, output_path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return counts
