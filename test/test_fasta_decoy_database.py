"""FASTA database construction and command-line integration."""
import hashlib
from importlib.resources import files
import subprocess
import sys

from Bio import SeqIO
import pytest
from msptools.FASTA.DecoyDatabase import create_decoy_database


def test_targets_decoys_and_headers(tmp_path):
    source, output = tmp_path/'input.fa', tmp_path/'output.fasta'
    source.write_text('>sp|P12345|TEST Protein description\nACDEKFGHR*\n>A short\nA\n')
    counts = create_decoy_database(source, output)
    records = list(SeqIO.parse(output, 'fasta'))
    assert [r.id for r in records] == ['sp|P12345|TEST', 'DECOY_sp|P12345|TEST',
                                      'A', 'DECOY_A']
    assert records[1].description == 'DECOY_sp|P12345|TEST Protein description'
    assert str(records[1].seq) == 'EDCAKRHGF*'
    assert counts == {'Targets': 2, 'Contaminants': 0, 'Decoys': 2,
                      'Decoys identical to their target': 1}


def test_bundled_contaminants_and_cli(tmp_path):
    source, output = tmp_path/'input.fa', tmp_path/'output.fa'
    source.write_text('>TEST example\nACDEKFGHR\n')
    result = subprocess.run([
        sys.executable, '-m', 'msptools.cli', 'FASTA', 'CreateDecoyDatabase',
        '--input-path', str(source), '--output-path', str(output),
        '--include-crap', '--decoy-prefix', 'REV_', '--contam-prefix', 'CRAP_'],
        check=True, capture_output=True, text=True)
    records = list(SeqIO.parse(output, 'fasta'))
    assert len(records) == 234
    assert records[2].id == 'CRAP_sp|ALBU_BOVIN|'
    assert records[3].id == 'REV_CRAP_sp|ALBU_BOVIN|'
    assert len({r.id for r in records}) == len(records)
    assert 'Contaminants: 116' in result.stdout
    assert 'Decoys: 117' in result.stdout
    data = files('msptools.data').joinpath('crap.fasta').read_bytes()
    assert hashlib.sha256(data).hexdigest() == (
        '4b0e6e97ab1d618baa38be612a787d884a781e07f627b919c4b29ac064db5382')


@pytest.mark.parametrize('content', [
    '', 'not FASTA\n', '>\nACDE\n', '>A\n', '>A\nAC-D\n', '>A\nAC*DE\n',
    '>DECOY_A\nACDE\n', '>A\nACDE\n>A\nFGHI\n',
    '>A\nACDE\n>DECOY_A\nFGHI\n',
])
def test_invalid_inputs_preserve_existing_output(tmp_path, content):
    source, output = tmp_path/'input.fa', tmp_path/'output.fa'
    source.write_text(content)
    output.write_text('keep existing output')
    with pytest.raises(ValueError):
        create_decoy_database(source, output)
    assert output.read_text() == 'keep existing output'
    assert not list(tmp_path.glob('.output.fa.*'))


@pytest.mark.parametrize('options', [
    {'decoy_prefix': ''}, {'contam_prefix': 'BAD PREFIX'},
    {'decoy_prefix': '>'}, {'decoy_prefix': 'CONTAM_'},
    {'contam_prefix': 'DECOY_CONTAM_'}, {'enzyme': 'not_an_enzyme'},
])
def test_invalid_options(tmp_path, options):
    source = tmp_path/'input.fa'
    source.write_text('>A\nACDE\n')
    with pytest.raises(ValueError):
        create_decoy_database(source, tmp_path/'output.fa', **options)


def test_input_output_aliases_rejected(tmp_path):
    source = tmp_path/'input.fa'
    source.write_text('>A\nACDE\n')
    alias = tmp_path/'alias.fa'
    alias.hardlink_to(source)
    for output in (source, alias):
        with pytest.raises(ValueError, match='paths must differ'):
            create_decoy_database(source, output)
    assert source.read_text() == '>A\nACDE\n'


def test_contaminant_accession_collision(tmp_path):
    source = tmp_path/'input.fa'
    source.write_text('>CONTAM_sp|ALBU_BOVIN|\nACDE\n')
    with pytest.raises(ValueError, match='Duplicate output accession'):
        create_decoy_database(source, tmp_path/'output.fa', include_crap=True)
    assert not (tmp_path/'output.fa').exists()


def test_enzyme_changes_decoy(tmp_path):
    source = tmp_path/'input.fa'
    source.write_text('>A\nACDEKFGHR\n')
    output = tmp_path/'output.fa'
    create_decoy_database(source, output, enzyme='no cleavage')
    records = list(SeqIO.parse(output, 'fasta'))
    assert str(records[1].seq) == 'RHGFKEDCA'
