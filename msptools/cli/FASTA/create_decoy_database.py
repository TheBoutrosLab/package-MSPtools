"""Create a combined target/decoy search database."""
import argparse
from msptools import logger
from msptools.FASTA.DecoyDatabase import create_decoy_database


def add_subparsers(subparsers):
    """Register the database builder."""
    parser = subparsers.add_parser(
        'CreateDecoyDatabase', help='Build a target/decoy protein FASTA database.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input-path', required=True, help='Input protein FASTA file.')
    parser.add_argument('--output-path', required=True, help='Output combined FASTA file.')
    parser.add_argument('--include-crap', action='store_true',
                        help='Add bundled cRAP contaminants and their decoys.')
    parser.add_argument('--decoy-prefix', default='DECOY_', help='Decoy accession prefix.')
    parser.add_argument('--contam-prefix', default='CONTAM_',
                        help='Added contaminant accession prefix.')
    parser.add_argument('--enzyme', default='Trypsin',
                        help='OpenMS enzyme name for peptide reversal.')
    parser.set_defaults(func=main)
    return parser


def main(args):
    """Build the database and print the record tally."""
    counts = create_decoy_database(
        args.input_path, args.output_path, args.include_crap,
        args.decoy_prefix, args.contam_prefix, args.enzyme)
    for label, count in counts.items():
        logger(f'{label}: {count}')
