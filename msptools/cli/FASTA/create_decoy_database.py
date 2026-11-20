"""Create a combined target/decoy search database."""
import argparse
from msptools import logger
from msptools.FASTA.DecoyDatabase import DecoyDatabase


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
    parser.add_argument('--enzyme',
                        help='OpenMS enzyme name for peptide reversal. If omitted, reverse'
                        ' whole proteins while preserving an initial M.')
    parser.set_defaults(func=main)
    return parser


def main(args):
    """Build the database and print the record tally."""
    database = DecoyDatabase(
        input_path=args.input_path, output_path=args.output_path,
        include_crap=args.include_crap, decoy_prefix=args.decoy_prefix,
        contam_prefix=args.contam_prefix, enzyme=args.enzyme)
    counts = database.build()
    for label, count in counts.items():
        logger(f'{label}: {count}')
