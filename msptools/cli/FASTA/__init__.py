"""FASTA commands."""
import argparse
from msptools.cli import common
from msptools.cli.FASTA import create_decoy_database


def parse_args(subparsers):
    """Register FASTA commands."""
    parser = subparsers.add_parser('FASTA', help='Manipulate protein FASTA files.',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    create_decoy_database.add_subparsers(parser.add_subparsers(dest='subcommand'))
    common.print_help_if_missing_args(parser)
