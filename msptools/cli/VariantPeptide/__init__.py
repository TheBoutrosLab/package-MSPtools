""" MSPtools VariantPeptide """
import argparse
from msptools.cli import common
from msptools.cli.VariantPeptide import annotate


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='VariantPeptide',
        help='Handle variant peptides.',
        description='Handle variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )
    annotate.parse_args(sub_subparsers)

    common.print_help_if_missing_args(parser)
