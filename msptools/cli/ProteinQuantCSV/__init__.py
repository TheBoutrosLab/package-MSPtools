""" MSPtools ProteinQuantCSV """
import argparse
from msptools.cli import common
from msptools.cli.ProteinQuantCSV import reheader, fix_zero_plex


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='ProteinQuantCSV',
        help='Manipulate protein quant CSV produced by ProteinQuantifier.',
        description='Manipulate protein quant CSV produced by ProteinQuantifier.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )
    reheader.add_subparsers(sub_subparsers)
    fix_zero_plex.add_subparsers(sub_subparsers)

    common.print_help_if_missing_args(parser)
