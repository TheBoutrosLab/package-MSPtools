""" MSPtools MzML """
import argparse
from msptools.cli import common
from msptools.cli.MzML import filter as mzml_filter


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='MzML',
        help='Manipulate mzML files.',
        description='Manipulate mzML files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )

    mzml_filter.parse_args(sub_subparsers)

    common.print_help_if_missing_args(parser)
