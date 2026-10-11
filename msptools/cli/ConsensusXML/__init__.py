""" MSPtools ConsensusXML """
import argparse
from msptools.cli import common
from msptools.cli.ConsensusXML import set_column_headers, mad_normalizer, remove_na


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='ConsensusXML',
        help='Manipulate ConsensusXML files.',
        description='Manipulate ConsensusXML files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )

    set_column_headers.add_subparsers(sub_subparsers)
    mad_normalizer.add_subparsers(sub_subparsers)
    remove_na.add_subparsers(sub_subparsers)

    common.print_help_if_missing_args(parser)
