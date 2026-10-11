""" MSPtools IdXML """
import argparse
from msptools.cli import common
from msptools.cli.IdXML import filter_protein_groups, set_database_path


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='IdXML',
        help='Manipulate idXML files.',
        description='Manipulate idXML files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )
    filter_protein_groups.parse_args(sub_subparsers)
    set_database_path.parse_args(sub_subparsers)

    common.print_help_if_missing_args(parser)
