""" MSPtools FeatureXML """
import argparse
from msptools.cli import common
from msptools.cli.FeatureXML import set_spectra_path, featurexml_filter


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='FeatureXML',
        help='Manipulate featureXML files.',
        description='Manipulate featureXML files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub_subparsers = parser.add_subparsers(
        dest='subcommand'
    )

    set_spectra_path.parse_args(sub_subparsers)
    featurexml_filter.parse_args(sub_subparsers)

    common.print_help_if_missing_args(parser)
