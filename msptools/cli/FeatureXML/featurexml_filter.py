""" Removes features from a featureXML file. """
import argparse
from pyopenms import FeatureXMLFile, DataProcessing
from msptools.cli import common
from msptools.FeatureXML.Filter import FeatureXMLFilter
from msptools.FileTypes import FileTypes


PROG_NAME = 'filter'

# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name=PROG_NAME,
        help='Filter features from featureXML.',
        description='Filter features from featureXML.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--input-path',
        type=str,
        help="Path to the input file. Valid format: ['.featureXML']"
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help="Path to the output file. Valid format: ['.featureXML']"
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ Main """
    common.validate_input_file(args.input_path, [FileTypes.FEATUREXML])
    common.validate_output_file(args.output_path, [FileTypes.FEATUREXML])

    filt = FeatureXMLFilter()

    FeatureXMLFile().load(args.input_path, filt.fmap)
    print(f"Number of spectra loaded: {filt.fmap.size()}", flush=True)

    filt.filter_feature_with_internal_id()

    process_actions = { DataProcessing.ProcessingAction.FILTERING }
    common.add_data_processing(filt.fmap_filtered, args, PROG_NAME, process_actions)

    FeatureXMLFile().store(args.output_path, filt.fmap_filtered)

    print(f"Number of spectra saved: {filt.fmap_filtered.size()}", flush=True)
