""" This script reads in a featureXML file and fixes the `spectra_file`
UserParam field. """
import argparse
from pyopenms import FeatureMap, FeatureXMLFile
from msptools.cli import common
from msptools.FileTypes import FileTypes


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser = argparse.ArgumentParser()
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='setSpectraPath',
        help='Set spectra file path of a featureXML file.',
        description='Set spectra file path of a featureXML file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input-path',
        type=str,
        help='Path to the input featureXML file.',
        required=True
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help='Path to the output featureXML file.',
        required=True
    )
    parser.add_argument(
        '--spectra-path',
        type=str,
        help='Path to the mzML file.',
        required=True
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args):
    """ File spectra file """
    input_path:str = args.input_path
    output_path:str = args.output_path

    common.validate_input_file(input_path, [FileTypes.FEATUREXML])
    common.validate_output_file(output_path, [FileTypes.FEATUREXML])

    fmap = FeatureMap()
    FeatureXMLFile().load(args.input_path, fmap)

    spectra_path = args.spectra_path.encode('utf-8')
    fmap.setPrimaryMSRunPath([spectra_path])

    FeatureXMLFile().store(args.output_path, fmap)
