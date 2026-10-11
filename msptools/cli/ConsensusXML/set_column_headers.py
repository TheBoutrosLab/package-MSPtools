""" This module takes a consensusXML file and set the column headers as
specified from command line. """
import argparse
from typing import Dict
from pyopenms import ConsensusMap, ConsensusXMLFile, ColumnHeader
from msptools.cli import common


INPUT_FILE_FORMATS = ['.consensusXML']
OUTPUT_FILE_FORMATS = ['.consensusXML']

# pylint: disable=W0212
def add_subparsers(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='setColumnHeaders',
        help='Set column headers.',
        description='Set column headers.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-path',
        type=str,
        help=f"Path to the input ConsensusXML file. Valid formats: {INPUT_FILE_FORMATS}",
        required=True
    )
    parser.add_argument(
        '-o', '--output-path',
        type=str,
        help=f"Path to the output idXML file. Valid formats: {OUTPUT_FILE_FORMATS}",
        required=True
    )
    parser.add_argument(
        '--column-headers',
        type=str,
        nargs='+',
        help="Column header names to set.",
        required=True
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ main """
    input_path:str = args.input_path
    output_path:str = args.output_path

    common.validate_file_format(input_path, INPUT_FILE_FORMATS, True)
    common.validate_file_format(output_path, OUTPUT_FILE_FORMATS)

    cmap = ConsensusMap()
    ConsensusXMLFile().load(input_path, cmap)

    if len(cmap.getColumnHeaders()) != len(args.column_headers):
        raise ValueError(
            'The size of --names does not match with the number of columns in'
            f" consensusXML. The consensusXML file has {len(cmap.getColumnHeaders())}"
            ' columns'
        )

    headers:Dict[int,ColumnHeader] = {}

    for i,(key, header) in enumerate(cmap.getColumnHeaders().items()):
        header.filename = args.column_headers[i]
        headers[key] = header

    cmap.setColumnHeaders(headers)

    ConsensusXMLFile().store(output_path, cmap)
