""" Set database path of an idXML file
"""
import argparse
from pyopenms import IdXMLFile
from msptools.cli import common


INPUT_FILE_FORMATS = ['.idXML']
DATABASE_FORMATS = ['.fasta', '.fa']
OUTPUT_FILE_FORMATS = ['.idXML']

# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='setDatabasePath',
        help='Set database file path of an idXML file.',
        description='Set database file path of an idXML file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input-path',
        type=str,
        help=f"Path to input idXML file. Valid formats: {INPUT_FILE_FORMATS}",
        required=True
    )
    parser.add_argument(
        '--database',
        type=str,
        help=f"Path to the database FASTA file. Valid formats: {DATABASE_FORMATS}",
        required=True
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help=f"Path to the output idXML file.  Valid formats: {OUTPUT_FILE_FORMATS}",
        required=True
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args):
    """ main """
    input_path:str = args.input_path
    output_path:str = args.output_path
    database:str = args.database

    common.validate_file_format(input_path, INPUT_FILE_FORMATS, True)
    common.validate_file_format(database, DATABASE_FORMATS, True)
    common.validate_file_format(output_path, OUTPUT_FILE_FORMATS)

    peptides = []
    proteins = []

    IdXMLFile().load(input_path, proteins, peptides)

    for prot_ident in proteins:
        params = prot_ident.getSearchParameters()
        params.db = database
        prot_ident.setSearchParameters(params)

    IdXMLFile().store(output_path, proteins, peptides)
