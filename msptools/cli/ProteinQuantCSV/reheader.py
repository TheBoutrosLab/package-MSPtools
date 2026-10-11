""" Replace the headers of the CSV/TSV file output by ProteinQuantifier with
the channel description. """
import argparse
from pyopenms import ConsensusXMLFile
from msptools.cli import common
from msptools.ProteinQuantCSV.Reheader import ProteinQuantCSVReheaderLFQ, ProteinQuantCSVReheaderTMT


INPUT_FILE_FORMATS = ['.csv']
CONSENSUS_XML_FORMATS = ['.consensusXML']
OUTPUT_FILE_FORMATS = ['.csv']

# pylint: disable=W0212
def add_subparsers(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='reheader',
        help='Replace the headers of the ProteinQuantifier output CSV with the'
        ' actual sample names.',
        description='Replace the headers of the ProteinQuantifier output CSV with the'
        ' actual sample names.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-path',
        type=str,
        help=f"Path to the protein quant CSV file. Valid formats: {INPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    parser.add_argument(
        '--consensusXML',
        type=str,
        help=f"Path to the consensusXML file. Valid formats: {CONSENSUS_XML_FORMATS}",
        metavar='<file>'
    )
    parser.add_argument(
        '-o', '--output-path',
        type=str,
        help=f"Path to the output CSV file. Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    parser.add_argument(
        '--header',
        type=str,
        nargs='*',
        help='Headers to replace.',
        metavar='<value>',
        default=None
    )
    parser.add_argument(
        '--experiment-type',
        type=str,
        choices=['LFQ', 'TMT'],
        help="Experiment type. Valid options: ['LFQ', 'TMT]",
        metavar='<choice>'
    )
    parser.add_argument(
        '--separator',
        type=str,
        default='\t',
        help='Character(s) used to separate fields.',
        metavar='<value>'
    )
    parser.add_argument(
        '--quoting',
        type=str,
        choices=['none', 'double', 'escape'],
        default='double',
        help='Method for quoting of strings',
        metavar='<choice>'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ Protein quantification CSV reheader """
    common.validate_file_format(args.input_path, INPUT_FILE_FORMATS, check_exist=True)
    common.validate_file_format(args.output_path, OUTPUT_FILE_FORMATS, check_exist=False)

    if args.experiment_type == 'LFQ':
        reheader = ProteinQuantCSVReheaderLFQ()
    else:
        reheader = ProteinQuantCSVReheaderTMT()
        common.validate_file_format(args.consensusXML, CONSENSUS_XML_FORMATS, check_exist=True)
        ConsensusXMLFile().load(args.consensusXML, reheader.consensus_map)
        reheader.get_channel_data()

    if args.quoting == 'double':
        reheader.quoting = '"'
    else:
        reheader.quoting = ''
    reheader.separator = args.separator

    with open(args.input_path, 'r') as handle:
        reheader.quant_handle = handle
        with open(args.output_path, 'w') as out_handle:
            reheader.rehead(out_handle, args.header)
