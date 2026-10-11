""" Fix zero plexes """
from __future__ import annotations
from typing import Dict, List, Tuple
import argparse
from contextlib import ExitStack
import sys
from pyopenms import ConsensusMap, ConsensusXMLFile, IdXMLFile
from msptools import logger
from msptools.cli import common
from msptools.FileTypes import FileTypes
from msptools.ProteinQuantCSV.ProteinQuantCSV import ProteinQuantCSV


# pylint: disable=W0212
def add_subparsers(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='fixZeroPlex',
        help='Set the intensity of all zero plexes to NA.',
        description='Set the intensity of all zero plexes to NA.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-path',
        type=str,
        help=f"Path to the protein quant CSV file. Valid formats: {FileTypes.CSV}",
        metavar='<file>',
        required=True
    )
    parser.add_argument(
        '--consensusXML',
        type=str,
        help=f"Path to the consensusXML file. Valid formats: {FileTypes.CONSENSUSXML}",
        metavar='<file>',
        default=None
    )
    parser.add_argument(
        '--protein-groups',
        type=str,
        help='Protein inference results for the identification runs that were '
        'used to annotate the input (e.g. via the ProteinInference tool). '
        'Information about indistinguishable proteins will be used for protein '
        f"quantification. Valid formats: {FileTypes.IDXML}",
        metavar='<file>',
        default=None
    )
    parser.add_argument(
        '--type',
        type=str,
        help="Quant CSV table for protein (default) or peptide.",
        default='protein',
        choices=['protein', 'peptide']
    )
    parser.add_argument(
        '-o', '--output-path',
        type=str,
        help=f"Path to the output CSV file. Valid formats: {FileTypes.CSV}",
        metavar='<file>',
        required=True
    )
    parser.add_argument(
        '--greedy-group-resolution',
        action='store_true',
        help='Pre-process identifications with greedy resolution of shared peptides'
        ' based on the protein group probabilities. (Only works with an idXML file'
        ' given as --protein-groups parameter).'
    )
    parser.add_argument(
        '--nocheck-intensity',
        action='store_true',
        help='Skip checking protein and peptide intensities. When set, peptides'
        ' and proteins not used for quantitation will be set to NA regardless of'
        ' their actual intensity from the input CSV file.'
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
    """ main """
    common.validate_input_file(args.input_path, [FileTypes.CSV])
    common.validate_output_file(args.output_path, [FileTypes.CSV])

    supported_formats = ['protein', 'peptide']
    if args.type not in supported_formats:
        logger(f"Unsupported format: '{format}'. Supported formats are: {supported_formats}")
        sys.exit(1)

    with ExitStack() as stack:
        in_handle = stack.enter_context(open(args.input_path, 'rt'))
        out_handle = stack.enter_context(open(args.output_path, 'wt'))

        quant_csv = ProteinQuantCSV.read(in_handle, args.type)
        summary:Dict[str, int] = {}
        logger("Setting all-zero plexes to nan..")
        quant_csv.remove_zero_plex(summary=summary)

        for plex, n in summary.items():
            logger(f"'{plex}': {n} features set to nan.")

        if args.consensusXML:
            logger("Setting unquantified features to nan..")
            cmap = ConsensusMap()
            ConsensusXMLFile().load(args.consensusXML, cmap)

            summary:List[Tuple[str,int]] = []

            if args.type == 'peptide':
                quant_csv.remove_unquantified_peptides(
                    cmap, summary=summary, check_intensity=(not args.nocheck_intensity)
                )
            else:
                protein_groups = None
                if args.protein_groups:
                    prot_ids, pep_ids = [], []
                    IdXMLFile().load(args.protein_groups, prot_ids, pep_ids)
                    if not prot_ids:
                        logger("Protein grouping information was not found.")
                    else:
                        protein_groups = prot_ids[0]
                quant_csv.removed_unquantified_proteins(
                    cmap=cmap,
                    protein_groups=protein_groups,
                    greedy_group_resolution=args.greedy_group_resolution,
                    summary=summary, check_intensity=(not args.nocheck_intensity)
                )

            for i, (sample, count) in enumerate(summary):
                logger(f"Column {i}, sample '{sample}': {count} features set to nan.")

        fieldnames = ProteinQuantCSV.get_fieldnames(in_handle, args.separator)
        quoting = ProteinQuantCSV.get_quoting(args.quoting)
        quant_csv.write(
            out_handle, fieldnames=fieldnames,
            delimiter=args.separator, quoting=quoting
        )
