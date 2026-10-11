""" MAD Normalizer """
import argparse
from pyopenms import ConsensusXMLFile, ConsensusMap, DataProcessing
from msptools.cli import common
from msptools.FileTypes import FileTypes
from msptools.algorithm.MADNormalizer import MADNormalizer


PROG_NAME = 'MADNormalizer'

# pylint: disable=W0212
def add_subparsers(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name=PROG_NAME,
        help='Cross-plex MAD normalization.',
        description='Cross-plex MAD normalization similar to TMT-integrator',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--input-path',
        type=str,
        help="Path to the input file. Valid format: ['.consensusXML']"
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help="Path to the output file. Valid format: ['.consensusXML']"
    )
    parser.add_argument(
        '--min-int',
        type=float,
        help='The minimal intensity to consider. Channels with intensity smaller'
        ' than this value will be set to zero. Defaults to 1e-08',
        default=1e-08
    )
    parser.add_argument(
        '--virtual-reference',
        action='store_true',
        help='Whether to treat the experiment using virtual reference. When true'
        ' the average (mean or median) intensity of each feature will be used as'
        ' the reference intensity during normalization.'
    )
    parser.add_argument(
        '--virtual-reference-method',
        choices=['median', 'mean'],
        default='median',
        help='The method to calculate virtual reference intensity.'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ main """
    common.validate_input_file(args.input_path, [FileTypes.CONSENSUSXML])
    common.validate_output_file(args.output_path, [FileTypes.CONSENSUSXML])

    cmap = ConsensusMap()
    ConsensusXMLFile().load(args.input_path, cmap)

    normalizer = MADNormalizer(
        cmap=cmap,
        virtual_reference=args.virtual_reference,
        virtual_reference_method=args.virtual_reference_method,
        min_int=args.min_int
    )
    new_cmap = normalizer.normalize()

    process_actions = { DataProcessing.ProcessingAction.NORMALIZATION }
    common.add_data_processing(new_cmap, args, PROG_NAME, process_actions)

    ConsensusXMLFile().store(args.output_path, new_cmap)
