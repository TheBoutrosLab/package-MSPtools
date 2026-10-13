""" MAD Normalizer """
import argparse
from pyopenms import ConsensusXMLFile, ConsensusMap, DataProcessing
from msptools.cli import common
from msptools.FileTypes import FileTypes
from msptools import logger
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
    parser.add_argument(
        '--quasi-intensity',
        action='store_true',
        help='Write reference_intensity * 2**normalized_log2_ratio to each output'
        ' channel instead of the normalized log2 ratio. Measured zero sample'
        ' intensities remain zero when the reference is valid; other undefined'
        ' values produce NaN.'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def print_summary(summary:dict[str, int]) -> None:
    """ Print an aligned tally using the CLI's timestamped logger. """
    width = max(len(label) for label in summary)
    logger('MAD normalization tally (features normalized = finite sample output):')
    logger(f"{'Metric':<{width}}  {'Count':>7}")
    logger(f"{'-' * width}  {'-' * 7}")
    for label, count in summary.items():
        logger(f'{label:<{width}}  {count:>7}')


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
        min_int=args.min_int,
        quasi_intensity=args.quasi_intensity
    )
    new_cmap = normalizer.normalize()

    process_actions = { DataProcessing.ProcessingAction.NORMALIZATION }
    common.add_data_processing(new_cmap, args, PROG_NAME, process_actions)

    ConsensusXMLFile().store(args.output_path, new_cmap)
    print_summary(normalizer.summary)
