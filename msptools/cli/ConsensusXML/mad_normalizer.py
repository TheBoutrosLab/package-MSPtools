""" MAD Normalizer """
import argparse
from pathlib import Path
from pyopenms import ConsensusXMLFile, ConsensusMap, DataProcessing
from msptools.cli import common
from msptools.FileTypes import FileTypes
from msptools import logger
from msptools.algorithm.MADNormalizer import MADNormalizer


PROG_NAME = 'MADNormalizer'


def positive_int(value:str) -> int:
    """ Parse a positive integer for a command-line option. """
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('expected a positive integer') from exc
    if number < 1:
        raise argparse.ArgumentTypeError('expected a positive integer')
    return number

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
        '--virtual-reference-min-channels', type=positive_int, default=1,
        help='Minimum number of finite reporter channels above --min-int needed'
        ' to calculate a virtual reference for a feature and plex.'
    )
    parser.add_argument(
        '--quasi-intensity',
        action='store_true',
        help='Write reference_intensity * 2**normalized_log2_ratio to each output'
        ' channel instead of the normalized log2 ratio. Measured zero sample'
        ' intensities remain zero when the reference is valid; other undefined'
        ' values produce NaN.'
    )
    parser.add_argument(
        '--layout', choices=['scan', 'linked'], default='scan',
        help='Input map layout: one plex per feature (scan) or allocated channel'
        ' areas from several plexes per feature (linked).'
    )
    parser.add_argument(
        '--bridge-scale', action='store_true',
        help='Scale linked quasi-intensities so each plex has the same median'
        ' across all its positive, finite reference areas.'
    )
    parser.add_argument(
        '--allow-virtual-bridge', action='store_true',
        help='Allow bridge scaling with virtual references. This assumes plexes'
        ' have comparable overall sample composition.'
    )
    parser.add_argument(
        '--scaling-report-path', type=str,
        help='Optional path for a per-plex tab-separated scaling report.'
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
    if Path(args.input_path).resolve() == Path(args.output_path).resolve():
        raise ValueError('MADNormalizer input and output paths must differ')
    if args.scaling_report_path:
        report_path = Path(args.scaling_report_path).resolve()
        if report_path in (Path(args.input_path).resolve(),
                           Path(args.output_path).resolve()):
            raise ValueError('Scaling report path must differ from input and output paths')

    cmap = ConsensusMap()
    ConsensusXMLFile().load(args.input_path, cmap)

    normalizer = MADNormalizer(
        cmap=cmap,
        virtual_reference=args.virtual_reference,
        virtual_reference_method=args.virtual_reference_method,
        virtual_reference_min_channels=args.virtual_reference_min_channels,
        min_int=args.min_int,
        quasi_intensity=args.quasi_intensity,
        layout=args.layout,
        bridge_scale=args.bridge_scale,
        allow_virtual_bridge=args.allow_virtual_bridge,
        scaling_report_path=args.scaling_report_path
    )
    new_cmap = normalizer.normalize()

    process_actions = { DataProcessing.ProcessingAction.NORMALIZATION }
    processing_args = argparse.Namespace(**{key: value for key, value in vars(args).items()
                                            if value is not None})
    common.add_data_processing(new_cmap, processing_args, PROG_NAME, process_actions)

    ConsensusXMLFile().store(args.output_path, new_cmap)
    print_summary(normalizer.summary)
    for row in normalizer.scaling_report:
        logger(f"Plex {row['plex_id']}: valid references={row['reference_feature_count']}, "
               f"median reference area={row['median_reference_area']}, "
               f"target={row['target_reference_area']}, "
               f"scale={row['plex_scale_factor']}")
    for plex, map_index in normalizer.fallback_channels:
        logger(f'Plex {plex}, map index {map_index}: zero or undefined MAD; '
               'used median centering without variance scaling.')
