""" MzMLFilter """
import argparse
from pyopenms import MzMLFile, ConsensusMap, ConsensusXMLFile, IdXMLFile
from msptools import warning
from msptools.cli import common
from msptools.MzML.Filter import MzMLFilter


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='filter',
        help='Filter spectra from mzML file.',
        description='Filter spectra from mzML file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input-path',
        type=str,
        help="Path to the input file. Valid format: ['.mzML']"
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help="Path to the output file. Valid format: ['.mzML']"
    )
    parser.add_argument(
        '--consensus-file',
        type=str,
        help="Path to the consensusXML file that the spectra with peptide "
        "identification data should be removed.",
        default=None
    )
    parser.add_argument(
        '--id-file',
        type=str,
        help="Path to the idXML file that the matched MS2 spectra will"
        " be removed.",
        default=None
    )
    parser.add_argument(
        '--rt-tolerance',
        type=float,
        help="RT tolerance in second.",
        default=0.1
    )
    parser.add_argument(
        '--mz-tolerance',
        type=float,
        help="MZ tolerance.",
        default=20
    )
    parser.add_argument(
        '--mz-tolerance-unit',
        type=str,
        help="MZ tolerance unit",
        choices=['ppm', 'Da'],
        default='ppm'
    )
    parser.add_argument(
        '--consensus-check-id',
        type=str,
        help='When true, only consensus features with peptide'
        ' identification are removed.',
        choices=['true', 'false'],
        default='true'
    )
    parser.add_argument(
        '--use-rt-mz',
        action='store_true',
        help='If the spectrum ID does not match with any feature from the '
        'IdXML or ConsensusXML, try to match by RT and MZ. Defaults to false.'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ main """
    if args.id_file and args.consensus_file:
        warning(
            '--id-file and --consensus-file are both specified but only --id-file'
            ' is used for filtering.'
        )

    if not args.id_file and not args.consensus_file:
        raise ValueError('Either of --id-file or --consensus-file must be given')

    mzml_filter = MzMLFilter()
    mzml_filter.input_path = args.input_path
    mzml_filter.output_path = args.output_path
    mzml_filter.id_file = args.id_file
    mzml_filter.consensus_file = args.consensus_file
    mzml_filter.rt_tolerance = args.rt_tolerance
    mzml_filter.mz_tolerance = args.mz_tolerance
    mzml_filter.mz_tolerance_unit = args.mz_tolerance_unit
    mzml_filter.consensus_check_id = args.consensus_check_id == 'true'
    mzml_filter.use_rt_mz = args.use_rt_mz

    MzMLFile().load(mzml_filter.input_path, mzml_filter.exp)

    mzml_filter.count_total_spectra()
    mzml_filter.exp.sortSpectra()
    if mzml_filter.use_rt_mz:
        mzml_filter.make_spectrum_dict()

    if mzml_filter.id_file:
        prot_ids = []
        pep_ids = []
        IdXMLFile().load(mzml_filter.id_file, prot_ids, pep_ids)
        mzml_filter.get_feature_pos_id(pep_ids)
    else:
        cmap = ConsensusMap()
        ConsensusXMLFile().load(mzml_filter.consensus_file, cmap)
        mzml_filter.get_feature_pos_consensus(cmap)

    mzml_filter.tally.n_psm = len(mzml_filter.feature_pos)
    exp2 = mzml_filter.filter_by_denylist()
    mzml_filter.tally.n_stored = exp2.size()
    mzml_filter.tally.log()
    MzMLFile().store(mzml_filter.output_path, exp2)
