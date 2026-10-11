""" Remove NA features """
import argparse
import numpy as np
from pyopenms import (
    ConsensusXMLFile, ConsensusMap, DataProcessing, ConsensusFeature,
    BaseFeature
)
from msptools.FileTypes import FileTypes
from msptools.cli import common
from msptools import logger


PROG_NAME = 'removeNA'

# pylint: disable=W0212
def add_subparsers(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name=PROG_NAME,
        help='Remove NA features',
        description='Remove NA features',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-path',
        type=str,
        help=f"Path to the input ConsensusXML file. Valid formats: {FileTypes.CONSENSUSXML}",
        required=True
    )
    parser.add_argument(
        '-o', '--output-path',
        type=str,
        help=f"Path to the output idXML file. Valid formats: {FileTypes.CONSENSUSXML}",
        required=True
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def main(args:argparse.Namespace):
    """ main """
    common.validate_input_file(args.input_path, [FileTypes.CONSENSUSXML])
    common.validate_output_file(args.output_path, [FileTypes.CONSENSUSXML])

    cmap1 = ConsensusMap()
    ConsensusXMLFile().load(args.input_path, cmap1)
    cmap2 = ConsensusMap(cmap1)
    cmap2.clear(False)

    summary = {i:0 for i in cmap1.getColumnHeaders()}

    for cf in cmap1:
        new_cf = ConsensusFeature(cf)
        new_cf.clear()

        for feature in cf.getFeatureList():
            intensity = feature.getIntensity()
            if np.isnan(intensity):
                summary[feature.getMapIndex()] += 1
                continue

            new_feature = BaseFeature()
            new_feature.setWidth(feature.getWidth())
            new_feature.setCharge(feature.getCharge())
            new_feature.setMZ(feature.getMZ())
            new_feature.setRT(feature.getRT())
            new_feature.setIntensity(intensity)
            map_index = feature.getMapIndex()
            new_cf.insert(map_index, new_feature)

        if len(new_cf.getFeatureList()) == 0:
            continue

        keys = []
        cf.getKeys(keys)
        for key in keys:
            val = new_cf.getMetaValue(key)
            new_cf.setMetaValue(key, val)
        cmap2.push_back(new_cf)
    cmap2.ensureUniqueId()

    process_actions = { DataProcessing.ProcessingAction.FILTERING }
    common.add_data_processing(cmap2, args, PROG_NAME, process_actions)

    ConsensusXMLFile().store(args.output_path, cmap2)

    for i,v in summary.items():
        column = cmap1.getColumnHeaders()[i]
        k = f"{column.filename} ({column.label})"
        logger(f"{k}: {v} features removed.")
