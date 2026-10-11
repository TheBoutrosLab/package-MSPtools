""" This module takes in the idXML file converted from a protXML output by
ProteinProphet and filters protein groups (indistinguishable proteins) with
a given cutoff of probability. The motivation is that ProteinProphet usually
assign a leader protein to a group of indistinguishable proteins. But it also
reports "unneeded" proteins with probability of zero. Including those protein
groups in the idXML file confuses ProteinQuantifier for not being able to use
shared peptides to quantitate the correct protein (group).
"""
from __future__ import annotations
import argparse
from pyopenms import (
    IdXMLFile, ProteinIdentification, ProteinGroup, ConsensusXMLFile, ConsensusMap,
    ProteinHit
)
from msptools.FileTypes import FileTypes
from msptools.cli import common


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction) -> argparse.Namespace:
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='filterProteinGroups',
        help='Filter protein groups.',
        description='Filter protein groups.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-path',
        type=str,
        help=f"Path to the input idXML file. Valid formats: {FileTypes.IDXML}",
        required=True
    )
    parser.add_argument(
        '--accept-list',
        type=str,
        default=None,
        help=f'Path to the consensusXML file that proteins are used as the accept'
        f' list. Valid formats: {FileTypes.CONSENSUSXML}'
    )
    parser.add_argument(
        '-o', '--output-path',
        type=str,
        help=f"Path to the output idXML file. Valid formats: {FileTypes.IDXML}",
        required=True
    )
    parser.add_argument(
        '--probability-cutoff',
        type=float,
        help='Probability cutoff for filtering. Protein groups with probability'
        'smaller or equal to this value will be removed.',
        required=True,
        default=0.0
    )
    parser.add_argument(
        '--grouping-algorithm',
        type=str,
        choices=['ProteinProphet', 'Fido'],
        help='Protein grouping algorithm used.',
        required=True
    )
    parser.add_argument(
        '--remove-decoy',
        action='store_true',
        help='Remove decoy proteins from protein groups.'
    )
    parser.add_argument(
        '--decoy-string',
        type=str,
        help='Decoy string (e.g. "DECOY_", "rev_")',
        default='DECOY_'
    )
    parser.add_argument(
        '--decoy-string-position',
        type=str,
        choices={'prefix', 'suffix'},
        help='Decoy string position.',
        default='prefix'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def create_protein_identification(proteins:ProteinIdentification) -> ProteinIdentification:
    """ Creates a new protein identification based on a given one. """
    new_proteins = ProteinIdentification()
    new_proteins.setHits(proteins.getHits())
    new_proteins.setScoreType(proteins.getScoreType())
    new_proteins.setSearchEngine(proteins.getSearchEngine())
    new_proteins.setSearchEngineVersion(proteins.getSearchEngineVersion())
    new_proteins.setSearchParameters(proteins.getSearchParameters())
    new_proteins.setSignificanceThreshold(proteins.getSignificanceThreshold())
    new_proteins.setIdentifier(proteins.getIdentifier())

    keys = []
    proteins.getKeys(keys)
    for key in keys:
        new_proteins.setMetaValue(key, proteins.getMetaValue(key))

    return new_proteins

def is_decoy_protein(accession:str, decoy_string:str, decoy_string_pos:str):
    """ Whether a given protein is decoy """
    supported_decoy_pos = {'prefix', 'suffix'}
    if decoy_string_pos not in supported_decoy_pos:
        raise ValueError(
            f"Unrecognized decoy string position of {decoy_string_pos}."
            f"Supported values: {supported_decoy_pos}"
        )

    if decoy_string_pos == 'prefix':
        return accession.startswith(decoy_string)
    return accession.endswith(decoy_string)

def remove_devoy_protein(protein_group:ProteinGroup, decoy_string:str,
        decoy_string_pos:str) -> ProteinGroup:
    """ Remove decoy protein from a protein group. """
    accessions = []
    for acc in protein_group.accessions:
        if is_decoy_protein(acc.decode('utf-8'), decoy_string, decoy_string_pos):
            continue
        accessions.append(acc)
    pg = ProteinGroup(protein_group)
    pg.accessions = accessions
    return pg

def load_accept_list(path:str) -> set[str]:
    """ Load protein accept list """
    if path is None:
        return None
    accept_list = set()
    cmap = ConsensusMap()
    ConsensusXMLFile().load(path, cmap)
    prot_id:ProteinIdentification
    protein:ProteinHit
    for prot_id in cmap.getProteinIdentifications():
        for protein in prot_id.getHits():
            acc = protein.getAccession()
            accept_list.add(acc)
    return accept_list

def remove_unaccepted_protein(protein_group:ProteinGroup, accept_list:set[str]):
    """ remove proteins abscent from the accept list """
    accessions = []
    for acc in protein_group.accessions:
        if acc.decode('utf-8') not in accept_list:
            continue
        accessions.append(acc)
    pg = ProteinGroup(protein_group)
    pg.accessions = accessions
    return pg

def main(args:argparse.Namespace):
    """ main """
    common.validate_input_file(args.input_path, [FileTypes.IDXML])
    common.validate_output_file(args.output_path, [FileTypes.IDXML])

    prot_idents:list[ProteinIdentification] = []
    pep_idents = []
    IdXMLFile().load(args.input_path, prot_idents, pep_idents)

    new_prot_idents = create_protein_identification(prot_idents[0])

    accept_list = load_accept_list(args.accept_list)

    if args.grouping_algorithm == 'ProteinProphet':
        protein_groups = prot_idents[0].getProteinGroups()
    else:
        protein_groups = prot_idents[0].getIndistinguishableProteins()

    for protein_group in protein_groups:
        if protein_group.probability <= args.probability_cutoff:
            continue
        if args.remove_decoy:
            protein_group = remove_devoy_protein(
                protein_group, args.decoy_string, args.decoy_string_position
            )
        if accept_list is not None:
            protein_group = remove_unaccepted_protein(
                protein_group, accept_list
            )
        if protein_group.accessions:
            new_prot_idents.insertProteinGroup(protein_group)
            new_prot_idents.insertIndistinguishableProteins(protein_group)

    IdXMLFile().store(args.output_path, [new_prot_idents], pep_idents)
