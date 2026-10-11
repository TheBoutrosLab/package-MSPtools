""" This script takes the peptide intensity matrix of the variant peptide
searching and quantitation results, and creates an annotation table of the
variant combination of each peptide sequence in each sample.

--quant-csv needs 3 columns and should not have column names:

1. plex_id: Plex ID (e.g., Prot_01, Prot_02).
2. database_identifier: Database ID (e.g., Coding, Noncoding, Fusion).
3. path: Path to the ProteinQuantifier peptide CSV file.

--database-csv should be the input to the search-ProteomicCustomDB
pipeline and must have 4 columns:

1. sample_name: Sample name (e.g., ACH-000004)
2. order: Order of database to search (starting from 1).
3. database_identifier: Database ID, which is mapped to the database_identifer
   column of --quant-csv.
4. path: Path to the merged and encoded database FASTA file.

--individual-fasta-csv is the input to the merge-ProteomicDB pipeline and must
have 4 columns below. Noted that this is not required for `--label-free`.

1. sample_id: Sample name (e.g., ACH-000004)
2. plex_id: The plex ID, which is mapped to the plex_id column of --quant-csv
3. database_identifier: Database ID, which is mapped to the database_identifer
   column of --quant-csv.
4. path: Path to the individual database FASTA file.

--gvf-csv is the input to the call-NonCanonicalPeptide pipeline that has 2
columns:

1. sample_id: Sample name, which is mapped to the sample_name column of
--database-csv and sample_id of --individual-fasta-csv
2. path: Path to the GVF file.
"""
import argparse
from pathlib import Path
from msptools.cli import common
from msptools.VariantPeptide.Annotator import VariantPeptideAnnotator


INPUT_FILE_FORMATS = ['.csv']

# pylint: disable=W0511, W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='annotate',
        help='Annotate ProteinQuantifier CSV file from custom DB searching.',
        description='Annotate ProteinQuantifier CSV file from custom DB searching.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--quant-csv',
        type=Path,
        help='CSV file of the ProteinQuantifier peptide CSV file of each plex'
        '  and database identifier',
        default=None
    )
    parser.add_argument(
        '--idxml-csv',
        type=Path,
        help='CSV file to the idXML files. Not compatible with --quant-csv.',
        default=None
    )
    parser.add_argument(
        '--database-csv',
        type=Path,
        help='CSV file of the merged database FASTA file of each plex and'
        ' database identifier.'
    )
    parser.add_argument(
        '--individual-fasta-csv',
        type=Path,
        help='CSV file of the individual FASTA files of each sample.',
        required=False
    )
    parser.add_argument(
        '--gvf-csv',
        type=Path,
        help='CSV file of the GVF files for each sample.',
        required=True
    )
    parser.add_argument(
        '--index-dir',
        type=Path,
        help="Annotation GTF file.",
        required=True
    )
    parser.add_argument(
        '--novel-orf-fasta',
        type=Path,
        help='Noncoding ORF FASTA file called by moPepGen callNoncoding',
        required=True
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory.',
        required=True
    )
    parser.add_argument(
        '--label-free',
        action='store_true',
        help='Whether the experiment is label free.'
    )
    parser.add_argument(
        '--order-source',
        type=str,
        help='Order of sources, separate by comma. E.g., SNP,SNV,Fusion',
        metavar='<value>'
    )
    parser.add_argument(
        '--group-source',
        type=str,
        help='Group sources. The peptides with sources grouped will be written'
        ' to the same FASTA file. E.g., "PointMutation:gSNP,sSNV'
        ' INDEL:gINDEL,sINDEL".',
        metavar='<value>',
        nargs='*'
    )
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser


def main(args:argparse.Namespace):
    """ main entrance """
    if args.quant_csv and args.idxml_csv:
        print('Both --quant-csv and --idxml-csv received. Ignoring --idxml-csv..', flush=True)
        args.idxml_csv = None

    if args.quant_csv:
        common.validate_file_format(args.quant_csv, INPUT_FILE_FORMATS, check_exist=True)
    elif args.idxml_csv:
        common.validate_file_format(args.idxml_csv, INPUT_FILE_FORMATS, check_exist=True)
    else:
        raise ValueError('Either --quant-csv or --idxml-csv should be provided.')

    common.validate_file_format(args.database_csv, INPUT_FILE_FORMATS, check_exist=True)
    common.validate_file_format(args.individual_fasta_csv, INPUT_FILE_FORMATS, check_exist=True)
    common.validate_file_format(args.gvf_csv, INPUT_FILE_FORMATS, check_exist=True)

    annotator = VariantPeptideAnnotator()
    annotator.label_free = args.label_free

    # source group
    raw_to_group_map = None
    group_to_raw_map = {}
    if args.group_source:
        raw_to_group_map = {}
        for it in args.group_source:
            key, val = it.split(':')
            group_to_raw_map[key] = val.split(',')
            for v in group_to_raw_map[key]:
                raw_to_group_map[v] = key
    if raw_to_group_map:
        annotator.group_map = raw_to_group_map

    # source order
    if args.order_source:
        source_order = {}
        i = 0
        for val in args.order_source.split(','):
            if '-' in val:
                source = frozenset(val.split('-'))
            else:
                source = val
            if source in source_order:
                raise ValueError(
                    f"Non-unique value found from `--group-source`: {source}"
                )
            source_order[source] = i
            i += 1
            if source in group_to_raw_map:
                for child_source in group_to_raw_map[source]:
                    source_order[child_source] = i
                    i += 1

        annotator.source_order = source_order
        annotator.static_source_order = True
        annotator.append_order_internal_sources()
        annotator.set_source_set_level()

    if args.quant_csv:
        print('Loading quantitate peptide CSV file..', flush=True)
        annotator.load_peptide_path_csv(args.quant_csv)
    else:
        print('Loading detected peptide CSV file..', flush=True)
        annotator.load_peptide_path_csv(args.idxml_csv)

    print('Loading merged DB CSV file..', flush=True)
    annotator.load_merged_db_csv(args.database_csv)

    if not annotator.label_free:
        print('Loading individual DB CSV file..', flush=True)
        annotator.load_individual_db_csv(args.individual_fasta_csv)

    print('Load GVF CSV file..', flush=True)
    annotator.load_gvf_csv(args.gvf_csv)

    print('Loading transcript annotation from GTF..', flush=True)
    annotator.load_transcript_annotation(args.index_dir)

    print('Loading noncoding ORF FASTA..', flush=True)
    annotator.load_orf_annotation(args.novel_orf_fasta)

    for plex_id in annotator.peptide_paths:
        for db_id, peptide_path in annotator.peptide_paths[plex_id].items():
            print(f"Annotating plex_id={plex_id}, db_id={db_id}..", flush=True)
            annotator.clear_tally()
            peptide_path = Path(peptide_path)
            if args.quant_csv:
                annotator.load_quantified_peptides(peptide_path)
            else:
                annotator.load_detected_peptides(peptide_path)
            annotator.load_annotation_data(plex_id, db_id)
            annotator.annotate()
            annotator.tally.log()
            output_file = f"{plex_id}_{db_id}.csv"
            output_path = args.output_dir/output_file
            annotator.write(output_path)
