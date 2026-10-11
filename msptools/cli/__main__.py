""" Console script main entrance """
import argparse
import sys
from msptools import cli, __version__


def main():
    """ The main entrance of this CLI software """
    parser = argparse.ArgumentParser(
        prog='msptools'
    )
    parser.add_argument(
        '-V', '--version',
        help='Show version number and exit.',
        action='store_true'
    )
    subparsers = parser.add_subparsers(
        dest='command'
    )

    cli.MzML.parse_args(subparsers)
    cli.IdXML.parse_args(subparsers)
    cli.FeatureXML.parse_args(subparsers)
    cli.ConsensusXML.parse_args(subparsers)
    cli.ProteinQuantCSV.parse_args(subparsers)
    cli.VariantPeptide.parse_args(subparsers)

    # allowing values to start with -, such as -100,3
    # https://stackoverflow.com/a/21446783/11081630
    for i, arg in enumerate(sys.argv):
        if arg[0] == '-' and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.version:
        print(f'MSPtools {__version__}', file=sys.stdout, flush=True)
        sys.exit()

    args.func(args)


if __name__ == '__main__':
    main()
