""" Common functions for CLI """
from __future__ import annotations
import argparse
from os import PathLike
import sys
from typing import List, Union, Set, TYPE_CHECKING
from pathlib import Path
from pyopenms import DataProcessing, DateTime
from msptools.FileTypes import FileTypes


if TYPE_CHECKING:
    from pyopenms import FeatureMap, ConsensusMap, MSExperiment


def print_help_if_missing_args(parser:argparse.ArgumentParser):
    """ If no args are provided, print help and exit """
    if len(sys.argv) in [2,3] and sys.argv[1:] == parser.prog.split(' ')[1:]:
        parser.print_help(sys.stderr)
        sys.exit(1)

def validate_file_format(file:PathLike, types:List[str], check_exist:bool=False):
    """ Validate the file type """
    suffixes = Path(file).suffixes
    actual_suffixes = [suffixes[-1]]
    if len(suffixes) > 1:
        actual_suffixes.append(''.join(suffixes[-2:]))
    if not any(suffix in types for suffix in actual_suffixes):
        raise ValueError(
            f'The file {file} is invalid. Valid file types are {types}.'
        )
    if check_exist:
        if not Path(file).exists():
            raise FileNotFoundError(f'File does not exist: {file}')

def get_file_type(path:str):
    """ Get file type from file path. """
    filename = Path(path).name
    fields = filename.split('.')
    if len(fields) == 1:
        return FileTypes('unknown')
    filetype = FileTypes(fields[-1])
    if filetype is FileTypes.UNKNOWN and len(fields) > 2:
        filetype = FileTypes('.'.join(fields[-2:]).upper())
    return filetype

def validate_input_file(path, types:List[FileTypes]):
    """ Validate input file """
    filetype = get_file_type(path)
    if filetype not in types:
        valid_types = [x.value for x in types]
        raise ValueError(f"Input file {path} is invalid. Valid types are: {valid_types}")

    if not Path(path).exists():
        raise ValueError(f"Input file {path} does not exist.")

def validate_output_file(path, types:List[FileTypes]):
    """ Validate output file. """
    filetype = get_file_type(path)
    if filetype not in types:
        valid_types = [x.value for x in types]
        raise ValueError(f"Output file {path} is invalid. Valid types are: {valid_types}")

def add_data_processing(x: Union[MSExperiment, FeatureMap, ConsensusMap],
        args:argparse.Namespace, prog_name:str, process_action:Set[int]):
    """ add DataProcessing """
    # pylint: disable=C0415
    from msptools import __name__ as msptools___name__
    from msptools import __version__ as msptools___version__

    dp = DataProcessing()
    # 11 is for 'FILTERING', from 'DataProcessing.h'
    dp.setProcessingActions(process_action)
    software = dp.getSoftware()
    software.setName(f"{msptools___name__}::{prog_name}")
    software.setVersion(msptools___version__)

    dp.setSoftware(software)
    dp.setCompletionTime(DateTime().now())

    # pylint: disable=W0212
    for k,v in args._get_kwargs():
        if k == 'func':
            continue
        if isinstance(v, Path):
            v = str(v)
        dp.setMetaValue(f"parameter: {k}", v)

    data_processing = x.getDataProcessing()
    data_processing.append(dp)

    x.setDataProcessing(data_processing)
