""" ProteinQuantCSV Metadata """
from __future__ import annotations
from typing import Dict, Tuple, Union, IO, List
import re

OPENMS_PROTEIN_QUANT_PREFIX = {
    'file_origin': '# Peptide abundances computed from file',
    'header_map': '# Files/samples associated with abundance values below',
    'parameters': '# Parameters (relevant only):'
}
FILE_ORIGIN_PREFIXES = (
    OPENMS_PROTEIN_QUANT_PREFIX['file_origin'],
    '# Protein abundances computed from file'
)

class ProteinQuantCSVMetadata:
    """ Metadata for ProteinQuantCSV """
    def __init__(self, file_origin:str=None,
        header_map:Dict[int,Tuple[str,str]]=None,
        parameters:Dict[str,Union[str,int]]=None,
        file_origin_prefix:str=OPENMS_PROTEIN_QUANT_PREFIX['file_origin']):
        """ constructor """
        self.file_origin = file_origin
        self.file_origin_prefix = file_origin_prefix
        self.header_map = header_map or {}
        self.parameters = parameters or {}

    def to_lines(self) -> List[str]:
        """ convert to lines """
        prefix = OPENMS_PROTEIN_QUANT_PREFIX
        parameters = ', '.join([f'{k}={v}' for k,v in self.parameters.items()])
        lines = [
            f"{self.file_origin_prefix} {self.file_origin}",
            f"{prefix['parameters']} {parameters}"
        ]
        if self.header_map:
            header_map_fields = []
            for k,v in self.header_map.items():
                field = f"{k}: '{v[0]}'"
                if v[1] is not None:
                    field += f" ('{v[1]}')"
                header_map_fields.append(field)
            header_map = ', '.join(header_map_fields)
            lines.append(f"{prefix['header_map']}: {header_map}")
        return lines

    @classmethod
    def parse(cls, handle:IO[str]):
        """ Parse metadata from ProteinQuantCSV file. """
        handle.seek(0)
        prefix = OPENMS_PROTEIN_QUANT_PREFIX
        file_origin = None
        file_origin_prefix = prefix['file_origin']
        header_map = {}
        parameters = {}
        for line in handle:
            if not line.startswith('#'):
                handle.seek(0)
                break
            if line.startswith(FILE_ORIGIN_PREFIXES):
                file_origin_prefix = next(p for p in FILE_ORIGIN_PREFIXES if line.startswith(p))
                file_origin = line.rstrip().removeprefix(file_origin_prefix)\
                    .strip().strip('"').strip("'")
            elif line.startswith(prefix['header_map']):
                for field in line.rstrip().split(':', 1)[-1].strip().split(', '):
                    i,v = field.split(': ')
                    p = re.search(r"'(.+?)'( \('(.+)'\))?", v)
                    if not p:
                        raise ValueError('Cannot parse CSV header.')
                    filename = p.group(1)
                    label = p.group(3)
                    header_map[int(i)] = (filename, label)
            elif line.startswith(prefix['parameters']):
                fields = line.rstrip().replace(prefix['parameters'], '').strip()
                if fields in ('(none)', ''):
                    continue
                for field in fields.split(','):
                    key, val = field.split('=')
                    parameters[key] = val
        return cls(file_origin, header_map, parameters, file_origin_prefix)

    def get_plex_2_column(self, offset:int=0) -> Dict[str, List[int]]:
        """ get plex to column map """
        mapper:Dict[str, List[int]] = {}
        for i, (plex, _) in self.header_map.items():
            if plex in mapper:
                mapper[plex].append(i + offset)
            else:
                mapper[plex] = [i + offset]
        return mapper

    def get_column_2_plex(self, offset:int=0) -> Dict[int, str]:
        """ get column to plex map """
        mapper:Dict[int, str] = {}
        for i, (plex, _) in self.header_map.items():
            mapper[i + offset] = plex
        return mapper
