""" ProteinQuantCSVReheader """
from abc import ABC, abstractmethod
from collections import deque
import re
from typing import IO, Dict, List, Tuple
from pyopenms import ConsensusMap
from .ProteinQuantCSVMetadata import ProteinQuantCSVMetadata


OPENMS_PROTEIN_QUANT_FILE_ORIGIN_PREFIX = '# Peptide abundances computed from file'
OPENMS_PROTEIN_QUANT_SAMPLE_INFO_PREFIX = '# Files/samples associated with abundance values below'


class ProteinQuantCSVReheader(ABC):
    """ Abstract classs for protein quant CSV reheader """
    def __init__(self, quant_handle:IO, separator:str='\t', quoting:str='"'):
        """ """
        self.quant_handle = quant_handle
        self.separator = separator
        self.quoting = quoting

    @staticmethod
    @abstractmethod
    def get_csv_header_data(line:str) -> List[Tuple[str,str]]:
        """  Get CSV header data from the metadata line in the CSV file. """

    @abstractmethod
    def get_fixed_header(self, line:str, metadata:ProteinQuantCSVMetadata) -> str:
        """ Get the fixed header """

    def get_fixed_header_manual(self, line:str, headers:List[str]) -> str:
        """ Get the fixed header manually """
        err_msg = f"Number of headers provided does not match with the file: {headers}"
        new_header_fields = []
        headers_queue = deque(headers)
        for field in line.rstrip().split(self.separator):
            raw_field = field
            field = field.strip(self.quoting)
            if not field.startswith('abundance'):
                new_header_fields.append(raw_field)
                continue
            if not headers_queue:
                raise ValueError(err_msg)
            new_header = headers_queue.popleft()
            new_header_fields.append(self.quoting + new_header + self.quoting)
        if headers_queue:
            raise ValueError(err_msg)
        return  self.separator.join(new_header_fields)

    def rehead(self, output:IO, headers:List[str]=None):
        """ Rehead and write to output. """
        metadata = ProteinQuantCSVMetadata.parse(self.quant_handle)
        for line in metadata.to_lines():
            output.write(line + '\n')

        self.quant_handle.seek(0)
        line:str = self.quant_handle.readline()
        while line.startswith('#'):
            line = self.quant_handle.readline()

        if not headers:
            new_header = self.get_fixed_header(line, metadata)
        else:
            new_header = self.get_fixed_header_manual(line, headers)
        output.write(new_header + '\n')

        line = self.quant_handle.readline()
        while line:
            output.write(line)
            line = self.quant_handle.readline()

class ProteinQuantCSVReheaderTMT(ProteinQuantCSVReheader):
    """ ProteinQuantCSVReheader """
    def __init__(self, consensus_map:ConsensusMap=None, quant_handle:IO=None,
            channel_id_map:Dict[str, Dict[str, int]]=None,
            channel_descriptions:List[str]=None,
            separator:str='\t', quoting:str='"'):
        """ Constructor """
        self.consensus_map = consensus_map or ConsensusMap()
        # mapping: filename -> channel -> index
        self.channel_id_map = channel_id_map or {}
        self.channel_descriptions = channel_descriptions or []
        self.separator = separator
        self.quoting = quoting
        super().__init__(quant_handle)

    def get_channel_data(self):
        """ Get channel data from the consensusXML file """
        for i, header in self.consensus_map.getColumnHeaders().items():
            filename = header.filename
            label = header.label
            if filename not in self.channel_id_map:
                self.channel_id_map[filename] = {}
            self.channel_id_map[filename][label] = i
            self.channel_descriptions.append(header.getMetaValue('channel_description'))

    @staticmethod
    def get_csv_header_data(line:str) -> List[Tuple[str,str]]:
        """ Get CSV header data from the metadata line in the CSV file. """
        data:List[Tuple[str,str]] = []
        for field in line.rstrip().split(':', 1)[-1].strip().split(', '):
            _,v = field.split(': ')
            p = re.search(r"'(.+)' \('(.+)'\)", v)
            if not p:
                raise ValueError('Cannot parse CSV header.')
            filename = p.group(1)
            label = p.group(2)
            data.append((filename, label))
        return data

    def get_fixed_header(self, line:str, metadata:ProteinQuantCSVMetadata) -> str:
        """ Get the fixed header """
        headers = []
        for field in line.rstrip().split(self.separator):
            raw_field = field
            field = field.strip(self.quoting)
            if not field.startswith('abundance'):
                headers.append(raw_field)
                continue
            if field == 'abundance':
                # When there is only 1 sample
                index = 0
            else:
                m = re.search('^abundance_([0-9]+)$', field)
                if not m:
                    raise ValueError(f"Failed to interpret the sample index from {raw_field}")
                index = int(m.group(1)) - 1
            filename, label = metadata.header_map[index + 1]
            column_index = self.channel_id_map[filename][label]
            new_header = self.channel_descriptions[column_index]
            headers.append(self.quoting + new_header + self.quoting)
        return  self.separator.join(headers)

class ProteinQuantCSVReheaderLFQ(ProteinQuantCSVReheader):
    """ Protein quant CSV reheader for LFQ """
    def __init__(self, quant_handle:IO=None, separator:str='\t', quoting:str='"'):
        """ """
        self.separator = separator
        self.quoting = quoting
        super().__init__(quant_handle)

    @staticmethod
    def get_csv_header_data(line:str) -> List[Tuple[str,str]]:
        """ Get CSV header data from the metadata line in the CSV file. """
        data:List[Tuple[str,str]] = []
        for field in line.rstrip().split(':', 1)[-1].strip().split(', '):
            _,v = field.split(': ')
            filename = v.strip("'").strip('"')
            data.append((filename, ''))
        return data

    def get_fixed_header(self, line:str, metadata:ProteinQuantCSVMetadata) -> str:
        """ Get the fixed header """
        headers = []
        for field in line.rstrip().split(self.separator):
            raw_field = field
            field = field.strip(self.quoting)
            if not field.startswith('abundance'):
                headers.append(raw_field)
                continue
            if field == 'abundance':
                # When there is only 1 sample
                index = 0
            else:
                m = re.search('^abundance_([0-9]+)$', field)
                if not m:
                    raise ValueError(f"Failed to interpret the sample index from {raw_field}")
                index = int(m.group(1)) - 1
            if metadata.header_map:
                filename,_ = metadata.header_map[index + 1]
            else:
                filename = metadata.file_origin
            headers.append(self.quoting + filename + self.quoting)
        return  self.separator.join(headers)
