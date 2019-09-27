#!/usr/bin/env python3

import struct
import autoprop
import xml.etree.ElementTree as etree
from pathlib import Path
from more_itertools import one

block_ids = {}
block_classes = {}
file_types = {
        0: 'unknown',
        1: 'dna',
}

def parse(path):
    dna = SnapGene()
    dna.parse(path)
    return dna

def write(path, dna):
    dna.write(path)


def blocks_from_file(path):
    bytes = Path(path).read_bytes()
    return blocks_from_bytes(bytes)

def blocks_from_bytes(bytes):
    try:
        i = 0
        blocks = []

        while i < len(bytes):
            j = i + 5
            id, size = struct.unpack('>BI', bytes[i:j])
            block_cls = block_classes.get(id, UndocumentedBlock)
            block = block_cls.from_bytes(bytes[j:j+size])
            block.id = id
            blocks.append(block)
            i = j + size

    except SnapGeneError as e:
        e.path = path
        raise e from None

    return blocks

def file_from_blocks(path, blocks):
    bytes = bytes_from_blocks(blocks)
    Path(path).write_bytes(bytes)

def bytes_from_blocks(blocks):
    bytes = b''
    for block in blocks:
        bytes += bytes_from_block(block)
    return bytes

def bytes_from_block(block):
    bytes = block.to_bytes()
    header = struct.pack('>BI', block.id, len(bytes))
    return header + bytes

def ztr_from_data(data):
    from subprocess import run, PIPE

    # The details of the ZTR fromat are described in this publication:
    #
    #   Bonfield and Staden.  ZTR: a new format for DNA sequence trace data.  
    #   Bioinformatics 18:1:3-10 (2002).  DOI: 10.1093/bioinformatics/18.1.3
    #
    # The tools for doing the conversion can be found on GitHub:
    #
    #   https://github.com/jkbonfield/io_lib.git
    # 
    # The same tools can also be found on AUR:
    #
    #   $ yay -S staden-io_lib
    #
    # The command for doing the conversion is:
    # 
    #   $ convert_trace < path/to/ab1 > path/to/ztr
    #
    # Note that this command automatically detects the format of the input 
    # trace, and a number of formats are supported.  The are a number of 
    # command-line flags available, see -h for more information.

    p = run(['convert_trace'], input=data, capture_output=True, check=True)
    return p.stdout

class SnapGene:

    def __init__(self, path=None):
        if path:
            self.parse(path)

        self.blocks = []
        self.input_path = None

    def parse(self, path):
        self.blocks = blocks_from_file(path)
        self.input_path = Path(path)

    def write(self, path=None):
        if path is None:
            path = self.input_path
        if path is None:
            raise ValueError("not originally parsed from *.dna file, output path required.")

        file_from_blocks(path, self.blocks)


    def find_blocks(self, cls):
        return [x for x in self.blocks if isinstance(x, cls)]

    def find_block(self, cls):
        return one(self.find_blocks(cls))

    def remove_blocks(self, cls):
        self.blocks = [x for x in self.blocks if not isinstance(x, cls)]

    def remove_block(self, cls):
        block = self.find_block(cls)
        self.blocks.remove(block)


    def add_trace(self, path):
        path = Path(path)

        # Convert the sequencing data to the ZTR format.
        data = path.read_bytes()
        ztr = ztr_from_data(data)

        # Figure out the next id from the AlignmentsBlock metadata.
        align_block = self.find_block(AlignmentsBlock)
        next_id = max((x.id for x in align_block.metadata), default=0) + 1
        next_order = max((x.sort_order for x in align_block.metadata), default=0) + 1

        # Make a new AlignedTraceBlock with the ZTR data.
        trace_block = AlignedTraceBlock()
        trace_block.bytes = ztr

        # Make a new AlignedSequenceBlock with the above id and trace block.
        seq_block = AlignedSequenceBlock()
        seq_block.seq_id = next_id
        seq_block.traces = [trace_block]

        # Add the sequence block to the file (the trace block will be added 
        # indirectly via the sequence block).
        self.blocks.append(seq_block)

        # Update the AlignmentsBlock metadata.
        meta = AlignmentMetadata()
        meta.id = next_id
        meta.name = path.stem
        meta.is_trace = True
        meta.sort_order = next_order
        align_block.metadata.append(meta)

    def clear_history(self):
        self.remove_blocks(HistoryBlock)
        self.remove_blocks(HistoryNodeBlock)


class Block:

    def __init_subclass__(cls):
        if hasattr(cls, 'id'):
            block_ids[cls.name] = cls.id
            block_classes[cls.id] = cls

    def __repr__(self):
        return f"<{self.__class__.__name__} id={self.id} {self.__repr_attrs__()}".strip() + ">"

    def __repr_attrs__(self):
        return ''

    @classmethod
    def from_bytes(cls, bytes):
        raise NotImplementedError

    def to_bytes(self):
        raise NotImplementedError


class UnparsedBlock(Block):

    def __init__(self, bytes=None):
        self.bytes = bytes

    def __repr_attrs__(self):
        bytes_repr = repr(self.bytes)[2:-1]
        if len(bytes_repr) > 32:
            bytes_repr = bytes_repr[:32] + '...'
        return f"UNPARSED bytes=b'{bytes_repr}'"

    @classmethod
    def from_bytes(cls, bytes):
        return cls(bytes)

    def to_bytes(self):
        return self.bytes


@autoprop
class HeaderBlock(Block):
    id = 9
    name = 'header'

    def __init__(self):
        self.type_id = None
        self.export_version = None
        self.import_version = None

    def __repr_attrs__(self):
        return f"type='{self.type}' export_version='{self.export_version}' import_version='{self.import_version}'"


    @classmethod
    def from_bytes(cls, bytes):
        magic_cookie = bytes[0:8].decode('ascii')

        if magic_cookie != "SnapGene":
            raise SnapGeneError("not a snapgene file")

        info = struct.unpack('>HHH', bytes[8:14])
        return cls.from_info(*info)

    @classmethod
    def from_info(cls, type, export_version, import_version):
        block = cls()
        block.type_id = type
        block.export_version = export_version
        block.import_version = import_version
        return block

    def to_bytes(self):
        return b'SnapGene' + struct.pack('>HHH',
                self.type_id, self.export_version, self.import_version)


    def get_type(self):
        return file_types[self.type_id]


class DnaBlock(Block):
    id = 0
    name = 'dna'

    def __init__(self):
        self.topology = None
        self.strandedness = None
        self.is_dam_methylated = None
        self.is_dcm_methylated = None
        self.is_ecoki_methylated = None
        self.sequence = None

    def __repr_attrs__(self):
        return f"{self.topology} {self.strandedness} sequence='{self.sequence[:32]}...'"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        props = bytes[0]
        block.topology = 'circular' if props & 0x01 else 'linear'
        block.strandedness = 'double' if props & 0x02 else 'single'
        block.is_dam_methylated = bool(props & 0x04)
        block.is_dcm_methylated = bool(props & 0x08)
        block.is_ecoki_methylated = bool(props & 0x10)

        block.sequence = bytes[1:].decode('ascii')

        return block

    def to_bytes(self):
        props = sum([
                0x01 if self.topology == 'circular' else 0x00,
                0x02 if self.strandedness == 'double' else 0x00,
                0x04 * self.is_dam_methylated,
                0x08 * self.is_dcm_methylated,
                0x10 * self.is_ecoki_methylated,
        ])
        return struct.pack('>B', props) + self.sequence.encode('ascii')

class CompressedDnaBlock(UnparsedBlock):
    id = 1
    name = 'compressed_dna'

class PrimerBlock(UnparsedBlock):
    id = 5
    name = 'primers'

class NotesBlock(UnparsedBlock):
    id = 6
    name = 'notes'

class HistoryBlock(UnparsedBlock):
    id = 7
    name = 'history'

class HistoryNodeBlock(UnparsedBlock):
    id = 11
    name = 'history_node'

class PropertiesBlock(UnparsedBlock):
    id = 8
    name = 'properties'

class FeaturesBlock(UnparsedBlock):
    id = 10
    name = 'features'

class AlignmentsBlock(Block):
    id = 17
    name = 'alignments'

    def __init__(self):
        self.trim_stringency = None
        self.metadata = []

    def __repr_attrs__(self):
        return f"meta={','.join(str(x.id) for x in self.metadata)}"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        xml = bytes.decode('utf8')
        root = etree.fromstring(xml)

        block.trim_stringency = root.attrib.get('trimStringency')
        block.metadata = [
                AlignmentMetadata.from_xml(child)
                for child in root
        ]
        return block

    def to_bytes(self):
        root = etree.Element('AlignableSequences', {
            'trimStringency': self.trim_stringency,
        })
        for meta in self.metadata:
            child = meta.to_xml()
            root.append(child)

        return etree.tostring(root)

class AlignmentMetadata:

    def __init__(self):
        self.id = None
        self.name = None
        self.is_trace = None
        self.sort_order = None

    def __repr__(self):
        return f"<AlignmentMetadata id={self.id} name='{self.name}' is_trace={int(self.is_trace)} sort_order={self.sort_order}>"

    @classmethod
    def from_xml(cls, element):
        seq = cls()
        seq.id = int(element.attrib.get('ID'))
        seq.name = element.attrib.get('name')
        seq.is_trace = bool(element.attrib.get('isTrace'))
        seq.sort_order = int(element.attrib.get('sortOrder'))
        return seq

    def to_xml(self):
        return etree.Element('Sequence', {
            'ID': str(self.id),
            'name': self.name,
            'isTrace': str(int(self.is_trace)),
            'sortOrder': str(self.sort_order),
        })


class AlignedSequenceBlock(Block):
    id = 16
    name = 'aligned_sequence'

    def __init__(self):
        self.seq_id = None
        self.traces = []

    def __repr_attrs__(self):
        return f"seq_id={self.seq_id}"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()
        block.seq_id, = struct.unpack('>I', bytes[:4])
        block.traces = blocks_from_bytes(bytes[4:])
        return block

    def to_bytes(self):
        bytes = struct.pack('>I', self.seq_id)
        bytes += bytes_from_blocks(self.traces)
        return bytes

class AlignedTraceBlock(UnparsedBlock):
    id = 18
    name = 'aligned_trace'

class UracilBlock(UnparsedBlock):
    id = 19
    name = 'uracils'

class DnaColorBlock(UnparsedBlock):
    id = 20
    name = 'dna_colors'

class UndocumentedBlock(UnparsedBlock):
    pass



class SnapGeneError(Exception):

    def __init__(self, message):
        self.message = message
        self.path = None

    def __str__(self):
        return self.message.format(path=self.path)
