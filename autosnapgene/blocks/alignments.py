#!/usr/bin/env python3

from ..parser import Block, Xml, Repr, UnparsedBlock
from ..parser import blocks_from_bytes, bytes_from_blocks

import struct
import xml.etree.ElementTree as etree

# AlignmentsBlock should inherit from Xml
# (then I can get rid of the etree import)

class AlignmentsBlock(Block):
    block_id = 17
    repr_attrs = ['metadata']

    def __init__(self):
        self.trim_stringency = None
        self.metadata = []

    def __repr_attr__(self, attr):
        if attr == 'metadata':
            return ','.join(str(x.id) for x in self.metadata) or 'none'

        else:
            return super().__repr_attr__(attr)

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
        block.metadata.sort(key=lambda x: x.sort_order)
        return block

    def to_bytes(self):
        root = etree.Element('AlignableSequences', {
            'trimStringency': self.trim_stringency,
        })
        for meta in self.metadata:
            child = meta.to_xml()
            root.append(child)

        return etree.tostring(root)

class AlignmentMetadata(Xml, Repr):

    class TrimmedRangeAttrib:

        @staticmethod
        def from_str(str):
            return tuple(int(x) for x in str.split('..'))

        @staticmethod
        def to_str(value):
            return '..'.join(str(x) for x in value)

    xml_tag = 'Sequence'
    xml_attrib_defs = [
            ('id', 'ID', Xml.IntAttrib),
            ('name', 'name', Xml.TextAttrib),
            ('is_visible', 'use', Xml.BoolAttrib),
            ('is_trace', 'isTrace', Xml.BoolAttrib),
            ('sort_order', 'sortOrder', Xml.IntAttrib),
            ('trimmed_range', 'trimmedRange', TrimmedRangeAttrib),
            ('is_manually_trimmed', 'manuallyTrimmed', Xml.BoolAttrib),
    ]

    def __repr__(self):
        return f"<AlignmentMetadata id={self.id} name='{self.name}' is_trace={int(self.is_trace)} sort_order={self.sort_order}>"

class AlignedSequenceBlock(Block):
    block_id = 16
    repr_attrs = ['id']

    def __init__(self):
        self.id = None
        self.traces = []

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()
        block.id, = struct.unpack('>I', bytes[:4])
        block.traces = blocks_from_bytes(bytes[4:])
        return block

    def to_bytes(self):
        bytes = struct.pack('>I', self.id)
        bytes += bytes_from_blocks(self.traces)
        return bytes

class AlignedTraceBlock(UnparsedBlock):
    block_id = 18


