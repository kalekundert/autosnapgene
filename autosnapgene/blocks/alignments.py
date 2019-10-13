#!/usr/bin/env python3

from ..parser import Block, Xml, Repr, UnparsedBlock
from ..parser import blocks_from_bytes, bytes_from_blocks

import struct

class AlignmentsBlock(Xml, Block):
    block_id = 17
    repr_attrs = ['metadata']

    class MetadataTag(Xml.AppendListTag):

        @staticmethod
        def from_xml(element):
            return AlignmentMetadata.from_xml(element)

        @staticmethod
        def to_xml(parent, tag, metadata):
            for meta in metadata:
                e = meta.to_xml()
                parent.append(e)

    xml_tag = 'AlignableSequences'
    xml_subtag_defs = [
            ('metadata', 'Sequence', MetadataTag, []),
    ]
    xml_attrib_defs = [
            ('trim_stringency', 'trimStringency', Xml.TextAttrib),
    ]

    def __repr_attr__(self, attr):
        if attr == 'metadata':
            return ','.join(str(x.id) for x in self.metadata) or 'none'

        else:
            return super().__repr_attr__(attr)

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


