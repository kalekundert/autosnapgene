#!/usr/bin/env python3

from ..parser import Block
import struct

class DnaBlock(Block):
    block_id = 0
    repr_attrs = 'sequence',

    def __init__(self):
        self.topology = 'linear'
        self.strandedness = 'double'
        self.is_dam_methylated = False
        self.is_dcm_methylated = False
        self.is_ecoki_methylated = False
        self.sequence = ''

    def __repr_attr__(self, attr):
        if attr == 'sequence':
            if len(self.sequence) > 32:
                return f"{self.sequence[:16]}...{self.sequence[-16:]}"
            else:
                return self.sequence

        else:
            return super().__repr_attr__(attr)

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

class ProteinBlock(Block):
    block_id = 21
    repr_attrs = 'sequence',

    # This block is undocumented.  I attempted to reverse-engineer the file 
    # format from some example files.  Some information might be unavailable.

    def __init__(self):
        self.props = None
        self.sequence = None

    def __repr_attr__(self, attr):
        if attr == 'sequence':
            if len(self.sequence) > 32:
                return f"sequence='{self.sequence[:16]}...{self.sequence[-16:]}'"
            else:
                return f"sequence='{self.sequence}'"

        else:
            return super().__repr_attr__(attr)

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        # I don't know what the 'props' byte means.  In DNA blocks, there is an 
        # corresponding byte that encodes a bit-field specifying a few 
        # properties of the DNA, but those same properties wouldn't apply to 
        # proteins.
        block.props = bytes[0]
        block.sequence = bytes[1:].decode('ascii')

        return block

    def to_bytes(self):
        return struct.pack('>B', self.props or 0) \
                + self.sequence.encode('ascii')

