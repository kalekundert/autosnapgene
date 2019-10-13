#!/usr/bin/env python3

from ..parser import Block
from ..errors import *

import autoprop
import struct

@autoprop
class HeaderBlock(Block):
    block_id = 9
    repr_attrs = 'type', 'export_version', 'import_version'

    file_types = {
            0: 'unknown',
            1: 'dna',
            2: 'protein',  # Just a guess...
    }

    def __init__(self):
        self.type_id = 1
        self.export_version = 14
        self.import_version = 14

    @classmethod
    def from_bytes(cls, bytes):
        magic_cookie = bytes[0:8].decode('ascii')

        if magic_cookie != "SnapGene":
            raise ParseError("not a snapgene file")

        info = struct.unpack('>HHH', bytes[8:14])
        return cls.from_info(*info)

    @classmethod
    def from_info(cls, type_id, export_version, import_version):
        block = cls()
        block.type_id = type_id
        block.export_version = export_version
        block.import_version = import_version
        return block

    def to_bytes(self):
        return b'SnapGene' + struct.pack('>HHH',
                self.type_id, self.export_version, self.import_version)

    def get_type(self):
        return self.file_types[self.type_id]


