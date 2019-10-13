#!/usr/bin/env python3

from ..parser import UndocumentedBlock
import struct

class RestrictionDigestBlock(UndocumentedBlock):
    block_id = 3
    repr_attrs = ['sites']

    # This is an undocumented block that seems to be necessary for SnapGene to 
    # recognize restriction digest sites.  Below is what I've learned about 
    # this block by trial-and-error:
    #
    # - If this block is completely missing from the file, SnapGene will behave 
    #   as if the file contains no restriction sites.  So no restriction sites 
    #   will be displayed, and the restriction-based cloning tools won't work.  
    #   This can be worked-around by making any edit to the sequence, but of 
    #   course this is inconvenient.
    #
    # - If the block is present in the file, restriction sites will work as 
    #   expected.  Interestingly, this is true even if the block is completely 
    #   empty (i.e. 0x030000).
    #
    # Below is what I know about the structure of the block:
    #
    # - The block begins with a single byte.  This byte has had a value of 0x01 
    #   in every file I've looked at.  I tried changing it to both 0x00 and 
    #   0xff; neither had any noticeable effect on SnapGene.
    #
    # - The next 4 bytes are an integer encoding the size of the ASCII string 
    #   to follow.  These bytes have had a value of 4721 in every file I've 
    #   looked at.
    #
    # - The ASCII string appears to be a comma-separated list of restriction 
    #   sites, encoded using degenerate codons.  In every file I've looked at, 
    #   there have been 469 such sites.  Note that this does not correspond to 
    #   the number of commercial sites (666) or the number of non-redundant 
    #   commercial sites (295).  Also note that neither removing this list (and 
    #   updating the preceding size) nor putting errors in the sequences has 
    #   any noticeable effect on SnapGene.
    #
    # - After the string are a few thousand bytes (2815 in the files I've 
    #   looked at) that are mostly zero, with a few ones sprinkled in there.  I 
    #   have no idea what this means.  Removing these bytes has no noticeable 
    #   effect on SnapGene.
    #
    # - At the end of the block, there is the following pattern:
    #
    #   - 0x000023c8 (int: 9160) repeated 4 times
    #   - 0x000023cc (int: 9164) repeated 3 times
    #   - 0x000023ce (int: 9166)
    #   - 0x000000000100
    #
    #   I don't know what any of that means.

    def __init__(self):
        self.unk_1 = 1
        self.sites = []
        self.unk_2 = b''

    def __repr_attr__(self, attr):
        if attr == 'sites':
            return str(len(self.sites))
        else:
            return super().__repr_attr__(attr)

    @classmethod
    def from_bytes(cls, bytes):
        block = super().from_bytes(bytes)

        _, n = struct.unpack('>BI', bytes[0:5])

        block.unk_1 = _
        block.sites = bytes[5:5+n].decode('ascii').split(',')
        block.unk_2 = bytes[5+n:]

        return block

    def to_bytes(self):
        sites_bytes = b','.join(x.encode('ascii') for x in self.sites)

        bytes = b''
        bytes += struct.pack('>BI', self.unk_1, len(sites_bytes))
        bytes += sites_bytes
        bytes += self.unk_2

        return bytes


