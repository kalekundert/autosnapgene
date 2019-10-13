#!/usr/bin/env python3

import autosnapgene as snap
from pprint import pprint

def test_blocks_from_file(examples):
    blocks = snap.parser.blocks_from_file(examples / 't7_promoter.dna')
    pprint(blocks)
    assert len(blocks) == 10

