#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path
from inspect import getmembers

EX = Path(__file__).parent / 'examples'

snap_getters = [
        k
        for k, v in getmembers(snap.SnapGene)
        if k.startswith('get_') and k not in {
            'get_references',
        }
]

@pytest.mark.parametrize('getter', snap_getters)
def test_getters_same_after_write(getter, parse_and_write):
    # Make sure no getters change value as a result of reading and writing the 
    # same file.
    in_paths = EX.glob('*.dna')

    for in_path in in_paths:
        in_dna, out_dna = parse_and_write(in_path.name)
        in_getter = getattr(in_dna, getter)
        out_getter = getattr(out_dna, getter)

        # Some getters will raise AttributeError if the file in question 
        # doesn't have the requested information.

        try:
            in_value = in_getter()
        except AttributeError:
            in_value = AttributeError

        try:
            out_value = out_getter()
        except AttributeError:
            out_value = AttributeError

        assert in_value == out_value


def test_blocks_from_file():
    blocks = snap.blocks_from_file(EX / 't7_promoter.dna')
    assert len(blocks) == 10

