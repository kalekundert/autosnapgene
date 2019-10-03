#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path

EX = Path(__file__).parent / 'examples'

@pytest.mark.skip
def test_parse_write(tmp_path):
    import inspect

    in_paths = EX.glob('*.dna')
    getters = [
            k
            for k, v in inspect.getmembers(snap.SnapGene)
            if k.startswith('get_')
    ]

    for in_path in in_paths:
        in_dna = snap.parse(in_path)
        out_path = tmp_path / in_path.name
        snap.write(out_path, in_dna)
        out_dna = snap.parse(out_path)

        # Make sure no getters change value as a result of reading and writing 
        # the same file.

        for getter in getters:
            in_getter = getattr(in_dna, getter)
            out_getter = getattr(out_dna, getter)
            assert in_getter() == out_getter()

    for dna in [in_dna, out_dna]:
        assert dna.sequence == 'TAATACGACTCACTATAGG'
        assert dna.topology == 'linear'
        assert dna.is_dam_methylated == False
        assert dna.is_dcm_methylated == False
        assert dna.is_ecoki_methylated == False


def test_blocks_from_file():
    blocks = snap.blocks_from_file(EX / 't7_promoter.dna')
    assert len(blocks) == 10

    pprint(blocks)

