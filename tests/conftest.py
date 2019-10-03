#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path

EX = Path(__file__).parent / 'examples'

@pytest.fixture
def parse_and_write(tmp_path):

    def helper(name):
        in_path = EX / name
        out_path = tmp_path / in_path.name

        in_dna = snap.parse(EX / 't7_promoter.dna')
        snap.write(out_path, in_dna)
        out_dna = snap.parse(out_path)

        return in_dna, out_dna

    return helper



