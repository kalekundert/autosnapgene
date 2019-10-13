#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path

@pytest.fixture
def examples():
    return Path(__file__).parent / 'examples'

@pytest.fixture
def parse_and_write(examples, tmp_path):

    def helper(name):
        in_path = examples / name
        out_path = tmp_path / in_path.name

        in_dna = snap.parse(examples / name)
        snap.write(out_path, in_dna)
        out_dna = snap.parse(out_path)

        print('input...');  yield in_dna
        print('output...'); yield out_dna

    return helper



