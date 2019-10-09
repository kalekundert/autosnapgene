#!/usr/bin/env python3

import pytest
import autosnapgene as snap

def test_getters(parse_and_write):
    for dna in parse_and_write('flag_tag.prot'):
        assert dna.sequence == 'DYKDDDDK'
        assert dna.protein_sequence == 'DYKDDDDK'

def test_setters():
    dna = snap.SnapGene()
    dna.protein_sequence = 'DYKDDDDK'
