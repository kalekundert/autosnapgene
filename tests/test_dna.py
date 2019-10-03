#!/usr/bin/env python3

import pytest
import autosnapgene as snap

bytes_params = [
            (b'\x00',     (    '',   'linear', 'single', 0, 0, 0)),
            (b'\x01',     (    '', 'circular', 'single', 0, 0, 0)),
            (b'\x02',     (    '',   'linear', 'double', 0, 0, 0)),
            (b'\x03',     (    '', 'circular', 'double', 0, 0, 0)),

            (b'\x04',     (    '',   'linear', 'single', 1, 0, 0)),
            (b'\x05',     (    '', 'circular', 'single', 1, 0, 0)),
            (b'\x06',     (    '',   'linear', 'double', 1, 0, 0)),
            (b'\x07',     (    '', 'circular', 'double', 1, 0, 0)),

            (b'\x08',     (    '',   'linear', 'single', 0, 1, 0)),
            (b'\x09',     (    '', 'circular', 'single', 0, 1, 0)),
            (b'\x0A',     (    '',   'linear', 'double', 0, 1, 0)),
            (b'\x0B',     (    '', 'circular', 'double', 0, 1, 0)),

            (b'\x10',     (    '',   'linear', 'single', 0, 0, 1)),
            (b'\x11',     (    '', 'circular', 'single', 0, 0, 1)),
            (b'\x12',     (    '',   'linear', 'double', 0, 0, 1)),
            (b'\x13',     (    '', 'circular', 'double', 0, 0, 1)),

            (b'\x00A',    (   'A',   'linear', 'single', 0, 0, 0)),
            (b'\x00AT',   (  'AT',   'linear', 'single', 0, 0, 0)),
            (b'\x00ATC',  ( 'ATC',   'linear', 'single', 0, 0, 0)),
            (b'\x00ATCG', ('ATCG',   'linear', 'single', 0, 0, 0)),
            (b'\x00N',    (   'N',   'linear', 'single', 0, 0, 0)),
]

def test_examples(parse_and_write):
    for dna in parse_and_write('t7_promoter.dna'):
        assert dna.sequence == 'TAATACGACTCACTATAGG'
        assert dna.topology == 'linear'
        assert dna.is_dam_methylated == False
        assert dna.is_dcm_methylated == False
        assert dna.is_ecoki_methylated == False

@pytest.mark.parametrize('bytes,params', bytes_params)
def test_from_bytes(bytes, params):
    dna = snap.DnaBlock.from_bytes(bytes)

    assert dna.sequence == params[0]
    assert dna.topology == params[1]
    assert dna.strandedness == params[2]
    assert dna.is_dam_methylated == params[3]
    assert dna.is_dcm_methylated == params[4]
    assert dna.is_ecoki_methylated == params[5]

@pytest.mark.parametrize('bytes,params', bytes_params)
def test_to_bytes(bytes, params):
    dna = snap.DnaBlock()

    dna.sequence = params[0]
    dna.topology = params[1]
    dna.strandedness = params[2]
    dna.is_dam_methylated = params[3]
    dna.is_dcm_methylated = params[4]
    dna.is_ecoki_methylated = params[5]

    assert dna.to_bytes() == bytes
