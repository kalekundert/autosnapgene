#!/usr/bin/env python3

import pytest
import autosnapgene.util as util

@pytest.mark.parametrize(
        'given, expected', [
            ('', ''),
            ('A', 'A'),
            ('AC', 'CA'),
            ('ACT', 'TCA'),
            ('ACTG', 'GTCA'),

            ('a', 'a'),
            ('ac', 'ca'),
            ('act', 'tca'),
            ('actg', 'gtca'),

            ('aC', 'Ca'),
            ('Ac', 'cA'),
        ],
)
def test_reverse(given, expected):
    assert util.reverse(given) == expected

@pytest.mark.parametrize(
        'given, expected', [
            ('', ''),
            ('A', 'T'),
            ('AC', 'TG'),
            ('ACT', 'TGA'),
            ('ACTG', 'TGAC'),
        ],
)
@pytest.mark.parametrize(
        'f', [str.upper, str.lower, str.title],
)
def test_complement(given, expected, f):
    assert util.complement(f(given)) == f(expected)

@pytest.mark.parametrize(
        'given, expected', [
            ('', ''),
            ('A', 'T'),
            ('AC', 'GT'),
            ('ACT', 'AGT'),
            ('ACTG', 'CAGT'),

            ('a', 't'),
            ('ac', 'gt'),
            ('act', 'agt'),
            ('actg', 'cagt'),

            ('aC', 'Gt'),
            ('Ac', 'gT'),
        ],
)
def test_reverse_complement(given, expected):
    assert util.reverse_complement(given) == expected

