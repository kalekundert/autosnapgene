#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path
from inspect import getmembers, signature
from pprint import pprint

def no_required_args(f):
    # Don't count the self parameter.
    non_self_params = list(signature(f).parameters.values())[1:]

    for p in non_self_params:
        if p.default == p.empty:
            return False

    return True

snap_getters = [
        k
        for k, v in getmembers(snap.SnapGene)
        if k.startswith('get_') and no_required_args(v)
]

@pytest.mark.parametrize('getter', snap_getters)
def test_getters_same_after_write(getter, examples, parse_and_write):
    # Make sure no getters change value as a result of reading and writing the 
    # same file.
    in_paths = examples.glob('*.dna')

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


