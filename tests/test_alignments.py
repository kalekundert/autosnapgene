#!/usr/bin/env python3

import pytest
import autosnapgene as snap
from pathlib import Path

def test_getters(parse_and_write):
    for dna in parse_and_write('puc19_bsai_abc.dna'):
        assert dna.count_traces() == count_seq_blocks(dna) == 3
        assert dna.trace_names == [
                'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_c']

@pytest.mark.parametrize(
        'path, count', [
            ('puc19_bsai.dna',     0),
            ('puc19_bsai_a.dna',   1),
            ('puc19_bsai_ab.dna',  2),
            ('puc19_bsai_abc.dna', 3),
])
def test_count_traces(examples, path, count):
    dna = snap.parse(examples / path)
    assert dna.count_traces() == count

@pytest.mark.parametrize(
        'path, name, has_trace', [
            ('puc19_bsai.dna',     'puc19_bsai_a', False),
            ('puc19_bsai.dna',     'puc19_bsai_b', False),
            ('puc19_bsai.dna',     'puc19_bsai_c', False),

            ('puc19_bsai_a.dna',   'puc19_bsai_a', True),
            ('puc19_bsai_a.dna',   'puc19_bsai_b', False),
            ('puc19_bsai_a.dna',   'puc19_bsai_c', False),

            ('puc19_bsai_ab.dna',  'puc19_bsai_a', True),
            ('puc19_bsai_ab.dna',  'puc19_bsai_b', True),
            ('puc19_bsai_ab.dna',  'puc19_bsai_c', False),

            ('puc19_bsai_abc.dna', 'puc19_bsai_a', True),
            ('puc19_bsai_abc.dna', 'puc19_bsai_b', True),
            ('puc19_bsai_abc.dna', 'puc19_bsai_c', True),
])
def test_have_trace(examples, path, name, has_trace):
    dna = snap.parse(examples / path)
    assert dna.has_trace(name) == has_trace

@pytest.mark.parametrize(
        'path, names', [
            ('puc19_bsai.dna', [
            ]),

            ('puc19_bsai_a.dna', [
                'puc19_bsai_a',
            ]),

            ('puc19_bsai_ab.dna', [
                'puc19_bsai_a',
                'puc19_bsai_b',
            ]),

            ('puc19_bsai_abc.dna', [
                'puc19_bsai_a',
                'puc19_bsai_b',
                'puc19_bsai_c',
            ]),
])
def test_trace_names(examples, path, names):
    dna = snap.parse(examples / path)
    assert dna.trace_names == names

def test_add_trace(examples):
    dna = snap.parse(examples / 'puc19_bsai.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

    dna.add_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == ['puc19_bsai_a']

    dna.add_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == ['puc19_bsai_a']

    dna.add_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == ['puc19_bsai_a', 'puc19_bsai_b']

    dna.add_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == ['puc19_bsai_a', 'puc19_bsai_b']

def test_append_trace(examples):
    dna = snap.parse(examples / 'puc19_bsai.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

    dna.append_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == [
            'puc19_bsai_a']

    dna.append_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_a']

    dna.append_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_a', 'puc19_bsai_b']

    dna.append_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 4
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_b']

def test_prepend_trace(examples):
    dna = snap.parse(examples / 'puc19_bsai.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

    dna.prepend_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == [
            'puc19_bsai_a']

    dna.prepend_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_a']

    dna.prepend_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_b', 'puc19_bsai_a', 'puc19_bsai_a']

    dna.prepend_trace(examples / 'puc19_bsai_b.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 4
    assert dna.trace_names == [
            'puc19_bsai_b', 'puc19_bsai_b', 'puc19_bsai_a', 'puc19_bsai_a']

def test_insert_trace(examples):
    dna = snap.parse(examples / 'puc19_bsai_ab.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b']

    dna.insert_trace(1, examples / 'puc19_bsai_c.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_c', 'puc19_bsai_b']

def test_remove_trace(examples):
    dna = snap.parse(examples / 'puc19_bsai_abc.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_c']

    dna.remove_trace('puc19_bsai_b')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_c']

    dna.remove_trace('puc19_bsai_c')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == [
            'puc19_bsai_a']

    dna.remove_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

    with pytest.raises(ValueError):
        dna.remove_trace('xxx')

def test_remove_trace_dups(examples):
    dna = snap.parse(examples / 'puc19_bsai_ab.dna')
    dna.append_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_a']

    dna.remove_trace('puc19_bsai_a')
    assert dna.count_traces() == count_seq_blocks(dna) == 1
    assert dna.trace_names == [
            'puc19_bsai_b']

def test_replace_target(examples):
    dna = snap.parse(examples / 'puc19_bsai_ab.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b']

    dna.replace_trace(examples / 'puc19_bsai_b', examples / 'puc19_bsai_c.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_c']

def test_replace_target_dups(examples):
    dna = snap.parse(examples / 'puc19_bsai_ab.dna')
    dna.append_trace(examples / 'puc19_bsai_a.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_a']

    dna.replace_trace(examples / 'puc19_bsai_a', examples / 'puc19_bsai_c.ab1')
    assert dna.count_traces() == count_seq_blocks(dna) == 2
    assert dna.trace_names == [
            'puc19_bsai_c', 'puc19_bsai_b']

def test_sort_traces(examples):
    dna = snap.parse(examples / 'puc19_bsai_abc.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_c']

    dna.sort_traces(reverse=True)
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_c', 'puc19_bsai_b', 'puc19_bsai_a']

def test_clear_traces(examples):
    dna = snap.parse(examples / 'puc19_bsai_abc.dna')
    assert dna.count_traces() == count_seq_blocks(dna) == 3
    assert dna.trace_names == [
            'puc19_bsai_a', 'puc19_bsai_b', 'puc19_bsai_c']

    dna.clear_traces()
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

    # Not an error to clar an empty sequence.
    dna.clear_traces()
    assert dna.count_traces() == count_seq_blocks(dna) == 0
    assert dna.trace_names == []

def test_extract_traces(examples, tmp_path):
    dna = snap.parse(examples / 'puc19_bsai_abc.dna')
    dna.extract_traces(tmp_path)

    assert (tmp_path / 'puc19_bsai_a.ztr').exists()
    assert (tmp_path / 'puc19_bsai_b.ztr').exists()
    assert (tmp_path / 'puc19_bsai_c.ztr').exists()


def count_seq_blocks(dna):
    return len(dna.find_blocks(snap.blocks.AlignedSequenceBlock))
