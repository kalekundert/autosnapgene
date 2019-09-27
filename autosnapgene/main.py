#!/usr/bin/env python3

from . import parser
from docopt import docopt
from functools import wraps
from textwrap import dedent
from pathlib import Path

def main():
    """\
    Perform batch operations on SnapGene files.

    Usage:
        autosnapgene add_traces <dna_path> <ab1_paths>...
        autosnapgene clear_history <dna_path>
    """
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in globals():
        return globals()[sys.argv[1]]()
    else:
        parse_cli()

def add_traces():
    """\
    Add sequencing traces to a SnapGene file.

    Usage:
        autosnapgene add_traces <dna_path> <ab1_paths>... [-o <dna_path>]

    Options:
        -o --out <dna_path>
    """
    args = parse_cli()
    dna = parser.parse(args['<dna_path>'])

    for ab1 in args['<ab1_paths>']:
        dna.add_trace(ab1)

    dna.write(args['--out'])

def clear_history():
    """\
    Remove all history from the given SnapGene file.

    Usage:
        autosnapgene clear_history <dna_path> [-o <dna_path>]

    Options:
        -o --out <dna_path>
    """
    args = parse_cli()
    dna = parser.parse(args['<dna_path>'])
    dna.clear_history()
    dna.write(args['--out'])


def parse_cli(main=None):
    import inspect

    if main is None:
        frame = inspect.stack()[1]
        main = globals()[frame.function]

    return docopt(dedent(main.__doc__))


