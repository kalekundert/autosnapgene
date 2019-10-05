#!/usr/bin/env python3

from . import parser
from pathlib import Path

def main():
    """\
    Perform batch operations on SnapGene files.

    Usage:
        autosnapgene <command> [<args>...] [options]

    Commands:
        trace
        history
    """
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in globals():
        return globals()[sys.argv[1]]()
    else:
        parse_cli()

def trace():
    """\
    Add, remove, or query sequencing traces within a SnapGene file.

    Usage:
        autosnapgene trace add <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace append <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace prepend <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace remove <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace sort <dna_path> [-o <dna_path>]
        autosnapgene trace clear <dna_path> [-o <dna_path>]
        autosnapgene trace extract <dna_path> <out_dir>

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.
    """
    args = parse_cli()
    dna = parser.parse(args['<dna_path>'])

    def apply_and_save(method):
        method()
        dna.write(args['--out'])
    def apply_ab1_and_save(method):
        for ab1 in args['<ab1_paths>']:
            method(Path(ab1))
        dna.write(args['--out'])


    if args['add']:
        apply_ab1_and_save(dna.add_trace)
    if args['append']:
        apply_ab1_and_save(dna.append_trace)
    if args['prepend']:
        apply_ab1_and_save(dna.prepend_trace)
    if args['remove']:
        apply_ab1_and_save(dna.remove_trace)
    if args['sort']:
        apply_and_save(dna.sort_traces)
    if args['clear']:
        apply_and_save(dna.clear_traces)
    if args['extract']:
        dna.extract_traces(args['<out_dir>'])

def history():
    """\
    Remove all history from the given SnapGene file.

    Usage:
        autosnapgene history clear <dna_path> [-o <dna_path>]

    Options:
        -o --out <dna_path>
    """
    args = parse_cli()
    dna = parser.parse(args['<dna_path>'])

    if args['clear']:
        dna.clear_history()
        dna.write(args['--out'])


def parse_cli(main=None):
    import inspect
    from docopt import docopt
    from textwrap import dedent

    if main is None:
        frame = inspect.stack()[1]
        main = globals()[frame.function]

    return docopt(dedent(main.__doc__))


