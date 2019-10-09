#!/usr/bin/env python3

from . import parser
from docopt import docopt
from pathlib import Path

commands = {}

def command(f):
    commands[f.__name__] = f
    return f

def parse_cli(main=None, **kwargs):
    import inspect
    from textwrap import dedent

    if main is None:
        frame = inspect.stack()[1]
        main = globals()[frame.function]

    return docopt(dedent(main.__doc__.format(**kwargs)))



def main():
    """\
    Perform batch operations on SnapGene files.

    Usage:
        autosnapgene <command> [<args>...] [options]

    Commands:
        {commands}
    """
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in commands:
        return commands[sys.argv[1]]()
    else:
        parse_cli(commands='\n'.join(commands))

@command
def seq():
    """\
    Query or edit the sequence represented by a SnapGene file.

    Usage:
        autosnapgene seq get <dna_path>
        autosnapgene seq set <dna_path> <seq> [-o <dna_path>]
        autosnapgene seq upper <dna_path> [-o <dna_path>]
        autosnapgene seq lower <dna_path> [-o <dna_path>]

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.

    """

    # - Need to think about CLI for get/set, topology/strandedness/methylation
    # - Upper/lower
    # - Also need to make setters in SnapGene class work is the block doesn't 
    #   already exist.

    args = parse_cli()
    dna = parser.parse(args['<dna_path>'])

    if args['get']:
        print(dna.sequence)

    if args['set']:
        dna.sequence = args['<seq>']
        dna.write(args['--out'])

    if args['upper']:
        dna.sequence = dna.sequence.upper()
        dna.write(args['--out'])

    if args['lower']:
        dna.sequence = dna.sequence.lower()
        dna.write(args['--out'])

@command
def info():
    pass

@command
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

@command
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


