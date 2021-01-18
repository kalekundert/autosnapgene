#!/usr/bin/env python3

from . import api

from docopt import docopt
from pathlib import Path
from nonstdlib import indices_from_str
import sys, os
import textwrap

commands = {}

def command(f):
    commands[f.__name__] = f
    return f

def parse_cli(main=None, argv=None, **kwargs):
    import inspect
    from textwrap import dedent

    if main is None:
        frame = inspect.stack()[1]
        main = globals()[frame.function]

    doc = dedent(main.__doc__.format(**kwargs))
    return docopt(doc, argv=argv)

def format_command_descriptions():
    # Note that the docstrings are all indented, and will be de-dented before 
    # rendering to the stdout.  That's why indent needs to be 8 and not 4.
    indent_len, pad_len, line_len = 8, 2, 79
    max_cmd_len = max(len(x) for x in commands)
    desc_len = line_len - indent_len - pad_len - max_cmd_len

    doc = ''
    indent = ' ' * indent_len

    for cmd, func in commands.items():
        pad = ' ' * (pad_len + (max_cmd_len - len(cmd)))
        desc = textwrap.shorten(
                (func.__doc__ or '').split('\n')[0],
                width=desc_len,
                placeholder='...',
        )
        doc += f'{indent}{cmd}{pad}{desc}\n'

    return doc.strip()


def main():
    """\
    Perform batch operations on SnapGene files.

    Usage:
        autosnapgene <command> [<args>...] [options]

    Commands:
        {commands}
    """
    import sys

    try:
        if len(sys.argv) > 1 and sys.argv[1] in commands:
            return commands[sys.argv[1]]()
        else:
            parse_cli(argv=['-h'], commands=format_command_descriptions())

    except FileNotFoundError as err:
        sys.exit(err)

@command
def seq():
    """\
    Query or edit the sequence represented by a SnapGene file.

    Usage:
        autosnapgene seq get <dna_path>
        autosnapgene seq len <dna_path>
        autosnapgene seq set <dna_path> <seq> [-o <dna_path>]
        autosnapgene seq upper <dna_path> [-o <dna_path>]
        autosnapgene seq lower <dna_path> [-o <dna_path>]

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.
    """
    args = parse_cli()
    dna = api.parse(args['<dna_path>'])

    if args['get']:
        print(dna.sequence)

    if args['len']:
        print(len(dna.sequence))

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
def feature():
    """\
    Add or remove features from a SnapGene file.

    Usage:
        autosnapgene feature add <dna_path> <name> <seq> [options]
        autosnapgene feature remove <dna_path> <name> [-o <dna_path>]
        autosnapgene feature clear <dna_path> <seq> [-o <dna_path>]

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.

        -c --color <color>
            The color of the feature being added in HTML format, e.g. '#99ccff'

        -t --type <type>
            The type of feature to add.  Some common options are "CDS", "RBS", 
            "promoter", "terminator", etc.

        -n --note <desc>
            A description of the feature being added.  If "-", the description 
            will be read from stdin.
            
        -T --translate
            Indicate that the feature being added should be translated.
    """
    args = parse_cli()
    dna = api.parse(args['<dna_path>'])

    if args['add']:
        feature = api.Feature.from_segment()
        feature.name = args['<name>']

        if args['--color']:
            feature.segment.color = args['--color']
        if args['--type']:
            feature.segment.type = args['--type']
        if args['--translate']:
            feature.segment.is_translated = True
        if args['--note']:
            if args['--note'] == '-':
                feature.qualifiers = {'note': sys.stdin.read()}
            else:
                feature.qualifiers = {'note': args['--note']}

        dna.write(args['--out'])

    if args['remove']:
        try: dna.remove_feature(args['<name>'])
        except api.FeatureNotFound: pass
        dna.write(args['--out'])

    if args['clear']:
        dna.clear_features()
        dna.write(args['--out'])

@command
def trace():
    """\
    Add, remove, or query sequencing traces within a SnapGene file.

    Usage:
        autosnapgene trace list <dna_path>
        autosnapgene trace add <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace append <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace prepend <dna_path> <ab1_paths>... [-o <dna_path>]
        autosnapgene trace remove <dna_path> <trace_name>... [-o <dna_path>]
        autosnapgene trace pick <dna_path> <trace_name> [-o <dna_path>]
        autosnapgene trace sort <dna_path> [-o <dna_path>]
        autosnapgene trace clear <dna_path> [-o <dna_path>]
        autosnapgene trace extract <dna_path> <out_dir>

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.
    """
    args = parse_cli()
    dna = api.parse(args['<dna_path>'])

    def apply_and_save(method):
        method()
        dna.write(args['--out'])

    def apply_ab1_and_save(method):
        for ab1 in args['<ab1_paths>']:
            method(Path(ab1))
        dna.write(args['--out'])

    def apply_name_and_save(method):
        for name in args['<trace_name>']:
            method(name)
        dna.write(args['--out'])

    if args['list']:
        for trace in sorted(dna.traces, key=lambda x: (x.sort_order, x.name)):
            print(trace.name)
    if args['add']:
        apply_ab1_and_save(dna.add_trace)
    if args['append']:
        apply_ab1_and_save(dna.append_trace)
    if args['prepend']:
        apply_ab1_and_save(dna.prepend_trace)
    if args['remove']:
        apply_name_and_save(dna.remove_trace)
    if args['pick']:
        apply_name_and_save(dna.pick_trace)
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
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.
    """
    args = parse_cli()
    dna = api.parse(args['<dna_path>'])

    if args['clear']:
        dna.clear_history()
        dna.write(args['--out'])

@command
def debug():
    """\
    Various utilities for inspecting and debugging SnapGene files.

    Usage:
        autosnapgene debug list-blocks <dna_path> (-I <id> | -i <i>)
        autosnapgene debug dump-blocks <dna_path> (-I <id> | -i <i>) [-bx]
        autosnapgene debug remove-blocks <dna_path> (-I <id> | -i <i>) [-o <dna_path>]
        autosnapgene debug parse <dna_path> [-o <dna_path>]

    Options:
        -o --out <dna_path>
            Save the modified file to the given path, and leave the input file 
            unmodified.  The default is to overwrite the input file.

        -I --id <id>
            Select blocks by id.  You can select multiple blocks by using 
            commas and dashes, e.g. 1-3,5

        -i --index <i>
            Select blocks by index, starting from 0.  You can select multiple 
            blocks by using commas and dashes, e.g. 1-3,5

        -b --bytes
            Dump block data in raw binary format, rather than in string format.  
            This will produce suitable input for a hex editor.

        -x --xml
            Dump the block data in pretty-printed XML format, rather than 
            string format.  This will only work for blocks that are actually 
            XML data.

    """
    args = parse_cli()
    dna = api.parse(args['<dna_path>'])

    def is_specified(i, block):
        if args['--id']:
            return block.block_id in indices_from_str(args['--id'])
        if args['--index']:
            return i in indices_from_str(args['--index'])

    if args['parse']:
        dna.write(args['--out'])

    if args['list-blocks']:
        for i, block in enumerate(dna.blocks):
            if is_specified(i, block):
                print(f"{i}: {block}")

    if args['dump-blocks']:
        for i, block in enumerate(dna.blocks):
            if is_specified(i, block):
                if args['--bytes']:
                    fp = os.fdopen(sys.stdout.fileno(), 'wb')
                    fp.write(block.bytes)
                if args['--xml']:
                    from xml.dom import minidom 
                    xml = minidom.parseString(block.bytes.decode('utf8'))
                    print(xml.toprettyxml())
                else:
                    print(block.bytes)
                    print()

    if args['remove-blocks']:
        dna.blocks = [
                b for i, b in enumerate(dna.blocks)
                if not is_specified(i, b)
        ]
        dna.write(args['--out'])
