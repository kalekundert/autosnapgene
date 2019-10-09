#!/usr/bin/env python3

import struct
import autoprop
import arrow
import xml.etree.ElementTree as etree
import html
from pathlib import Path
from more_itertools import one

block_ids = {}
block_classes = {}

def parse(path):
    dna = SnapGene()
    dna.parse(path)
    return dna

def write(path, dna):
    dna.write(path)


def blocks_from_file(path):
    bytes = Path(path).read_bytes()
    return blocks_from_bytes(bytes)

def blocks_from_bytes(bytes):
    try:
        i = 0
        blocks = []

        while i < len(bytes):
            j = i + 5
            id, size = struct.unpack('>BI', bytes[i:j])
            cls = block_classes.get(id, UndocumentedBlock)
            block = cls.from_bytes(bytes[j:j+size])
            block.block_id = id
            blocks.append(block)
            i = j + size

    except SnapGeneError as e:
        e.path = path
        raise e from None

    return blocks

def file_from_blocks(path, blocks):
    bytes = bytes_from_blocks(blocks)
    Path(path).write_bytes(bytes)

def bytes_from_blocks(blocks):
    bytes = b''
    for block in blocks:
        bytes += bytes_from_block(block)
    return bytes

def bytes_from_block(block):
    bytes = block.to_bytes()
    header = struct.pack('>BI', block.block_id, len(bytes))
    return header + bytes

def ztr_from_data(data):
    from subprocess import run, PIPE

    # The details of the ZTR fromat are described in this publication:
    #
    #   Bonfield and Staden.  ZTR: a new format for DNA sequence trace data.  
    #   Bioinformatics 18:1:3-10 (2002).  DOI: 10.1093/bioinformatics/18.1.3
    #
    # The tools for doing the conversion can be found on GitHub:
    #
    #   https://github.com/jkbonfield/io_lib.git
    # 
    # The same tools can also be found on AUR:
    #
    #   $ yay -S staden-io_lib
    #
    # The command for doing the conversion is:
    # 
    #   $ convert_trace < path/to/ab1 > path/to/ztr
    #
    # Note that this command automatically detects the format of the input 
    # trace, and a number of formats are supported.  The are a number of 
    # command-line flags available, see -h for more information.

    p = run(['convert_trace'], input=data, capture_output=True, check=True)
    return p.stdout

@autoprop
class SnapGene:

    def __init__(self, path=None):
        if path:
            self.parse(path)

        self.blocks = []
        self.input_path = None

    def parse(self, path):
        self.blocks = blocks_from_file(path)
        self.input_path = Path(path)

    def write(self, path=None):
        if path is None:
            path = self.input_path
        if path is None:
            raise ValueError("not originally parsed from *.dna file, output path required.")

        file_from_blocks(path, self.blocks)


    def make_block(self, cls):
        block = cls()
        self.blocks.append(block)
        return block

    def find_blocks(self, cls):
        return [x for x in self.blocks if isinstance(x, cls)]

    def find_block(self, cls):
        return one(
                self.find_blocks(cls),
                too_short=BlockNotFound(f"{self._this_seq} doesn't have any {cls.__name__} blocks."),
                too_long=AssertionError,
        )

    def find_or_make_block(self, cls):
        try:
            return self.find_block(cls)
        except BlockNotFound:
            return self.make_block()

    def remove_blocks(self, cls):
        self.blocks = [x for x in self.blocks if not isinstance(x, cls)]

    def remove_block(self, cls):
        try:
            block = self.find_block(cls)
            self.blocks.remove(block)
        except BlockNotFound:
            pass

    # DNA

    def get_sequence(self):
        """
        Get the DNA or protein sequence stored in this file.

        This method will not indicate whether it returned a protein or DNA 
        sequence, so if that is not clear from the context, use 
        `get_dna_sequence()` or `get_protein_sequence()` to be sure.  If the 
        file doesn't contain either type of sequence, a `BlockNotFound` 
        exception will be raised.
        """
        try:
            return self.dna_sequence
        except BlockNotFound:
            return self.protein_sequence

    def set_sequence(self, value):
        """
        Set the DNA or protein sequence stored in this file.

        If the file does not yet contain a sequence, the DNA sequence will be 
        set.
        """
        try:
            self.find_block(DnaBlock).sequence = value
        except BlockNotFound:
            try:
                self.find_block(ProteinBlock).sequence = value
            except BlockNotFound:
                self.make_block(DnaBlock).sequence = value

    def get_topology(self):
        return self.find_block(DnaBlock).topology

    def set_topology(self, value):
        self.find_or_make_block(DnaBlock).topology = value

    def get_strandedness(self):
        return self.find_block(DnaBlock).strandedness

    def set_strandedness(self, value):
        self.find_or_make_block(DnaBlock).strandedness = value

    def get_is_dam_methylated(self):
        return self.find_block(DnaBlock).is_dam_methylated

    def set_is_dam_methylated(self, value):
        self.find_or_make_block(DnaBlock).is_dam_methylated = value

    def get_is_dcm_methylated(self):
        return self.find_block(DnaBlock).is_dcm_methylated

    def set_is_dcm_methylated(self, value):
        self.find_or_make_block(DnaBlock).is_dcm_methylated = value

    def get_is_ecoki_methylated(self):
        return self.find_block(DnaBlock).is_ecoki_methylated

    def set_is_ecoki_methylated(self, value):
        self.find_or_make_block(DnaBlock).is_ecoki_methylated = value

    def get_dna_sequence(self):
        return self.find_block(DnaBlock).sequence

    def set_dna_sequence(self, value):
        self.find_or_make_block(DnaBlock).sequence = value

    # Protein

    def get_protein_sequence(self):
        return self.find_block(ProteinBlock).sequence

    def set_protein_sequence(self, value):
        self.find_or_make_block(ProteinBlock).sequence = value

    # Header
    
    def get_file_type(self):
        return self.find_block(HeaderBlock).file_type

    def set_file_type(self, value):
        self.find_block(HeaderBlock).file_type = value

    def get_import_version(self):
        return self.find_block(HeaderBlock).import_version

    def set_import_version(self, value):
        self.find_block(HeaderBlock).import_version = value

    def get_export_version(self):
        return self.find_block(HeaderBlock).export_version

    def set_export_version(self, value):
        self.find_block(HeaderBlock).export_version = value

    # Notes

    def get_plasmid_type(self):
        """
        "Natural" or "Synthetic".
        """
        return self.find_block(NotesBlock).type

    def set_plasmid_type(self, value):
        self.find_block(NotesBlock).type = value

    def get_custom_map_label(self):
        return self.find_block(NotesBlock).custom_map_label

    def set_custom_map_label(self, value):
        self.find_block(NotesBlock).custom_map_label = value

    def get_use_custom_map_label(self):
        return self.find_block(NotesBlock).use_custom_map_label

    def set_use_custom_map_label(self, value):
        self.find_block(NotesBlock).use_custom_map_label = value

    def get_is_confirmed_experimentally(self):
        """
        True if this sequence has been experimentally confirmed.
        """
        return self.find_block(NotesBlock).is_confirmed_experimentally

    def set_is_confirmed_experimentally(self, value):
        self.find_block(NotesBlock).is_confirmed_experimentally = value

    def get_description(self):
        """
        A description of the sequence.
        """
        return self.find_block(NotesBlock).description

    def set_description(self, value):
        self.find_block(NotesBlock).description = value

    def get_date_created(self):
        """
        The date the sequence was created.
        """
        return self.find_block(NotesBlock).date_created

    def set_date_created(self, value):
        self.find_block(NotesBlock).date_created = value

    def get_date_last_modified(self):
        """
        The date the sequence was last modified.
        """
        return self.find_block(NotesBlock).date_last_modified

    def set_date_last_modified(self, value):
        self.find_block(NotesBlock).date_last_modified = value

    def get_accession_number(self):
        return self.find_block(NotesBlock).accession_number

    def set_accession_number(self, value):
        self.find_block(NotesBlock).accession_number = value

    def get_code_number(self):
        return self.find_block(NotesBlock).code_number

    def set_code_number(self, value):
        self.find_block(NotesBlock).code_number = value

    def get_author(self):
        """
        The creator of this sequence.
        """
        return self.find_block(NotesBlock).author

    def set_author(self, value):
        self.find_block(NotesBlock).author = value

    def get_organism(self):
        """
        The organism this sequence derives from.
        """
        return self.find_block(NotesBlock).organism

    def set_organism(self, value):
        self.find_block(NotesBlock).organism = value

    def get_sequence_class(self):
        return self.find_block(NotesBlock).sequence_class

    def set_sequence_class(self, value):
        self.find_block(NotesBlock).sequence_class = value

    def get_transformed_into(self):
        """
        The organism/strain being used to propagate this sequence in the lab.
        """
        return self.find_block(NotesBlock).transformed_into

    def set_transformed_into(self, value):
        self.find_block(NotesBlock).transformed_into = value

    def get_comments(self):
        """
        Miscellaneous comments on this sequence.
        """
        return self.find_block(NotesBlock).comments

    def set_comments(self, value):
        self.find_block(NotesBlock).comments = value

    def get_references(self):
        return self.find_block(NotesBlock).references

    def set_references(self, value):
        self.find_block(NotesBlock).references = value

    # Alignment

    def get_traces(self, name=None):
        """
        Return information about all of the alignments/traces associated with 
        the sequence.

        If a name is provided, only traces with that name will be returned.  If 
        a Pathlib.Path() is given as the name,  the stem of that path will be 
        taken as the name.
        """
        try:
            metadata = self.find_block(AlignmentsBlock).metadata
        except BlockNotFound:
            return []

        if name is None:
            return metadata
        else:
            if isinstance(name, Path):
                name = name.stem
            return [x for x in metadata if x.name == name]

    def get_trace_names(self):
        return [x.name for x in self.get_traces()]

    def has_trace(self, name):
        """
        Return True is the sequence contains a trace with the given name.

        If a Pathlib.Path() is given as the name,  the stem of that path will 
        be taken as the name.
        """
        return bool(self.get_traces(name))

    def add_trace(self, path, name=None):
        """
        Add the given trace to the sequence, if the sequence doesn't already 
        have a trace of the same name (or the given name).

        The trace is always added after any existing traces.
        """

        if not self.has_trace(name or path):
            self.append_trace(path, name=name)
        else:
            self.replace_trace(name or path, path, new_name=name)

    def append_trace(self, path, name=None):
        """
        Add the given trace to this sequence after any existing traces.

        Unlike add_trace(), this function adds the trace unconditionally, which 
        may result in duplicates.
        """
        self.insert_trace(self.count_traces(), path, name=name)

    def prepend_trace(self, path, name=None):
        """
        Add the given trace to this sequence before any existing traces.

        Unlike add_trace(), this function adds the trace unconditionally, which 
        may result in duplicates.
        """
        self.insert_trace(0, path, name=name)

    def insert_trace(self, i, path, name=None):
        """
        Add the given trace to this sequence at the given index.

        Unlike add_trace(), this function adds the trace unconditionally, which 
        may result in duplicates.
        """
        path = Path(path)

        # Convert the sequencing data to the ZTR format.
        data = path.read_bytes()
        ztr = ztr_from_data(data)

        # Figure out the next id from the AlignmentsBlock metadata.
        align_block = self.find_block(AlignmentsBlock)
        next_id = max((x.id for x in align_block.metadata), default=0) + 1
        next_order = max((x.sort_order for x in align_block.metadata), default=0) + 1

        # Make a new AlignedTraceBlock with the ZTR data.
        trace_block = AlignedTraceBlock()
        trace_block.bytes = ztr

        # Make a new AlignedSequenceBlock with the above id and trace block.
        seq_block = AlignedSequenceBlock()
        seq_block.id = next_id
        seq_block.traces = [trace_block]

        # Add the sequence block to the file (the trace block will be added 
        # indirectly via the sequence block).
        self.blocks.append(seq_block)

        # Update the AlignmentsBlock metadata.
        meta = AlignmentMetadata()
        meta.id = next_id
        meta.name = name or path.stem
        meta.is_trace = True
        align_block.metadata.insert(i, meta)

        for i, meta in enumerate(align_block.metadata):
            meta.sort_order = i

    def remove_trace(self, name):
        """
        Remove the trace with the given name.

        If a Pathlib.Path() is given as the name,  the stem of that path will 
        be taken as the name.  If there are multiple traces with the same name, 
        all will be removed.  If there are no traces with the given name, a 
        ValueError will be raised.
        """

        if isinstance(name, Path):
            name = name.stem

        # Remove the blocks/metadata corresponding to the given name.

        align_block = self.find_block(AlignmentsBlock)
        seq_blocks = {
                block.id: block
                for block in self.find_blocks(AlignedSequenceBlock)
        }

        found_name = False
        for meta in align_block.metadata[:]:
            if meta.name == name:
                align_block.metadata.remove(meta)
                self.blocks.remove(seq_blocks[meta.id])
                found_name = True

        if not found_name:
            raise ValueError(f"no trace named '{name}'")

        # Make the id numbers contiguous.  I don't think this is necessary, but 
        # it seems like the right thing to do.

        for new_id, meta in enumerate(align_block.metadata):
            old_id = meta.id
            meta.id = seq_blocks[old_id].id = new_id

    def rename_trace(self, old_name, new_name):
        """
        Rename the given trace.

        If multiple traces have the same name, they will all be renamed.
        """
        traces = self.get_trace(old_name)
        if not traces:
            raise ValueError("no trace named '{old_name}'")

        for meta in traces:
            meta.name = new_name

    def replace_trace(self, old_name, path, new_name=None):
        """
        Replace the trace with the given name with the given path.

        If there are multiple traces with the given name, the first will be 
        replaced and the rest will be removed.
        """

        if isinstance(old_name, Path):
            old_name = old_name.stem

        traces = self.get_traces()
        if not traces:
            raise ValueError("no trace named '{old_name}'")

        for i, meta in enumerate(traces):
            if meta.name == old_name:
                break

        self.remove_trace(old_name)
        self.insert_trace(i, path, name=new_name)

    def count_traces(self):
        return len(self.get_traces())

    def sort_traces(self, key=lambda x: x.name, reverse=False):
        """
        Rearrange the traces according to the given key function.

        If no key function is given, the traces will be sorted alphabetically.  
        The argument to the key function will be an AlignmentMetadata object.
        """
        traces = self.get_traces()
        traces.sort(key=key, reverse=reverse)
        for i, meta in enumerate(traces):
            meta.sort_order = i

    def clear_traces(self):
        """
        Remove all traces from the sequence.
        """
        self.remove_block(AlignmentsBlock)
        self.remove_blocks(AlignedSequenceBlock)

    def extract_traces(self, dir):
        """
        Save any traces associated with this sequence as separate files in the 
        given directory.

        The traces will be saved in the ZTR format, which is the format used 
        internally by SnapGene.
        """
        dir = Path(dir)
        dir.mkdir(parents=True, exist_ok=True)

        seq_blocks = {
                block.id: block
                for block in self.find_blocks(AlignedSequenceBlock)
        }

        for meta in self.get_traces():
            seq_block = seq_blocks[meta.id]
            trace_blocks = seq_block.traces

            for i, trace_block in enumerate(trace_blocks):
                suffix = f'_{i+1}' if len(trace_blocks) > 1 else ''
                path = dir / f'{meta.name}{suffix}.ztr'
                path.write_bytes(trace_block.bytes)

    # History

    def clear_history(self):
        self.remove_blocks(HistoryBlock)
        self.remove_blocks(HistoryNodeBlock)

    # add_primer
    # remove_primer
    # clear_primers
    # count_primers
    # extract_primers
    
    # add_feature
    # remove_feature
    # clear_features
    # count_features
    # extract_features

    @property
    def _this_seq(self):
        return f"'{self.input_path}'" if self.input_path else "this sequence"

class Block:

    def __init_subclass__(cls):
        if hasattr(cls, 'block_id'):
            block_ids[cls.block_name] = cls.block_id
            block_classes[cls.block_id] = cls

    def __repr__(self):
        return f"<{self.__class__.__name__} block_id={self.block_id} {self.__repr_attrs__()}".strip() + ">"

    def __repr_attrs__(self):
        return ''

    @classmethod
    def from_bytes(cls, bytes):
        raise NotImplementedError

    def to_bytes(self):
        raise NotImplementedError


class UnparsedBlock(Block):

    def __init__(self, bytes=None):
        self.bytes = bytes

    def __repr_attrs__(self):
        bytes_repr = repr(self.bytes)[2:-1]
        if len(bytes_repr) > 32:
            bytes_repr = bytes_repr[:32] + '...'
        return f"UNPARSED bytes=b'{bytes_repr}'"

    @classmethod
    def from_bytes(cls, bytes):
        return cls(bytes)

    def to_bytes(self):
        return self.bytes


@autoprop
class HeaderBlock(Block):
    block_id = 9
    block_name = 'header'

    file_types = {
            0: 'unknown',
            1: 'dna',
            2: 'protein',  # Just a guess...
    }

    def __init__(self):
        self.type_id = None
        self.export_version = None
        self.import_version = None

    def __repr_attrs__(self):
        return f"type='{self.type}' export_version='{self.export_version}' import_version='{self.import_version}'"

    @classmethod
    def from_bytes(cls, bytes):
        magic_cookie = bytes[0:8].decode('ascii')

        if magic_cookie != "SnapGene":
            raise SnapGeneError("not a snapgene file")

        info = struct.unpack('>HHH', bytes[8:14])
        return cls.from_info(*info)

    @classmethod
    def from_info(cls, type_id, export_version, import_version):
        block = cls()
        block.type_id = type_id
        block.export_version = export_version
        block.import_version = import_version
        return block

    def to_bytes(self):
        return b'SnapGene' + struct.pack('>HHH',
                self.type_id, self.export_version, self.import_version)

    def get_type(self):
        return self.file_types[self.type_id]


class DnaBlock(Block):
    block_id = 0
    block_name = 'dna'

    def __init__(self):
        self.topology = None
        self.strandedness = None
        self.is_dam_methylated = None
        self.is_dcm_methylated = None
        self.is_ecoki_methylated = None
        self.sequence = None

    def __repr_attrs__(self):
        if len(self.sequence) > 32:
            truncated_seq = f"{self.sequence[:16]}...{self.sequence[-16:]}"
        else:
            truncated_seq = self.sequence

        return f"{self.topology} {self.strandedness} sequence='{truncated_seq}'"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        props = bytes[0]
        block.topology = 'circular' if props & 0x01 else 'linear'
        block.strandedness = 'double' if props & 0x02 else 'single'
        block.is_dam_methylated = bool(props & 0x04)
        block.is_dcm_methylated = bool(props & 0x08)
        block.is_ecoki_methylated = bool(props & 0x10)

        block.sequence = bytes[1:].decode('ascii')

        return block

    def to_bytes(self):
        props = sum([
                0x01 if self.topology == 'circular' else 0x00,
                0x02 if self.strandedness == 'double' else 0x00,
                0x04 * self.is_dam_methylated,
                0x08 * self.is_dcm_methylated,
                0x10 * self.is_ecoki_methylated,
        ])
        return struct.pack('>B', props) + self.sequence.encode('ascii')

class ProteinBlock(Block):
    block_id = 21
    block_name = 'protein'

    # This block is undocumented.  I attempted to reverse-engineer the file 
    # format from some example files.  Some information might be unavailable.

    def __init__(self):
        self.props = None
        self.sequence = None

    def __repr_attrs__(self):
        if len(self.sequence) > 32:
            return f"sequence='{self.sequence[:16]}...{self.sequence[-16:]}'"
        else:
            return f"sequence='{self.sequence}'"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        # I don't know what the 'props' byte means.  In DNA blocks, there is an 
        # corresponding byte that encodes a bit-field specifying a few 
        # properties of the DNA, but those same properties wouldn't apply to 
        # proteins.
        block.props = bytes[0]
        block.sequence = bytes[1:].decode('ascii')

        return block

    def to_bytes(self):
        return struct.pack('>B', self.props or 0) \
                + self.sequence.encode('ascii')

class CompressedDnaBlock(UnparsedBlock):
    block_id = 1
    block_name = 'compressed_dna'

class PrimerBlock(UnparsedBlock):
    block_id = 5
    block_name = 'primers'

class NotesBlock(Block):
    block_id = 6
    block_name = 'notes'

    class TextTag:

        @staticmethod
        def from_xml(element):
            return element.text

        @staticmethod
        def to_xml(element, value):
            element.text = value

    class BoolTag:

        @staticmethod
        def from_xml(element):
            return {'0': False, '1': True}[element.text]

        @staticmethod
        def to_xml(element, value):
            element.text = str(int(value))

    class DateTag:

        @staticmethod
        def from_xml(element):
            return arrow.get(element.text, 'YYYY.M.D')

        @staticmethod
        def to_xml(element, value):
            element.text = value.format('YYYY.M.D')

    class HtmlTag(TextTag):

        # I made this class because I thought that I'd have to specially escape 
        # and unescape HTML content, but it turns out that the XML library 
        # takes care of that automatically.  So this class ends up being 
        # functionally the same as TextTag.  I'm keeping it because it's still 
        # semantic.

        pass

    class ReferencesTag:

        @staticmethod
        def from_xml(element):
            return [
                    Reference.from_xml(child)
                    for child in element
            ]

        @staticmethod
        def to_xml(element, value):
            for ref in value:
                child = ref.to_xml()
                element.append(child)

    attr_defs = [
            ('uuid', 'UUID', TextTag),
            ('type', 'Type', TextTag),
            ('custom_map_label', 'CustomMapLabel', TextTag),
            ('use_custom_map_label', 'UseCustomMapLabel', BoolTag),
            ('is_confirmed_experimentally', 'ConfirmedExperimentally', BoolTag),
            ('description', 'Description', HtmlTag),
            ('date_created', 'Created', DateTag),
            ('date_last_modified', 'LastModified', DateTag),
            ('accession_number', 'AccessionNumber', TextTag),
            ('code_number', 'CodeNumber', TextTag),
            ('author', 'CreatedBy', TextTag),
            ('organism', 'Organism', TextTag),
            ('sequence_class', 'SequenceClass', TextTag),
            ('transformed_into', 'TransformedInto', TextTag),
            ('comments', 'Comments', HtmlTag),
            ('references', 'References', ReferencesTag),
    ]
    parsers_by_attr = {
            attr: (tag, parser)
            for attr, tag, parser in attr_defs
    }
    parsers_by_tag = {
            tag: (attr, parser)
            for attr, tag, parser in attr_defs
    }

    def __init__(self):
        self.attrs = {}

    def __repr_attrs__(self):
        attrs = 'type', 'created_by', 'last_modified'
        return ' '.join(
                f'{k}="{getattr(self, k)}"'
                for k in attrs
                if k in self.attrs
        )

    def __getattr__(self, attr):
        if attr in self.attrs:
            return self.attrs[attr]

        elif attr in self.parsers_by_attr:
            raise AttributeError(f"note '{attr}' not defined for this sequence.")

        else:
            did_you_mean = '\n    '.join(self.parsers_by_attr)
            raise AttributeError(f"unknown note '{attr}', did you mean:\n    {did_you_mean}")

    def __setattr__(self, attr, value):
        if attr in self.parsers_by_attr:
            self.attrs[attr] = value
        else:
            super().__setattr__(attr, value)

    def __delattr__(self, attr):
        if attr in self.attrs:
            del self.attrs[attr]

        elif attr in self.parsers_by_attr:
            raise AttributeError(f"note '{attr}' not defined for this sequence.")

        else:
            did_you_mean = '\n    '.join(self.parsers_by_attr)
            raise AttributeError(f"unknown note '{attr}', did you mean:\n    {did_you_mean}")



    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        xml = bytes.decode('utf8')
        root = etree.fromstring(xml)

        for element in root:
            attr, parser = cls.parsers_by_tag[element.tag]
            block.attrs[attr] = parser.from_xml(element)

        return block

    def to_bytes(self):
        root = etree.Element('Notes')

        for attr in self.attrs:
            tag, parser = self.parsers_by_attr[attr]
            element = etree.SubElement(root, tag)
            parser.to_xml(element, self.attrs[attr])

        return etree.tostring(root)

class Reference:

    def __init__(self):
        self.title = None
        self.pubmed_id = None
        self.journal = None
        self.authors = None

    def __repr__(self):
        return f"<Reference id={self.id} name='{self.name}' is_trace={int(self.is_trace)} sort_order={self.sort_order}>"

    @classmethod
    def from_xml(cls, element):
        unescape = lambda x: html.unescape(x) if x is not None else x

        ref = cls()
        ref.title = unescape(element.attrib.get('title'))
        ref.pubmed_id = element.attrib.get('pubMedID')
        ref.journal = unescape(element.attrib.get('journal'))
        ref.authors = unescape(element.attrib.get('authors'))

        return ref

    def to_xml(self):
        def attr(k, v, f=lambda x: x):
            if v is not None:
                attrs[k] = f(v)

        attrs = {}
        attr('title', self.title, html.escape)
        attr('pubMedID', self.pubmed_id)
        attr('journal', self.journal, html.escape)
        attr('authors', self.authors, html.escape)

        return etree.Element('Reference', attrs)


class HistoryBlock(UnparsedBlock):
    block_id = 7
    block_name = 'history'

class HistoryNodeBlock(UnparsedBlock):
    block_id = 11
    block_name = 'history_node'

class PropertiesBlock(UnparsedBlock):
    block_id = 8
    block_name = 'properties'

class FeaturesBlock(UnparsedBlock):
    block_id = 10
    block_name = 'features'

class AlignmentsBlock(Block):
    block_id = 17
    block_name = 'alignments'

    def __init__(self):
        self.trim_stringency = None
        self.metadata = []

    def __repr_attrs__(self):
        return f"meta={','.join(str(x.id) for x in self.metadata)}"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()

        xml = bytes.decode('utf8')
        root = etree.fromstring(xml)

        block.trim_stringency = root.attrib.get('trimStringency')
        block.metadata = [
                AlignmentMetadata.from_xml(child)
                for child in root
        ]
        block.metadata.sort(key=lambda x: x.sort_order)
        return block

    def to_bytes(self):
        root = etree.Element('AlignableSequences', {
            'trimStringency': self.trim_stringency,
        })
        for meta in self.metadata:
            child = meta.to_xml()
            root.append(child)

        return etree.tostring(root)

class AlignmentMetadata:

    def __init__(self):
        self.id = None
        self.name = None
        self.is_trace = None
        self.sort_order = None

    def __repr__(self):
        return f"<AlignmentMetadata id={self.id} name='{self.name}' is_trace={int(self.is_trace)} sort_order={self.sort_order}>"

    def __eq__(self, other):
        return all([
            self.id == other.id,
            self.name == other.name,
            self.is_trace == other.is_trace,
            self.sort_order == other.sort_order,
        ])

    @classmethod
    def from_xml(cls, element):
        seq = cls()
        seq.id = int(element.attrib.get('ID'))
        seq.name = element.attrib.get('name')
        seq.is_trace = bool(element.attrib.get('isTrace'))
        seq.sort_order = int(element.attrib.get('sortOrder'))
        return seq

    def to_xml(self):
        return etree.Element('Sequence', {
            'ID': str(self.id),
            'name': self.name,
            'isTrace': str(int(self.is_trace)),
            'sortOrder': str(self.sort_order),
        })


class AlignedSequenceBlock(Block):
    block_id = 16
    block_name = 'aligned_sequence'

    def __init__(self):
        self.id = None
        self.traces = []

    def __repr_attrs__(self):
        return f"id={self.id}"

    @classmethod
    def from_bytes(cls, bytes):
        block = cls()
        block.id, = struct.unpack('>I', bytes[:4])
        block.traces = blocks_from_bytes(bytes[4:])
        return block

    def to_bytes(self):
        bytes = struct.pack('>I', self.id)
        bytes += bytes_from_blocks(self.traces)
        return bytes

class AlignedTraceBlock(UnparsedBlock):
    block_id = 18
    block_name = 'aligned_trace'

class UracilBlock(UnparsedBlock):
    block_id = 19
    block_name = 'uracils'

class DnaColorBlock(UnparsedBlock):
    block_id = 20
    block_name = 'dna_colors'

class UndocumentedBlock(UnparsedBlock):
    pass



class SnapGeneError(Exception):

    def __init__(self, message):
        self.message = message
        self.path = None

    def __str__(self):
        return self.message.format(path=self.path)
class BlockNotFound(AttributeError):
    pass
