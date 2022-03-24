#!/usr/bin/env python3

import autoprop
from pathlib import Path
from copy import deepcopy
from more_itertools import one
from . import parser, blocks

# Import classes that should be part of the public API:
from .blocks import AlignmentMetadata
from .blocks import Feature, FeatureSegment
from .blocks import Reference
from .errors import *

def parse(path, block_classes=None):
    """
    Parse the given file and return a `SnapGene` object.

    This is an alias for ``Snapgene(path)``.
    """
    return SnapGene(path)

def write(path, dna):
    """
    Write the given `SnapGene` object (dna) to the given path.

    This is an alias for ``dna.write(path)``.
    """
    dna.write(path)

@autoprop
class SnapGene:

    def __init__(self, path=None):
        if path:
            self.parse(path)
        else:
            self.reset()

    def reset(self):
        self.blocks = [
                blocks.HeaderBlock(),
                blocks.RestrictionDigestBlock(),
        ]
        self.input_path = None

    def parse(self, path, block_classes=None):
        self.blocks = parser.blocks_from_file(path, block_classes)
        self.input_path = Path(path)

    def write(self, path=None):
        if path is None:
            path = self.input_path
        if path is None:
            raise ValueError("not originally parsed from *.dna file, output path required.")

        # Make sure there's a header.
        assert isinstance(self.blocks[0], blocks.HeaderBlock)

        parser.file_from_blocks(path, self.blocks)


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
            return self.make_block(cls)

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
            self.find_block(blocks.DnaBlock).sequence = value
        except BlockNotFound:
            try:
                self.find_block(blocks.ProteinBlock).sequence = value
            except BlockNotFound:
                self.make_block(blocks.DnaBlock).sequence = value

    def get_topology(self):
        return self.find_block(blocks.DnaBlock).topology

    def set_topology(self, value):
        self.find_or_make_block(blocks.DnaBlock).topology = value

    def get_strandedness(self):
        return self.find_block(blocks.DnaBlock).strandedness

    def set_strandedness(self, value):
        self.find_or_make_block(blocks.DnaBlock).strandedness = value

    def get_is_dam_methylated(self):
        return self.find_block(blocks.DnaBlock).is_dam_methylated

    def set_is_dam_methylated(self, value):
        self.find_or_make_block(blocks.DnaBlock).is_dam_methylated = value

    def get_is_dcm_methylated(self):
        return self.find_block(blocks.DnaBlock).is_dcm_methylated

    def set_is_dcm_methylated(self, value):
        self.find_or_make_block(blocks.DnaBlock).is_dcm_methylated = value

    def get_is_ecoki_methylated(self):
        return self.find_block(blocks.DnaBlock).is_ecoki_methylated

    def set_is_ecoki_methylated(self, value):
        self.find_or_make_block(blocks.DnaBlock).is_ecoki_methylated = value

    def get_dna_sequence(self):
        return self.find_block(blocks.DnaBlock).sequence

    def set_dna_sequence(self, value):
        self.find_or_make_block(blocks.DnaBlock).sequence = value

    # Protein

    def get_protein_sequence(self):
        return self.find_block(blocks.ProteinBlock).sequence

    def set_protein_sequence(self, value):
        self.find_or_make_block(blocks.ProteinBlock).sequence = value

    # Features

    def get_features(self):
        """
        Return a list of all the features in this sequence.
        """
        try:
            return self.find_block(blocks.FeaturesBlock).features
        except BlockNotFound:
            return []

    def get_feature(self, name):
        """
        Return the feature with the given name.

        If multiple features have the given name, only the first one 
        encountered will be returned.  A `FeatureNotFound` exception will be 
        raised if no feature with the given name can be found.
        """
        for feat in self.features:
            if feat.name == name:
                return feat

        raise FeatureNotFound(f"no feature named '{name}'")

    def count_features(self):
        """
        Return the number of features in this sequence.
        """
        return len(self.features)

    def add_feature(self, feature, seq=None):
        """
        Add the given feature to the sequence.

        The argument should be a Feature instance.  Colors and positions are 
        specified as segments, which you can provide by filling in the 
        segments attribute of the Feature instance with FeatureSegments 
        instances.

        If you specify the optional *seq* argument, the position of the feature 
        will automatically be set to the position of that subsequence in the 
        full sequence.  If the given subsequence appears multiple times, 
        multiple features will be created.  If the subsequence doesn't appear, 
        a `SequenceNotFound` exception will be raised.  The given feature must 
        have either 0 or 1 segments, otherwise it isn't clear how the position 
        should be set.
        """
        block = self.find_or_make_block(blocks.FeaturesBlock)
        new_features = []

        if not seq:
            feature.id = block.next_id
            new_features.append(feature)

        else:
            import re, copy

            # Find all occurrences of the given sequence.

            positions = [
                    m.start()
                    # Use lookahead to allow overlapping matches.
                    for m in re.finditer(f'(?={seq.upper()})', self.sequence.upper())
            ]
            if not positions:
                raise SequenceNotFound(f"'{seq}' not found in sequence")

            # Add a copy of the given feature positioned over each occurrence 
            # of the given sequence.

            segments = getattr(feature, 'segments', [])
            if len(segments) > 1:
                raise ValueError(f"{feature} has multiple segments, unclear how to set position from sequence.")

            for i in positions:
                feat = deepcopy(feature)
                if not segments:
                    feat.segment = FeatureSegment()

                # Indexing starts at 1, per the spec.
                feat.segment.range = (i + 1, i + len(seq))
                feat.id = block.next_id
                new_features.append(feat)

        block.features += new_features
        return new_features

    def remove_feature(self, name):
        """
        Remove the feature with the given name from this sequence.

        Only the feature annotation is removed; the sequence corresponding to 
        the feature remains.  If no feature with the given name can be found, a 
        `FeatureNotFound` exception is raised.
        """
        try:
            block = self.find_block(blocks.FeaturesBlock)
        except BlockNotFound:
            raise FeatureNotFound(f"no feature named '{name}'")

        found_name = False
        for feat in block.features[:]:
            if feat.name == name:
                block.features.remove(feat)
                found_name = True

        if not found_name:
            raise FeatureNotFound(f"no feature named '{name}'")

        # Make the id numbers contiguous.  I don't think this is necessary, but 
        # it seems like the right thing to do.

        for new_id, feat in enumerate(block.features):
            feat.id = new_id

    def clear_features(self):
        """
        Remove all features from the sequence.
        """
        self.remove_block(blocks.FeaturesBlock)

    def extract_features(self):
        raise NotImplementedError

    # Alignment

    def get_traces(self, name=None):
        """
        Return information about all of the alignments/traces associated with 
        the sequence.

        If a name is provided, only traces with that name will be returned.  If 
        a `pathlib.Path` is given as the name,  the stem of that path will be 
        taken as the name.
        """
        try:
            metadata = self.find_block(blocks.AlignmentsBlock).metadata
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

        If a `pathlib.Path` is given as the name,  the stem of that path will 
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
        ztr = parser.ztr_from_data(data)

        # Figure out the next id from the AlignmentsBlock metadata.
        align_block = self.find_or_make_block(blocks.AlignmentsBlock)
        next_id = max((x.id for x in align_block.metadata), default=0) + 1
        next_order = max((x.sort_order for x in align_block.metadata), default=0) + 1

        # Make a new AlignedTraceBlock with the ZTR data.
        trace_block = blocks.AlignedTraceBlock()
        trace_block.bytes = ztr

        # Make a new AlignedSequenceBlock with the above id and trace block.
        seq_block = blocks.AlignedSequenceBlock()
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

        self.sync_trace_metadata()

    def remove_trace(self, name):
        """
        Remove the trace with the given name.

        If a `pathlib.Path` is given as the name,  the stem of that path will 
        be taken as the name.  If there are multiple traces with the same name, 
        all will be removed.  If there are no traces with the given name, a 
        ValueError will be raised.
        """

        if isinstance(name, Path):
            name = name.stem

        # Remove the blocks/metadata corresponding to the given name.

        align_block = self.find_block(blocks.AlignmentsBlock)
        seq_blocks = {
                block.id: block
                for block in self.find_blocks(blocks.AlignedSequenceBlock)
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

        self.sync_trace_metadata()

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

    def sync_trace_metadata(self):
        align_block = self.find_block(blocks.AlignmentsBlock)

        # Remove old entries from the metadata.  I don't know how these entries 
        # appear.  I initially wrote the code to get rid of these entries 
        # because I thought they might be causing problems, but further 
        # experimentation made that seem like the wrong hypothesis.  I'm 
        # keeping the code, though, because it doesn't seem like it could hurt.  
        existing_ids = set()
        for block in self.find_blocks(blocks.AlignedSequenceBlock):
            existing_ids.add(block.id)

        align_block.metadata = [
                meta
                for meta in align_block.metadata
                if meta.id in existing_ids
        ]

        # Make sure the sort order matches the actual list order.
        for i, meta in enumerate(align_block.metadata):
            meta.sort_order = i

    def pick_trace(self, name):
        """
        Hide every trace except for the one specified.
        """
        if isinstance(name, Path):
            name = name.stem

        for trace in self.traces:
            trace.is_visible = (trace.name == name)
            debug(trace)

    def count_traces(self):
        """
        Return the number of traces associated with the sequence.
        """
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
        self.remove_block(blocks.AlignmentsBlock)
        self.remove_blocks(blocks.AlignedSequenceBlock)

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
                for block in self.find_blocks(blocks.AlignedSequenceBlock)
        }

        for meta in self.get_traces():
            seq_block = seq_blocks[meta.id]
            trace_blocks = seq_block.traces

            for i, trace_block in enumerate(trace_blocks):
                suffix = f'_{i+1}' if len(trace_blocks) > 1 else ''
                path = dir / f'{meta.name}{suffix}.ztr'
                path.write_bytes(trace_block.bytes)

    # Primers

    # add_primer
    # remove_primer
    # clear_primers
    # count_primers
    # extract_primers
    
    # History

    def clear_history(self):
        """
        Remove all history from the sequence.
        """
        self.remove_blocks(blocks.HistoryBlock)
        self.remove_blocks(blocks.HistoryNodeBlock)

    # Notes

    def get_plasmid_type(self):
        """
        "Natural" or "Synthetic".
        """
        return self.find_block(blocks.NotesBlock).type

    def set_plasmid_type(self, value):
        self.find_block(blocks.NotesBlock).type = value

    def get_custom_map_label(self):
        return self.find_block(blocks.NotesBlock).custom_map_label

    def set_custom_map_label(self, value):
        self.find_block(blocks.NotesBlock).custom_map_label = value

    def get_use_custom_map_label(self):
        return self.find_block(blocks.NotesBlock).use_custom_map_label

    def set_use_custom_map_label(self, value):
        self.find_block(blocks.NotesBlock).use_custom_map_label = value

    def get_is_confirmed_experimentally(self):
        """
        True if this sequence has been experimentally confirmed.
        """
        return self.find_block(blocks.NotesBlock).is_confirmed_experimentally

    def set_is_confirmed_experimentally(self, value):
        self.find_block(blocks.NotesBlock).is_confirmed_experimentally = value

    def get_description(self):
        """
        A description of the sequence.
        """
        return self.find_block(blocks.NotesBlock).description

    def set_description(self, value):
        self.find_block(blocks.NotesBlock).description = value

    def get_date_created(self):
        """
        The date the sequence was created.
        """
        return self.find_block(blocks.NotesBlock).date_created

    def set_date_created(self, value):
        self.find_block(blocks.NotesBlock).date_created = value

    def get_date_last_modified(self):
        """
        The date the sequence was last modified.
        """
        return self.find_block(blocks.NotesBlock).date_last_modified

    def set_date_last_modified(self, value):
        self.find_block(blocks.NotesBlock).date_last_modified = value

    def get_accession_number(self):
        return self.find_block(blocks.NotesBlock).accession_number

    def set_accession_number(self, value):
        self.find_block(blocks.NotesBlock).accession_number = value

    def get_code_number(self):
        return self.find_block(blocks.NotesBlock).code_number

    def set_code_number(self, value):
        self.find_block(blocks.NotesBlock).code_number = value

    def get_author(self):
        """
        The creator of this sequence.
        """
        return self.find_block(blocks.NotesBlock).author

    def set_author(self, value):
        self.find_block(blocks.NotesBlock).author = value

    def get_organism(self):
        """
        The organism this sequence derives from.
        """
        return self.find_block(blocks.NotesBlock).organism

    def set_organism(self, value):
        self.find_block(blocks.NotesBlock).organism = value

    def get_sequence_class(self):
        return self.find_block(blocks.NotesBlock).sequence_class

    def set_sequence_class(self, value):
        self.find_block(blocks.NotesBlock).sequence_class = value

    def get_transformed_into(self):
        """
        The organism/strain being used to propagate this sequence in the lab.
        """
        return self.find_block(blocks.NotesBlock).transformed_into

    def set_transformed_into(self, value):
        self.find_block(blocks.NotesBlock).transformed_into = value

    def get_comments(self):
        """
        Miscellaneous comments on this sequence.
        """
        return self.find_block(blocks.NotesBlock).comments

    def set_comments(self, value):
        self.find_block(blocks.NotesBlock).comments = value

    def get_references(self):
        return self.find_block(blocks.NotesBlock).references

    def set_references(self, value):
        self.find_block(blocks.NotesBlock).references = value

    # File format
    
    def get_file_type(self):
        return self.find_block(blocks.HeaderBlock).file_type

    def set_file_type(self, value):
        self.find_block(blocks.HeaderBlock).file_type = value

    def get_import_version(self):
        return self.find_block(blocks.HeaderBlock).import_version

    def set_import_version(self, value):
        self.find_block(blocks.HeaderBlock).import_version = value

    def get_export_version(self):
        return self.find_block(blocks.HeaderBlock).export_version

    def set_export_version(self, value):
        self.find_block(blocks.HeaderBlock).export_version = value

    @property
    def _this_seq(self):
        return f"'{self.input_path}'" if self.input_path else "this sequence"



