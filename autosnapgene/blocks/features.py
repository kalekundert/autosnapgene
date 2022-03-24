#!/usr/bin/env python3

from ..parser import Block, Xml, Repr
from ..util import reverse_complement

import autoprop
import xml.etree.ElementTree as etree
from more_itertools import one, always_iterable

@autoprop
class FeaturesBlock(Xml, Block):
    block_id = 10
    repr_attrs = ['features']

    # I'm not totally sure whats up with the id numbers in this block.  They 
    # aren't mentioned in the 2015 spec...

    class FeatureTag(Xml.AppendListTag):

        @staticmethod
        def from_xml(element):
            return Feature.from_xml(element)

        @staticmethod
        def to_xml(parent, tag, features):
            for feature in features:
                e = feature.to_xml()
                parent.append(e)

    xml_tag = 'Features'
    xml_subtag_defs = [
            ('features', 'Feature', FeatureTag, []),
    ]
    xml_attrib_defs = [
            ('_next_id', 'nextValidID', Xml.IntAttrib),

            # I don't know what this is at all...
            ('_recycled_ids', 'recycledIDs', Xml.TextAttrib),
    ]

    def __repr_attr__(self, attr):
        if attr == 'features':
            from textwrap import shorten
            return shorten(
                    ', '.join(x.name for x in self.features),
                    width=32,
                    placeholder='...',
            )

        else:
            return super().__repr_attr__(attr)

    def to_xml(self):
        self._next_id = len(self.features)
        return super().to_xml()

    def get_next_id(self):
        return len(self.features)


@autoprop
class Feature(Xml, Repr):
    repr_attrs = 'name', 'type'

    class SegmentTag(Xml.AppendListTag):

        @staticmethod
        def from_xml(element):
            return FeatureSegment.from_xml(element)

        @staticmethod
        def to_xml(parent, tag, segments):
            for segment in segments:
                e = segment.to_xml()
                parent.append(e)

    class QualifierTag(Xml.UpdateDictTag):
        data_types = {
                'int': int,
                'text': str,
        }
        data_formats = {v: k for k, v in data_types.items()}

        @classmethod
        def from_xml(cls, element):
            name = element.attrib['name']

            def get_value(sub):
                # This isn't in the spec, but I assume I'll never get multiple 
                # attributes in a <V> tag.
                # 
                # However, it is possible to get no attributes.  I've only 
                # encountered this is one example: a translated feature that 
                # was too short to actually have a translation (e.g. 1 bp).  In 
                # this case, value probably should be "", since "translation" 
                # is normally text.  Assuming that integer data will never just 
                # be left out like this, I'm going to interpret no attributes 
                # as "".
                if not sub.attrib:
                    return ""

                assert len(sub.attrib) == 1, etree.tostring(element)

                key, value = sub.attrib.popitem()
                return cls.data_types.get(key, str)(value)

            if len(element) == 1:
                value = get_value(element.find('V'))
            else:
                value = [get_value(x) for x in element.findall('V')]

            return {name: value}

        @classmethod
        def to_xml(cls, parent, tag, values):
            for name, value in values.items():
                q = etree.SubElement(parent, tag)
                q.attrib['name'] = str(name)

                for value_i in always_iterable(value):
                    v = etree.SubElement(q, 'V')

                    # Special-case empty strings as described in from_xml().
                    if value_i == "":
                        continue

                    data_format = cls.data_formats[type(value_i)]
                    v.attrib[data_format] = str(value_i)

    class DirectionalityAttrib(Xml.EnumAttrib):
        value_from_str = {
                '0': 'none',
                '1': 'forward',
                '2': 'backward',
                '3': 'bidirectional',
        }

    class CleavageArrowsAttrib:

        @staticmethod
        def from_str(str):
            return [int(x) for x in str.split(',')]

        @staticmethod
        def to_str(value):
            return ','.join(str(x) for x in value)

    xml_tag = 'Feature'
    xml_subtag_defs = [
            ('segments', 'Segment', SegmentTag, []),
            ('qualifiers', 'Q', QualifierTag, {}),
    ]
    xml_attrib_defs = [
            ('id', 'recentID', Xml.IntAttrib),
            ('name', 'name', Xml.TextAttrib),
            ('type', 'type', Xml.TextAttrib),
            ('directionality', 'directionality', DirectionalityAttrib),
            ('reading_frame', 'readingFrame', Xml.IntAttrib),
            ('cleavage_arrows', 'cleavageArrows', CleavageArrowsAttrib),
            ('allow_segment_overlaps', 'allowSegmentOverlaps', Xml.BoolAttrib),
            ('swapped_segment_numbering', 'swappedSegmentNumbering', Xml.BoolAttrib),
            ('max_run_on', 'maxRunOn', Xml.IntAttrib),
            ('max_fused_run_on', 'maxFusedRunOn', Xml.IntAttrib),
            ('detection_mode', 'detectionMode', Xml.TextAttrib),

            # Translation-related attributes:
            ('genetic_code_id', 'geneticCode', Xml.TextAttrib),
            ('first_codon_met', 'translateFirstCodonAsMet', Xml.BoolAttrib),
            ('consecutive_translation_numbering', 'consecutiveTranslationNumbering', Xml.BoolAttrib),
            ('consecutive_numbering_start', 'consecutiveNumberingStartsFrom', Xml.IntAttrib),
            ('translated_mw', 'translationMW', Xml.FloatAttrib),
            ('hits_stop_codon', 'hitsStopCodon', Xml.BoolAttrib),
    ]

    @classmethod
    def from_segment(cls, **kwargs):
        """
        Instantiate a feature with a single segment.

        Keyword arguments corresponding to any attribute of either the Feature 
        or FeatureSegment classes are accepted.
        """
        feature_kwargs = {
                k: v
                for k, v in kwargs.items() 
                if k in cls._defined_names
        }
        segment_kwargs = {
                k: v
                for k, v in kwargs.items() 
                if k not in cls._defined_names
        }

        feature = cls(**feature_kwargs)
        feature.segments = [FeatureSegment(**segment_kwargs)]
        return feature

    def get_segment(self):
        """
        If this feature has only one segment, return it.  Otherwise, raise an 
        exception.
        """
        return one(self.segments)

    def set_segment(self, segment):
        """
        Remove any segments associated with this feature, and replace them with 
        the given segment.

        After calling this setter (or assigning to this attribute), the feature 
        will have exactly one segment.  It will still be possible to add or 
        remove segments later, though.
        """
        self.segments = [segment]

    def get_range(self):
        return self.begin, self.end

    def get_begin(self):
        return min(x.range[0] for x in self.segments)

    def get_end(self):
        return max(x.range[1] for x in self.segments)

    def get_sequence(self, full_seq):
        # It'd be nice to get rid of the *full_seq* argument, but that would 
        # require this block to have access to the top-level `SnapGene` object.  
        # This might be ok, but it's a bigger change than I want to make right 
        # now.
        start, end = self.range
        seq = full_seq[start:end]

        if self.directionality == 'backward':
            return reverse_complement(seq)
        else:
            return seq

@autoprop
class FeatureSegment(Xml, Repr):
    repr_attrs = 'type', 'color', 'range'

    class RangeAttrib:

        @staticmethod
        def from_str(str):
            # Make the range indices compatible with the conventions for python 
            # slicing.
            i, j = tuple(int(x) for x in str.split('-'))
            return i - 1, j

        @staticmethod
        def to_str(value):
            i, j = value
            return f'{i+1}-{j}'

    xml_tag = 'Segment'
    xml_attrib_defs = [
            ('name', 'name', Xml.TextAttrib),
            ('range', 'range', RangeAttrib),
            ('display', 'type', Xml.TextAttrib),
            ('color', 'color', Xml.TextAttrib),
            ('is_translated', 'translated', Xml.BoolAttrib),
            ('translation_start_number', 'translationNumberingStartsFrom', Xml.IntAttrib),
    ]

    def get_begin(self):
        return self.range[0]

    def get_end(self):
        return self.range[1]


