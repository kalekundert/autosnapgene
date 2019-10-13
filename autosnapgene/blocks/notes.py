#!/usr/bin/env python3

from ..parser import Block, Xml, Repr

class NotesBlock(Xml, Block):
    block_id = 6
    repr_attrs = 'type', 'author', 'last_modified'

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

    xml_tag = 'Notes'
    xml_subtag_defs = [
            ('uuid', 'UUID', Xml.TextTag),
            ('type', 'Type', Xml.TextTag),
            ('author', 'CreatedBy', Xml.TextTag),
            ('description', 'Description', Xml.HtmlTag),
            ('comments', 'Comments', Xml.HtmlTag),
            ('references', 'References', ReferencesTag),
            ('transformed_into', 'TransformedInto', Xml.TextTag),
            ('is_confirmed_experimentally', 'ConfirmedExperimentally', Xml.BoolTag),
            ('organism', 'Organism', Xml.TextTag),
            ('accession_number', 'AccessionNumber', Xml.TextTag),
            ('code_number', 'CodeNumber', Xml.TextTag),
            ('sequence_class', 'SequenceClass', Xml.TextTag),
            ('custom_map_label', 'CustomMapLabel', Xml.TextTag),
            ('use_custom_map_label', 'UseCustomMapLabel', Xml.BoolTag),
            ('date_created', 'Created', Xml.DateTag),
            ('date_last_modified', 'LastModified', Xml.DateTag),
    ]

class Reference(Xml, Repr):
    repr_attrs = 'pubmed_id',
    xml_tag = 'Reference'
    xml_attrib_defs = [
            ('title', 'title', Xml.TextAttrib),
            ('pubmed_id', 'pubMedID', Xml.TextAttrib),
            ('journal', 'journal', Xml.TextAttrib),
            ('authors', 'authors', Xml.TextAttrib),
    ]


