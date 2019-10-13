#!/usr/bin/env python3

import pytest
import autosnapgene as snap
import xml.etree.ElementTree as etree
from arrow import get as datetime
from pprint import pprint

Xml = snap.parser.Xml

class DummyXml(Xml):
    xml_tag = 'Dummy'

    class MyEnumAttrib(Xml.EnumAttrib):
        value_from_str = {
                'left': 'LEFT',
                'right': 'RIGHT',
        }

    class MyAppendListTag(Xml.AppendListTag):

        @staticmethod
        def from_xml(element):
            return element.text

        @staticmethod
        def to_xml(parent, tag, values):
            for value in values:
                e = etree.SubElement(parent, tag)
                e.text = value

    class MyUpdateDictTag(Xml.UpdateDictTag):

        @staticmethod
        def from_xml(element):
            k, v = element.text.split(':')
            return {k: int(v)}

        @staticmethod
        def to_xml(parent, tag, values):
            for k, v in values.items():
                e = etree.SubElement(parent, tag)
                e.text = f'{k}:{v}'

    xml_attrib_defs = [
            ('text_attrib', 'text', Xml.TextAttrib),
            ('bool_attrib', 'bool', Xml.BoolAttrib),
            ('enum_attrib', 'enum', MyEnumAttrib),
            ('int_attrib', 'int', Xml.IntAttrib),
            ('float_attrib', 'float', Xml.FloatAttrib),
            ('default_attrib', 'default', Xml.TextAttrib, ''),
    ]
    xml_subtag_defs = [
            ('text_tag', 'Text', Xml.TextTag),
            ('bool_tag', 'Bool', Xml.BoolTag),
            ('date_tag', 'Date', Xml.DateTag),
            ('list_tag', 'List', MyAppendListTag),
            ('dict_tag', 'Dict', MyUpdateDictTag),
            ('default_tag', 'Default', Xml.TextTag, ''),
    ]

@pytest.mark.parametrize(
        'bytes, expected', [(
                b'', [
            ]), (
                b'\x01\x00\x00\x00\x00', [
                    (1, b''),
            ]), (
                b'\x01\x00\x00\x00\x05Hello', [
                    (1, b'Hello'),
            ]), (
                b'\x01\x00\x00\x00\x00\x02\x00\x00\x00\x00', [
                    (1, b''),
                    (2, b''),
            ]), (
                b'\x01\x00\x00\x00\x05Hello\x02\x00\x00\x00\x06world!', [
                    (1, b'Hello'),
                    (2, b'world!'),
            ])
])
def test_blocks_from_bytes(bytes, expected):
    blocks = snap.parser.blocks_from_bytes(bytes, {})

    assert len(blocks) == len(expected)

    for block, (id, bytes) in zip(blocks, expected):
        assert block.block_id == id
        assert block.bytes == bytes

@pytest.mark.parametrize(
        'bytes', [
            # Not enough bytes for header.
            b'\x01',
            b'\x01\x00',
            b'\x01\x00\x00',
            b'\x01\x00\x00\x00',

            # Not enough bytes for content.
            b'\x01\x00\x00\x00\x01',
])
def test_blocks_from_bytes_errors(bytes):
    with pytest.raises(snap.ParseError):
        snap.parser.blocks_from_bytes(bytes, {})

def test_blocks_from_file(examples):
    blocks = snap.parser.blocks_from_file(examples / 't7_promoter.dna')
    pprint(blocks)
    assert len(blocks) == 10

@pytest.mark.parametrize(
        'xml, raw_expected', [(
                b'<Dummy />', {
            }),

            # Attributes
            (
                b'<Dummy text="hello" />', {
                    'text_attrib': 'hello',
            }), (
                b'<Dummy bool="1" />', {
                    'bool_attrib': True,
            }), (
                b'<Dummy enum="left" />', {
                    'enum_attrib': 'LEFT',
            }), (
                b'<Dummy int="13" />', {
                    'int_attrib': 13,
            }), (
                b'<Dummy float="3.14" />', {
                    'float_attrib': pytest.approx(3.14),
            }), (
                b'<Dummy default="set" />', {
                    'default_attrib': 'set',
            }), (
                # Can't parse, but should still write correctly.
                b'<Dummy unknown="!@#$" />', {
            }),

            # Sub tags
            (
                b'<Dummy><Text>Hello</Text></Dummy>', {
                    'text_tag': 'Hello',
            }), (
                b'<Dummy><Bool>1</Bool></Dummy>', {
                    'bool_tag': True,
            }), (
                b'<Dummy><Date>2019.10.13</Date></Dummy>', {
                    'date_tag': datetime(2019, 10, 13),
            }), (
                b'<Dummy><List>A</List></Dummy>', {
                    'list_tag': ['A'],
            }), (
                b'<Dummy><List>A</List><List>B</List></Dummy>', {
                    'list_tag': ['A', 'B'],
            }), (
                b'<Dummy><Dict>A:1</Dict></Dummy>', {
                    'dict_tag': {'A': 1},
            }), (
                b'<Dummy><Dict>A:1</Dict><Dict>B:2</Dict></Dummy>', {
                    'dict_tag': {'A': 1, 'B': 2},
            }), (
                b'<Dummy><Default>Set</Default></Dummy>', {
                    'default_tag': 'Set',
            }), (
                # Can't parse, but should still write correctly.
                b'<Dummy><Unknown>!@#$</Unknown></Dummy>', {
            })
])
def test_xml_parse_and_write(xml, raw_expected):
    expected = {
            'default_attrib': '',
            'default_tag': '',
            **raw_expected
    }

    x = DummyXml.from_bytes(xml)
    for attr, value in expected.items():
        assert getattr(x, attr) == value

    unset_attrs = [x for x in x._defined_names if x not in expected]
    for attr in unset_attrs:
        with pytest.raises(AttributeError):
            getattr(x, attr)

    assert x.to_bytes() == xml


