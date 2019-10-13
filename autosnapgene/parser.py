#!/usr/bin/env python3

import struct
import arrow
import xml.etree.ElementTree as etree

from pathlib import Path
from copy import deepcopy
from .errors import *

def blocks_from_file(path, block_classes=None):
    try:
        bytes = Path(path).read_bytes()
        return blocks_from_bytes(bytes)

    except ParseError as e:
        e.path = path
        raise e from None

def blocks_from_bytes(bytes, block_classes=None):
    i = 0
    blocks = []

    if block_classes is None:
        # Make sure all of the subclasses have been loaded.
        from . import blocks as _
        block_classes = Block.block_classes

    while i < len(bytes):
        j = i + 5
        if len(bytes) < j:
            raise ParseError("unexpected EOF")

        id, size = struct.unpack('>BI', bytes[i:j])
        if len(bytes) < j + size:
            raise ParseError("unexpected EOF")

        content = bytes[j:j+size]
        cls = block_classes.get(id, UndocumentedBlock)
        block = cls.from_bytes(content)
        block.block_id = id
        block.bytes = content
        blocks.append(block)
        i = j + size

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

class Repr:

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.__repr_attrs__()}".strip() + ">"

    def __repr_attrs__(self):
        if not self.repr_attrs:
            raise NotImplementedError

        return ' '.join(
                f'{k}="{self.__repr_attr__(k)}"'
                for k in self.repr_attrs
                if hasattr(self, k)
        )

    def __repr_attr__(self, attr):
        return getattr(self, attr)


class Xml:
    """
    A class that can read/write its attributes to/from XML.
    """
    xml_tag = None
    xml_subtag_defs = []
    xml_attrib_defs = []

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


    class AppendListTag:
        # For when the same tag may appear multiple times, and each appearance 
        # should add the value to a growing list.

        @staticmethod
        def setattr(obj, name, value):
            if not hasattr(obj, name):
                setattr(obj, name, [value])
            else:
                getattr(obj, name).append(value)

    class UpdateDictTag:

        @staticmethod
        def setattr(obj, name, value):
            if not hasattr(obj, name):
                setattr(obj, name, value)
            else:
                getattr(obj, name).update(value)



    class TextAttrib:

        @staticmethod
        def from_str(str):
            return str

        @staticmethod
        def to_str(value):
            return value

    class BoolAttrib:

        @staticmethod
        def from_str(str):
            return {'0': False, '1': True}[str]

        @staticmethod
        def to_str(value):
            return str(int(value))

    class EnumAttrib:
        value_from_str = {}

        def __init_subclass__(cls):
            cls.str_from_value = {v: k for k, v in cls.value_from_str.items()}

        @classmethod
        def from_str(cls, str):
            return cls.value_from_str[str]

        @classmethod
        def to_str(cls, value):
            return cls.str_from_value[value]

    class IntAttrib:

        @staticmethod
        def from_str(str):
            return int(str)

        @staticmethod
        def to_str(value):
            return str(value)

    class FloatAttrib:

        @staticmethod
        def from_str(str):
            return float(str)

        @staticmethod
        def to_str(value):
            return str(value)

    def __init__(self, **kwargs):
        for name, value in self._defaults.items():
            setattr(self, name, deepcopy(value))

        for name, value in kwargs.items():
            if name in self._defined_names:
                setattr(self, name, value)
            else:
                did_you_mean = '\n    '.join(self._defined_names)
                raise AttributeError(f"'{name}' is not a valid attribute of {self.__class__.__name__}, did you mean:\n    {did_you_mean}")

        # Keep track of unexpected attributes and subtags, so that we can write 
        # everything we read, even if we don't understand it all.
        self._unparsed_attribs = {}
        self._unparsed_subtags = []

    def __init_subclass__(cls):
        super().__init_subclass__()

        def find_dups(xs):
            seen = set()
            dups = set()

            for x in xs:
                if x in seen:
                    dups.add(x)
                seen.add(x)

            return dups

        def check_dups(xs):
            dups = find_dups(xs)
            if dups:
                dups_str = '\n    '.join(dups)
                raise ValueError("The following attributes are defined more than once:\n    {dups_str}")

        cls._defaults = {
                defs[0]: defs[3]
                for defs in cls.xml_attrib_defs + cls.xml_subtag_defs
                if len(defs) == 4
        }

        cls._attrib_parsers_by_name = {
                name: (attrib, parser)
                for name, attrib, parser, *_ in cls.xml_attrib_defs
        }
        cls._attrib_parsers_by_attrib = {
                attrib: (name, parser)
                for name, attrib, parser, *_ in cls.xml_attrib_defs
        }
        cls._subtag_parsers_by_name = {
                name: (tag, parser)
                for name, tag, parser, *_ in cls.xml_subtag_defs
        }
        cls._subtag_parsers_by_tag = {
                tag: (name, parser)
                for name, tag, parser, *_ in cls.xml_subtag_defs
        }

        cls._attrib_names = [
                name for name, *_ in cls.xml_attrib_defs
        ]
        cls._subtag_names = [
                name for name, *_ in cls.xml_subtag_defs
        ]
        cls._defined_names = cls._attrib_names + cls._subtag_names

        check_dups(cls._defined_names)

    def __getattr__(self, name):
        if name in self._defined_names:
            raise AttributeError(f"'{name}' not defined for {self.__class__.__name__}.")

        else:
            did_you_mean = '\n    '.join(self._defined_names)
            raise AttributeError(f"'{name}' is not a valid attribute of {self.__class__.__name__}, did you mean:\n    {did_you_mean}")

    def __delattr__(self, name):
        try:
            super().__delattr__(name)

        except AttributeError:
            if name in self._defined_names:
                raise AttributeError(f"'{name}' not defined for {self.__class__.__name__}.")

            else:
                did_you_mean = '\n    '.join(self._defined_names)
                raise AttributeError(f"'{name}' is not a valid attribute of {self.__class__.__name__}, did you mean:\n    {did_you_mean}")

    def __eq__(self, other):
        undef = object()
        return all([
            getattr(self, x, undef) == getattr(other, x, undef)
            for x in self._defined_names
        ])


    @classmethod
    def from_bytes(cls, bytes):
        xml = bytes.decode('utf8')
        root = etree.fromstring(xml)
        return cls.from_xml(root)

    @classmethod
    def from_xml(cls, root):
        if not cls.xml_tag:
            raise NotImplementedError(f"'{cls.__qualname__}.xml_tag' not defined.")
        if root.tag != cls.xml_tag:
            raise ValueError(f"expected <{self.xml_tag}>, but got {root}")

        self = cls()

        for attrib in root.attrib:
            try:
                name, parser = cls._attrib_parsers_by_attrib[attrib]
            except KeyError:
                self._unparsed_attribs[attrib] = root.attrib[attrib]
            else:
                value = parser.from_str(root.attrib[attrib])
                setattr(self, name, value)

        for element in root:
            try:
                name, parser = cls._subtag_parsers_by_tag[element.tag]
            except KeyError:
                self._unparsed_subtags.append(element)
            else:
                value = parser.from_xml(element)
                getattr(parser, 'setattr', setattr)(self, name, value)

        return self

    def to_bytes(self):
        root = self.to_xml()
        return etree.tostring(root)

    def to_xml(self):
        from inspect import signature

        if not self.xml_tag:
            raise NotImplementedError("'{self.__class__.__qualname__}.xml_tag' not defined.")

        root = etree.Element(self.xml_tag)

        def has_default_value(name):
            # In some cases, the absence of an attribute/subtag implies some 
            # default value.  This function is used to avoid writing 
            # attributes/subtags with such default values.
            if name in self._defaults:
                return getattr(self, name) == self._defaults[name]
            else:
                return False

        for attrib, value in self._unparsed_attribs.items():
            root.attrib[attrib] = value

        for name in self._attrib_names:
            if not hasattr(self, name): continue
            if has_default_value(name): continue
            attrib, parser = self._attrib_parsers_by_name[name]
            value = getattr(self, name)
            root.attrib[attrib] = parser.to_str(value)

        for element in self._unparsed_subtags:
            root.append(element)

        for name in self._subtag_names:
            if not hasattr(self, name): continue
            if has_default_value(name): continue
            tag, parser = self._subtag_parsers_by_name[name]
            sig = signature(parser.to_xml)

            # For most parsers, it's convenient if we take care of making the 
            # element.  But a few of the more complex parsers need to customize 
            # this process.  So we basically overload the 'to_xml()' method and 
            # inspect the signature to see which behavior the parser wants.
            # 
            # Note that we could've gotten similar behavior by having a 
            # superclass method that creates the element and calls a 
            # overload-able method to customize it.  But for aesthetic reasons, 
            # I don't want to require the parsers to inherit from anything.

            if len(sig.parameters) == 2:
                element = etree.SubElement(root, tag)
                parser.to_xml(element, getattr(self, name))
            else:
                parser.to_xml(root, tag, getattr(self, name))

        return root

class Block(Repr):
    block_classes = {}

    def __init_subclass__(cls):
        super().__init_subclass__()

        if hasattr(cls, 'block_id'):
            Block.block_classes[cls.block_id] = cls

        if not hasattr(cls, 'repr_attrs'):
            cls.repr_attrs = ['block_id']
        else:
            if isinstance(cls.repr_attrs, str):
                raise ValueError(f"{cls.__qualname__}.repr_attrs should be a list, not a str")
            cls.repr_attrs = ['block_id', *[
                    x for x in cls.repr_attrs
                    if x != 'block_id'
            ]]

    @classmethod
    def from_bytes(cls, bytes):
        raise NotImplementedError

    def to_bytes(self):
        raise NotImplementedError


class UnparsedBlock(Block):
    repr_attrs = ['bytes']

    def __repr_attr__(self, attr):
        if attr == 'bytes':
            bytes_repr = repr(self.bytes)[2:-1]
            if len(bytes_repr) > 32:
                bytes_repr = bytes_repr[:32] + '...'
            return bytes_repr

        else:
            return super().__repr_attr__(attr)

    @classmethod
    def from_bytes(cls, bytes):
        return cls()

    def to_bytes(self):
        return self.bytes


class UndocumentedBlock(UnparsedBlock):
    pass



