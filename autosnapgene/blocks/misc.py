#!/usr/bin/env python3

from ..parser import UnparsedBlock

class CompressedDnaBlock(UnparsedBlock):
    block_id = 1

class PropertiesBlock(UnparsedBlock):
    block_id = 8

class UracilBlock(UnparsedBlock):
    block_id = 19

class DnaColorBlock(UnparsedBlock):
    block_id = 20

