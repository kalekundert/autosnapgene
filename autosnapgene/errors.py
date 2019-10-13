#!/usr/bin/env python3

class ParseError(Exception):

    def __init__(self, message):
        self.message = message
        self.path = None

    def __str__(self):
        return self.message.format(path=self.path)

class BlockNotFound(AttributeError):
    pass

class SequenceNotFound(Exception):
    pass

class FeatureNotFound(Exception):
    pass
