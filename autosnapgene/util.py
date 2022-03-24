#!/usr/bin/env python3

def reverse(sequence):
    return sequence[::-1]

def complement(sequence):
    complements = str.maketrans('ACTGactg', 'TGACtgac')
    return sequence.translate(complements)

def reverse_complement(sequence):
    return reverse(complement(sequence))


