************
AutoSnapGene
************
.. image:: https://img.shields.io/pypi/v/autosnapgene.svg
   :target: https://pypi.python.org/pypi/autosnapgene

.. image:: https://img.shields.io/pypi/pyversions/autosnapgene.svg
   :target: https://pypi.python.org/pypi/autosnapgene

.. image:: https://img.shields.io/travis/kalekundert/autosnapgene.svg
   :target: https://travis-ci.org/kalekundert/autosnapgene

.. image:: https://img.shields.io/coveralls/kalekundert/autosnapgene.svg
   :target: https://coveralls.io/github/kalekundert/autosnapgene?branch=master

.. image:: https://readthedocs.org/projects/autosnapgene/badge/?version=latest
   :target: https://autosnapgene.readthedocs.io/en/latest/?badge=latest

`SnapGene <https://www.snapgene.com/>`_ is an excellent tool for viewing and 
editing plasmid maps, but it's difficult to programmatically work with the 
binary files (``*.dna``) that SnapGene reads and writes.  This is what 
AutoSnapGene helps with.  AutoSnapGene is both a python library for parsing 
SnapGene files, and a set of command-line tools for automating simple tasks.  
You can use these tools to:

- Automatically add sequencing data to your plasmids.

- Synchronize features and primers over all of your plasmids at once.

- Seamlessly use the sequences stored in your SnapGene file in other scripts.

- Programmatically build sequences from combinations of smaller parts.

