AutoSnapGene
============

.. image:: https://img.shields.io/pypi/v/autosnapgene.svg
   :target: https://pypi.python.org/pypi/autosnapgene

.. image:: https://img.shields.io/pypi/pyversions/autosnapgene.svg
   :target: https://pypi.python.org/pypi/autosnapgene

.. image:: https://img.shields.io/travis/kalekundert/autosnapgene.svg
   :target: https://travis-ci.org/kalekundert/autosnapgene

.. image:: https://img.shields.io/coveralls/kalekundert/autosnapgene.svg
   :target: https://coveralls.io/github/kalekundert/autosnapgene?branch=master

AutoSnapGene is a set of tools for automatically processing SnapGene files.  It includes a parser capable of reading and writing SnapGene files, in addition to a set of command-line utilities that implement common operations.

Installation
------------
AutoSnapGene is not yet available on PyPI.  In the meantime::

   $ git clone git@github.com:kalekundert/autosnapgene.git
   $ cd autosnapgene
   $ flit install
   
Usage
-----
Use ``--help`` for brief descriptions of the available commands::

   $ autosnapgene --help
   Perform batch operations on SnapGene files.

   Usage:
       autosnapgene add_traces <dna_path> <ab1_paths>...
       autosnapgene clear_history <dna_path>
