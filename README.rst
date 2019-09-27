AutoSnapGene
============

.. image:: https://img.shields.io/pypi/v/snapgene_tools.svg
   :target: https://pypi.python.org/pypi/snapgene_tools

.. image:: https://img.shields.io/pypi/pyversions/snapgene_tools.svg
   :target: https://pypi.python.org/pypi/snapgene_tools

.. image:: https://img.shields.io/travis/kalekundert/snapgene_tools.svg
   :target: https://travis-ci.org/kalekundert/snapgene_tools

.. image:: https://img.shields.io/coveralls/kalekundert/snapgene_tools.svg
   :target: https://coveralls.io/github/kalekundert/snapgene_tools?branch=master

AutoSnapGene is a set of tools for automatically processing SnapGene files.  It includes a parser capable of reading and writing SnapGene files, in addition to a set of command-line utilities that implement common operations.

Installation
------------
AutoSnapGene is not yet available on PyPI.  In the meantime::

   $ git clone git@github.com:kalekundert/autosnapgene.git
   $ cd autosnapgene
   $ flit install
   
Usage
-----
Use --help for brief descriptions of the available commands::

   $ autosnapgene --help
   Perform batch operations on SnapGene files.

   Usage:
       autosnapgene add_traces <dna_path> <ab1_paths>...
       autosnapgene clear_history <dna_path>
