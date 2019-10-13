************************
Working with the history
************************

Although SnapGene's history feature is often quite useful, it can also cause 
your files to grow ridiculously big.  You can use AutoSnapGene to clear history 
in files that don't need it.

.. note::

   Support for history is currently quite rudimentary, in no small part because  
   history is probably the most complicated aspect of a SnapGene file to parse.  
   If you would like more control over the history, submit a `bug report
   <https://github.com/kalekundert/autosnapgene/issues>`_ and I'll try to help 
   you out.

Command-line
============
Clear the history from a SnapGene file::

   $ autosnapgene history clear pKBK076.dna

Python API
==========
Clear the history from a SnapGene file::

   >>> import autosnapgene as snap
   >>> dna = snap.parse('pKBK076.dna')
   >>> dna.clear_history()

