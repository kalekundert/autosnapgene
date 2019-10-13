**********************
Working with sequences
**********************

The most fundamental part of a SnapGene file is the sequence.  Below are the 
utilities for extracting and modifying that information.

Command-line
============
Get the sequence from a SnapGene file::

   $ autosnapgene seq get t7_promoter.dna
   TAATACGACTCACTATAGG

Change the sequence of a SnapGene file::

   $ autosnapgene seq set t7_promoter.dna TAATACGACTCACTATAGG

Make a sequence uppercase or lowercase::

   $ autosnapgene seq upper t7_promoter.dna
   $ autosnapgene seq lower t7_promoter.dna

Python API
==========
Below are some examples of how the python API can be used to get 
sequence-related information.  See the :doc:`api` page for more details.

Get the sequence from a SnapGene file::

   >>> import autosnapgene as snap
   >>> t7 = snap.parse('t7_promoter.dna')
   >>> t7.sequence
   'TAATACGACTCACTATAGG'

Get various other properties of the sequence::

   >>> t7.topology
   'linear'
   >>> t7.strandedness
   'double'
   >>> t7.is_dam_methylated
   False
   >>> t7.is_dcm_methylated
   False
   >>> t7.is_ecoki_methylated
   False

All of the above attributes can be manipulated however you like.  For example, 
we can add another "G" to the end of the sequence::

   >>> t7.sequence += 'G'
   >>> t7.sequence
   TAATACGACTCACTATAGGG

You can also make new sequences from scratch::

   tac = snap.SnapGene()
   tac.sequence = 'TTGACAATTAATCATCGGCTCGTATAATG'

Save the changes::

   >>> t7.write()                     # Overwrite the original file.
   >>> tac.write('tac_promoter.dna')  # Make a new file.
