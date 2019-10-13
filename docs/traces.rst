******************************
Working with sequencing traces
******************************

SnapGene makes it easy to check sequencing data for your plasmids by having the 
ability to show that data aligned with the plasmid sequence.  You can make it 
even easier to check sequencing data by automatically adding sequencing data to 
your plasmids.

.. note::
   
   In order to add traces using AutoSnapGene, you must install Staden io_lib as 
   described on the :doc:`installation page <installation>`.

Command-line
============
Below are some examples of how the python API can be used to manipulate 
sequencing traces.  See the :doc:`api` page for more details.

Add a sequencing trace to a plasmid file::

   $ autosnapgene trace add pKBK076.dna 76{A,B,C}.ab1

Remove a sequencing trace from a file::

   $ autosnapgene trace add pKBK076.dna 76A

Remove all sequencing traces from a file::

   $ autosnapgene trace clear pKBK076.dna

Alphabetically sort the sequencing traces in a file::

   $ autosnapgene trace sort pKBK076.dna

Copy any sequencing traces into separate files.  The files will be in the ZTR 
format, which is the format used internally by SnapGene::

   $ mkdir ztr
   $ autosnapgene trace extract pKBK076.dna ztr

Python API
==========
Add a sequencing trace to a file.  Note that the name of the trace will be its 
file name without the extension::

   >>> import autosnapgene as snap
   >>> dna = snap.parse('pKBK076.dna')
   >>> dna.add_trace('76A.ab1')

Count the number of sequencing traces in a file::

   >>> dna.count_traces()
   1

Determine is a file contains a certain trace::

   >>> dna.has_trace('76A')
   True
   >>> dna.has_trace('76B')
   False

Remove a sequencing trace from a file::

   >>> dna.remove_trace('76A')

Remove all sequencing traces from a file::

   >>> dna.clear_traces()

Sort the sequencing traces in a file.  You can provide a key function to sort 
the traces by any criteria you want.  The key function should accept a single 
argument, which will be an AlignmentMetaData instance with attributes such as 
``name``, ``id``, ``is_visible``, ``is_trace``, and some others.  For example, 
the following snippet will sort the traces first be visibility, and second by 
name (alphabetically)::

   >>> dna.sort_traces(key=lambda x: (x.is_visible, x.name))
   $ autosnapgene trace sort pKBK076.dna

Copy any sequencing traces into separate files.  The files will be in the ZTR 
format, which is the format used internally by SnapGene::

   >>> dna.extract_traces('ztr')
