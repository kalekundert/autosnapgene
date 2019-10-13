************************
Working with annotations
************************

SnapGene files contain a lot of metadata describing the sequence in question.  
You can access all of this information using AutoSnapGene.

Command-line
============
The command-line interface for interacting with annotations has not been 
implemented yet.

Python API
==========
Below are some examples of how the python API can be used to get sequence 
annotations. See the :doc:`api` page for more details.

The `pUC19 plasmid map`__ available for download from the SnapGene website is 
well-annotated::

   >>> import autosnapgene as snap
   >>> puc = snap.parse('puc19.dna')
   >>> puc.author
   'New England Biolabs'
   >>> puc.organism
   'Escherichia coli'
   >>> puc.plasmid_type
   'Synthetic'
   >>> puc.is_confirmed_experimentally
   False

__ https://www.snapgene.com/resources/plasmid-files/?set=basic_cloning_vectors&plasmid=pUC19

Get a description of the plasmid.  Note that the description (like many of the 
text-based metadata fields) is an HTML string::

   >>> puc.description
   '<html><body>Standard <i>E. coli</i> vector with a multiple cloning site (MCS) for DNA cloning. The MCS is reversed in pUC18.</body></html>'
   >>> puc.comments
   'See also GenBank accession L09137.'

Get modification dates for the plasmid file.  AutoSnapGene uses `arrow.Arrow 
<arrow.arrow.Arrow>` to represent dates.  These objects mimic the 
`datetime.datetime` API, but are more convenient in a number of ways::

   >>> puc.date_created
   <Arrow [2012-03-08T00:00:00+00:00]>
   >>> puc.date_last_modified
   <Arrow [2013-11-25T00:00:00+00:00]>

Get any references associated with the plasmid.  Detailed information is 
available for each reference::

   >>> puc.references
   [<Reference pubmed_id="2985470">]
   >>> puc.references[0].title
   '<html><body>Improved M13 phage cloning vectors and host strains: nucleotide sequences of the M13mp18 and pUC19 vectors.</body></html>'
   >>> puc.references[0].pubmed_id
   '2985470'
   >>> puc.references[0].journal
   '<html><body>Gene 1985;33:103-19.</body></html>'
   >>> puc.references[0].authors
   '<html><body>Yanisch-Perron C, Vieira J, Messing J.</body></html>'

Any of the above attributes can be modified::

   >>> puc.is_confirmed_experimentally = True
   >>> puc.write()
