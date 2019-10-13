*********************
Working with features
*********************

A feature is a sequence annotated with information like a name, a function, a 
description, etc.  Maintaining a complete set of features is critical to making 
a plasmid map understandable.  AutoSnapGene can help automate this maintenance, 
giving you more confidence that all of your plasmid maps have the most 
up-to-date annotations.

Command-line
============
Add a feature to a sequence::

   $ autosnapgene feature add puc19.dna "lac operator" TTGTTATCCGCTCACAA

Remove a feature from a sequence::

   $ autosnapgene feature remove puc19.dna "lac operator"

Remove all features from a sequence::

   $ autosnapgene feature clear puc19.dna

.. note::

   I intend to extend this API at some point, so that it is possible to apply 
   an external list of features (e.g. in an Excel file) to a large number of 
   SnapGene files.  If you would find this feature (or any other feature) 
   useful, make a `bug report 
   <https://github.com/kalekundert/autosnapgene/issues>`_ or `pull request 
   <https://github.com/kalekundert/autosnapgene/pulls>`_ to encourage me to 
   work on it faster!

Python API
==========
To add a feature to a sequence, you first need to create and fill in an 
`autosnapgene.Feature` object::

   >>> import autosnapgene as snap
   >>> puc = snap.parse('puc19.dna')
   >>>
   >>> feat = snap.Feature()
   >>> feat.name = "lac operator"
   >>> feat.type = 'protein_bind'

A description of the feature can be added as a "note" qualifier.  Qualifiers 
are a concept from `GenBank files 
<http://www.insdc.org/files/feature_table.html#3.3>`_ that are copied in 
SnapGene, so read the link for more information on the many other qualifiers 
available for you to specify.  Also note that the "note" qualifier, like many 
text fields in SnapGene, is HTML-formatted to allow for simple markup.  Unicode 
characters (e.g. β) are also supported::

   >>> feat.qualifiers = {
   ...         'note': """\
   ... <html><body>The <i>lac</i> repressor binds to the <i>lac</i> operator to
   ... inhibit transcription in <i>E. coli</i>. This inhibition can be relieved
   ... by adding lactose or isopropyl-β-D-thiogalactopyranoside
   ... (IPTG).</body></html>""",
   ... }

Features must also contain one or more `autosnapgene.FeatureSegment` objects, 
so it's important to consider what information goes in the feature and what 
goes in the segments.  Information in the feature applies to the whole feature.  
This includes things like the name and type of the feature, what direction (if 
any) the feature goes in, how the feature should be translated, etc.  
Information in the segments applies only to that segment.  This includes 
another name and a color, among other things::

   >>> # In the common case that there is only one segment, access it using
   >>> # the `feat.segment` shorthand.  In the more general case where there
   >>> # may be multiple segments, use `feat.segments` to work with a list of
   >>> # segments instead.
   >>>
   >>> feat.segment = snap.FeatureSegment()
   >>> feat.segment.color = '#31849b'

The snippets above show how to instantiate a feature from scratch, but the same 
thing could also be done more succinctly using the 
`autosnapgene.Feature.from_segment` factory function.  This function can only 
make single-segment features, but it takes care of properly delegating each 
argument to either the feature or the segment::

   >>> feat = snap.Feature.from_segment(
   ...         name="lac operator",
   ...         type='protein_bind',
   ...         qualifiers={
   ...            'note': """\
   ... <html><body>The <i>lac</i> repressor binds to the <i>lac</i> operator to
   ... inhibit transcription in <i>E. coli</i>. This inhibition can be relieved
   ... by adding lactose or isopropyl-β-D-thiogalactopyranoside
   ... (IPTG).</body></html>""",
   ...         },
   ...         color='#31849b',
   ... )

At this point, we still have not specified where the feature will go.  The 
indices of the feature (relative to the start of the sequence) are specified in 
the segments.  It may be counter-intuitive that the feature itself has no 
position information: the segments are totally responsible for identifying 
where the feature is located in the plasmid.  It may also be counter-intuitive 
that neither the feature nor its segments contain any sequence information.  
Instead, the sequence of the feature is inferred from the location of the 
feature and the sequence of the whole file.

Since it's more natural to think about features as being associated with 
particular sequences, the `SnapGene.add_feature()` method takes a sequence as 
an argument and uses it to set the position of the segments in the given 
feature.  If the given sequence occurs multiple times, it will add a copy of 
the given feature for each occurrence::

   >>> puc.add_feature(feat, 'TTGTTATCCGCTCACAA')
   >>> puc.write()
   
You can also get information about the features in an existing file::

   >>> pprint(puc.features)
   [<Feature name="lac operator" type="protein_bind">,
    <Feature name="M13 rev" type="primer_bind">,
    <Feature name="M13 fwd" type="primer_bind">,
    <Feature name="lac promoter" type="promoter">,
    <Feature name="MCS" type="misc_feature">,
    <Feature name="AmpR promoter" type="promoter">,
    <Feature name="lacZα" type="CDS">,
    <Feature name="ori" type="rep_origin">,
    <Feature name="AmpR" type="CDS">]

::

   >>> lac = puc.get_feature("lac promoter")
   >>> lac.type = 'promoter'
   >>> lac.qualifiers
   {'note': '<html><body>promoter for the <i>E. coli lac</i> 
   operon</body></html>'}
   >>> pprint(lac.segments)
   [<FeatureSegment color="#ffffff" range="(513, 519)">,
    <FeatureSegment color="#ffffff" range="(520, 537)">,
    <FeatureSegment color="#ffffff" range="(538, 543)">]

It's informative to look at some of the attributes of a translated feature::

   >>> ampr = puc.get_feature("AmpR")
   >>> ampr.type
   'CDS'
   >>> ampr.directionality
   'backward'
   >>> ampr.reading_frame
   -1
   >>> ampr.cleavage_arrows
   [2417]
   >>> ampr.consecutive_translation_numbering
   True
   >>> ampr.translated_mw
   31558.16
   >>> ampr.hits_stop_codon
   True
   >>> pprint(ampr.segments)
   [<FeatureSegment color="#ccffcc" range="(1626, 2417)">,
   <FeatureSegment color="#ccffcc" range="(2418, 2486)">]
   >>> ampr.segments[0].display
   'standard'
   >>> ampr.segments[0].is_translated
   True

Remove a feature from a sequence::

   >>> # Primers shouldn't be features...
   >>> puc.remove_feature("M13 fwd")
   >>> puc.write()

Remove all features from a sequence::

   >>> puc.clear_features()
   >>> puc.write()

