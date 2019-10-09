#!/usr/bin/env python3

import pytest
import autosnapgene as snap

def test_getters_t7(parse_and_write):
    for dna in parse_and_write('t7_promoter.dna'):
        assert len(dna.features) == 1

        feat = dna.features[0]
        assert feat.id == 0
        assert feat.name == "T7 promoter"
        assert feat.type == 'promoter'
        assert feat.directionality == 'forward'
        assert feat.allow_segment_overlaps == False
        assert feat.consecutive_translation_numbering == True
        assert feat.qualifiers == {
                'note': "<html><body>promoter for bacteriophage T7 RNA polymerase</body></html>",
        }
        assert len(feat.segments) == 1

        seg = feat.segments[0]
        assert seg.begin == 1
        assert seg.end == 19
        assert seg.range == (1, 19)
        assert seg.color == '#ffffff'
        assert seg.display == 'standard'

def test_getters_flag(parse_and_write):
    for dna in parse_and_write('flag_tag.dna'):
        assert len(dna.features) == 1

        feat = dna.features[0]
        assert feat.id == 0
        assert feat.name == "FLAG"
        assert feat.type == 'CDS'
        assert feat.directionality == 'forward'
        assert feat.translated_mw == pytest.approx(1012.98)
        assert feat.allow_segment_overlaps == False
        assert feat.consecutive_translation_numbering == True
        assert feat.detection_mode == 'exactProteinMatch'
        assert feat.qualifiers == {
                'codon_start': 1,
                'product': "<html><body>FLAG® epitope tag, followed by an enterokinase cleavage site</body></html>",
                'transl_table': 1,
                'translation': "DYKDDDDK",
        }
        assert len(feat.segments) == 1

        seg = feat.segments[0]
        assert seg.begin == 1
        assert seg.end == 24
        assert seg.range == (1, 24)
        assert seg.color == '#cc99b2'
        assert seg.display == 'standard'
        assert seg.is_translated == True

xml_t7 = '''\
<Features nextValidID="1">
 <Feature recentID="0" name="T7 promoter" directionality="1" type="promoter" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1">
  <Segment range="1-19" color="#ffffff" type="standard"/>
  <Q name="note">
   <V text="&lt;html&gt;&lt;body&gt;promoter for bacteriophage T7 RNA polymerase&lt;/body&gt;&lt;/html&gt;"/>
  </Q>
 </Feature>
</Features>
'''
xml_flag = '''\
<Features nextValidID="1">
 <Feature recentID="0" name="FLAG" directionality="1" translationMW="1012.98" type="CDS" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1" detectionMode="exactProteinMatch">
  <Segment range="1-24" color="#cc99b2" type="standard" translated="1"/>
  <Q name="codon_start">
   <V int="1"/>
  </Q>
  <Q name="product">
   <V text="&lt;html&gt;&lt;body&gt;FLAG\xc2\xae epitope tag, followed by an enterokinase cleavage site&lt;/body&gt;&lt;/html&gt;"/>
  </Q>
  <Q name="transl_table">
   <V int="1"/>
  </Q>
  <Q name="translation">
   <V text="DYKDDDDK"/>
  </Q>
 </Feature>
</Features>
'''
