#!/usr/bin/env python3

import autosnapgene as snap
from arrow import get as datetime

def test_from_bytes():
    bytes = b'''\
<Notes>
<UUID>d528ff2f-6579-48ab-9301-ab40df3f2505</UUID>
<Type>Natural</Type>
<CustomMapLabel>T7 promoter</CustomMapLabel>
<UseCustomMapLabel>1</UseCustomMapLabel>
<Description>&lt;html>&lt;body>Promoter for bacteriophage T7 RNA polymerase.&lt;/body>&lt;/html></Description>
<Created UTC="5:0:0">2012.6.28</Created>
<LastModified UTC="5:0:0">2012.9.6</LastModified>
<Organism>bacteriophage T7</Organism>
<SequenceClass>PHG</SequenceClass>
<TransformedInto>unspecified</TransformedInto>
</Notes>'''

    notes = snap.NotesBlock.from_bytes(bytes)

    assert notes.uuid == 'd528ff2f-6579-48ab-9301-ab40df3f2505'
    assert notes.type == 'Natural'
    assert notes.custom_map_label == 'T7 promoter'
    assert notes.use_custom_map_label == True
    assert notes.description == '<html><body>Promoter for bacteriophage T7 RNA polymerase.</body></html>'
    assert notes.date_created == datetime(2012, 6, 28)
    assert notes.date_last_modified == datetime(2012, 9, 6)
    assert notes.organism == 'bacteriophage T7'
    assert notes.sequence_class == 'PHG'
    assert notes.transformed_into == 'unspecified'

def test_to_bytes():
    notes = snap.NotesBlock()
    xml = lambda x: x.replace('\n', '').encode('utf8')

    assert notes.to_bytes() == xml('''\
<Notes />''')

    notes.type = 'Natural'
    assert notes.to_bytes() == xml('''\
<Notes>
<Type>Natural</Type>
</Notes>''')

    # Note that python escapes "&", "<", and ">", while snapgene just seems to 
    # escape "<".  I think python is doing the right thing here, so I'm going 
    # to see if it works.
    notes.description = '<html>'
    assert notes.to_bytes() == xml('''\
<Notes>
<Type>Natural</Type>
<Description>&lt;html&gt;</Description>
</Notes>''')

    notes.date_created = datetime(2012, 6, 28)
    assert notes.to_bytes() == xml('''\
<Notes>
<Type>Natural</Type>
<Description>&lt;html&gt;</Description>
<Created>2012.6.28</Created>
</Notes>''')

    notes.is_confirmed_experimentally = True
    assert notes.to_bytes() == xml('''\
<Notes>
<Type>Natural</Type>
<Description>&lt;html&gt;</Description>
<Created>2012.6.28</Created>
<ConfirmedExperimentally>1</ConfirmedExperimentally>
</Notes>''')

    del notes.type
    del notes.description
    del notes.date_created
    del notes.is_confirmed_experimentally
    assert notes.to_bytes() == xml('''\
<Notes />''')



