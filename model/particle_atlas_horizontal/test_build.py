"""Scope, complete coverage and reference/rendering regressions for the spread."""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import pytest
from . import build, drawing


@pytest.fixture(scope='module')
def data():return build.load()


def test_every_family_and_each_owed_bridge_survives(data):
    build.validate_content(data)
    assert len(data['families']['lanes'])==9
    assert len(data['families']['calibration_paths'])==10
    assert len(data['theorems']['roles'])==5
    assert data['v1']['physical_particle_map']==[]


@pytest.mark.parametrize('kind',['missing_target','duplicate_target','promoted_identity','promoted_edge','lost_strategy'])
def test_visual_shortcuts_cannot_silently_promote_the_science(data,kind):
    bad=deepcopy(data)
    if kind=='missing_target':bad['families']['lanes'][0]['target_ids'].pop()
    elif kind=='duplicate_target':bad['families']['lanes'][0]['target_ids'][0]='muon'
    elif kind=='promoted_identity':bad['families']['lanes'][0]['targets'][0]['source_physical_identification']=True
    elif kind=='promoted_edge':bad['story']['edges'][2]['status']='exact'
    else:bad['families']['calibration_paths'][0]['removal_strategy']=''
    with pytest.raises(ValueError):build.validate_content(bad)


def test_render_has_real_measurement_links_and_all_derivation_paths(data):
    svg=drawing.render(data)
    root=ET.fromstring(svg)
    assert float(root.get('width'))>float(root.get('height'))
    assert (float(root.get('width')),float(root.get('height'))) == (48000,47065)
    hrefs=[n.get('href') for n in root.iter('{http://www.w3.org/2000/svg}a')]
    assert any('OperationalBornBridge.lean#L214' in h for h in hrefs)
    assert any('PrimitivePortFrameQuotient.lean#L631' in h for h in hrefs)
    assert any('ModularClockComparison.lean#L40' in h for h in hrefs)
    assert sum('.pdf#page=' in h for h in hrefs)>20
    text=' '.join(''.join(root.itertext()).split())
    assert 'no single sharp value' in text
    assert 'one history: 1/8' in text
    assert 'opposite winding does not identify an antiparticle' in text
    assert text.count('Pass:')==10
    for lane in data['families']['lanes']:
        assert lane['title'] in text


def test_render_uses_certified_complete_probabilities(data):
    value=deepcopy(data)
    value['measurement']['position_probabilities']=['1/3','1/3','1/3','0']
    svg=drawing.render(value)
    # A display-only mutation must reach actual bars; it is never a source certificate.
    assert '>1/3</text>' in svg


def test_reader_embeds_the_vector_and_preserves_every_reference_anchor(data):
    text=build.reader(drawing.render(data),data)
    assert '<svg ' in text and 'fetch(' not in text
    assert 'calibrations:[32000,6150]' in text
    assert 'pointermove' in text and 'ArrowRight' in text
    assert 'ModularClockComparison.lean#L40' in text
    assert 'FiniteBuschGleason.lean#L1133' in text
    assert '>Nine dependencies</button>' in text
    assert 'model/particle_primary_locality/PROOF.md' in text
    assert 'model/particle_joint_reconstruction/DEPENDENCIES.md' in text
    assert 'model/particle_joint_reconstruction/CANDIDATE.md' in text


def test_generated_outputs_match_their_exact_source_inventory():
    svg,html,audit=build.build()
    assert build.SVG.read_text()==svg
    assert build.HTML.read_text()==html
    assert build.AUDIT.read_text()==audit
    assert json.loads(audit)['source_physical_identifications']==[]
