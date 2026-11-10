"""Scope, certified inputs and layout regressions for the explainer atlas."""
from copy import deepcopy
import hashlib
import json
import struct
import xml.etree.ElementTree as ET

import pytest

from . import build
from .canvas import Canvas, text_width


@pytest.fixture(scope='module')
def data():
    return build.load()


def test_render_is_dark_compact_and_complete(data):
    svg = build.render(data)
    root = ET.fromstring(svg)
    assert root.get('width') == '1600'
    # Includes the mobile-object, common-clock and complete charged-response panels.
    assert 12000 < int(root.get('height'))
    text = ' '.join(''.join(root.itertext()).split())
    for phrase in ('From a twist to a particle', 'Topological defect:', 'Why nothing local can erase it',
                   'Position is where a record is written', 'Velocity is the slope of the records',
                   'Mass is the price of standing still', 'Observers compress the records into fields'):
        assert phrase in text
    assert text.count('Computed in the model'.upper()) >= 9
    assert '—' not in text


def test_every_label_meets_a_minimum_size(data):
    svg = build.render(data)
    sizes = [float(n.get('font-size')) for n in ET.fromstring(svg).iter('{http://www.w3.org/2000/svg}text')]
    assert min(sizes) >= 11.5


def test_detector_inset_reads_the_certified_sector_law(data):
    bad = deepcopy(data)
    bad['attachment']['same_time_scalar_records']['-1']['occupied'] = '1/2'
    assert '1 of 2 runs, twist present' in build.render(bad)
    assert '3 of 4 runs, twist present' in build.render(data)


def test_dependency_strip_follows_the_register(data):
    done = [n['status'] == 'finite_constructed' for n in data['nodes']]
    assert done[:3] == [True, True, True] and not any(done[3:])


def test_evidence_panel_uses_receipts_and_keeps_domains_separate(data):
    def visible_text(value):
        return ' '.join(' '.join(ET.fromstring(value).itertext()).split())

    svg = build.render(data)
    text = visible_text(svg)
    assert '22,531 native pulses' in text
    assert 'Source-derived GeV mass table: empty' in text
    assert 'rare empty-matter sector' in text
    assert 'separate supplied contact action' in text
    assert 'Static fusion is not a moving collision' in text
    assert '18 named targets; no identified species' in text
    assert 'Six source identities; seven records' in text
    assert 'A dressed charged energy channel' in text
    assert 'pair-energy transfers in [77/8, 271/8] source units' in text
    assert 'probability above 1/40' in text
    assert len(data['identity_tests']['target_rows']) == 18
    changed = deepcopy(data)
    changed['candidates']['collision_status']['separate_contact_extension']['native_pulses'] = 12345
    assert '12,345 native pulses' in build.render(changed)
    changed['identity_tests']['calculated_source_results']['charged']['original_compensated_dressed_energy_window']['original_dressed_recorded_simple_strict_lower'] = '1/41'
    assert 'probability above 1/41' in visible_text(build.render(changed))
    assert 'Repair treats every seam alike, so momentum is conserved' not in text
    assert 'that standing cost is the defect’s rest energy' not in text


def test_candidate_mass_promotion_is_rejected(monkeypatch):
    original = build.receipt

    def promoted(path):
        value = original(path)
        if path == build.CANDIDATES:
            value['mass_predictions_GeV'] = [{'electron': '0.000511'}]
        return value
    monkeypatch.setattr(build, 'receipt', promoted)
    with pytest.raises(ValueError, match='candidate map changed'):
        build.load()


def test_broken_producer_custody_is_rejected(monkeypatch):
    monkeypatch.setattr(build, 'sha', lambda _: '0' * 64)
    with pytest.raises(ValueError, match='stale source custody'):
        build.receipt(build.CANDIDATES)


def test_a_physical_promotion_is_rejected(data, monkeypatch):
    original = build.json.loads

    def promoted(text):
        value = original(text)
        if isinstance(value, dict) and 'holonomy_transport' in value:
            value['holonomy_transport']['physical_particle_assignments'] = ['electron']
        return value
    monkeypatch.setattr(build.json, 'loads', promoted)
    with pytest.raises(ValueError):
        build.load()


def test_overflowing_text_fails_the_build():
    c = Canvas()
    c.text(0, 0, 'a label that is far too long for its slot', 20, maxw=40)
    assert c.overflow and text_width('abc', 20) > 0


def test_generated_outputs_match_their_inputs():
    svg, audit = build.build()
    assert build.SVG.read_text() == svg
    assert build.AUDIT.read_text() == audit
    assert json.loads(audit)['source_physical_identifications'] == []


def test_png_binds_the_same_svg_and_dimensions():
    raw = build.PNG.read_bytes()
    assert raw[:8] == b'\x89PNG\r\n\x1a\n'
    root = ET.fromstring(build.SVG.read_bytes())
    assert struct.unpack('>II', raw[16:24]) == build.png_dimensions(build.SVG.read_text())
    assert max(struct.unpack('>II', raw[16:24])) <= 32000
    provenance = b'Source SVG SHA256\0' + hashlib.sha256(build.SVG.read_bytes()).hexdigest().encode()
    assert provenance in raw


def test_grouped_catalogue_and_current_response_bounds(data):
    text = ' '.join(' '.join(ET.fromstring(build.render(data)).itertext()).split())
    catalogue = data['identity_tests']['labelled_catalogue']
    for family in catalogue['protected_families']:
        for row in family['rows']:
            assert row['label'] in text and row['source_id'] in text
    assert 'Probability > 3/80' in text
    assert 'Total probability ≥ 3481/160000' in text
    assert '40 component witnesses' in text
    assert 'Final energy is not transfer energy or rest mass' in text
    assert 'A catalogue is not a cosmological birth order' in text
    assert 'Source ID: bare_native_source' in text
    assert 'Six tests connect a source response to a particle' in text
    assert 'OPH_PARTICLE_ATLAS_V2' not in build.render(data)
