"""The repair diagram keeps its receipt ancestry and physical claim boundary."""
from copy import deepcopy
import xml.etree.ElementTree as ET

import pytest

from . import repair_panels
from .canvas import Canvas


@pytest.fixture(scope='module')
def repair():
    return repair_panels.load()


def text_and_canvas(value):
    c = Canvas()
    repair_panels.draw(c, {'repair': value}, 0)
    svg = '<svg xmlns="http://www.w3.org/2000/svg">' + ''.join(c.items) + '</svg>'
    return ' '.join(' '.join(ET.fromstring(svg).itertext()).split()), c


def test_one_program_readout_and_boundaries_are_visible(repair):
    text, c = text_and_canvas(repair)
    assert not c.overflow
    for phrase in (
        'Repair can move disagreement into memory',
        'Commit a classical outcome at every event',
        'Reuse coherent memory during propagation',
        'Success: 131/810 per attempt',
        '1,440 source-and-memory states',
        '2,073,600 Hermitian directions',
        '7,200 states with one unitary clock',
        'Return-amplitude change: -1/6750',
        'The five code states are not five species.',
        'The wave experiment’s link background is supplied and static.',
        'source-derived GeV mass follows from this static-link construction.',
    ):
        assert phrase in text
    assert 'No dynamical gauge backreaction' in text


def test_numerical_labels_follow_the_bound_receipt(repair):
    changed = deepcopy(repair)
    changed['reference']['success_probability'] = '1/7'
    changed['gauged_walk']['static_flux_witness']['amplitude_difference'] = '-1/6000'
    text, _ = text_and_canvas(changed)
    assert 'Success: 1/7 per attempt' in text
    assert 'Return-amplitude change: -1/6000' in text


@pytest.mark.parametrize('claim', ['physical_particle_identification', 'continuum_limit', 'GeV_mass_prediction'])
def test_physical_promotion_is_rejected(monkeypatch, repair, claim):
    changed = deepcopy(repair)
    changed['scope'][claim] = True
    monkeypatch.setattr(repair_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='physical scope promotion'):
        repair_panels.load()


def test_stale_repair_source_is_rejected(monkeypatch, repair):
    changed = deepcopy(repair)
    path = next(iter(changed['input_sha256']))
    changed['input_sha256'][path] = '0' * 64
    monkeypatch.setattr(repair_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='stale coherent repair input'):
        repair_panels.load()


def test_boundary_and_observer_comparison_remain_separate_programs(repair):
    c = Canvas()
    repair_panels.draw(c, {'repair': repair, 'protection': repair_panels.load_protection(),
                          'agreement': repair_panels.load_agreement()}, 0)
    svg = '<svg xmlns="http://www.w3.org/2000/svg">' + ''.join(c.items) + '</svg>'
    text = ' '.join(' '.join(ET.fromstring(svg).itertext()).split())
    assert not c.overflow
    assert 'probability 1/60' in text and 'old detector gives 1/120' in text
    assert 'does not certify every possible observer response' in text
    assert '30-bit link register' in text and 'all 5 region faces' in text
    assert '20 checked moves and 30 retained holds' in text
    assert 'no proved joint evolution or stationary reference with the coherent wave' in text


def test_boundary_protection_cannot_be_promoted_to_all_link_moves(monkeypatch):
    changed = repair_panels.load_protection()
    changed['protection']['arbitrary_link_updates_protect_initial_pair'] = True
    monkeypatch.setattr(repair_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='protection or physical scope'):
        repair_panels.load_protection()


def test_finite_chart_agreement_cannot_be_promoted_to_full_response(monkeypatch):
    changed = repair_panels.load_agreement()
    changed['A1_complete_response_certified'] = True
    monkeypatch.setattr(repair_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='scope promotion'):
        repair_panels.load_agreement()
