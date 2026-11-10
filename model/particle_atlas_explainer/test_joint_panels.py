"""The joined interaction graphic uses exact receipts without promoting scope."""
from copy import deepcopy
import xml.etree.ElementTree as ET

import pytest

from . import joint_panels
from .canvas import Canvas


@pytest.fixture(scope='module')
def joint():
    return joint_panels.load()


def visible(value):
    c = Canvas()
    joint_panels.draw(c, {'joint': value}, 0)
    svg = '<svg xmlns="http://www.w3.org/2000/svg">' + ''.join(c.items) + '</svg>'
    return ' '.join(' '.join(ET.fromstring(svg).itertext()).split()), c


def test_both_reciprocal_readouts_and_preparation_domains_are_visible(joint):
    text, c = visible(joint)
    assert not c.overflow
    for phrase in ('72,000 basis states', 'Click probability: 1/120 versus 0',
                   'Face detector: 1 versus 0', 'Return: 1 unread; 1/2 after readout',
                   '281/1080', '36,000 basis states', '281/2160',
                   'the retained link frame to compare amplitudes',
                   'Every failure is retained.', 'rest masses and an interacting relativistic continuum are not identified.'):
        assert phrase in text


def test_frequency_and_species_boundaries_are_visible(joint):
    text, c = visible(joint)
    assert not c.overflow
    for phrase in ('Fresh preparation success: 1/10', '283/7200',
                   '40 face/filter/clock records', 'quadratic gapless spectator branch',
                   'Equal triplet spectra cannot be three unequal lepton masses',
                   'Fermionic hopping and tensor-wave exchange are distinct spatial laws'):
        assert phrase in text


def test_numeric_preparation_label_is_receipt_driven(joint):
    changed = deepcopy(joint)
    changed['reference']['success_probability_from_uniform_basis_ensemble'] = '2/7'
    assert 'probability 2/7' in visible(changed)[0]


def test_unconstrained_return_cannot_be_advertised_as_invariant(monkeypatch, joint):
    changed = deepcopy(joint)
    changed['interaction']['constrained_coherence']['return_projector_commutes_with_Gamma'] = False
    monkeypatch.setattr(joint_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='admitted gauge sector'):
        joint_panels.load()


@pytest.mark.parametrize('claim', ['physical_particle_species', 'physical_energy_or_GeV_mass',
                                 'interacting_3plus1_continuum'])
def test_physical_promotion_rejected(monkeypatch, joint, claim):
    changed = deepcopy(joint)
    changed['scope'][claim] = True
    monkeypatch.setattr(joint_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='scope promotion'):
        joint_panels.load()
