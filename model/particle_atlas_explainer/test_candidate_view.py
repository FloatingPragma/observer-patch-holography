"""Current candidate view preserves the numerical base and source scope."""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET

import pytest

from . import candidate_view as view


@pytest.fixture(scope='module')
def inputs():
    return view.load()


def test_retained_numerical_base_and_full_source_catalogue(inputs):
    comparison, original, confirmation = inputs
    result = view.render(original, confirmation)
    root = ET.fromstring(result)
    ns = '{http://www.w3.org/2000/svg}'
    panel = root.find(ns+'g[@id="candidate-confirmation-panel"]')
    assert panel is not None
    old_root = ET.fromstring(original)
    for old, retained in zip(list(old_root), list(root)):
        assert ET.tostring(old) == ET.tostring(retained)
    text = ' '.join(''.join(panel.itertext()).split())
    for family in confirmation['labelled_catalogue']['protected_families']:
        for row in family['rows']:
            assert row['label'] in text and row['source_id'] in text
    assert 'A catalogue is not a cosmological birth order' in text
    assert 'Six tests connect a source response to a particle' in text
    assert comparison['physical_mass_predictions'] == []
    assert int(root.get('height')) > int(old_root.get('height'))


def test_confirmed_names_and_masses_are_rejected(monkeypatch, inputs):
    comparison, original, confirmation = inputs
    monkeypatch.setattr(view.base, 'build', lambda: (comparison, original))
    forged = deepcopy(confirmation)
    forged['physical_particle_map'] = [{'source': 'loop', 'particle': 'photon'}]
    monkeypatch.setattr(view.build, 'receipt', lambda _: forged)
    with pytest.raises(ValueError):
        view.load()


def test_generated_current_view_binds_inputs():
    svg, audit = view.produce()
    assert view.OUTPUT.read_text() == svg
    assert view.AUDIT.read_text() == audit
    receipt = json.loads(audit)
    assert receipt['protected_catalogue_count'] == 11
    assert receipt['source_excitation_count'] == 16
    assert receipt['physical_particle_map'] == receipt['mass_predictions_GeV'] == []
