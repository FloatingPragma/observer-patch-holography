from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet._produce()


def test_exact_dispersion_and_all_sources_bound(receipt):
    packet.validate(receipt)
    pins = receipt['source_bindings']
    assert 'model/mass_generation_audit/results.json' in pins
    assert 'model/periodic_momentum/momentum.py' in pins
    assert 'paper/tex_fragments/HUBBARD_RELATIVISTIC_SECTOR.tex' in pins


@pytest.mark.parametrize('key', ['running_family_derived_from_fixed_bare_point',
    'native_ground_preparation_executed', 'native_local_CAR_pole_proved',
    'three_dimensional_massive_species_identified', 'physical_SI_mass_calculated'])
def test_scope_promotions_rejected(receipt, key):
    value = deepcopy(receipt)
    value['scope'][key] = True
    with pytest.raises(ValueError): packet.validate(value)


def test_same_interpreter_source_edit_rejected(monkeypatch):
    changed = packet.source_bindings()
    changed['model/hubbard_relativistic_sector/dispersion.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='source changed since import'):
        packet._produce()
