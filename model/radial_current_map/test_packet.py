"""Reject deleted complement evidence, altered currents and physical promotions."""
from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet._produce()


def test_fresh_scope_custody_and_replay(receipt):
    packet.validate(receipt, receipt)
    assert receipt['source_bindings'] == packet.source_bindings()
    assert not any(p in receipt['source_bindings'] for p in
                   ['model/verify.py', 'model/completion/program_class_v2.json'])


@pytest.mark.parametrize('field', ['physical_defect_to_particle_map_complete',
    'eta_identified_with_electroweak_or_Lorentz_spin',
    'native_repair_preserves_eta_proved', 'all_onsite_even_Hermitian_CAR_currents_classified'])
def test_scope_changes_rejected(receipt, field):
    changed = deepcopy(receipt); changed['scope'][field] = not changed['scope'][field]
    with pytest.raises(ValueError): packet.validate(changed, receipt)


def test_false_ground_doublet_and_deleted_quartic_equations_rejected(receipt):
    changed = deepcopy(receipt)
    changed['actual_low_current_map_and_leakage']['actual_interacting_low_representation']['eta_Casimir'][1] = '3/4'
    with pytest.raises(ValueError): packet.validate(changed, receipt)
    changed = deepcopy(receipt)
    del changed['exact_source_current_classification']['all_three_edge_classifications']
    with pytest.raises(ValueError): packet.validate(changed, receipt)


def test_boolean_physical_count_and_stale_import_rejected(receipt, monkeypatch):
    changed = deepcopy(receipt); changed['physical_comparison_count_added'] = False
    with pytest.raises(ValueError): packet.validate(changed, receipt)
    monkeypatch.setattr(packet, '_IMPORT_BINDINGS', {})
    with pytest.raises(ValueError, match='after.*import'): packet._produce()
