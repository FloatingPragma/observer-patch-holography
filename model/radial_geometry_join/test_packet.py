from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_cold_receipt_matches_retained_full_source(fresh):
    packet.validate(strict_json.load(packet.OUTPUT), fresh)
    assert fresh['native_histories']['complete_native_histories_checked'] == 79872
    assert 'model/radial_geometry_join/PROOF.md' in fresh['source_bindings']


@pytest.mark.parametrize('damage', ['gap', 'contrast', 'failure', 'work', 'proof', 'parent', 'count'])
def test_changed_join_or_erased_custody_fails(fresh, damage):
    value = deepcopy(fresh)
    if damage == 'gap':
        value['same_source_bounds']['actual_charged_gap_lower'] = '1'
    elif damage == 'contrast':
        value['same_source_bounds']['simple_geometry_contrast_strict_lower'] = '1'
    elif damage == 'failure':
        value['quantum_instrument']['native_first_failure_slots'] = 0
    elif damage == 'work':
        value['quantum_instrument']['field_energy_injection'] = '0'
    elif damage == 'proof':
        value['source_bindings'].pop('model/radial_geometry_join/PROOF.md')
    elif damage == 'parent':
        value['inherited_custody'].pop('radial_interacting_core')
    else:
        value['physical_comparison_count_added'] = False
    with pytest.raises(ValueError):
        packet.validate(value, fresh)


@pytest.mark.parametrize('key', list(packet.scope()))
def test_every_state_operation_or_physical_boundary_is_retained(fresh, key):
    value = deepcopy(fresh)
    value['scope'][key] = not value['scope'][key]
    with pytest.raises(ValueError, match='scope'):
        packet.validate(value, fresh)


def test_cold_replay_uses_source_not_mutated_imported_bounds(monkeypatch, fresh):
    monkeypatch.setattr(packet.bounds, 'joined_bounds', lambda: {'counterfeit': True})
    packet.validate(packet.report(), fresh)
