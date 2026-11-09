import json

import pytest

from . import packet


def test_replay_retains_scope_and_rejects_physical_promotion(tmp_path):
    value = packet.build()
    path = tmp_path/'receipt.json'
    path.write_text(packet.canonical(value))
    packet.verify(path)
    value['scope']['physical_particle_assignments'] = ['electron']
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='receipt differs'):
        packet.verify(path)


def test_replay_rejects_changed_source_coefficients(tmp_path):
    value = packet.build()
    value['exact_source_algebra']['cubic_hessian_through_second_order'] = '3*J times I3'
    path = tmp_path/'receipt.json'
    path.write_text(packet.canonical(value))
    with pytest.raises(ValueError, match='receipt differs'):
        packet.verify(path)


def test_native_preparation_dependencies_and_finite_scope_bound():
    value = packet.build()
    assert {
        'model/native_spin_observer/precision.py',
        'model/native_spin_observer/algebra.py',
        'model/spin_twist_closure/__init__.py',
        'code/source_selection_model/response.json',
        'evidence/native_geometric_source_20260925/spec.json',
    } <= value['source_bindings'].keys()
    prep = value['short_complete_native_preparation']
    assert 'two-state spins' in prep['finite_product_ground_bound']['ground_Hilbert_space']
    assert not prep['infinite_volume_Yarotsky_threshold_certified_for_this_coupling']
