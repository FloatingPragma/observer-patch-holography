"""Full-source intertwining, physical seed, complete weights and custody errors."""
import numpy as np
import pytest
import sympy as s

from . import model
from .symmetries import certificate as symmetry_certificate
from ..periodic_family.model import local_yukawa


def test_full_native_integer_action_matches_quantum_Yukawa_columns():
    masks = (0, 1 << 11, 1 << 10, (1 << 30) - 1,
             sum(1 << j for j in range(0, 30, 2)), 123456789, 987654321)
    for mask in masks:
        for scalar in range(5):
            expected = dict(local_yukawa(mask, scalar))
            if scalar:
                expected[mask, scalar] = expected.get((mask, scalar), 0) + 1
            expected = {key: value for key, value in expected.items() if value}
            assert model.column((mask, scalar)) == expected


def test_exact_cyclic_operator_and_complete_energy_records():
    result = model.spectrum()
    assert result['full_seed_component_dimension'] == 113948
    assert result['equitable_dimension'] == 403
    assert result['squared_physical_cyclic_dimension'] == 25
    assert result['full_physical_cyclic_dimension'] == 49
    assert len(result['complete_energy_records']) == 49
    assert sum(row['energy_interval'] == ['1', '1'] for row in result['complete_energy_records']) == 1
    assert s.Rational(result['odd_minus_vacuum_component_ground_interval'][1]) < 0
    assert not result['global_ground_across_all_local_sectors_proved']
    assert not result['seed_is_Ddagger_on_interacting_ground']


def test_actual_Gauss_spin_and_centered_charge_conjugation():
    result = symmetry_certificate()
    assert result['actual_local_spin_Casimir'] == '3/4'
    assert result['centered_charge_L_minus_one'] == '0'
    assert result['charge_conjugate_ground_is_opposite_spin_partner']
    assert result['wrong_relative_seed_sign_rejected']
    assert not result['raw_L_is_identified_with_measured_lepton_charge']


def test_partition_rejects_wrong_source_and_false_intertwiner(tmp_path, monkeypatch):
    original = model.PARTITION
    with np.load(original, allow_pickle=False) as archive:
        colors = archive['colors'].copy(); digest = str(archive['source_digest'])
    target = tmp_path / 'partition.npz'
    monkeypatch.setattr(model, 'PARTITION', target)
    model.quotient.cache_clear()
    np.savez_compressed(target, colors=colors, source_digest=np.array('wrong source'))
    with pytest.raises(ValueError, match='different full source'):
        model.quotient()
    bad = colors.copy(); bad[100] = (int(bad[100]) + 1) % 403
    np.savez_compressed(target, colors=bad, source_digest=np.array(digest))
    with pytest.raises(ValueError):
        model.quotient()
    model.quotient.cache_clear()


def test_missing_declared_and_transitive_inputs_fail_closed(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-odd-local-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-odd-local-import.py'})
    with pytest.raises(FileNotFoundError, match='missing neutral local import'):
        custody.input_hashes()
