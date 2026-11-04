"""Complete correlated-reference repair identity and hopping distinction."""
import numpy as np
import pytest

from model.defect_field_bridge import full_domain
from model.native_repair_spin_limit.algebra import branch_matrices, protected
from . import repair


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_coherent_archive_preserves_all_protected_reference_columns(law):
    matrices = branch_matrices(2, law, (0, 1))
    P = protected(2)
    # Compare the entire protected isometry, which implies the equality
    # after tensoring any untouched reference, not just diagonal inputs.
    actual = np.vstack([K[:, P] for K in matrices])
    expected = np.zeros_like(actual)
    expected[:16] = np.eye(16, dtype=np.int64)[:, P]
    expected[16*16:17*16] = np.eye(16, dtype=np.int64)[:, P]
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(actual.T@actual, 2*np.eye(4, dtype=np.int64))
    assert full_domain.RAW_SLOTS == 16


@pytest.mark.parametrize('sites,law', [(2, 'C'), (2, 'E'), (3, 'C'), (3, 'E')])
def test_hopping_is_retained_as_a_nonidentity_countercontrol(sites, law):
    value = repair.audit(sites, law)
    assert value['held_collision_columns'] > 0
    assert value['hopping_outside_protected_Frobenius_squared_for_T2'] > 0
    assert value['retained_coin_raw_slots'] == 32
    assert not value['interleaved_Hamiltonian_evolution_is_identity']
