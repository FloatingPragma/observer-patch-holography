"""Independent complete-carrier and retained-record controls."""
import numpy as np
import pytest

from . import algebra as a
from ..defect_field_bridge import full_domain


@pytest.mark.parametrize('sites', [2, 3, 4])
@pytest.mark.parametrize('law', ['C', 'E'])
def test_every_native_branch_preserves_defect_number_and_protected_spin(sites, law):
    value = a.source_audit(sites, law, (0, sites-1))
    assert value['full_Fock_dimension'] == 4**sites
    assert value['coin_raw_slots'] == 32


@pytest.mark.parametrize('law', ['C', 'E'])
def test_executed_source_dilation_has_matching_coherent_protected_coin(law):
    actual, target, _ = full_domain.source_isometry(2, 2, law)
    P = a.protected(2)
    expected = np.zeros_like(target[:, P])
    masks = np.flatnonzero(P)
    for j, mask in enumerate(masks):
        expected[mask, j] = 1/np.sqrt(2)
        expected[mask+16, j] = 1/np.sqrt(2)
    assert np.linalg.norm(actual[:, P]-expected) < full_domain.SOURCE_BOUND
    assert np.array_equal(target[:, P] != 0, expected != 0)
    state = actual[:, masks[0]]
    assert abs(state[masks[0]]*state[masks[0]+16].conjugate()) > .49


def test_unital_spin_embedding_preserves_products_and_all_onsite_charge_projections():
    rng = np.random.default_rng(118)
    first, second = (rng.integers(-3, 4, size=(8, 8)) for _ in range(2))
    J, K = a.spin_embedding(3, first), a.spin_embedding(3, second)
    assert np.array_equal(J@K, a.spin_embedding(3, first@second))
    assert np.array_equal(J.T, a.spin_embedding(3, first.T))
    assert np.array_equal(a.spin_embedding(3, np.eye(8, dtype=int)), np.eye(64))
    parity = np.array([(-1)**mask.bit_count() for mask in range(64)])
    assert not np.any((parity[:, None]-parity[None, :])*J)
    for x in range(3):
        q = ~a.protected(3, (x,))
        assert not np.any((q[:, None].astype(int)-q[None, :])*J)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_charge_diagonal_hopping_and_repair_energy_increment_have_required_support(law):
    a.hopping_audit(law)
    T = a.hopping_twice(3, ((0, 1, -1), (1, 2, 1)))
    operators = [a.annihilator(3, m) for m in range(6)]
    for j, first in enumerate(operators):
        for k, second in enumerate(operators):
            assert np.array_equal(first@second+second@first, np.zeros((64, 64)))
            assert np.array_equal(first@second.T+second.T@first, np.eye(64)*(j == k))
    assert np.array_equal(T, T.T)
    assert np.linalg.norm(T@np.diag(a.defects(3))-np.diag(a.defects(3))@T) > 1


@pytest.mark.parametrize('law', ['C', 'E'])
def test_full_coherent_record_correction_is_charge_supported_for_joint_readouts(law):
    rng = np.random.default_rng(118)
    spin, empty, tape = [rng.integers(-2, 3, size=(n, n)) for n in (4, 4, 32)]
    correction = a.record_correction_twice(spin, tape, empty, law)
    P, N = a.protected(2), a.defects(2)
    assert not np.any(correction[:, P])
    assert not np.any(correction[P, :])
    assert not np.any((N[:, None]-N[None, :])*correction)
    assert np.any(correction)
    # A classical coin substituted on a coherent tape fails on protected input.
    cross_coin = np.zeros((32, 32), dtype=int)
    cross_coin[0, 16] = cross_coin[16, 0] = 1
    identity = np.eye(4, dtype=int)
    good = a.record_correction_twice(identity, cross_coin, np.zeros((4, 4), dtype=int), law)
    assert not np.any(good[np.ix_(P, P)])
    wrong = good+2*a.spin_embedding(2, identity)
    assert np.linalg.norm(wrong[np.ix_(P, P)], 2) == 2


def test_dropping_collision_branches_breaks_complete_law():
    matrices = a.branch_matrices(2, 'C', (0, 1))
    effects = sum(K.T@K for j, K in enumerate(matrices) if j % 16 != 8)/2
    assert not np.array_equal(effects, np.eye(16))
    assert effects[3, 3] == 0
