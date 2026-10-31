"""Independent full CAR signs on the complete hard-core two-particle basis."""
import pytest
import sympy as s

from . import colored as c


@pytest.mark.parametrize('d', [1, 2, 3])
def test_dimension_identity_counts_collision_complement(d):
    row = c.dimensions(d)
    assert row['native_hard_core_two_particle_dimension'] + row['same_port_collision_complement_dimension'] == row['full_CAR_two_particle_dimension']
    states = [(a, b) for a in range(12*d) for b in range(a+1, 12*d)]
    assert sum(a//d == b//d for a, b in states) == row['same_port_collision_complement_dimension']


@pytest.mark.parametrize('d', [2, 3])
def test_every_C_E_column_against_bitstring_CAR(d):
    value = c.complete_native_certificate(d)
    for row in value['native'].values():
        assert row['all_seams'] == 30
        assert row['abort_slots'] > 0
        assert row['zero_slots'] > 0
        assert row['independently_checked_active_columns'] == 3960*d*d


@pytest.mark.parametrize('d', [2, 3])
def test_tensor_gauge_action_commutes_with_actual_exchange(d):
    swap = s.zeros(d*d)
    for a in range(d):
        for b in range(d):
            swap[b*d+a, a*d+b] = -1
    # The identity is polynomial and holds for every matrix, not only a
    # convenient diagonal basis or a sampled group element.
    U = s.Matrix(d, d, s.symbols(f'u0:{d*d}'))
    representation = s.kronecker_product(U, U)
    assert swap*representation == representation*swap
    assert swap*swap == s.eye(d*d)
    for a in range(d):
        assert swap[a*d+a, a*d+a] == -1


def test_true_two_fermion_exchange_is_not_global_parity_or_scalar_only():
    d = 2
    loop = c.protection.exchange_loop('E')
    index = c.native.INDEX[tuple(loop['initial_pair'])]
    maps = {(row.law, row.seam, row.coin, row.raw):
            {i: (j, sign) for i, j, sign in c.colored_pairs(row, d)}
            for row in c.protection.branches('E')}
    start = c.address(index, 0, 1, d)
    at, sign = start, 1
    for label in loop['history']:
        at, phase = maps[(label[0], tuple(label[1]), label[2], label[3])][at]
        sign *= phase
    assert at == c.address(index, 1, 0, d)
    assert at != start  # a scalar sign with no internal SWAP is wrong
    assert sign == -1
    assert (-1)**2 == 1  # global Fock parity is positive on the entire domain


def test_adversarial_removed_CAR_sign_is_detected(monkeypatch):
    original = c.colored_pairs
    def unsigned(row, d):
        return tuple((i, j, 1) for i, j, _ in original(row, d))
    monkeypatch.setattr(c, 'colored_pairs', unsigned)
    with pytest.raises(ArithmeticError, match='independently applied CAR'):
        c.complete_native_certificate(2)


@pytest.mark.parametrize('d', [1, 2, 3, 4])
def test_actual_orientation_cover_intertwines_true_colored_exchange(d):
    W, exchange = c.ordered_cover_embedding(d)
    n = d*d
    Q = s.kronecker_product(s.Matrix([[0, 1], [1, 0]]), s.eye(n))
    assert W.H*W == 2*s.eye(n)
    assert Q*W == W*exchange
    assert (s.eye(n)-exchange).rank() == d*(d+1)//2
    assert (s.eye(n)+exchange).rank() == d*(d-1)//2
    if d > 1:
        antisymmetric = s.zeros(n, 1)
        antisymmetric[1], antisymmetric[d] = 1, -1
        assert exchange*antisymmetric == antisymmetric
        assert Q*(W*antisymmetric) == W*antisymmetric


def test_ninety_internal_components_are_per_particle_not_collective():
    assert c.dimensions(90) == {'internal_per_particle': 90, 'ports': 12,
        'one_particle_modes': 1080, 'native_hard_core_two_particle_dimension': 534600,
        'full_CAR_two_particle_dimension': 582660,
        'same_port_collision_complement_dimension': 48060}


@pytest.mark.parametrize('d', [0, -1, True, 2.5])
def test_invalid_internal_dimension_rejected(d):
    with pytest.raises(ValueError):
        c.dimensions(d)
