"""Independent full-CAR action, eta dynamics and bounded exact-admission tests."""
import math
import warnings

import numpy as np
import pytest
from scipy.linalg import expm

from model.full_fock_field.field import FieldGeometry
from model.interacting_field import lattice
from model.field_symmetry import symmetry as s


def jw(modes):
    result = []
    for mode in range(modes):
        value = np.ones((1, 1), complex)
        for bit in reversed(range(modes)):
            value = np.kron(value, np.array([[0, 1], [0, 0]]) if bit == mode else
                            np.diag([1, -1]) if bit < mode else np.eye(2))
        result.append(value)
    return result


@pytest.fixture(scope='module')
def square():
    geometry = FieldGeometry((2, 2, 1))
    return geometry, s.operators(geometry), jw(8)


def test_independent_jordan_wigner_action_and_full_exterior_frame(square):
    geometry, ops, a = square
    n = [op.conj().T@op for op in a]
    # Site order is (0,0),(0,1),(1,0),(1,1). Use literal coordinate bonds,
    # not FieldGeometry.groups or the production routing helper.
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]))
    h = sum(n[2*x]@n[2*x+1] for x in range(4))
    for x, y, axis in ((0, 1, 1), (0, 2, 0), (1, 3, 0), (2, 3, 1)):
        directed = sum(1j*sigma[axis][p, q]/16*a[2*x+p].conj().T@a[2*y+q]
                       for p in range(2) for q in range(2))
        h += directed+directed.conj().T
    np.testing.assert_array_equal(ops.H.matrix(ops.scale).toarray(), h)
    frame = ops.frame.matrix().toarray()
    # Independently lift the one-body frame by every exterior-power minor.
    one = s.one_body_frame(geometry).toarray()
    for count in (0, 1, 2, 3, 4, 5, 6, 7, 8):
        addresses = [mask for mask in range(256) if mask.bit_count() == count]
        for mask in addresses:
            occupied = [j for j in range(8) if mask >> j & 1]
            target_modes = [int(np.flatnonzero(one[:, j])[0]) for j in occupied]
            target = sum(1 << j for j in target_modes)
            determinant = np.linalg.det(one[np.ix_(sorted(target_modes), occupied)]) if count else 1
            assert frame[target, mask] == determinant
            assert np.count_nonzero(frame[:, mask]) == 1


@pytest.mark.parametrize('theta,phase,duration', [(math.pi/2, 0., 1.), (.37, -.8, -.6), (0., .1, 2.)])
def test_complete_eta_product_evolves_at_U_under_H_and_is_stationary_under_K(square, theta, phase, duration):
    geometry, ops, _ = square
    initial = s.eta_product_state(geometry, theta, phase)
    h, k = ops.H.matrix(ops.scale).toarray(), ops.K.matrix(ops.scale).toarray()
    target = s.eta_product_state(geometry, theta, phase-geometry.interaction*duration)
    np.testing.assert_allclose(expm(-1j*duration*h)@initial, target, rtol=0, atol=2e-15)
    np.testing.assert_allclose(expm(-1j*duration*k)@initial, initial, rtol=0, atol=2e-15)
    plus = ops.eta_plus.matrix().toarray()
    expected = 2*math.sin(theta)*np.exp(-1j*phase+1j*duration)
    assert abs(np.vdot(target, plus@target)-expected) < 2e-15
    # The exact stationary K tower is not the ground reference.
    assert np.linalg.eigvalsh(k)[0] < -1
    assert abs(np.vdot(initial, k@initial)) < 1e-15


def test_every_fixed_number_eta_tower_column_and_every_complete_bond(square):
    geometry, ops, _ = square
    state = np.eye(256, dtype=complex)[:, 0]
    for pairs in range(geometry.volume+1):
        np.testing.assert_allclose(ops.H.matrix(ops.scale)@state, pairs*state, atol=1e-15, rtol=0)
        np.testing.assert_allclose(ops.K.matrix(ops.scale)@state, 0, atol=1e-15, rtol=0)
        np.testing.assert_allclose(ops.N.matrix()@state, 2*pairs*state, atol=1e-15, rtol=0)
        for bond in ops.bonds:
            np.testing.assert_allclose(bond.matrix()@state, 0, atol=1e-15, rtol=0)
        if pairs < geometry.volume:
            state = ops.eta_plus.matrix()@state/math.sqrt((pairs+1)*(geometry.volume-pairs))
    assert abs(np.vdot(state, state)-1) < 2e-15


def test_component_hop_is_not_the_complete_eta_preserving_bond(square):
    geometry, ops, a = square
    initial = s.eta_product_state(geometry)
    # First y-bond: the two disjoint component hops commute; neither alone
    # annihilates the eta state. This is the relevant source grouping premise.
    first = (a[0].conj().T@a[3]+a[3].conj().T@a[0])/16
    second = -(a[1].conj().T@a[2]+a[2].conj().T@a[1])/16
    np.testing.assert_array_equal(first@second, second@first)
    assert np.linalg.norm(first@initial) > .02
    assert np.linalg.norm(second@initial) > .02
    np.testing.assert_allclose((first+second)@initial, 0, atol=1e-16, rtol=0)
    np.testing.assert_allclose(expm(-.7j*first)@expm(-.7j*second)@initial,
                               initial, atol=2e-16, rtol=0)


@pytest.mark.parametrize('sides', [(2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1), (2, 2, 2)])
def test_exact_complete_carrier_certificate(sides):
    result = s.certify(FieldGeometry(sides))
    volume = math.prod(sides)
    assert result['complete_Fock_dimension'] == 4**volume
    assert sum(result['number_sector_dimensions']) == 4**volume
    assert result['complete_CAR_Ward_identities_checked'] == 8*volume
    assert result['all_spatial_bonds_annihilate_entire_eta_tower'] is True


@pytest.mark.parametrize('speed,interaction', [(0., 1.), (.125, 0.), (.25, 2.), (.5, .375)])
def test_same_identities_hold_for_admitted_dyadic_action_parameters(speed, interaction):
    assert s.certify(FieldGeometry((2, 1, 1), speed=speed, interaction=interaction))['spin_eta_algebra_exact']


def test_literal_action_binding_rejects_poisoned_spin_chart(monkeypatch):
    monkeypatch.setattr(lattice, 'SIGMA', tuple(-sigma for sigma in lattice.SIGMA))
    # field.py imported the same tuple by name, so change its consumed chart.
    from model.full_fock_field import field
    monkeypatch.setattr(field, 'SIGMA', lattice.SIGMA)
    with pytest.raises(ArithmeticError, match='actual field action'):
        s.operators(FieldGeometry((2, 2, 1)))


@pytest.mark.parametrize('speed,interaction', [(1e308, 1.), (.125, 1e308), (.1, 1.)])
def test_exact_admission_rejects_unbounded_or_nonbinary_small_denominator_coefficients(speed, interaction):
    with pytest.raises(ValueError, match='(numerator|denominator)'):
        s.operators(FieldGeometry((1, 1, 1), speed=speed, interaction=interaction))


@pytest.mark.parametrize('field,value', [('speed', True), ('speed', float('nan')),
    ('interaction', float('inf')), ('sides', (2., 1, 1)), ('sides', (0, 1, 1))])
def test_mutated_geometry_is_revalidated(field, value):
    geometry = FieldGeometry((2, 1, 1))
    object.__setattr__(geometry, field, value)
    with pytest.raises(ValueError):
        s.operators(geometry)


@pytest.mark.parametrize('site,axis', [((True, 0, 0), 0), ((0, 0, 0), False),
    ((0, 0, 0), 1.), ((-1, 0, 0), 0), ([0, 0, 0], 0), ((4097, 0, 0), 0)])
def test_coordinate_and_axis_domains_are_strict(site, axis):
    with pytest.raises(ValueError):
        s.link_sign(site, axis)


def test_explicit_volume_and_pair_state_domains_are_bounded():
    with pytest.raises(ValueError, match='eight'):
        s.operators(FieldGeometry((3, 3, 1)))
    for angle in (True, float('nan'), float('inf')):
        with pytest.raises(ValueError, match='finite real'):
            s.eta_product_state(theta=angle)
    result = s.eta_product_state()
    with pytest.raises(ValueError):
        result[0] = 1


def test_no_cached_array_metadata_can_poison_a_fresh_frame():
    frame = s.onsite_frame((1, 2, 3))
    # NumPy 2.5 warns about this deliberately hostile metadata mutation.
    # Scope its allowance to the mutation, never production evaluation.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Setting the shape on a NumPy array', category=DeprecationWarning)
        frame.shape = (4,)
    assert s.onsite_frame((1, 2, 3)).shape == (2, 2)
    assert s.frame_certificate()['exact_link_identities'] == 192
