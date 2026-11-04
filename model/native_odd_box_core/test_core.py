"""Independent source matrices and full-Fock comparisons of the finite theorem."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import eigh

from ..full_fock_field.field import FieldGeometry, FullFock
from . import core


def full_chain(coupling):
    geometry = FieldGeometry((3, 1, 1), speed=1., interaction=float(coupling))
    carrier = FullFock(6)
    h = geometry.single_h().toarray()
    H = carrier.hamiltonian(geometry).toarray()
    numbers = np.array([mask.bit_count() for mask in range(64)])
    H += np.diag(float(coupling)*(-numbers/2+3/4))
    W = np.diag([float(coupling)*sum((((mask >> (2*x)) & 1)-.5)*
                                    (((mask >> (2*x+1)) & 1)-.5)
                                    for x in range(3)) for mask in range(64)])
    values, vectors = eigh(h)
    negative = values < -1e-10
    sea = vectors[:, negative]
    H -= values[negative].sum()*np.eye(64)
    zero = np.kron(np.array([[1.], [0.], [1.]])/np.sqrt(2), np.eye(2))
    columns = []
    for occupation in ((), (0,), (1,), (0, 1)):
        orbitals = np.column_stack([sea, zero[:, list(occupation)]]) if occupation else sea
        vector = np.zeros(64, complex)
        for occupied in combinations(range(6), orbitals.shape[1]):
            vector[sum(1 << j for j in occupied)] = np.linalg.det(orbitals[list(occupied), :])
        columns.append(vector)
    return H, W, np.column_stack(columns), carrier


@pytest.mark.parametrize('length', [3, 5])
def test_entire_native_cube_one_body_spectrum_and_weighted_sea(length):
    h = FieldGeometry((length,)*3, speed=1., interaction=0.).single_h().toarray()
    values, vectors = eigh(h)
    momenta = np.cos(np.pi*np.arange(1, length+1)/(length+1))
    expected = []
    for a in momenta:
        for b in momenta:
            for c in momenta:
                radius = np.sqrt(a*a+b*b+c*c)
                expected.extend([-radius, radius])
    assert np.max(abs(values-np.sort(expected))) < 1e-12
    negative = values < -1e-10
    assert negative.sum() == length**3-1 and sum(abs(values) < 1e-10) == 2
    z = np.array([1. if j % 2 == 0 else 0. for j in range(length)])
    z /= np.linalg.norm(z)
    phi = np.kron(np.kron(z, z), z)
    assert np.linalg.norm(h@np.kron(phi[:, None], np.eye(2))) < 1e-12
    R = vectors[:, negative]@vectors[:, negative].conj().T
    weighted = sum(phi[x]**2*R[2*x:2*x+2, 2*x:2*x+2] for x in range(length**3))
    assert np.linalg.norm(weighted-(1-(2/(length+1))**3)*np.eye(2)/2) < 1e-12


def test_original_full_Fock_centering_and_sea_Wick_compression():
    H, W, P, _ = full_chain(1.)
    assert np.linalg.norm(P.conj().T@P-np.eye(4)) < 1e-12
    assert np.linalg.norm((H-W)@P) < 1e-12
    compression = P.conj().T@W@P
    constant = compression[0, 0]
    assert np.linalg.norm(compression-np.diag([constant, constant-.25, constant-.25, constant])) < 1e-12
    assert np.linalg.norm(compression-np.diag([constant, constant, constant, constant+.5])) > .5
    assert abs(np.linalg.norm(W, 2)-.75) < 1e-12


@pytest.mark.parametrize('coupling', [F(1, 256), F(1, 1024), F(1, 4096)])
def test_all_Fock_levels_sector_order_and_ground_projector_bound(coupling):
    H, W, P, _ = full_chain(coupling)
    bounds = core.core_bounds(coupling, dimensions=1)
    values, vectors = eigh(H)
    compressed = np.linalg.eigvalsh(P.conj().T@W@P)
    epsilon = float(bounds['ordered_low_energy_error'])
    assert max(abs(values[:4]-compressed)) < epsilon+1e-12
    assert values[4] >= float(bounds['rest_of_Fock_spectrum_above_shifted_zero'])-1e-12
    sectors = {}
    for number in (2, 3, 4):
        rows = [mask for mask in range(64) if mask.bit_count() == number]
        sectors[number] = np.linalg.eigvalsh(H[np.ix_(rows, rows)])
    assert abs(sectors[3][0]-sectors[3][1]) < 1e-12
    assert abs(sectors[2][0]-sectors[4][0]) < 1e-12
    gap = sectors[2][0]-sectors[3][0]
    assert float(bounds['actual_charged_gap_lower']) <= gap <= float(bounds['actual_charged_gap_upper'])
    actual = vectors[:, :2]@vectors[:, :2].conj().T
    trial = P[:, 1:3]@P[:, 1:3].conj().T
    assert np.linalg.norm(actual-trial, 2) <= float(bounds['projection_angle_bound'])+1e-12
    assert np.linalg.norm((actual-trial)/2, 'nuc') <= 2*float(bounds['projection_angle_bound'])+1e-12


def test_exact_full_Fock_particle_hole_includes_staggered_sign():
    H, W, _, carrier = full_chain(F(1, 8))
    ladders = [carrier.annihilator(j).toarray() for j in range(6)]
    transformations = []
    for staggered in (True, False):
        ph = np.zeros((64, 64), complex)
        for mask in range(64):
            vector = np.zeros(64, complex)
            vector[63] = 1
            for mode in reversed([j for j in range(6) if mask & (1 << j)]):
                sign = (-1)**(mode % 2)
                if staggered:
                    sign *= (-1)**(mode//2)
                vector = sign*ladders[mode ^ 1]@vector
            ph[:, mask] = vector
        assert np.array_equal(ph.T@ph, np.eye(64))
        transformations.append(ph)
    ph, wrong = transformations
    assert np.max(abs(ph@H@ph.T-H)) < 1e-12
    assert np.linalg.norm(wrong@H@wrong.T-H) > 1
    assert np.array_equal(ph@W@ph.T, W)


def test_rational_cube_and_controlled_family_keep_boundary_scope():
    value = core.certificate()
    cube = value['controlled_27_site_cube']
    assert cube['full_Fock_dimension'] == 2**54
    assert cube['actual_charged_gap_lower'] == '253315/274792972288'
    assert cube['projection_angle_bound'] == '81/524126'
    assert cube['ground_particle_number'] == 27
    assert value['source_anchor']['zero_orbital_central_probability'] == '0'
    assert not value['volume_uniform_nonzero_charged_gap_proved']
    gaps = [F(row['actual_charged_gap_upper']) for row in value['controlled_odd_length_family']]
    assert all(a > b > 0 for a, b in zip(gaps, gaps[1:]))


def test_unresolved_ordering_is_not_declared_ground():
    value = core.core_bounds(F(1, 100))
    assert not value['charged_gap_certified_positive']
    assert value['exact_low_Kramers_ground_multiplicity'] is None
    assert value['ground_particle_number'] is None


@pytest.mark.parametrize('args', [dict(coupling=0), dict(coupling=-1), dict(coupling=.001),
                                 dict(coupling=True), dict(coupling=1), dict(length=4),
                                 dict(length=True), dict(length=17), dict(dimensions=2),
                                 dict(speed=0), dict(speed=1.)])
def test_invalid_or_uncontrolled_source_parameters(args):
    with pytest.raises(ValueError):
        core.core_bounds(**args)


def test_actual_source_link_mutation_is_rejected(monkeypatch):
    original = FieldGeometry.single_h
    def changed(geometry):
        result = original(geometry)
        result.data[0] *= 2
        return result
    monkeypatch.setattr(FieldGeometry, 'single_h', changed)
    with pytest.raises(ArithmeticError, match='links'):
        core.exact_source_anchor()
