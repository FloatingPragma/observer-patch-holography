"""Independent full-Fock checks of the controlled contact-interacting core."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import eigh

from . import core, lattice


def annihilator(modes, mode):
    out = np.zeros((2**modes, 2**modes))
    for mask in range(2**modes):
        if mask & (1 << mode):
            out[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return out


def full_chain(coupling):
    B = np.array(lattice.jacobi(1), float)
    h = np.kron(B, np.array([[0., 1], [1, 0]]))
    modes = [annihilator(6, k) for k in range(6)]
    H = sum(h[j, k]*modes[j].T@modes[k] for j in range(6) for k in range(6))
    diagonal = [sum((((mask >> (2*x)) & 1)-.5)*(((mask >> (2*x+1)) & 1)-.5)
                    for x in range(3)) for mask in range(64)]
    W = coupling*np.diag(diagonal)
    values, vectors = eigh(h)
    sea = vectors[:, values < -1e-10]
    ground_energy = values[values < -1e-10].sum()
    phi = np.array([1., 2., 1.])/np.sqrt(6)
    zero = np.kron(phi[:, None], np.eye(2))
    columns = []
    for occupation in ((), (0,), (1,), (0, 1)):
        orbitals = np.column_stack([sea, zero[:, list(occupation)]]) if occupation else sea
        vector = np.zeros(64, complex)
        for occupied in combinations(range(6), orbitals.shape[1]):
            vector[sum(1 << j for j in occupied)] = np.linalg.det(orbitals[list(occupied), :])
        columns.append(vector)
    P = np.column_stack(columns)
    return h, H-ground_energy*np.eye(64), W, P, sea


def test_complete_27_site_cube_zero_pair_gap_and_weighted_sea_identity():
    h = np.array(lattice.transformed_box(1), complex)
    values, vectors = eigh(h)
    negative = values < -1e-10
    zeros = np.abs(values) < 1e-10
    assert negative.sum() == 26 and zeros.sum() == 2
    assert np.min(abs(values[~zeros])) > 6/5
    phi = np.array([1., 2., 1.])/np.sqrt(6)
    orbital = np.kron(np.kron(phi, phi), phi)
    pair = np.kron(orbital[:, None], np.eye(2))
    assert np.linalg.norm(h@pair) < 1e-12
    assert np.linalg.norm(pair.conj().T@pair-np.eye(2)) < 1e-12
    ipr = np.sum(orbital**4)
    assert abs(ipr-1/8) < 1e-14
    R = vectors[:, negative]@vectors[:, negative].conj().T
    weighted = sum(orbital[x]**2*R[2*x:2*x+2, 2*x:2*x+2] for x in range(27))
    assert np.linalg.norm(weighted-(1-ipr)*np.eye(2)/2) < 1e-12


def test_full_Fock_compression_includes_the_filled_sea_and_centering_term():
    h, H, W, P, sea = full_chain(1.)
    assert np.linalg.norm(P.conj().T@P-np.eye(4)) < 1e-12
    assert np.linalg.norm(H@P) < 1e-12
    actual = P.conj().T@W@P
    C = actual[0, 0].real
    expected = np.diag([C, C-1/4, C-1/4, C])
    assert np.linalg.norm(actual-expected) < 1e-12
    # A quartic-only projection discards the sea-induced centered core term.
    wrong = np.diag([C, C, C, C+1/2])
    assert np.linalg.norm(actual-wrong) > .5
    assert np.linalg.eigvalsh(H)[4] > 6/5


@pytest.mark.parametrize('coupling', [F(1, 256), F(1, 1024), F(1, 4096)])
def test_actual_all_sector_interacting_spectrum_obeys_bounds_and_exact_symmetries(coupling):
    _, H, W, P, _ = full_chain(float(coupling))
    values = np.linalg.eigvalsh(H+W)
    compressed = np.linalg.eigvalsh(P.conj().T@W@P)
    estimate = core.core_bounds(coupling, 1)
    epsilon = float(estimate['ordered_low_energy_error'])
    assert np.max(abs(values[:4]-compressed)) <= epsilon+1e-12
    assert values[4] >= float(estimate['rest_of_Fock_spectrum_above_shifted_zero'])-1e-12
    sectors = {}
    for number in (2, 3, 4):
        masks = [mask for mask in range(64) if mask.bit_count() == number]
        sectors[number] = np.linalg.eigvalsh((H+W)[np.ix_(masks, masks)])
    assert abs(sectors[3][0]-sectors[3][1]) < 1e-12
    assert abs(sectors[2][0]-sectors[4][0]) < 1e-12
    gap = sectors[2][0]-sectors[3][0]
    assert float(estimate['actual_charged_gap_lower']) <= gap <= float(estimate['actual_charged_gap_upper'])
    assert gap > 0


def test_exact_source_algebra_and_small_cube_interaction_certificate():
    value = core.certificate()
    cube = value['controlled_27_site_cube']
    assert cube['full_Fock_dimension'] == 2**54
    assert F(cube['actual_charged_gap_lower']) > 0
    assert cube['exact_low_Kramers_ground_multiplicity'] == 2
    assert not value['two_to_54_Fock_matrix_executed']
    assert value['supplied_linear_mass_texture_required']


def test_unresolved_level_ordering_does_not_promote_the_Kramers_pair_to_ground():
    value = core.core_bounds(F(1, 100), 3)
    assert not value['charged_gap_certified_positive']
    assert value['exact_low_Kramers_ground_multiplicity'] is None


@pytest.mark.parametrize('coupling,dimension', [(0, 3), (-1, 3), (1, 3), (1., 3), (True, 3), (F(1, 1024), 2)])
def test_invalid_or_uncontrolled_core_parameters_are_rejected(coupling, dimension):
    with pytest.raises(ValueError):
        core.core_bounds(coupling, dimension)
