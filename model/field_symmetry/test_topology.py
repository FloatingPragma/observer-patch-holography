"""Independent exact controls for the field symmetry's topology and domain.

The small CAR oracles below use only Gaussian integers.  They do not use the
producer's frame, generator construction, spectral vectors or residual gates.
The periodic graphs are explicit countercontrols, not admitted open boxes.
"""
from collections import Counter
from itertools import product

import numpy as np
import pytest
from scipy import sparse
import sympy as sp

from model.full_fock_field.field import FieldGeometry
from model.field_symmetry import symmetry


PAULI = (np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.array([[1, 0], [0, -1]], complex))
IDENTITY = np.eye(2, dtype=complex)


def _zero(matrix):
    value = sparse.csr_matrix(matrix).copy()
    value.eliminate_zeros()
    return value.nnz == 0


def _commutator(left, right):
    return left @ right - right @ left


def _annihilators(volume):
    dimension = 1 << (2*volume)
    result = []
    for mode in range(2*volume):
        columns = [mask for mask in range(dimension) if mask >> mode & 1]
        rows = [mask ^ (1 << mode) for mask in columns]
        values = [(-1)**((mask & ((1 << mode)-1)).bit_count()) for mask in columns]
        result.append(sparse.csr_matrix((values, (rows, columns)),
                                       shape=(dimension, dimension), dtype=complex))
    return tuple(result)


def _oracle(sites, edges):
    """Unscaled hopping, U=1; edges carry their literal 2x2 forward matrix."""
    volume = len(sites)
    ladders = _annihilators(volume)
    dimension = 1 << (2*volume)
    zero = sparse.csr_matrix((dimension, dimension), dtype=complex)
    number = zero.copy()
    contact = zero.copy()
    pair_creators = []
    spin = [zero.copy() for _ in range(3)]
    for site, coordinates in enumerate(sites):
        up, down = ladders[2*site:2*site+2]
        nu, nd = up.getH() @ up, down.getH() @ down
        number += nu+nd
        contact += nu @ nd
        pair_creators.append(up.getH() @ down.getH())
        for axis in range(3):
            sign = (-1)**(sum(coordinates)-coordinates[axis])
            for row, column in product(range(2), repeat=2):
                spin[axis] += sign*PAULI[axis][row, column]*(
                    ladders[2*site+row].getH() @ ladders[2*site+column])
    hopping = zero.copy()
    for left, right, link in edges:
        for row, column in product(range(2), repeat=2):
            term = link[row, column]*(
                ladders[2*left+row].getH() @ ladders[2*right+column])
            hopping += term+term.getH()
    eta = sum(((-1)**sum(site)*pair for site, pair in
               zip(sites, pair_creators, strict=True)), zero.copy())
    return hopping, contact, number, tuple(spin), eta, tuple(pair_creators)


def _open_edges(sites):
    lookup = {site: index for index, site in enumerate(sites)}
    result = []
    for left, site in enumerate(sites):
        for axis in range(3):
            following = list(site)
            following[axis] += 1
            right = lookup.get(tuple(following))
            if right is not None:
                result.append((left, right, 1j*PAULI[axis]))
    return tuple(result)


@pytest.mark.parametrize('sides', [(1, 1, 1), (2, 1, 1), (3, 1, 1), (2, 2, 1)])
def test_full_carrier_centers_distinguish_even_and_odd_site_boxes(sides):
    geometry = FieldGeometry(sides, speed=2., interaction=1.)
    operators = symmetry.operators(geometry)
    sites = geometry.sites
    _, _, number, spin, eta, _ = _oracle(sites, _open_edges(sites))
    assert _zero(operators.N.matrix()-number)
    for actual, expected in zip(operators.spin_twice, spin, strict=True):
        assert _zero(actual.matrix()-expected)
    assert _zero(operators.eta_plus.matrix()-eta)
    n = np.asarray(number.diagonal().real, dtype=int)
    spin_weights = np.asarray(spin[2].diagonal().real, dtype=int)
    eta_weights = np.asarray(operators.eta_z_twice.matrix().diagonal().real, dtype=int)
    spin_center = np.where(spin_weights % 2, -1, 1)
    eta_center = np.where(eta_weights % 2, -1, 1)
    parity = np.where(n % 2, -1, 1)
    assert np.array_equal(spin_center, parity)
    assert np.array_equal(eta_center, (-1)**geometry.volume*parity)
    assert np.array_equal(spin_center*eta_center,
                          np.full(1 << geometry.modes, (-1)**geometry.volume))


@pytest.mark.parametrize('volume', range(1, 6))
def test_irrep_character_matches_every_original_occupation_weight(volume):
    # Enumerate the original CAR basis rather than reusing the producer's
    # binomial highest-weight formula.  The chain's hidden Sz has no stagger.
    original = Counter()
    for mask in range(1 << (2*volume)):
        spin = sum(((mask >> (2*site)) & 1)-((mask >> (2*site+1)) & 1)
                   for site in range(volume))
        original[spin, mask.bit_count()-volume] += 1
    reconstructed = Counter()
    for row in symmetry.irrep_multiplicities(volume):
        spin, eta, multiplicity = (row[key] for key in
                                  ('twice_spin', 'twice_eta', 'multiplicity'))
        assert (spin+eta-volume) % 2 == 0
        for spin_weight in range(-spin, spin+1, 2):
            for eta_weight in range(-eta, eta+1, 2):
                reconstructed[spin_weight, eta_weight] += multiplicity
    assert reconstructed == original


def test_background_pi_flux_survives_endpoint_gauge_changes():
    ux, uy, uz = (1j*matrix for matrix in PAULI)
    links = (ux, uy, ux.conj().T, uy.conj().T)
    assert np.array_equal(np.linalg.multi_dot(links), -IDENTITY)
    # Different endpoint frames telescope around the loop; they cannot change
    # this background into four +I links, regardless of the occupied state.
    frames = (uz, ux, uy, IDENTITY)
    transformed = tuple(frames[index] @ link @ frames[(index+1) % 4].conj().T
                        for index, link in enumerate(links))
    assert np.array_equal(np.linalg.multi_dot(transformed), -IDENTITY)
    assert any(not np.array_equal(link, IDENTITY) for link in transformed)
    sites = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0))
    declared = tuple(symmetry.onsite_frame(site) for site in sites)
    scalar_links = tuple(declared[index] @ link @ declared[(index+1) % 4].conj().T
                         for index, link in enumerate(links))
    assert all(np.array_equal(link, IDENTITY) or np.array_equal(link, -IDENTITY)
               for link in scalar_links)
    assert np.array_equal(np.linalg.multi_dot(scalar_links), -IDENTITY)


@pytest.mark.parametrize('axis', range(3))
def test_odd_periodic_seam_breaks_open_box_generators(axis):
    sites = tuple(tuple(position if coordinate == axis else 0 for coordinate in range(3))
                  for position in range(3))
    link = 1j*PAULI[axis]
    open_edges = ((0, 1, link), (1, 2, link))
    hopping, _, _, spin, eta, _ = _oracle(sites, open_edges)
    assert all(_zero(_commutator(hopping, generator)) for generator in spin)
    assert _zero(_commutator(hopping, eta))
    periodic, _, _, _, _, _ = _oracle(sites, open_edges+((2, 0, link),))
    assert _zero(_commutator(periodic, spin[axis]))
    assert all(not _zero(_commutator(periodic, spin[other]))
               for other in range(3) if other != axis)
    assert not _zero(_commutator(periodic, eta))
    # This is a boundary-condition obstruction, not a ban on all twisted odd
    # rings.  A separately supplied U(2) boundary link restores both here.
    twisted, _, _, _, _, _ = _oracle(sites, open_edges+((2, 0, -1j*IDENTITY),))
    assert all(_zero(_commutator(twisted, generator)) for generator in spin)
    assert _zero(_commutator(twisted, eta))


def test_eta_survives_noncentral_su2_flux_while_the_hidden_spin_does_not():
    sites = tuple(product(range(2), range(2), range(1)))
    edges = list(_open_edges(sites))
    original, _, _, spin, eta, _ = _oracle(sites, edges)
    assert all(_zero(_commutator(original, generator)) for generator in spin)
    left, right, _ = edges[0]
    edges[0] = (left, right, IDENTITY)
    altered, _, _, _, _, _ = _oracle(sites, edges)
    assert _zero(_commutator(altered, eta))
    assert sum(not _zero(_commutator(altered, generator)) for generator in spin) == 2
    ux, uy = (1j*matrix for matrix in PAULI[:2])
    wilson = IDENTITY @ uy @ ux.conj().T @ uy.conj().T
    assert not (np.array_equal(wilson, IDENTITY) or np.array_equal(wilson, -IDENTITY))


def test_only_global_staggered_onsite_pair_profile_commutes_with_hopping():
    sites = tuple(product(range(2), range(2), range(1)))
    hopping, contact, number, spin, eta, pairs = _oracle(sites, _open_edges(sites))
    commutators = [_commutator(hopping, pair).tocoo() for pair in pairs]
    entries = sorted({(int(row), int(column)) for matrix in commutators
                      for row, column in zip(matrix.row, matrix.col, strict=True)})
    dictionaries = [dict(zip(zip(matrix.row, matrix.col, strict=True), matrix.data,
                             strict=True)) for matrix in commutators]
    # Exact integer linear equations, including both real and imaginary parts;
    # no numerical rank threshold supplies the one-dimensional conclusion.
    equations = []
    for entry in entries:
        values = [mapping.get(entry, 0j) for mapping in dictionaries]
        equations.extend(([int(value.real) for value in values],
                          [int(value.imag) for value in values]))
    kernel = sp.Matrix(equations).nullspace()
    staggered = sp.Matrix([(-1)**sum(site) for site in sites])
    assert len(kernel) == 1
    assert kernel[0] == staggered
    assert all(not _zero(matrix) for matrix in commutators)
    assert _zero(_commutator(hopping, eta))
    assert _zero(_commutator(hopping+contact, eta)-eta)
    assert _zero(_commutator(2*(hopping+contact)-number, eta))
    # A local density perturbation preserves charge and spin, but not eta.
    a = _annihilators(len(sites))
    local_number = a[0].getH() @ a[0]+a[1].getH() @ a[1]
    assert _zero(_commutator(local_number, number))
    assert all(_zero(_commutator(local_number, generator)) for generator in spin)
    assert not _zero(_commutator(local_number, eta))
