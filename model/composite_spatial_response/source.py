"""Exact signed spatial covariance, source commutators and multiplet maps."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import numpy as np
from scipy import sparse

from ..charged_composite_search import source as previous
from ..coarse_energy_ward import ward

parent = previous.parent


def probe(site):
    if type(site) is not int or not 0 <= site < 8:
        raise ValueError('one of the eight source sites is required')
    return _probe(site)


@lru_cache(maxsize=8)
def _probe(site):
    adjacency = parent.ground.band_data()[0]
    neighbors = tuple(int(y) for y in np.flatnonzero(adjacency[site]))
    remove = parent.car(4, 4, 2*site+1)
    result = sum(parent.car(4, 2, 2*y, True) @
                 parent.car(4, 3, 2*y+1) @ remove for y in neighbors)
    if site == 0 and (result-previous.probe()).nnz:
        raise ArithmeticError('origin probe differs from the certified source')
    return result.tocsr()


def fock_matrix(up, down, axes, shift):
    rows, signs = parent.fock_action(up, down, axes, shift)
    return sparse.csr_matrix((signs, (rows, np.arange(len(rows)))),
                             shape=(len(rows), len(rows)), dtype=np.int64)


def covariance_certificate():
    count = 0
    for axes, shift, permutation, signs in parent.automorphisms():
        rg = fock_matrix(4, 4, axes, shift)
        rc = fock_matrix(5, 2, axes, shift)
        for x in range(8):
            if (rc @ probe(x)-int(signs[x])*probe(int(permutation[x])) @ rg).nnz:
                raise ArithmeticError('full-occupation composite spatial covariance failed')
            count += 1
    parent.shift_gram_signs()
    return {'complete_operator_covariance_checks': count,
            'input_occupation_rows': len(parent.masks(4, 4)),
            'output_occupation_rows': len(parent.masks(5, 2)),
            'normalized_signed_spatial_lifts': 48,
            'projective_translation_signs': parent.shift_gram_signs().tolist(),
            'ordinary_continuous_momentum_inferred': False}


def spatial_commutant_certificate():
    """The complete signed site action has commutant span_C{I,A}."""
    identity = np.eye(8, dtype=np.int64)
    adjacency = parent.ground.band_data()[0]
    matrices = []
    for _, _, permutation, signs in parent.automorphisms():
        r = np.zeros((8, 8), dtype=np.int64)
        r[permutation, np.arange(8)] = signs
        if not np.array_equal(r @ adjacency, adjacency @ r):
            raise ArithmeticError('source hopping does not commute with its site action')
        matrices.append(r)
    equations = np.vstack([np.kron(identity, r)-np.kron(r.T, identity) for r in matrices])
    kernel = np.column_stack((identity.ravel(order='F'), adjacency.ravel(order='F')))
    if np.any(equations @ kernel) or not np.array_equal(adjacency @ adjacency, 3*identity):
        raise ArithmeticError('source kinetic commutant basis identity failed')
    prime = 1009
    if any(prime % d == 0 for d in range(2, isqrt(prime)+1)):
        raise ArithmeticError('commutant rank modulus is not prime')
    rank = parent.modular_rank(equations, prime)
    independent = parent.modular_rank(kernel, prime)
    if rank != 62 or independent != 2:
        raise ArithmeticError('source site commutant not completely exhausted')
    return {'complex_matrix_unknowns': 64, 'commutator_equations': int(equations.shape[0]),
            'prime': prime, 'exact_commutator_rank': rank, 'independent_exact_kernel': independent,
            'kernel_basis': ['I8', 'actual signed hopping matrix A'],
            'A_squared': '3I8', 'source_band_projectors': '(I8 +/- A/sqrt(3))/2',
            'source_band_ranks': [4, 4], 'actual_A': adjacency.tolist()}


def multiplet_certificate(block):
    full = block['full_basis']
    current = full
    up, down = 5, 2
    cumulative = 1
    rows = []
    for factor in (3, 4, 3):
        lower = parent.ladder(up-1, down+1, 'spin_plus').T.tocsr()
        if (parent.centered(up-1, down+1) @ lower-
                lower @ parent.centered(up, down)).nnz:
            raise ArithmeticError('spin-ladder source intertwiner failed')
        if np.any(lower.T @ (lower @ current)-factor*current):
            raise ArithmeticError('spin-three-halves ladder normalization failed')
        current = lower @ current
        up -= 1
        down += 1
        cumulative *= factor
        rows.append({'up': up, 'down': down, 'step_norm_squared': factor,
                     'cumulative_norm_squared': cumulative})
    raising = parent.ladder(6, 3, 'eta_minus').T.tocsr()
    if (parent.centered(6, 3) @ raising-raising @ parent.centered(5, 2)).nnz:
        raise ArithmeticError('eta-raising source intertwiner failed')
    if np.any(raising.T @ (raising @ full)-full):
        raise ArithmeticError('eta-doublet raising is not isometric')
    return {'spin_dimension': 4, 'spin_lowering': rows,
            'eta_dimension': 2, 'eta_raise_N7_to_N9_norm_squared': 1,
            'same_centered_energy_on_every_multiplet_member': True,
            'eta_raised_probe_is_bare_C_adjoint': False}


@lru_cache(maxsize=1)
def clifford_source():
    a = tuple(ward.annihilator(p) for p in range(16))
    adjacency = parent.ground.band_data()[0]
    kinetic = ward.total((a[2*x+s].adjoint() @ a[2*y+s]).scale(F(int(adjacency[x,y]), 16))
                         for x in range(8) for y in range(8) for s in (0, 1)
                         if adjacency[x,y])
    number = ward.total(z.adjoint() @ z for z in a)
    contact = ward.total((a[2*x].adjoint() @ a[2*x]) @
                         (a[2*x+1].adjoint() @ a[2*x+1]) for x in range(8))
    return a, kinetic, contact, number


def commutator_certificate():
    a, kinetic, contact, number = clifford_source()
    rows = []
    adjacency = parent.ground.band_data()[0]
    for x in range(8):
        c = ward.total((a[2*int(y)].adjoint() @ a[2*int(y)+1]) @ a[2*x+1]
                       for y in np.flatnonzero(adjacency[x]))
        nup = a[2*x].adjoint() @ a[2*x]
        forced = (ward.Clifford.scalar(F(1, 2))-nup) @ c
        tk = ward.commutator(kinetic, c)
        if ward.commutator(number, c) != -c:
            raise ArithmeticError('composite charge commutator failed')
        if ward.commutator(contact, c) != -(nup @ c):
            raise ArithmeticError('contact force differs from its full-CAR formula')
        if forced.adjoint() @ forced != (c.adjoint() @ c).scale(F(1, 4)):
            raise ArithmeticError('centered contact rate norm identity failed')
        if ward.commutator(kinetic+contact-number.scale(F(1, 2)), c) != tk+forced:
            raise ArithmeticError('kinetic/contact/chemical rate decomposition failed')
        if tk.is_zero() or forced.is_zero():
            raise ArithmeticError('a nonzero rate component was dropped')
        rows.append({'site': x, 'probe_Clifford_terms': len(c.terms),
                     'kinetic_rate_Clifford_terms': len(tk.terms),
                     'centered_contact_rate_Clifford_terms': len(forced.terms),
                     'kinetic_rate_norm_bound': str(tk.norm_bound()),
                     'centered_contact_rate_norm_bound': str(forced.norm_bound()),
                     'all_full_Fock_operator_coefficients_checked': True})
    return {'dimension': 65536, 'site_checks': rows,
            'equation': '[K,C_x]=[T,C_x]+(1/2-n_xup)C_x',
            'chemical_number_change': '-1',
            'current_is_charged_composite_rate_not_conserved_flux': True,
            'centered_contact_rate_norm_squared_is_probe_norm_squared_over_four': True,
            'contact_force_identically_zero': False}


def kinetic(up, down):
    ids = parent.masks(up, down)
    doublons = np.array([sum(((int(m) >> (2*x)) & 3) == 3 for x in range(8))
                         for m in ids], dtype=np.int64)
    diagonal = 16*doublons-8*(up+down)
    return parent.centered(up, down)-sparse.diags(diagonal, format='csr', dtype=np.int64)


def kinetic_source_anchor():
    """Rebuild the same sector kinetic matrices from their one-body CAR law."""
    adjacency = parent.ground.band_data()[0]
    rows = []
    for up, down in ((4, 4), (5, 2)):
        n = len(parent.masks(up, down))
        reconstructed = sparse.csr_matrix((n, n), dtype=np.int64)
        for x in range(8):
            for y in np.flatnonzero(adjacency[x]):
                for spin in (0, 1):
                    target = (up-(spin == 0), down-(spin == 1))
                    term = parent.car(*target, 2*x+spin, True) @ parent.car(up, down, 2*int(y)+spin)
                    reconstructed += int(adjacency[x,y])*term
        if (reconstructed-kinetic(up, down)).nnz:
            raise ArithmeticError('CAR kinetic law differs from the actual complete sector source')
        rows.append({'up': up, 'down': down, 'complete_source_rows': n,
                     'nonzero_kinetic_entries': reconstructed.nnz})
    return rows


def direct_rate_maps(site=0):
    c = probe(site)
    kg, kc = parent.centered(4, 4), parent.centered(5, 2)
    full = kc @ c-c @ kg
    tg, tc = kinetic(4, 4), kinetic(5, 2)
    free = tc @ c-c @ tg
    nup = np.array([(int(mask) >> (2*site)) & 1 for mask in parent.masks(5, 2)], dtype=np.int64)
    force = sparse.diags(8-16*nup, format='csr', dtype=np.int64) @ c
    if (full-free-force).nnz:
        raise ArithmeticError('direct full-occupation rate maps disagree with Clifford identity')
    return {'probe': c, 'full_rate16': full, 'kinetic_rate16': free,
            'centered_contact_rate16': force}
