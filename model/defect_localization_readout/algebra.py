"""Exact source-branch symmetry and certified local cube excitation."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import numpy as np
import sympy as s

from model.defect_field_bridge import full_domain
from model.centered_source_ground import dimer
from model.native_ordered_ground import source, certificate
from model.centered_particle_spectrum.certificate import decimal_bounds
from model.protected_spin_limit.experiment import full_hamiltonian
from model.strong_coupling_refinement.algebra import geometry, exchange_four

ROOT = Path(__file__).resolve().parents[2]


def rotation(mask, sites):
    """Second quantization of the common internal Pauli X, with CAR signs."""
    if type(sites) is not int or not 1 <= sites <= 8:
        raise ValueError('one through eight sites required')
    if type(mask) is not int or not 0 <= mask < 1 << (2*sites):
        raise ValueError('complete Fock address required')
    out = 0
    sign = 1
    for x in range(sites):
        local = (mask >> (2*x)) & 3
        out |= ((local & 1) << 1 | (local >> 1)) << (2*x)
        if local == 3:
            sign = -sign
    return out, sign


def spin_mask(word, sites):
    return sum(1 << (2*x + ((word >> x) & 1)) for x in range(sites))


def branch_certificate():
    comparisons = 0
    for sites in (2, 3, 4):
        for seam in combinations(range(sites), 2):
            for law in ('C', 'E'):
                buckets = {(coin, raw): [] for coin in (0, 1)
                           for raw in range(full_domain.RAW_SLOTS)}
                for coin in (0, 1):
                    for mask in range(1 << (2*sites)):
                        raw, out, sign, current = full_domain.branch_action(
                            mask, law, coin, sites=sites, flavors=2, seam=seam)
                        other, phase = rotation(mask, sites)
                        raw2, out2, sign2, current2 = full_domain.branch_action(
                            other, law, coin, sites=sites, flavors=2, seam=seam)
                        rotated_out, out_phase = rotation(out, sites)
                        if (raw2, out2, phase*sign2, current2) != (
                                raw, rotated_out, sign*out_phase, current):
                            raise ArithmeticError('a complete raw branch violates spin parity')
                        if mask.bit_count() != out.bit_count():
                            raise ArithmeticError('occupation charge changed')
                        buckets[coin, raw].append(out)
                        comparisons += 1
                if any(len(row) != len(set(row)) for row in buckets.values()):
                    raise ArithmeticError('raw branch ceased to be a partial isometry')
                if sum(map(len, buckets.values())) != 2*(1 << (2*sites)):
                    raise ArithmeticError('complete coin/raw record lost an input')
    protected = 0
    for seam in source.cube()[0]:
        for law in ('C', 'E'):
            for coin in (0, 1):
                for word in range(256):
                    mask = spin_mask(word, 8)
                    if full_domain.branch_action(mask, law, coin, sites=8,
                            flavors=2, seam=seam) != (0, mask, 1, 0):
                        raise ArithmeticError('a cube repair changes protected spin data')
                    protected += 1
    # Anchor the symmetry in the actual centered source Hamiltonian, including
    # all charges. Omitting the minus sign on a doubly occupied site fails.
    r = s.zeros(16)
    for mask in range(16):
        out, sign = rotation(mask, 2)
        r[out, mask] = sign
    actual = dimer.data()['K']
    if r*r != s.eye(16) or r*actual != actual*r:
        raise ArithmeticError('actual full-Fock source parity is not conserved')
    _, signed_edges = geometry((2, 2, 1))
    full = full_hamiltonian(4, signed_edges, 2, 8)
    if not np.array_equal(full, full.astype(np.int64)):
        raise ArithmeticError('integer source witness was rounded')
    rows = [rotation(mask, 4) for mask in range(256)]
    permutation = [row[0] for row in rows]
    signs = np.array([row[1] for row in rows], dtype=np.int64)
    if not np.array_equal(full[np.ix_(permutation, permutation)]*signs[:, None]*signs[None, :], full):
        raise ArithmeticError('the actual signed Hubbard/repair source violates parity')
    return {'all_seam_full_Fock_columns': comparisons,
            'protected_cube_columns': protected,
            'coin_raw_slots_per_seam_and_law': 32,
            'full_Fock_parity_commutes_with_each_record_branch': True,
            'source_dimer_Hamiltonian_commutes_with_parity': True,
            'signed_four_site_full_source_Hamiltonian_commutes_with_parity': True,
            'protected_repairs_are_identity_with_fair_coin_records': True,
            'occupation_charge_changed_by_local_spin_creation': False}


def cube_certificate():
    with np.load(ROOT/'model/native_ordered_ground/witness.npz', allow_pickle=False) as archive:
        witness = {key: archive[key] for key in archive.files}
    checked = certificate.verify(witness)
    h, scale = source.matrix()
    h = h.astype(object)
    vector = witness['vectors'][:, 0].astype(object)
    error = F(checked['normalized_ground_vector_error'])
    size = len(vector)
    flip = np.arange(size) ^ (size-1)
    z = np.array([1-2*(i & 1) for i in range(size)], dtype=object)
    x = np.arange(size) ^ 1
    if not np.array_equal(h[np.ix_(flip, flip)], h):
        raise ArithmeticError('same Hamiltonian broke global spin parity')
    if any(z[flip[i]] != -z[i] for i in range(size)):
        raise ArithmeticError('local creation did not flip parity')
    # All positive Perron coordinates imply that the unique ground is even.
    if not all(v > 0 for v in vector):
        raise ArithmeticError('Perron sign control lost')
    norm = sum(int(v)**2 for v in vector)
    raw_x = F(sum(int(vector[i])*int(vector[x[i]]) for i in range(size)), norm)
    contrast = (raw_x-2*error, raw_x+2*error)
    difference = z[:, None]*h*z[None, :]-h
    energy = F(sum(int(vector[i])*int(difference[i, j])*int(vector[j])
                   for i in range(size) for j in range(size)), scale*norm)
    delta_energy = (energy-2*F(13, 4)*error, energy+2*F(13, 4)*error)
    gap = tuple(F(x) for x in checked['finite_field_gap_interval'])
    if contrast[0] <= F(67, 100) or delta_energy[0] <= gap[0]:
        raise ArithmeticError('strict local signal or positive odd response energy lost')
    # Z cancels from every partial-trace coefficient on the outside algebra.
    exterior_terms = 0
    for a in range(size):
        for b in range(size):
            if (a & 1) == (b & 1):
                if z[a]*z[b] != 1:
                    raise ArithmeticError('creation changes the initial exterior state')
                exterior_terms += 1
    # Removing the selection term preserves parity and leaves the same state.
    field = sum((np.eye(size, dtype=object)[np.arange(size) ^ (1 << site)]
                 for site in range(8)), np.zeros((size, size), dtype=object))
    h0 = h + (scale//8)*field
    if scale % 8 or not np.array_equal(h0[np.ix_(flip, flip)], h0):
        raise ArithmeticError('zero-field source parity or scale changed')
    _, signed_edges = geometry((2, 2, 2))
    physical = exchange_four(8, signed_edges).toarray().astype(object)
    marshall = np.array([(-1)**sum((mask >> x) & 1 for x in source.cube()[1])
                         for mask in range(size)], dtype=object)
    if not np.array_equal(4*h0, scale*marshall[:, None]*physical*marshall[None, :]):
        raise ArithmeticError('the zero-field cube differs from the actual source spin reduction')
    difference0 = z[:, None]*h0*z[None, :]-h0
    e0 = F(sum(int(vector[i])*int(difference0[i, j])*int(vector[j])
               for i in range(size) for j in range(size)), scale*norm)
    return {'sites': 8, 'edges': 12, 'coupling': '1', 'selection_field': '1/4',
            'reference': 'unique actual H_B ground; native preparation has its stated error',
            'ground_vector_error': str(error),
            'creation': 'one-site logical Z at cube vertex 0',
            'reference_global_spin_parity': 1, 'created_global_spin_parity': -1,
            'initial_exterior_partial_trace_coefficients_checked': exterior_terms,
            'initial_exterior_state_unchanged': True,
            'local_X_pointer_probability_contrast': decimal_bounds(*contrast, 12),
            'same_H_B_excitation_energy': decimal_bounds(*delta_energy, 12),
            'same_H_B_odd_sector_gap_lower': str(gap[0]),
            'zero_field_energy_difference_on_same_input': decimal_bounds(e0-6*error, e0+6*error, 12),
            'zero_field_input_claimed_stationary': False,
            'strict_contrast_lower_used': '67/100'}
