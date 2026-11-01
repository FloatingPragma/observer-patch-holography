"""Exact interacting singlet spectrum and a fixed-number local density probe."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np
from scipy import sparse
import sympy as sy

from . import cube


Z = sy.Symbol('z')
Q2 = Z**2-16*Z-12
C64 = Z**3-16*Z**2-12*Z+64
C128 = Z**3-16*Z**2-12*Z+128


@lru_cache(maxsize=1)
def singlet_data():
    ops = cube.operators()
    h = cube.restrict(ops.hubbard_H.real, 2).toarray()
    if cube.restrict(ops.hubbard_H.imag, 2).nnz:
        raise ArithmeticError('the actual Hubbard-frame cube action is not real')
    indices = {int(mask): j for j, mask in enumerate(cube.weights(2))}
    basis, locations = [], []
    for x in range(8):
        for y in range(x, 8):
            col = np.zeros(120, dtype=np.int64)
            col[indices[(1 << (2*x)) | (1 << (2*y+1))]] = 1
            if x != y: col[indices[(1 << (2*x+1)) | (1 << (2*y))]] = -1
            basis.append(col); locations.append((x, y))
    v = np.column_stack(basis); gram = np.diag(v.T@v)
    compressed = v.T@h@v
    if np.any(compressed % gram[:, None]):
        raise ArithmeticError('unnormalized singlet basis lost integer action coefficients')
    compressed //= gram[:, None]
    if np.any(h@v-v@compressed):
        raise ArithmeticError('the exact36-dimensional singlet sector is not invariant')
    eta = np.zeros(36, dtype=np.int64)
    for j, (x, y) in enumerate(locations):
        if x == y: eta[j] = (-1)**sum(cube.GEOMETRY.sites[x])
    if np.any(compressed@eta-16*eta):
        raise ArithmeticError('the eta-four parent is not the exact U eigenvector')
    q = np.zeros(36, dtype=np.int64); q[0] = 8; q -= eta
    return h, v, gram, compressed, eta, q


def polynomial_certificate():
    h, v, gram, small, eta, _ = singlet_data()
    actual = sy.Matrix(small).charpoly(Z).as_expr()
    expected = Z**9*(Z-16)*(Z**2-12)**3*Q2*C64**3*C128**3
    if sy.expand(actual-expected) != 0:
        raise ArithmeticError('exact interacting cube singlet polynomial changed')
    highest = [j for j, mask in enumerate(cube.weights(2)) if int(mask) & 0xAAAA == 0]
    triplet = h[np.ix_(highest, highest)]
    if sy.Matrix(triplet).charpoly(Z).as_expr().expand() != (Z**16*(Z**2-12)**6).expand():
        raise ArithmeticError('exact cube triplet polynomial changed')
    factors = (Z, Z**2-12, Q2, C64, C128)
    if any(sy.gcd(p, sy.diff(p, Z)) != 1 for p in factors):
        raise ArithmeticError('a claimed distinct-energy factor is not squarefree')
    if any(sy.gcd(a, b) != 1 for j, a in enumerate(factors) for b in factors[j+1:]):
        raise ArithmeticError('claimed spectral factors have a common energy')
    full = Z**57*(Z**2-12)**21*Q2*C64**3*C128**3
    return {'spectral_variable': 'z=16 E_parent=16(E_N6-2)',
            'complete_spin_singlet_parent_polynomial': str(sy.factor(actual)),
            'eta3_spin_singlet_polynomial': str(sy.factor(actual/(Z-16))),
            'highest_spin_triplet_polynomial': str(sy.factor(Z**16*(Z**2-12)**6)),
            'complete_eta3_119_state_polynomial': str(sy.factor(full)),
            'eta3_distinct_energies': 11,
            'zero_z_multiplicity': 57, 'each_z_plus_minus_sqrt12_multiplicity': 21,
            'each_quadratic_root_multiplicity': 1, 'each_cubic_root_multiplicity': 3,
            'genuinely_contact_shifted_dimensions': 20,
            'spin_alignment_used_for_spectrum': 'H has hidden SU2; repair transparency instead follows eta orthogonality',
            'physical_mass_or_Lorentz_spin_identified': False}


@lru_cache(maxsize=1)
def probe_resolvent():
    _, _, gram, small, _, q = singlet_data()
    denominator = sy.Poly(Q2*C64*C128, Z)
    powers, power = [], q.copy()
    for _ in range(9):
        powers.append(power); power = small@power
    coefficients = [int(x) for x in denominator.all_coeffs()]
    annihilated = sum(coefficient*powers[8-j] for j, coefficient in enumerate(coefficients))
    if np.any(annihilated):
        raise ArithmeticError('claimed probe minimal polynomial does not annihilate its source column')
    krylov = np.column_stack(powers[:8])
    # h^7 q fits int64; its Gram products need arbitrary-precision integers.
    exact_krylov = sy.Matrix(krylov)
    krylov_gram = exact_krylov.T*sy.diag(*map(int, gram))*exact_krylov
    determinant = krylov_gram.det()
    if determinant <= 0:
        raise ArithmeticError('local density probe is not cyclic on all eight interacting energies')
    moments = [sy.Rational(int(q@(gram*power)), 56) for power in powers[:8]]
    numerator = sum(sum(coefficients[k]*moments[j-k] for k in range(j+1))*Z**(7-j) for j in range(8))
    if moments[:3] != [1, 16, sy.Rational(1840, 7)]:
        raise ArithmeticError('the local interacting probe mean/variance changed')
    return sy.factor(numerator), denominator.as_expr(), int(determinant)


def local_probe_certificate():
    ops = cube.operators(); indices = {int(mask): j for j, mask in enumerate(cube.weights(6))}
    reference = np.zeros(8008, dtype=np.int64)
    for sites in combinations(range(8), 3):
        mask = sum(3 << (2*site) for site in sites)
        reference[indices[mask]] = (-1)**sum(sum(cube.GEOMETRY.sites[site]) for site in sites)
    occupation = np.array([(int(mask) & 3).bit_count() for mask in cube.weights(6)], dtype=np.int64)
    probe = (4*occupation-3)*reference
    if int(reference@reference) != 56 or int(probe@probe) != 840 or reference@probe:
        raise ArithmeticError('exact fixed-number reference/probe normalization failed')
    for part, eigenvalue in ((ops.H.real, 48), (ops.H.imag, 0)):
        if np.any(cube.restrict(part, 6)@reference-eigenvalue*reference):
            raise ArithmeticError('the eta-four N6 reference is not an exact energy eigenstate')
    t, q8, _ = cube.parent_lift()
    parent_index = next(j for j, mask in enumerate(cube.weights(2)) if mask == 3)
    if np.any(t@q8[:, parent_index].toarray().ravel()-2*probe):
        raise ArithmeticError('local number probe differs from the eta-raised interacting parent')
    numerator, denominator, determinant = probe_resolvent()
    rows = []
    for factor, multiplicity in ((Q2, 1), (C64, 3), (C128, 3)):
        for (left, right), count in sy.Poly(factor, Z).intervals(eps=sy.Rational(1, 10**14)):
            if count != 1: raise ArithmeticError('probe root is not simple')
            root = sy.N((left+right)/2, 50)
            weight = float(sy.N(numerator.subs(Z, root)/sy.diff(denominator, Z).subs(Z, root), 45))
            if weight <= 0: raise ArithmeticError('positive spectral probe weight was lost numerically')
            rows.append({'factor': str(factor), 'parent_eigenspace_multiplicity': multiplicity,
                         'z_interval': [str(left), str(right)],
                         'N6_energy_interval': [str(2+left/16), str(2+right/16)],
                         'local_probe_frequency_interval': [str(left/16-1), str(right/16-1)],
                         'normalized_probe_weight': weight})
    if abs(math.fsum(row['normalized_probe_weight'] for row in rows)-1) > 2e-12:
        raise ArithmeticError('the complete local probe spectral mass was not retained')
    # Rows retain both local empty/doublon states; columns retain all128
    # exterior pair masks. Unpaired-mode columns vanish on both sides.
    pair_reference = np.zeros(256, dtype=np.int64)
    for mask in range(256):
        if mask.bit_count() == 3:
            pair_reference[mask] = (-1)**sum(sum(cube.GEOMETRY.sites[j]) for j in range(8) if mask >> j & 1)
    coefficient = pair_reference.reshape(128, 2).T
    phased = np.array([1, -1])[:, None]*coefficient
    if np.any(coefficient.T@coefficient-phased.T@phased):
        raise ArithmeticError('the local phase changed the full exterior reduced state')
    if int(np.sum(coefficient*phased)) != 14:
        raise ArithmeticError('the exact local disturbance overlap changed')
    return {'reference_CAR_number': 6, 'reference_eta': 4,
            'reference': '(eta_plus)^3 vacuum / sqrt2016',
            'reference_energy_H': '3', 'reference_support_occupation_columns': 56,
            'local_operator': 'n_up(origin)+n_down(origin)',
            'local_number_mean': '3/4', 'local_number_variance': '15/16',
            'normalized_probe': '(4 n_origin-3)|R>/sqrt15',
            'probe_eta': 3, 'probe_spin': 0,
            'exact_probe_resolvent_numerator': str(numerator),
            'exact_probe_resolvent_denominator': str(denominator),
            'eight_column_Krylov_Gram_determinant': str(determinant),
            'eight_interacting_energy_lines': rows,
            'normalized_probe_mean_H': '3', 'normalized_probe_H_variance': '3/112',
            'number_preserving_phase_probe': 'F=exp(i*pi*n_origin/2); F|R>=(|R>-sqrt15 |probe>)/4',
            'phase_probe_elastic_weight': '1/16', 'phase_probe_inelastic_weight': '15/16',
            'phase_probe_energy_variance': '45/1792',
            'exact_full_exterior_reduced_state_unchanged': True,
            'local_disturbance_pure_state_trace_distance': 'sqrt15/4',
            'fixed_number_Ramsey_formula': 'C(t)=1/16+(15/16) sum_j w_j exp[-it(z_j/16-1)]',
            'charged_vacuum_superposition_required': False,
            'source_preparation_or_detector_executed_here': False,
            'reference_is_physical_ground_state': False}
