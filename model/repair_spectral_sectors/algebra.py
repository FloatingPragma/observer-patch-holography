"""Algebraic eta refinement, finite-prime cyclicity and dark exact spectra.

All exact rank certificates use the actual rational Hubbard-frame generator
and analytic dimer columns in Z[sqrt(17)].  Reduction modulo101 with
sqrt(17)->44 proves a lower rank bound; the exact SU2xSU2 count supplies the
upper bound.  No modular rank deficiency is used as an upper bound.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
import math

import numpy as np
from scipy import sparse
import sympy as sy

from ..common_spin_protection import spectrum, spin
from ..field_symmetry import symmetry
from ..field_stress.energy import GaussianInteger as GI, commutator
from ..observer import native
from ..exchange_sector import statistics
from ..charged_defect import dictionary


PRIME = 101
ROOT17 = 44
PAIRS = tuple(combinations(range(3), 2))
WHERE = spectrum.weight_indices(4, 2)


@lru_cache(maxsize=1)
def exact_data():
    """Return integer16K, raising-spin, lowering-eta and unnormalized seeds."""
    ops = spin.operators()
    if not commutator(ops.K, ops.eta_minus).is_zero():
        raise ArithmeticError('the actual centered field breaks eta lowering covariance')
    if any(not commutator(ops.K, value).is_zero() for value in ops.spin_twice):
        raise ArithmeticError('the actual centered field breaks spin covariance')
    frame = ops.frame
    h = frame@ops.K@frame.adjoint()
    if ops.scale != 16 or h.imag.nnz:
        raise ArithmeticError('the consumed action is not the fixed real16K')
    eta = frame@ops.eta_minus@frame.adjoint()
    if eta.imag.nnz:
        raise ArithmeticError('the consumed eta lowering map is not real integer')
    n2 = spectrum.weight_indices(2, 2)
    n4m2 = spectrum.weight_indices(4, 4)
    lowering = eta.real[n2][:, WHERE].toarray()
    raising = spectrum.raising()[n4m2][:, WHERE].toarray().astype(np.int64)
    g0, g1, v = np.zeros(16, dtype=np.int64), np.zeros(16, dtype=np.int64), np.zeros(16, dtype=np.int64)
    g0[6], g0[9], g0[3], g0[12] = -1, 1, 4, 4
    g1[3] = g1[12] = -1
    v[1], v[4] = 1, -1
    constants, radicals = [], []
    for pair in PAIRS:
        factors0 = [v if j in pair else g0 for j in range(3)]
        factors1 = [v if j in pair else g1 for j in range(3)]
        constants.append(np.kron(factors0[2], np.kron(factors0[1], factors0[0]))[WHERE])
        radicals.append(np.kron(factors1[2], np.kron(factors1[1], factors1[0]))[WHERE])
    c, d = np.column_stack(constants), np.column_stack(radicals)
    if any(np.any(operator@seeds) for operator in (lowering, raising) for seeds in (c, d)):
        raise ArithmeticError('the actual repaired columns fail highest-spin/lowest-eta conditions')
    return h.real[WHERE][:, WHERE].toarray(), raising, lowering, c, d


def modular_rank(columns, prime=PRIME):
    """Exact streaming elimination with a recorded row/column pivot witness."""
    if type(prime) is not int or prime != PRIME:
        raise ValueError('the fixed independently checked prime101 is required')
    value = np.asarray(columns)
    if value.ndim != 2 or value.dtype.kind not in 'iu':
        raise ValueError('integer matrix required')
    pivots, witness = {}, []
    for column, vector in enumerate(value.T):
        reduced = vector.copy() % prime
        for row, basis in pivots.items():
            reduced = (reduced-reduced[row]*basis) % prime
        support = np.flatnonzero(reduced)
        if len(support):
            row = int(support[0]); coefficient = int(reduced[row])
            pivots[row] = (reduced*pow(coefficient, -1, prime)) % prime
            witness.append([row, column, coefficient])
    return len(pivots), witness


def cyclic_certificate():
    h, raising, lowering, c, d = exact_data()
    if not sy.isprime(PRIME) or ROOT17**2 % PRIME != 17 or (4*17) % PRIME == 0:
        raise ArithmeticError('invalid separable prime specialization')
    # Characteristic-zero rank is at least the specialized rank.  The
    # independent spin/eta representation count is exactly90, not inferred
    # from failed pivots or floating eigenvalues.
    expected = [row for row in symmetry.irrep_multiplicities(6)
                if row['twice_spin'] == row['twice_eta'] == 2]
    if expected != [{'twice_spin': 2, 'twice_eta': 2, 'multiplicity': 90}]:
        raise ArithmeticError('the exact joint-irrep count changed')
    constraints = np.vstack((raising, lowering))
    rank, _ = modular_rank(constraints.T)
    if rank != 30:
        raise ArithmeticError('highest-spin/lowest-eta constraints lost rank')
    # 16K and seeds have integer coefficients.  No inversion or normalized
    # irrational column is reduced; only pivot operations are in the field.
    current = (c+ROOT17*d) % PRIME
    layers = []
    for _ in range(43):
        layers.append(current)
        current = (h@current) % PRIME
    krylov = np.column_stack(layers)
    rank, witness = modular_rank(krylov)
    if rank != 90:
        raise ArithmeticError('the old repair support is not cyclic on the90-dimensional active sector')
    # All norm denominators are nonzero at this reduction; normalization is
    # not needed for rank, but this also excludes a vanished projector chart.
    norm = int(np.dot((c+ROOT17*d)[:, 0], (c+ROOT17*d)[:, 0])) % PRIME
    if norm == 0:
        raise ArithmeticError('specialized seed norm vanished')
    # Source ground_column uses b/a=-4t/(sqrt(U²+16t²)+U).
    # At t=1/16,U=1 this is exactly4-sqrt17, so its frame column is
    # proportional to g0+sqrt17*g1.  The actual global frame adds the
    # harmless pair-dependent signs(-,+,-), checked separately below.
    root = sy.sqrt(17)
    if sy.simplify(-sy.Rational(1, 4)/(root/4+1)-(4-root)) != 0:
        raise ArithmeticError('analytic consumed ground ratio is not4-sqrt17')
    normalized = (c+math.sqrt(17)*d)/math.sqrt(272-64*math.sqrt(17))
    actual = spin.operators().frame.matrix()@dictionary.embedded_columns(3, 'hole')[:, [0, 4, 8]]
    target = np.zeros_like(actual)
    target[WHERE] = normalized*np.array([-1, 1, -1])
    seed_residual = float(np.linalg.norm(actual-target))
    if seed_residual > 2e-12:
        raise ArithmeticError('cyclic seeds differ from the actual repair support columns')
    return {'prime': PRIME, 'sqrt17_image': ROOT17,
            'minimal_polynomial': 't^2-17', 'discriminant': 68,
            'discriminant_nonzero_mod_prime': True,
            'generator_integer_scale': 16,
            'rank_inputs_have_no_denominators': True,
            'specialized_unnormalized_seed_norm': norm,
            'exact_unnormalized_seed_norm_squared': '272-64*sqrt(17)',
            'exact_source_ground_doublon_to_singlet_ratio': '4-sqrt(17)',
            'global_frame_seed_signs': [-1, 1, -1],
            'actual_repair_seed_binding_residual': seed_residual,
            'seed_binding_numerical_admission': 2e-12,
            'weight_space_dimension': 120,
            'exact_spin_eta_constraint_rank': 30,
            'exact_joint_irrep_multiplicity': 90,
            'krylov_degrees_included': [0, 42],
            'krylov_columns': 129, 'specialized_rank': rank,
            'pivot_witness': witness,
            'rank_conclusion_over_Q_sqrt17': 90,
            'active_joint_star_algebra': 'M90 tensor I3',
            'active_joint_commutant': 'I90 tensor M3',
            'all105_multiplicities_declared_active': False}


@lru_cache(maxsize=1)
def dark_columns():
    """2D=eta+ applied to every all-up two-particle occupation column."""
    ops = spin.operators()
    plus = (ops.frame@ops.eta_plus@ops.frame.adjoint()).real
    pairs = tuple(combinations(range(6), 2))
    addresses = [sum(1 << (2*x) for x in pair) for pair in pairs]
    return plus[:, addresses].tocsr(), pairs


def dark_certificate():
    numerator, pairs = dark_columns()
    if (numerator.T@numerator-4*sparse.eye(15, dtype=np.int64)).nnz:
        raise ArithmeticError('eta-raised dark columns are not normalized by2')
    adjacency = np.zeros((15, 15), dtype=np.int64)
    index = {pair: j for j, pair in enumerate(pairs)}
    for j, pair in enumerate(pairs):
        for origin in pair:
            for destination in (origin-1, origin+1):
                if 0 <= destination < 6 and destination not in pair:
                    output = tuple(sorted((destination, next(x for x in pair if x != origin))))
                    adjacency[index[output], j] += 1
    ops = spin.operators()
    h = (ops.frame@ops.K@ops.frame.adjoint()).real
    if (h@numerator-numerator@sparse.csr_matrix(adjacency-16*np.eye(15, dtype=np.int64))).nnz:
        raise ArithmeticError('eta-raised exact free-pair energy intertwiner failed')
    _, _, _, c, d = exact_data()
    restricted = numerator[WHERE].toarray()
    if np.any(restricted.T@c) or np.any(restricted.T@d):
        raise ArithmeticError('a charged repair support overlaps the dark eta sector')
    variable = sy.Symbol('y')
    polynomial = sy.Matrix(adjacency).charpoly(variable).as_expr()
    cubics = [variable**3-7*variable-7, variable**3-7*variable+7,
              variable**3-2*variable**2-variable+1, variable**3+2*variable**2-variable-1]
    if sy.expand(polynomial-variable**3*math.prod(cubics)) != 0:
        raise ArithmeticError('exact dark spectral polynomial changed')
    if any(sy.gcd(f, sy.diff(f, variable)) != 1 or f.subs(variable, 0) == 0 for f in cubics):
        raise ArithmeticError('dark nonzero bands are not simple')
    if any(sy.gcd(a, b) != 1 for j, a in enumerate(cubics) for b in cubics[j+1:]):
        raise ArithmeticError('dark polynomial factors share a spectral root')
    energies = sorted((math.cos(p*math.pi/7)+math.cos(q*math.pi/7))/8-1
                      for p, q in combinations(range(1, 7), 2))
    return {'number': 4, 'total_spin': 1, 'total_eta': 2,
            'multiplicity': 15, 'spin_dimension': 3, 'complete_sector_dimension': 45,
            'exact_isometry_denominator': 2,
            'energy_formula': '(cos(p*pi/7)+cos(q*pi/7))/8-1, 1<=p<q<=6',
            'scaled_spectral_variable': 'y=16*(K+1)',
            'integer_characteristic_polynomial': str(sy.factor(polynomial)),
            'distinct_energies': 13, 'multiplicity_of_K_minus_one': 3,
            'other_energy_multiplicities': [1]*12,
            'numerical_display_energies': energies,
            'repair_native_outcomes_identically_zero': True,
            'repair_outside_outcome_is_identity': True,
            'dark_joint_commutant_complex_dimension': 189,
            'active_plus_dark_joint_commutant_complex_dimension': 198,
            'physical_particle_or_mass_identified': False}


def repair_corner_certificate():
    """Exact M3 repair corner, bound to every actual complete field outcome."""
    ports = spin.BLOCK_PORTS
    basis = dictionary.embedded_columns(3, 'hole')[:, [0, 4, 8]]
    signatures = [[] for _ in PAIRS]
    neighbors = {j: set() for j in range(3)}
    count, residual = 0, 0.
    for law in ('C', 'E'):
        for seam in combinations(range(3), 2):
            actual = dict(dictionary.full_extension(3, 'hole', 0, seam, law))
            native_seam = tuple(sorted(ports[j] for j in seam))
            if native_seam not in statistics.SEAMS:
                raise ArithmeticError('the local block edge is not an actual native seam')
            for coin in (0, 1):
                for raw in range(4 if law == 'C' else 8):
                    matrix = np.zeros((3, 3), dtype=np.int64)
                    for source, pair in enumerate(PAIRS):
                        state = tuple(int(j in (ports[pair[0]], ports[pair[1]])) for j in range(12))
                        output, (_, current) = native.native(state, native_seam, coin)
                        record = {0: 0, 1: 1, -1: 2}[current]
                        if law == 'E':
                            record += 4*((state[native_seam[0]]+state[native_seam[1]]) % 2)
                        if raw != record:
                            continue
                        target_pair = tuple(j for j, p in enumerate(ports) if output[p])
                        labels = [seam[1] if j == seam[0] else seam[0] if j == seam[1] else j
                                  for j in pair] if current else list(pair)
                        sign = -1 if labels[0] > labels[1] else 1
                        target = PAIRS.index(target_pair)
                        matrix[target, source] = sign
                        neighbors[source].add(target); neighbors[target].add(source)
                    for j in range(3):
                        signatures[j].append(int(np.vdot(matrix[:, j], matrix[:, j])))
                    operator = actual[('native', coin, raw)]
                    residual = max(residual, float(np.linalg.norm(operator@basis-basis@matrix/math.sqrt(2))))
                    count += 1
            residual = max(residual, float(np.linalg.norm(actual[('outside',)]@basis)))
    if len(set(map(tuple, signatures))) != 3 or any(len(neighbors[j]-{j}) != 2 for j in range(3)):
        raise ArithmeticError('actual complete repair does not generate all three corner matrix units')
    if residual > 2e-12:
        raise ArithmeticError('the exact corner dictionary does not match full-Fock repair')
    return {'actual_native_triangle': list(ports),
            'raw_slots_including_zero': count,
            'effect_signatures': signatures,
            'corner_dimension': 3, 'corner_star_algebra_dimension': 9,
            'all_corner_positions_separated': True, 'corner_transition_graph_complete': True,
            'full_Fock_dictionary_max_residual': residual,
            'numerical_binding_admission': 2e-12,
            'outside_outcome_annihilates_corner': True}
