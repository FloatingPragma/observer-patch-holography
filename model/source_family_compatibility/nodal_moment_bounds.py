"""Full-CAR local constraints on the actual cubic nodal response moment."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import isqrt

import sympy as s

from . import momentum_response as momentum
from .local_matter import annihilator


def _rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational source input required')
    value = F(value)
    if value < 0 or (positive and value == 0):
        raise ValueError('nonnegative input, positive for a source scale, required')
    return value


def sqrt_interval(value, bits=40):
    value = _rational(value)
    if type(bits) is not int or not 1 <= bits <= 128:
        raise ValueError('integer square-root precision in [1,128] required')
    scale = 1 << bits
    low = F(isqrt(value.numerator*scale*scale//value.denominator), scale)
    high = low if low*low == value else low+F(1, scale)
    if not low*low <= value <= high*high:
        raise ArithmeticError('outward square-root enclosure failed')
    return low, high


def source_orbitals():
    """Eight actual orthogonal orbitals supported on thirteen cubic sites."""
    origin = (0, 0, 0)
    near = tuple(tuple(sign if k == axis else 0 for k in range(3))
                 for axis in range(3) for sign in (-1, 1))
    axial = tuple(tuple(2*x for x in y) for y in near)
    sites = (origin,)+near+axial
    # Column1 is the integer numerator of the normalized neighbour orbital.
    numerators = s.zeros(13, 8)
    numerators[0, 0] = 1
    for j in range(1, 7):
        numerators[j, 1] = 1
    for j in range(6):
        numerators[j+7, j+2] = 1
    norm_squares = (1, 6, 1, 1, 1, 1, 1, 1)
    eta = (1, -1, 1, 1, 1, 1, 1, 1)
    if numerators.T*numerators != s.diag(*norm_squares):
        raise ArithmeticError('source orbitals are not orthogonal canonical modes')
    site_eta = s.diag(*((-1)**(sum(x) % 2) for x in sites))
    if numerators.T*site_eta*numerators != s.diag(*(a*b for a, b in zip(norm_squares, eta, strict=True))):
        raise ArithmeticError('coherent neighbour orbital has the wrong eta sign')
    checked = 0
    for perm, signs in momentum.cubic_rotations():
        moved = tuple(momentum.rotate(x, perm, signs) for x in sites)
        if moved[0] != origin or set(moved[1:7]) != set(near) or set(moved[7:]) != set(axial):
            raise ArithmeticError('actual cubic lift does not preserve the source star')
        if any(momentum.rotation_phase(x, perm) != 1 for x in sites):
            raise ArithmeticError('source axial orbitals acquire an unaccounted phase')
        checked += 1
    return {'physical_sites': sites, 'canonical_orbital_count': 8,
            'physical_site_count': 13, 'orbital_integer_norm_squares': norm_squares,
            'orbital_eta_signs': eta, 'proper_cubic_lifts_checked': checked,
            'neighbour_orbital_is_supplied_extra_mode': False,
            'all_six_axial_correlations_equal_requires_stated_cubic_average': True}


def _generators(a, site, epsilon, identity):
    up, down = a[2*site:2*site+2]
    nu, nd = up.H*up, down.H*down
    return ((up.H*down, (nu-nd)/2),
            (epsilon*up.H*down.H, (nu+nd-identity)/2))


def _dot(left, right, dimension):
    return sum((z*w+(p*q.H+p.H*q)/2
                for (p, z), (q, w) in zip(left, right, strict=True)), s.zeros(dimension))


def local_doublet_audit():
    """Exact spin/eta decomposition on every occupation of one orbital."""
    a = tuple(annihilator(1, j) for j in range(2))
    identity = s.eye(4)
    nu, nd = a[0].H*a[0], a[1].H*a[1]
    spin = nu+nd-2*nu*nd
    eta = identity-spin
    if spin*spin != spin or eta*eta != eta or spin*eta != s.zeros(4):
        raise ArithmeticError('local parity projectors are not complementary')
    for sign in (-1, 1):
        groups = _generators(a, 0, sign, identity)
        for (plus, z), projection in zip(groups, (spin, eta), strict=True):
            minus = plus.H
            if plus*minus-minus*plus != 2*z or z*plus-plus*z != plus:
                raise ArithmeticError('source local doublet fails the su2 algebra')
            if z*z+(plus*minus+minus*plus)/2 != s.Rational(3, 4)*projection:
                raise ArithmeticError('spin or eta Casimir does not retain every local sector')
            if plus*projection != plus or z*projection != z:
                raise ArithmeticError('local generator leaks across parity type')
        for x in (groups[0][0], groups[0][0].H, groups[0][1]):
            for y in (groups[1][0], groups[1][0].H, groups[1][1]):
                if x*y != s.zeros(4) or y*x != s.zeros(4):
                    raise ArithmeticError('spin and eta generators are not disjoint doublets')
    return {'local_Fock_dimension': 4, 'spin_doublet_rank': 2, 'eta_doublet_rank': 2,
            'both_staggered_eta_signs_checked': True,
            'occupation_sectors_discarded': False}


def spin_star(weights, up_count=None):
    """Exact matrix4H for H=sum_j weights_j S0.Sj, in a literal spin basis."""
    weights = tuple(weights)
    if any(type(x) not in (int, F) for x in weights):
        raise ValueError('exact rational spin coupling required')
    weights = tuple(F(x) for x in weights)
    n = len(weights)+1
    if up_count is not None and (type(up_count) is not int or not 0 <= up_count <= n):
        raise ValueError('spin magnetization sector is outside its range')
    basis = tuple(x for x in range(1 << n) if up_count is None or x.bit_count() == up_count)
    index = {x: j for j, x in enumerate(basis)}
    matrix = s.zeros(len(basis))
    for x in basis:
        col = index[x]
        for j, weight in enumerate(weights, 1):
            opposite = (x & 1) != ((x >> j) & 1)
            matrix[col, col] += s.Rational(weight.numerator, weight.denominator)*(-1 if opposite else 1)
            if opposite:
                matrix[index[x ^ 1 ^ (1 << j)], col] += 2*s.Rational(weight.numerator, weight.denominator)
    return matrix


def full_CAR_three_orbital_audit(weak_sign=1):
    """All64 CAR states check active/inactive spin/eta blocks with actual signs."""
    if type(weak_sign) is not int or weak_sign not in (-1, 1):
        raise ValueError('weak coupling sign must be minus or plus one')
    a = tuple(annihilator(3, j) for j in range(6))
    identity = s.eye(64)
    groups = [_generators(a, j, sign, identity) for j, sign in enumerate((1, -1, 1))]
    H = 2*_dot(groups[0], groups[1], 64)+weak_sign*_dot(groups[0], groups[2], 64)
    if H != H.H:
        raise ArithmeticError('full-CAR invariant star is not self-adjoint')
    covered = set()
    for types in product(range(2), repeat=3):
        basis = [x for x in range(64)
                 if tuple(((x >> (2*j)) & 3).bit_count() % 2 for j in range(3)) == types]
        if len(basis) != 8 or covered.intersection(basis):
            raise ArithmeticError('full-CAR type blocks do not partition the source')
        covered.update(basis)
        other = [x for x in range(64) if x not in basis]
        if H.extract(basis, other) != s.zeros(8, 56):
            raise ArithmeticError('source star mixes local spin/eta types')
        active = (2 if types[1] == types[0] else 0, weak_sign if types[2] == types[0] else 0)
        block = 4*H.extract(basis, basis)
        if block.charpoly().all_coeffs() != spin_star(active).charpoly().all_coeffs():
            raise ArithmeticError('CAR sign or inactive-doublet spectrum differs')
    if covered != set(range(64)):
        raise ArithmeticError('full-CAR audit omitted an occupation sector')
    full = spin_star((2, weak_sign))
    x = s.Symbol('lambda')
    expected = ((x-3)**4*(x*x+6*x-3)**2 if weak_sign == 1
                else (x-1)**4*(x*x+2*x-27)**2)
    if s.Poly(full.charpoly().as_expr(), x) != s.Poly(expected, x):
        raise ArithmeticError('two-leaf universal bound changed')
    return {'full_Fock_dimension': 64, 'weak_coupling_sign': weak_sign,
            'complete_local_type_blocks_checked': 8,
            'inactive_leaf_multiplicities_retained': True,
            'two_leaf_minimum_exact': ('-(3+2sqrt(3))/4' if weak_sign == 1
                                       else '-(1+2sqrt(7))/4')}


def eight_spin_star_audit(weak_sign=1):
    """Exact256-state magnetization charpolys verify every angular momentum block."""
    if type(weak_sign) is not int or weak_sign not in (-1, 1):
        raise ValueError('weak coupling sign must be minus or plus one')
    factors = Counter()
    dimensions = []
    for up_count in range(9):
        matrix = spin_star((12,)+(weak_sign,)*6, up_count)
        dimensions.append(matrix.rows)
        coefficient, rows = s.factor_list(matrix.charpoly().as_expr())
        if coefficient != 1:
            raise ArithmeticError('spin-star polynomial is not monic')
        for polynomial, multiplicity in rows:
            factors[tuple(int(x) for x in s.Poly(polynomial).all_coeffs())] += int(multiplicity)
    expected = {(1, 36): 5, (1, -12): 15,
                (1, 26, -368): 27, (1, -14): 45, (1, -8): 9,
                (1, 26, -384): 25, (1, -16): 35, (1, -6): 15,
                (1, 26, -408): 7, (1, -18): 9, (1, -4): 5}
    if weak_sign == -1:
        expected = {(1, 36): 5, (1, -12): 15,
                    (1, 22, -512): 27, (1, -10): 45, (1, -16): 9,
                    (1, 22, -528): 25, (1, -8): 35, (1, -18): 15,
                    (1, 22, -552): 7, (1, -6): 9, (1, -20): 5}
    if dict(factors) != expected or sum(dimensions) != 256:
        raise ArithmeticError('literal spin-star spectrum differs from its complete reduction')
    reduced = []
    for j, multiplicity in ((1, 9), (2, 5), (3, 1)):
        off = weak_sign*s.sqrt(j*(j+1))/2
        block = s.Matrix([[-9, off], [off, 3-s.Rational(weak_sign, 2)]])
        polynomial = ([1, 26, -360-4*j*(j+1)] if weak_sign == 1
                      else [1, 22, -504-4*j*(j+1)])
        if (4*block).charpoly().all_coeffs() != polynomial:
            raise ArithmeticError('singlet/triplet mixing has incorrect normalization')
        reduced.append({'weak_total_spin': j, 'weak_spin_multiplicity': multiplicity,
                        'F_equals_j_scaled4_characteristic_coefficients': polynomial,
                        'F_equals_j_plus1_energy': F(3)+F(weak_sign*j, 2),
                        'F_equals_j_minus1_energy': F(3)-F(weak_sign*(j+1), 2)})
    return {'all_active_spin_dimension': 256, 'weak_coupling_sign': weak_sign,
            'magnetization_block_dimensions': dimensions,
            'exact_characteristic_factors_for_4H': [
                {'coefficients': list(k), 'multiplicity': v} for k, v in sorted(expected.items())],
            'weak_spin_zero_multiplicity': 5, 'weak_spin_zero_actual_energies': [F(-9), F(3)],
            'coupled_blocks': reduced,
            'minimum_energy_exact': ('(-13-sqrt(577))/4' if weak_sign == 1
                                     else '(-11-sqrt(673))/4'),
            'all_active_spin_minimum_energy_multiplicity': 7,
            'full_eight_orbital_CAR_dimension': 65536,
            'full_CAR_bound_uses_type_blocks_and_maximally_mixed_extension': True,
            'all_65536_states_diagonalized_as_one_dense_matrix': False}


def source_observable(weak_sign=1):
    """Hstar=12Q0r+sign*sum_axial Q0z in the thirteen-site CAR algebra."""
    if type(weak_sign) is not int or weak_sign not in (-1, 1):
        raise ValueError('weak coupling sign must be minus or plus one')
    source_orbitals()
    qstar = momentum.invariant_observables()['qstar']
    return momentum.add(momentum.scale(qstar, 12),
                        *(momentum.scale(momentum.spin_eta_dot(0, j, 1, 1), weak_sign)
                          for j in range(7, 13)))


def atom_budget(c=F(1), u=F(8), bits=40):
    """Source-derived degree-two upper bound; no positive atom is presumed."""
    c, u = _rational(c, positive=True), _rational(u, positive=True)
    lo, hi = sqrt_interval(577, bits)
    lower, upper = (F(31)+lo)/4, (F(31)+hi)/4
    r_upper = c*c*upper
    alpha_squared = u*u/4
    atom_upper = r_upper/(alpha_squared+r_upper)
    if upper >= 16:
        raise ArithmeticError('admitted source nodal Gram enclosure must be below16')
    return {'c': c, 'U': u, 'sqrt577_interval': (lo, hi),
            'nodal_residual_Gram_over_c_squared_upper_exact': '(31+sqrt(577))/4',
            'nodal_residual_Gram_upper_constant_enclosure': (lower, upper),
            'nodal_residual_Gram_upper': r_upper,
            'all_momentum_Y_Gram_upper': r_upper,
            'all_momentum_residual_Gram_upper': r_upper,
            'original_zero_atom_upper': atom_upper,
            'default_ballistic_all_packet_fixed_energy_limsup_upper': (
                atom_upper if u == 8*c else None),
            'default_ballistic_all_packet_escaping_weight_liminf_lower': (
                1-atom_upper if u == 8*c else None),
            'default_U_over_c_8_atom_upper_exact': '(31+sqrt(577))/(95+sqrt(577))',
            'degree_two_Christoffel_bound_used': True,
            'actual_zero_atom_positive_proved': False,
            'source_eighth_moment_computed': False,
            'zero_atom_to_shrinking_nodal_window_uses_fixed_reference_and_homogeneity': True,
            'physical_species_identified': [], 'physical_mass_prediction': None}


def certificate():
    orbitals = source_orbitals()
    doublets = local_doublet_audit()
    full = [full_CAR_three_orbital_audit(sign) for sign in (1, -1)]
    stars = [eight_spin_star_audit(sign) for sign in (1, -1)]
    readouts = [momentum.score_plan(source_observable(sign)) for sign in (1, -1)]
    if any(row['score_absolute_bound'] != F(117, 2) for row in readouts):
        raise ArithmeticError('actual source star readout normalization changed')
    plus_low, plus_high = sqrt_interval(577)
    minus_low, minus_high = sqrt_interval(673)
    gamma_plus = (F(31)+plus_low)/4
    gamma_minus = (F(29)+minus_high)/4
    if gamma_minus >= gamma_plus or plus_high >= 25:
        raise ArithmeticError('paired endpoint bounds do not imply the admitted global ceiling')
    return {'same_centered_cubic_full_CAR_source': True,
            'source_orbitals': orbitals, 'local_spin_eta_doublets': doublets,
            'full_CAR_sign_and_inactive_leaf_controls': full,
            'exact_complete_spin_stars': stars,
            'source_Hplus_Hminus_complete_Pauli_score_plans': readouts,
            'source_moment_budget': atom_budget(),
            'source_endpoint_expectation_identities': [
                'G(S=0)/c_squared=9/2-omega(Hplus)',
                'G(S=3)/c_squared=9/2-omega(Hminus)'],
            'plus_endpoint_upper_exact': '(31+sqrt(577))/4',
            'minus_endpoint_upper_exact': '(29+sqrt(673))/4',
            'source_global_Y_Gram_over_c_squared_strict_upper': F(14),
            'global_default_fixed_window_limsup_strict_upper': F(7, 15),
            'global_default_escape_liminf_strict_lower': F(8, 15),
            'source_cubic_average_consumed': True,
            'spin_only_ground_or_antiferromagnetic_order_assumed': False,
            'universal_local_bound_attained_in_the_actual_infinite_ground': False,
            'native_large_tape_executed': False,
            'physical_species_identified': [], 'physical_mass_prediction': None}
