"""A biased native interaction selects an approximate Bell reference.

The Hamiltonian is the actual four-carrier exchange Gram interaction plus
the actual native near-fivefold class sum on each geometric object. Its
finite-bias ground is derived, not a prescribed pure-triplet reference.
An exact Schur bound proves a uniform projected gap for every rational
bias >=1. A separate finite penalty bound retains all1296 color states.

The phase bias, geometric protection and cooling remain added architecture.
These statements concern the specified effective color Hamiltonian; they
do not replace it by recorded C/E, certify an autonomous bath, or assert a
complete48-port finite-pulse implementation of the virtual interaction.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as s

from . import native_persistent_response as native

ROOT = Path(__file__).resolve().parents[2]
DOMAIN = s.QQ.algebraic_field(s.sqrt(5))


def _dm(matrix):
    return s.Matrix(matrix).to_DM(extension=True).convert_to(DOMAIN)


def _zero(matrix):
    return _dm(matrix).is_zero_matrix


def _entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


def _rational(value, name):
    """Detach and normalize original exact integers before any comparison.

    Fraction slots are writable. Its comparison and copying methods assume
    a positive integer denominator, so neither is a validation boundary.
    """
    if type(value) is int:
        return F(value)
    if type(value) is not F:
        raise ValueError('an original exact rational '+name+' is required')
    numerator, denominator = value.numerator, value.denominator
    if type(numerator) is not int or type(denominator) is not int or denominator <= 0:
        raise ValueError('valid integer Fraction slots and a positive denominator are required')
    return F(numerator, denominator)


def _bias(value):
    value = _rational(value, 'native bias')
    if value < 1:
        raise ValueError('an original exact rational native bias at least one is required')
    return value


@lru_cache(maxsize=1)
def _native_data():
    leakage, B = native.composite_interaction()
    coordinates = native.geometric_coordinates()
    metric = s.kronecker_product(coordinates['metric'], coordinates['metric'])
    if not _zero(B.T*metric-metric*B):
        raise ArithmeticError('actual four-carrier Gram lost metric self-adjointness')
    p = (_dm(coordinates['left_inverse'])*_dm(native.objects.geometric_projectors()['three'])
         *_dm(coordinates['frame'])).to_Matrix()
    q = s.eye(6)-p
    P = _dm(s.kronecker_product(q, q)); Bd = _dm(B)
    A = P*Bd*P
    eigenvalues = (s.Rational(92, 5), s.Rational(1364, 75), s.Rational(160, 9))
    polynomial = P
    for value in eigenvalues:
        polynomial = polynomial*(A-P*DOMAIN.convert(value))
    if not polynomial.is_zero_matrix:
        raise ArithmeticError('actual BB block has a different spectral polynomial')
    ranks = []
    for value in eigenvalues:
        projector = P
        for other in eigenvalues:
            if other != value:
                projector = projector*(A-P*DOMAIN.convert(other))*DOMAIN.convert(1/(value-other))
        if not (projector*projector-projector).is_zero_matrix:
            raise ArithmeticError('BB spectral projector is not idempotent')
        ranks.append(s.expand(s.trace(projector.to_Matrix())))
    if ranks != [1, 5, 3]:
        raise ArithmeticError('BB singlet, quintet and triplet multiplicities differ')
    off_square = P*Bd*Bd*P-A*A
    off_frobenius_squared = s.expand(s.trace(off_square.to_Matrix()))
    if off_frobenius_squared != s.Rational(1519, 3375):
        raise ArithmeticError('complete off-BB response norm differs')
    x = s.Symbol('x')
    polynomial = B.charpoly(x).as_expr()
    roots = s.roots(polynomial, x)
    beta = s.Rational(166, 9)
    if sum(roots.values()) != 36 or any(s.simplify(root >= 0) is not s.true
            or s.simplify(root <= beta) is not s.true for root in roots):
        raise ArithmeticError('complete projected36 Gram norm is not bounded by166/9')
    # Reconstruct the bias on the full36 color carrier, before restriction.
    class_sum = s.zeros(36)
    class_count = 0
    for row in native.objects.geometric_projectors()['rotation_rows']:
        if row['order'] == 5 and s.sympify(row['character3']) > 0:
            rotation = native.objects.permutation_matrix(row['axis_action'])
            class_sum += s.kronecker_product(rotation, rotation)
            class_count += 1
    geometric = native.objects.geometric_projectors()['six']
    reduced = (_dm(coordinates['left_inverse'])*_dm(class_sum)*_dm(coordinates['frame'])).to_Matrix()
    expected = (2+2*s.sqrt(5))*p+(2-2*s.sqrt(5))*q
    if (class_count != 12 or class_sum != class_sum.T or not _zero(reduced-expected)
            or not _zero(class_sum*geometric-geometric*class_sum)):
        raise ArithmeticError('the bias is not the actual native near-fivefold class sum')
    inv_metric = coordinates['metric'].inv()
    singlets = s.Matrix.hstack((p*inv_metric).reshape(36, 1), (q*inv_metric).reshape(36, 1))
    bell = s.Matrix([[s.Rational(92, 5), -s.Rational(2, 45)],
                     [-s.Rational(2, 45), s.Rational(92, 5)]])
    if (not _zero(singlets.T*metric*singlets-3*s.eye(2))
            or not _zero(B*singlets-singlets*bell)):
        raise ArithmeticError('the biased Bell block is not an invariant actual source subspace')
    bounds = []
    for colors in ((0, 0), (0, 1), (1, 0), (1, 1)):
        interaction = native.native_interaction(*colors)
        if interaction != interaction.T:
            raise ArithmeticError('a native cross-pair color Gram is not self-adjoint')
        bounds.append(max(sum(abs(value) for value in interaction.row(i)) for i in range(36)))
    if sum(bounds) != s.Rational(260, 9):
        raise ArithmeticError('full1296 color interaction bound changed')
    # Remove the invariant two-singlet plane. On its BB part the energy
    # is at least -a+g; outside BB the bias costs delta and B<=beta.
    # The Schur complement above energy -a+g/2 is positive if
    # delta > beta-a+g/2+2*ell/g. For bias>=1, delta=4sqrt5*bias>8.
    a, g, ell = s.Rational(92, 5), s.Rational(16, 75), off_frobenius_squared
    threshold = beta-a+g/2+2*ell/g
    if not threshold < 8:
        raise ArithmeticError('uniform finite-bias gap comparison failed')
    return {'B': s.ImmutableMatrix(B), 'metric': s.ImmutableMatrix(metric),
        'first_projector': s.ImmutableMatrix(p), 'second_projector': s.ImmutableMatrix(q),
        'BB_projector': s.ImmutableMatrix(P.to_Matrix()),
        'class_sum': s.ImmutableMatrix(class_sum), 'reduced_class_sum': s.ImmutableMatrix(reduced),
        'singlets': s.ImmutableMatrix(singlets), 'BB_spectrum': tuple(zip(eigenvalues, ranks, strict=True)),
        'B_norm_upper': beta, 'off_BB_Frobenius_squared': ell, 'schur_threshold': threshold,
        'full_color_norm_upper': sum(bounds), 'cross_pair_norm_bounds': tuple(bounds),
        'native_geometric_leakage_witness': leakage['leakage_norm_squared_exact'],
        'characteristic_factors': tuple((str(factor), multiplicity) for factor, multiplicity
                                         in s.factor_list(polynomial)[1])}


@lru_cache(maxsize=32)
def _exact_data(bias):
    source = _native_data()
    lam = s.Rational(bias.numerator, bias.denominator)
    b, delta = s.Rational(2, 45), 4*s.sqrt(5)*lam
    radius = s.sqrt(delta**2+b*b)
    norm = s.sqrt(b*b+(delta+radius)**2)
    first, second = -b/norm, (delta+radius)/norm
    # Stable equivalent form, including for very large rational biases.
    p = b*b/(2*radius*(radius+delta))
    q = 1-p
    block = s.Matrix([[delta, b], [b, -delta]])
    vector = s.Matrix([-b, delta+radius])
    if any(s.simplify(value) != 0 for value in block*vector+radius*vector):
        raise ArithmeticError('finite-bias ground vector is not an exact Bell eigenvector')
    geometric = s.zeros(36, 1); padded = s.zeros(64, 1)
    for j in range(6):
        amplitude = (first if j < 3 else second)/s.sqrt(3)
        geometric[6*j+j] = amplitude
        padded[8*j+j] = amplitude
    C = source['reduced_class_sum']
    H = -source['B']+lam*(s.kronecker_product(C, s.eye(6))+s.kronecker_product(s.eye(6), C))
    ground_energy = 4*lam-s.Rational(92, 5)-radius
    return {'bias': bias, 'first_triplet_weight': p, 'second_triplet_weight': q,
        'ground_geometric': s.ImmutableMatrix(geometric), 'ground_padded': s.ImmutableMatrix(padded),
        'ground_native_coordinates': s.ImmutableMatrix(source['singlets']*s.Matrix([first, second])/s.sqrt(3)),
        'ground_energy': ground_energy, 'projected_gap_lower_bound': s.Rational(8, 75),
        'hamiltonian36': s.ImmutableMatrix(H), 'metric36': source['metric'],
        'Bell_gap': 2*radius, 'bias_delta': delta, 'mixing_b': b}


def exact_data(bias=F(1)):
    """Detached mapping of immutable exact matrices on stated native charts."""
    result = dict(_exact_data(_bias(bias)))
    # Do not expose the mutable Fraction stored in the cached record/key.
    value = result['bias']
    result['bias'] = F(value.numerator, value.denominator)
    return result


def ground64(bias=F(1)):
    exact = exact_data(bias)['ground_padded']
    values = np.asarray(exact.evalf(32), float).reshape(64)
    if not np.isfinite(values).all() or abs(np.linalg.norm(values)-1) > 2e-14:
        raise ArithmeticError('finite source-derived padded ground vector required')
    if any(value == 0 and source != 0 for source, value in zip(exact, values, strict=True)):
        raise ArithmeticError('binary64 export cannot retain a nonzero exact ground component')
    return np.frombuffer(values.tobytes(), dtype=np.float64)


def protection_bound(bias=F(1), penalty=2**30):
    """A full1296 color-space bound, with no geometric postselection.

    H_U=U(I-P)+H_base, P=P_geometric tensor P_geometric. The bias commutes
    with P, so only the actual exchange couples P and its complement.
    A Feshbach/Schur bound shifts each of the36 low eigenvalues by at most
    eps=beta_full^2/(U-2h), h=beta_full+24*bias. The full ground has outside
    component <=beta_full/(U-2h); its normalized inside component differs
    from the projected ground by at most eps/gamma in sine distance.

    The damage-Gram comparison requires Hermitian collective extensions
    L_X J = J L_X_geom for the geometric isometric embedding J, with
    ||L_X||<=2 for every unit Hilbert--Schmidt current. This intertwining
    preserves the geometric subspace and its variance; equality of the
    compression J^* L_X J alone does not preserve the second moment.
    """
    bias = _bias(bias)
    U = _rational(penalty, 'protection penalty')
    if U <= 0:
        raise ValueError('an original positive exact rational protection penalty required')
    beta, gamma = F(260, 9), F(8, 75)
    h = beta+24*bias
    admitted = U > 2*h
    if admitted:
        denominator = U-2*h
        epsilon = beta*beta/denominator
        outside = beta/denominator
        distance = outside+epsilon/gamma
        gap = gamma-epsilon
        admitted = gap > 0
    else:
        epsilon = outside = distance = gap = None
    return {'bias': str(bias), 'penalty': str(U), 'admitted': admitted,
        'full_color_dimension': 1296, 'retained_low_band_dimension': 36,
        'full_exchange_operator_norm_upper': str(beta), 'base_Hamiltonian_operator_norm_upper': str(h),
        'projected_gap_lower_bound': str(gamma),
        'low_energy_error_upper': str(epsilon) if epsilon is not None else None,
        'outside_geometric_ground_amplitude_upper': str(outside) if outside is not None else None,
        'full_ground_trace_distance_to_projected_ground_upper': str(distance) if distance is not None else None,
        'full_ground_gap_lower_bound': str(gap) if gap is not None else None,
        'uniform_time_trace_distance_from_prepared_projected_reference_upper': str(2*distance) if distance is not None else None,
        'full_ground_damage_Gram_operator_error_upper': str(24*distance) if distance is not None else None,
        'uniform_time_damage_Gram_operator_error_upper': str(48*distance) if distance is not None else None,
        'stability_argument': 'The full ground is stationary. Unitary invariance and the triangle inequality give D(exp(-itH_U)psi,psi)<=2d for all realt when D(psi,ground)<=d; no phase comparison or postselection is used.',
        'damage_Gram_argument': 'For HS-normalized currentX and a stated full-space collective extension with norm<=2, |Var_rho(L_X)-Var_sigma(L_X)|<=6 D(rho,sigma)||L_X||^2<=24 D(rho,sigma); taking the quadratic-form supremum bounds the Gram operator.',
        'damage_Gram_extension_scope': 'Hermitian collective extensions satisfy L_X J=J L_X_geom for the geometric isometric embedding J, with norm<=2 per unit Hilbert--Schmidt current; compression equality alone is insufficient. This bound does not itself execute an extension or a protective Hamiltonian.',
        'projection_or_postselection_applied': False,
        'geometric_penalty_is_an_added_service': True,
        'domain': 'specified full1296 corrected-color Hamiltonian; each pair Gram used full10626 native virtual outputs',
        'full48_port_finite_source_pulse_error_certified': False}


def report():
    source = _native_data()
    examples = []
    for bias in (F(1), F(16)):
        data = exact_data(bias)
        examples.append({'bias': str(bias), 'first_triplet_weight_exact': str(data['first_triplet_weight']),
            'second_triplet_weight_exact': str(data['second_triplet_weight']),
            'first_triplet_weight_diagnostic': float(data['first_triplet_weight']),
            'ground_energy_exact': str(data['ground_energy']),
            'projected_gap_lower_bound': str(data['projected_gap_lower_bound']),
            'ground_geometric_exact': [str(x) for x in data['ground_geometric']],
            'protection': protection_bound(bias)})
    return {'schema': 'oph.native_response_biased.v1',
        'Hamiltonian': '-B_four_carrier + bias*(C_native tensor I + I tensor C_native)',
        'class_sum_source': 'sum of the twelve actual native near-fivefold port rotations on each36-color object',
        'geometric_class_sum_values': ['2+2*sqrt(5)', '2-2*sqrt(5)'],
        'BB_spectral_projector_certificate': [{'eigenvalue': str(value), 'rank': int(rank)}
            for value, rank in source['BB_spectrum']],
        'complete_projected36_Gram_characteristic_factors': [list(row) for row in source['characteristic_factors']],
        'complete_projected36_Gram_operator_norm_upper': str(source['B_norm_upper']),
        'off_BB_operator_norm_squared_upper': str(source['off_BB_Frobenius_squared']),
        'sufficient_Schur_delta_threshold': str(source['schur_threshold']),
        'uniform_parameter_domain': 'original exact rational bias>=1',
        'uniform_projected_ground_gap_lower_bound': '8/75',
        'ground_formula': '-sqrt(p)*Omega_A+sqrt(1-p)*Omega_B; p=b^2/[2*r*(r+delta)], b=2/45, delta=4*sqrt(5)*bias, r=sqrt(delta^2+b^2)',
        'exact_pure_B_ground_at_any_finite_bias': False,
        'full_color_cross_pair_norm_bounds': [str(value) for value in source['cross_pair_norm_bounds']],
        'actual_unprotected_geometric_leakage_witness_squared': source['native_geometric_leakage_witness'],
        'examples': examples,
        'scope': {'ground_derived_from_actual_source_interaction_and_declared_bias': True,
            'empirical_mass_or_SM_target_input_used': False,
            'bias_protection_or_bath_selected_by_CE_repair': False,
            'finite_bias_reference_has_exact_twelve_observable_stabilizer': False,
            'physical_gauge_current_selection_or_laboratory_time_derived': False}}


def source_bindings():
    paths = set(native.source_bindings())|{'model/common_source/native_response_biased.py',
        'model/common_source/test_native_response_biased.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
