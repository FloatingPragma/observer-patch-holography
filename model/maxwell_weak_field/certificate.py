"""Exact carrier checks and conditional weak-field error certificates.

Derivative-bound inputs are caller-supplied uniform bounds on the whole
stated ball. Arithmetic validation does not certify an arbitrary field law.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from math import isqrt
import json
from pathlib import Path

import sympy as s

from ..maxwell_local_selection import certificate as local

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def rational(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, Q, s.Rational)):
        raise ValueError(f'{name} requires an exact rational, not a rounded value')
    result = Q(value)
    if result < 0 or (positive and result == 0):
        raise ValueError(f'{name} must be '+('positive' if positive else 'nonnegative'))
    return result


def signed_rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, Q, s.Rational)):
        raise ValueError(f'{name} requires an exact rational')
    return Q(value)


def sqrt_upper(value, *, bits=48):
    """Outward dyadic enclosure with a direct squared rational certificate."""
    value = rational(value, 'squared bound')
    if type(bits) is not int or bits < 0:
        raise ValueError('nonnegative integer precision required')
    scale = 1 << bits
    floor = isqrt(value.numerator*scale*scale//value.denominator)
    if Q(floor*floor, scale*scale) < value:
        floor += 1
    result = Q(floor, scale)
    local.require(result*result >= value, 'square-root enclosure failed')
    return result


def finite_time_bound(third_energy_bound, second_flow_bound, radius,
                      initial_radius, time):
    """Energy-confined C3/C2 theorem in the positive energy-Hessian norm.

    The caller certifies all analytic hypotheses in PROOF.md, including
    uniform derivative bounds on an open neighborhood of the closed ball.
    """
    M = rational(third_energy_bound, 'third energy derivative bound')
    L = rational(second_flow_bound, 'second flow derivative bound')
    R = rational(radius, 'domain radius', positive=True)
    r = rational(initial_radius, 'initial radius')
    t = signed_rational(time, 'time')
    a = M*R/3
    local.require(a < 1, 'positive quadratic energy enclosure required: M R/3 < 1')
    amplification = (1+a)/(1-a)
    trajectory_squared = amplification*r*r
    local.require(trajectory_squared < R*R, 'strict interior energy confinement required')
    return {
        'energy_relative_radius': str(a),
        'energy_norm_squared_amplification': str(amplification),
        'trajectory_norm_squared_upper': str(trajectory_squared),
        'error_norm_upper': str(L*trajectory_squared*abs(t)/2),
        'norm': 'positive vacuum energy-Hessian norm on R^30 x closed face flux',
        'caller_uniform_analytic_bounds_verified_by_this_API': False,
        'scope': 'all exact trajectories covered by the stated analytic hypotheses and initial norm bound',
    }


def quartic_bound(electric, magnetic, quartic, initial_radius, time):
    """Global energy and Maxwell-comparison bounds for the explicit quartic law."""
    epsilon = rational(electric, 'electric susceptibility', positive=True)
    beta = rational(magnetic, 'magnetic susceptibility', positive=True)
    lam = rational(quartic, 'quartic coefficient')
    r = rational(initial_radius, 'initial radius')
    t = signed_rational(time, 'time')
    amplitude_squared = r*r+lam*r**4/(2*beta*beta)
    error_squared = 6*lam*lam*t*t*amplitude_squared**3/(epsilon*beta**3)
    return {'trajectory_norm_squared_upper': str(amplitude_squared),
            'error_norm_squared_upper': str(error_squared),
            'error_norm_upper': str(sqrt_upper(error_squared)),
            'all_real_times_in_supplied_clock': True,
            'absolute_amplitude_order_at_fixed_time': 3,
            'relative_amplitude_order_at_fixed_time': 2}


def quartic_derivative_bounds(electric, magnetic, quartic, radius):
    epsilon = rational(electric, 'electric susceptibility', positive=True)
    beta = rational(magnetic, 'magnetic susceptibility', positive=True)
    lam = rational(quartic, 'quartic coefficient')
    R = rational(radius, 'domain radius', positive=True)
    # ||C||_2^2 <= ||C||_1 ||C||_infinity = 6 on the actual tables.
    coefficient = lam*sqrt_upper(Q(6)/(epsilon*beta**3))
    return {'third_energy_derivative_upper': str(6*lam*R/beta**2),
            'second_flow_derivative_upper': str(6*coefficient*R),
            'nonlinear_cubic_norm_coefficient_upper': str(coefficient)}


@lru_cache(maxsize=1)
def vacuum_certificate():
    D, C, actions = local.carrier()
    magnetic = local.classify()
    # For diagonal susceptibility, conjugation by signed U permutes weights.
    equations = []
    orbit = set()
    for U, _, _ in actions:
        for e in range(30):
            target = next(i for i in range(30) if U[i, e])
            if e == 0:
                orbit.add(target)
            row = [0]*30
            row[e] += 1
            row[target] -= 1
            equations.append(row)
    rank = s.Matrix(equations).rank()
    local.require(rank == 29 and len(orbit) == 30, 'electric seam transitivity failed')
    pivots = C.rref()[1]
    selector = s.zeros(30, len(pivots))
    for j, p in enumerate(pivots):
        selector[p, j] = 1
    basis = C[:, list(pivots)]
    local.require(len(pivots) == basis.rank() == 19 and C*selector == basis,
                  'closed-flux preimage certificate failed')
    local.require(s.ones(1, 20)*basis == s.zeros(1, 19), 'basis not closed flux')
    # Constant-flux Hessian extensions are invisible to every curl potential.
    invisible = s.ones(20, 20)
    local.require(C.T*invisible*C == s.zeros(30), 'ambient extension control failed')
    curl_norm_1 = max(sum(abs(C[i, j]) for i in range(20)) for j in range(30))
    curl_norm_inf = max(sum(abs(C[i, j]) for j in range(30)) for i in range(20))
    local.require((curl_norm_1, curl_norm_inf) == (2, 3), 'curl norm bound changed')
    return {'magnetic_vacuum_classification': magnetic,
            'diagonal_electric_weights': 30,
            'electric_rotation_constraint_rank_exact': rank,
            'electric_invariant_diagonal_dimension': 30-rank,
            'closed_flux_dimension': basis.rank(),
            'closed_flux_basis_curl_columns': list(pivots),
            'curl_preimages_verified_exactly': True,
            'restricted_magnetic_hessian': 'beta I on B0, by curl surjectivity and polarization',
            'ambient_constant_flux_extension_unique': False,
            'curl_operator_norm_squared_upper': int(curl_norm_1*curl_norm_inf),
            'energy_hessian': 'diag(epsilon I_30,beta I_B0)',
            'nonzero_susceptibilities_are_premises': True}


def fields(E, B):
    E, B = map(local.exact_matrix, (E, B))
    local.require(E.shape == (30, 1) and B.shape == (20, 1), 'wrong field dimensions')
    local.require(all(value.is_real is True for value in (*E, *B)), 'real field entries required')
    local.require(sum(B) == 0, 'closed magnetic flux required')
    return E, B


def quartic_field(E, B, electric=1, magnetic=1, quartic=1):
    epsilon = s.Rational(rational(electric, 'electric susceptibility', positive=True))
    beta = s.Rational(rational(magnetic, 'magnetic susceptibility'))
    lam = s.Rational(rational(quartic, 'quartic coefficient'))
    E, B = fields(E, B)
    _, C, _ = local.carrier()
    energy = epsilon*(E.dot(E))/2+beta*(B.dot(B))/2+lam*sum(b**4 for b in B)/4
    gradient_B = beta*B+lam*B.applyfunc(lambda b: b**3)
    derivative = (C.T*gradient_B/epsilon).col_join(-C*E)
    return energy, derivative


def gyroscopic_field(E, B):
    E, B = fields(E, B)
    D, C, _ = local.carrier()
    z = C.applyfunc(abs).T*B
    longitudinal = D*D.T
    skew = longitudinal*s.diag(*z)-s.diag(*z)*longitudinal
    return (C.T*B+skew*E).col_join(-C*E), skew


@lru_cache(maxsize=1)
def nonlinear_certificate():
    D, C, actions = local.carrier()
    E = s.Matrix(s.symbols('E0:30'))
    B = s.Matrix(s.symbols('B0:20'))
    epsilon, beta, lam = s.symbols('epsilon beta lambda', positive=True)
    magnetic_energy = beta*B.dot(B)/2+lam*sum(b**4 for b in B)/4
    gradient_B = s.Matrix([s.diff(magnetic_energy, b) for b in B])
    local.require(gradient_B == beta*B+lam*B.applyfunc(lambda b: b**3),
                  'quartic energy gradient failed')
    dE, dB = C.T*gradient_B/epsilon, -C*E
    balance = s.expand(epsilon*E.dot(dE)+gradient_B.dot(dB))
    local.require(balance == 0 and (D.T*dE).applyfunc(s.expand) == s.zeros(12, 1),
                  'quartic conservation identity failed')
    local.require(s.ones(1, 20)*dB == s.zeros(1, 1), 'quartic flux conservation failed')
    cube = B.applyfunc(lambda b: b**3)
    for _, V, _ in actions:
        local.require((V*B).applyfunc(lambda b: b**3) == V*cube,
                      'quartic rotation covariance failed')
    # Independent exact field with nonzero correction, and vacuum jet.
    E0, B0 = C.row(0).T, C[:, 0]
    t = s.symbols('a', real=True)
    energy, nonlinear = quartic_field(E0, B0, s.Rational(7, 3), s.Rational(11, 5), s.Rational(2, 7))
    G, _, reversal = local.evolution(s.Rational(7, 3), s.Rational(11, 5))
    x = E0.col_join(B0)
    linear = G*x
    _, scaled = quartic_field(t*E0, t*B0, s.Rational(7, 3), s.Rational(11, 5), s.Rational(2, 7))
    local.require(nonlinear != linear, 'nonlinear countercontrol has vanished')
    local.require((scaled-t*linear-t**3*(nonlinear-linear)).applyfunc(s.expand) == s.zeros(50, 1),
                  'cubic weak-field scaling failed')
    _, reverse = quartic_field(E0, -B0, s.Rational(7, 3), s.Rational(11, 5), s.Rational(2, 7))
    local.require(reverse == -reversal*nonlinear, 'quartic reversal failed')
    # A positive quartic action has a zero quadratic susceptibility.
    _, degenerate = quartic_field(t*E0, t*B0, 1, 0, 1)
    local.require(degenerate[:30, :].diff(t).subs(t, 0) == s.zeros(30, 1),
                  'degenerate magnetic susceptibility control failed')
    # Reversible, covariant conservation alone does not preserve Gauss.
    gyro, skew = gyroscopic_field(E0, B0)
    local.require(skew.T == -skew and x.dot(gyro) == 0, 'gyro energy identity failed')
    local.require(D.T*E0 == s.zeros(12, 1) and sum(B0) == 0, 'gyro input constraints failed')
    gauss_derivative = D.T*gyro[:30, :]
    local.require(gauss_derivative != s.zeros(12, 1), 'nonlinear Gauss counterexample failed')
    gyro_reverse, _ = gyroscopic_field(E0, -B0)
    local.require(gyro_reverse == -s.diag(*([1]*30+[-1]*20))*gyro, 'gyro reversal failed')
    longitudinal = D*D.T
    absolute_curl = C.applyfunc(abs)
    for U, V, _ in actions:
        unsigned = U.applyfunc(abs)
        local.require(U*longitudinal*U.T == longitudinal, 'longitudinal covariance failed')
        local.require(absolute_curl.T*V == unsigned*absolute_curl.T, 'absolute incidence covariance failed')
        g_rot, _ = gyroscopic_field(U*E0, V*B0)
        local.require(g_rot == s.diag(U, V)*gyro, 'gyro field covariance failed')
    return {'quartic_gradient_energy_Gauss_and_flux_identities_exact': True,
            'quartic_componentwise_covariance_checks': len(actions),
            'quartic_nonzero_cubic_correction_exact': True,
            'quartic_susceptibility_counterexample': 'sum(B_f^4)/4 has zero magnetic Hessian',
            'nonlinear_Gauss_counterexample': {
                'rotation_covariance_checks': len(actions),
                'energy_Faraday_reversal_preserved': True,
                'input_electric': 'first oriented face cycle C[0,:]^T',
                'input_magnetic': 'curl of first seam potential C[:,0]',
                'initial_Gauss_charge_zero': True,
                'Gauss_derivative': [str(v) for v in gauss_derivative]},
            'full_nonlinear_Maxwell_uniqueness_claimed': False}


def build():
    radius = Q(1, 4)
    estimates = quartic_derivative_bounds(Q(7, 3), Q(11, 5), Q(2, 7), radius)
    sample_general = finite_time_bound(Q(estimates['third_energy_derivative_upper']),
                                      Q(estimates['second_flow_derivative_upper']),
                                      radius, Q(1, 100), Q(3, 2))
    paths = local.SOURCE_PATHS+(
        'model/maxwell_local_selection/certificate.py',
        'model/maxwell_local_selection/PROOF.md',
        'model/maxwell_weak_field/__init__.py',
        'model/maxwell_weak_field/certificate.py',
        'model/maxwell_weak_field/PROOF.md',
        'model/maxwell_weak_field/test_certificate.py')
    return {'schema': 'oph.maxwell_weak_field.v1',
            'vacuum': vacuum_certificate(),
            'smooth_evolution': {
                'vacuum_equilibrium_derived': True,
                'linear_generator': 'G(E,B)=((beta/epsilon) C^T B,-C E) on closed flux',
                'quadratic_action_assumed': False,
                'linear_evolution_assumed': False,
                'generic_remainder_absolute_order': 2,
                'full_nonlinear_Gauss_conservation_derived_for_every_flow': False},
            'nonlinear': nonlinear_certificate(),
            'bound_example': {
                'parameters': {'epsilon': '7/3', 'beta': '11/5', 'lambda': '2/7',
                               'domain_radius': '1/4', 'initial_radius': '1/100', 'time': '3/2'},
                'quartic_uniform_derivative_bounds': estimates,
                'general_theorem': sample_general,
                'quartic_theorem': quartic_bound(Q(7, 3), Q(11, 5), Q(2, 7), Q(1, 100), Q(3, 2))},
            'scope': {'explicit_nonzero_susceptibility_and_locality_hypotheses': True,
                      'native_repair_selects_energy_or_flow': False,
                      'physical_clock_or_coefficients_selected': False,
                      'continuum_or_particle_pole_identification': False,
                      'Lean_kernel_checked': False},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(build(), sort_keys=True, indent=2)+'\n'
    path = HERE/'results.json'
    if args.check:
        local.require(path.read_text() == text, 'weak-field certificate or source custody changed')
        print('Exact weak-field Maxwell checks and rational bounds passed.')
    else:
        path.write_text(text)


if __name__ == '__main__':
    main()
