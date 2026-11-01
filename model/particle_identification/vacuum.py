"""Exact signs of the cube's positive greater spectral measure.

This tests a proposed vacuum, not the existence of the source theory. An
excited reference may perfectly consistently have this positive measure
on both signs of energy transfer. Retarded commutators are not being used
as a positivity test. The lowering unitary is an algebraic witness on the
admitted finite field algebra, not an executed source extraction protocol.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..cube_eta_protection import spectrum


FACTORS = (
    ('Q2', (1, -16, -12), ((-1, 0), (16, 17))),
    ('C64', (1, -16, -12, 64), ((-3, -2), (1, 2), (16, 17))),
    ('C128', (1, -16, -12, 128), ((-3, -2), (2, 3), (16, 17))),
)


def evaluate(coefficients, x):
    value = F()
    for c in coefficients:
        value = value*x+c
    return value


def root_interval(coefficients, bracket, steps=80):
    """Exact bisection of a proven unique sign-changing real root."""
    if type(steps) is not int or not 1 <= steps <= 256:
        raise ValueError('one through256 exact bisections required')
    left, right = map(F, bracket)
    fl, fr = evaluate(coefficients, left), evaluate(coefficients, right)
    if left >= right or fl*fr >= 0:
        raise ValueError('a proper strict sign-changing root bracket is required')
    for _ in range(steps):
        middle = (left+right)/2
        fm = evaluate(coefficients, middle)
        if fm == 0:
            return middle, middle
        if fl*fm < 0:
            right, fr = middle, fm
        else:
            left, fl = middle, fm
    return left, right


def interval_polynomial(coefficients, interval):
    low = high = F()
    left, right = interval
    for coefficient in coefficients:
        products = (low*left, low*right, high*left, high*right)
        low, high = min(products)+coefficient, max(products)+coefficient
    return low, high


def ratio_interval(numerator, denominator):
    if denominator[0] <= 0 <= denominator[1]:
        raise ValueError('denominator interval must exclude zero')
    values = tuple(a/b for a in numerator for b in denominator)
    return min(values), max(values)


def outward_decimal(interval, places=12):
    scale = 10**places
    lower = interval[0]*scale
    upper = interval[1]*scale
    return F(lower.numerator//lower.denominator, scale), F(-((-upper.numerator)//upper.denominator), scale)


def exact_coefficients(polynomial):
    return tuple(F(int(x.p), int(x.q)) for x in s.Poly(polynomial, spectrum.Z).all_coeffs())


@lru_cache(maxsize=1)
def certificate():
    numerator, denominator, determinant = spectrum.probe_resolvent()
    z = spectrum.Z
    expected_denominator = s.prod(sum(c*z**(len(coeffs)-i-1) for i, c in enumerate(coeffs))
                                  for _, coeffs, _ in FACTORS)
    if s.expand(denominator-expected_denominator) != 0 or determinant <= 0:
        raise ArithmeticError('the exact eight-line source resolvent changed')
    num = exact_coefficients(numerator)
    derivative = exact_coefficients(s.diff(denominator, z))
    rows = []
    negative_low = negative_high = F()
    for name, coefficients, brackets in FACTORS:
        # Degree-many disjoint strict sign changes prove that these are all
        # roots and that each is real and simple, without floating root fits.
        if len(brackets) != len(coefficients)-1:
            raise ArithmeticError('the brackets must exhaust the polynomial degree')
        if any(a[1] > b[0] for a, b in zip(brackets, brackets[1:])):
            raise ArithmeticError('the within-factor brackets overlap')
        for bracket in brackets:
            interval = root_interval(coefficients, bracket)
            omega = tuple(x/16-1 for x in interval)
            sign = -1 if omega[1] < 0 else 1 if omega[0] > 0 else 0
            if sign == 0:
                raise ArithmeticError('a probe frequency has not been separated from zero')
            residue = ratio_interval(interval_polynomial(num, interval),
                                     interval_polynomial(derivative, interval))
            if residue[0] <= 0:
                raise ArithmeticError('the positive spectral residue was not certified')
            if sign < 0:
                negative_low += residue[0]; negative_high += residue[1]
            lo, hi = outward_decimal(residue)
            rows.append({'factor': name, 'coarse_z_bracket': list(bracket),
                         'endpoint_polynomial_values': [str(evaluate(coefficients, F(x))) for x in bracket],
                         'z_interval': [str(x) for x in interval],
                         'frequency_sign': sign,
                         'omega_interval': [str(x) for x in omega],
                         'positive_normalized_probe_weight_interval': [str(lo), str(hi)]})
    if sum(r['frequency_sign'] < 0 for r in rows) != 5 or sum(r['frequency_sign'] > 0 for r in rows) != 3:
        raise ArithmeticError('the exact source frequency count changed')
    mass = outward_decimal((negative_low, negative_high))
    if not F(1, 40) < mass[0] < mass[1] < F(13, 500):
        raise ArithmeticError('negative positive-measure mass lost its stated rational enclosure')
    # Exact polynomial congruence proves the residue at BOTH Q2 roots:
    # w(z)=(2z+3)/266; at z-=8-2sqrt19 it is (19-4sqrt19)/266.
    residue_q2 = (2*z+3)/266
    if s.rem(numerator-s.diff(denominator, z)*residue_q2, spectrum.Q2, z) != 0:
        raise ArithmeticError('quadratic-root residue identity failed')
    if not F(87, 20)**2 < 19 < F(109, 25)**2:
        raise ArithmeticError('exact square-root bounds changed')
    return {
        'reference': 'the actually prepared fixed-N6 cube eta reference R',
        'reference_energy_H': '3',
        'normalized_centered_probe': 'q=(4 n_origin-3)R/sqrt15',
        'greater_correlation_convention': '<R|delta_n(t) delta_n(0)|R>=(15/16) sum_j w_j exp(-i omega_j t)',
        'positive_measure_convention': 'Fourier transform with exp(+i omega t); omega_j=E_j-E_R=z_j/16-1',
        'measure_is_retarded_commutator': False,
        'retarded_caution': 'a ground-state retarded commutator also contains both signs; that fact is not this obstruction',
        'exact_eight_column_cyclic_Gram_determinant': str(determinant),
        'root_count_proof': 'degree-many disjoint strict rational sign changes for each exact factor',
        'negative_frequency_lines': 5, 'positive_frequency_lines': 3, 'zero_frequency_lines': 0,
        'spectral_rows': rows,
        'normalized_probe_negative_mass_interval': [str(x) for x in mass],
        'normalized_probe_negative_mass_simple_bounds': ['1/40', '13/500'],
        'actual_density_negative_mass_interval': [str(F(15, 16)*x) for x in mass],
        'quadratic_lower_line': {
            'z': '8-2 sqrt19', 'H_energy': '5/2-sqrt19/8',
            'omega': '-1/2-sqrt19/8', 'omega_strictly_less_than': '-1',
            'normalized_probe_weight': '(1-4/sqrt19)/14 = 1/14-2sqrt19/133',
            'normalized_probe_weight_rational_bounds': ['1/174', '9/1526'],
            'actual_density_weight': '(15/16)*(1-4/sqrt19)/14',
        },
        'lowering_witness': {
            'state': 'psi_minus=P_(H=5/2-sqrt19/8) q / sqrt(w_minus)',
            'orthogonal_to_R': True, 'same_CAR_number': 6,
            'unitary': 'U=I-|R><R|-|psi_minus><psi_minus|+|R><psi_minus|+|psi_minus><R|',
            'energy_change_from_R': '-1/2-sqrt19/8 < -1',
            'passivity_scope': 'full finite unitary algebra, already within the number-preserving protected sector',
            'native_extraction_protocol_executed': False,
            'executed_local_F_mean_energy_change': '0',
        },
        'positive_affine_energy_conversion_removes_obstruction': False,
        'same_number_chemical_potential_removes_obstruction': False,
        'global_time_reversal_removes_mixed_signs': False,
        'R_is_a_ground_state_of_the_admitted_field_algebra': False,
        'R_is_passive_under_all_number_preserving_unitaries': False,
        'excited_reference_or_thermal_response_inconsistent': False,
        'other_vacuum_candidates_excluded': False,
        'physical_particle_identification_established': False,
    }
