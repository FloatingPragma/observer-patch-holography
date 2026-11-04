from fractions import Fraction as F
from math import comb, factorial

import pytest
import sympy as s

from . import bounds


def test_profiles_and_exact_derivative_margins():
    x = s.symbols('x', real=True)
    polynomial = (1-x*x)**3
    assert s.integrate(polynomial**2, (x, -1, 1)) == s.Rational(2048, 3003)
    for degree in (0, 1, 2):
        derivative = s.diff(polynomial, x, degree)
        assert derivative.subs(x, -1) == derivative.subs(x, 1) == 0
    assert bounds.profile_certificate()['preparation_uses_interaction_phase'] is False
    assert bounds.analytic_certificate()['density_bounds'] == ['a^3/100', 'a^3/2']
    # Two derivatives of sqrt(rho)*q, including both distinct Cook packets.
    incoming = 2000000*2+2*500*48+192
    first = 500*2+48
    both = 2*incoming+2*3220*first+11000000*2
    assert both == 36845504 < 40000000


def test_finite_schedule_satisfies_all_rational_constraints_without_giant_tapes():
    p = bounds.plan()
    a, error = p['spacing'], p['component_error']
    T, R, L = (p[k] for k in ('absolute_echo_segment_duration', 'normalized_packet_support_radius', 'box_padding'))
    assert F(160000000, 1)/(a**3*T) <= error
    assert F(40000000, 1)/(a**5*R) <= error
    assert L >= 288*T/a**2
    # The bit-count implication avoids allocating 2**L for the enormous box.
    assert L >= 2+bounds.log2_ceiling(1/error)
    V, M, n = (p[k] for k in ('box_sites', 'field_bits', 'first_order_steps_per_echo_segment'))
    assert 6*T*T*(10*V/a**2)**2/n <= error
    assert p['native_pulse_log2_upper'] > M
    assert p['dyadic_integrated_port_fraction_bits']-p['native_pulse_log2_upper'] >= bounds.log2_ceiling(58/error)
    assert p['continuum_phase_error_upper'] == F(2101, 32768)
    assert p['valid_Y_plus_probability_upper'] < F(1, 8)
    assert p['complete_law_event_separation_lower'] > F(3, 10)
    assert p['fine_spacing_native_experiment_executed'] is False


@pytest.mark.parametrize('a,error', [(F(1, 65535), F(1, 100)), (0, F(1, 100)),
                                    (True, F(1, 100)), (F(1, 65536), F(1, 99)),
                                    (F(1, 65536), 0)])
def test_invalid_schedule_domains_fail(a, error):
    with pytest.raises(ValueError):
        bounds.plan(a, error)


def test_logarithmic_budget_matches_small_exact_pulse_counts():
    for M in (1, 2, 4, 10, 100):
        V, n = M, 3
        count = 11*M*(1 << M)+768*n*V*(M+1)+10
        logarithm = max(M+(11*M).bit_length(), (768*n*V*(M+1)+10).bit_length())+1
        assert count <= 1 << logarithm
    for value in (F(1, 2), F(1), F(2), F(9, 4), F(17), F(257, 2)):
        exponent = bounds.log2_ceiling(value)
        assert 2**exponent >= value
        assert exponent == 0 or 2**(exponent-1) < value


def test_complete_history_sampling_bound_is_an_exact_positive_series():
    result = bounds.sampling_certificate()
    series = sum((F(10**n, factorial(n)) for n in range(9)), F())
    assert series > 4000
    assert result['failure_upper'] == 4/series < F(1, 1000)
    assert result['rejected_histories_counted'] is True
