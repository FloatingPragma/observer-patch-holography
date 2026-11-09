"""Exact native record rates with a common clock factor removed from ratios.

These are rates of the complete-attempt binary carrier. Their identification
with physical particle energies is not assumed.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt

import sympy as s

from model.native_quantum_reconstruction import certificate as source
from model.observer import quantum_mixing


def require(test, message):
    if not test:
        raise ValueError(message)


@lru_cache(maxsize=1)
def source_matrices():
    counts, denominator = source.source_counts()
    p = s.ImmutableMatrix(counts)/denominator
    b = s.ImmutableMatrix(source.STATES)
    a = s.zeros(12)
    for i, j in quantum_mixing.seams():
        a[i, j] = a[j, i] = 1
    laplacian = 5*s.eye(12)-a
    return p, b, s.ImmutableMatrix(laplacian)


def verify_source(p, b, laplacian):
    """All 66 inputs and 60 seam/coin outcomes enter this identity."""
    require(p == p.T and p*s.ones(66, 1) == s.ones(66, 1),
            'complete symmetric attempt kernel required')
    require(b.T*b == 10*s.eye(12)+s.ones(12), 'all binary total-two states required')
    require(p*b == b*(s.eye(12)-laplacian/60),
            'record clock or the source occupation dynamics changed')
    z = s.Symbol('z')
    expected = z*(z-6)**5*((z-5)**2-5)**3
    require(s.expand(laplacian.charpoly(z).as_expr()-expected) == 0,
            'actual source Laplacian spectrum changed')


def neg_log_bounds(p_lo, p_hi, terms=24):
    """Rational bounds on -log p for 0<p_lo<=p<=p_hi<=1."""
    require(type(p_lo) is Q and type(p_hi) is Q and 0 < p_lo <= p_hi <= 1,
            'positive exact rational contraction interval required')
    require(type(terms) is int and terms >= 1, 'positive integer series degree required')

    def series(p):
        y = (1-p)/(1+p)
        total = 2*sum((y**(2*k+1)/Q(2*k+1) for k in range(terms)), Q())
        tail = 2*y**(2*terms+1)/(Q(2*terms+1)*(1-y*y))
        return total, total+tail

    return series(p_hi)[0], series(p_lo)[1]


def mode_intervals():
    unit = 10**42
    low = isqrt(5*unit*unit)
    root_lo, root_hi = Q(low, unit), Q(low+1, unit)
    require(root_lo**2 < 5 < root_hi**2, 'sqrt(5) interval failed')
    return {'slow': ((55+root_lo)/60, (55+root_hi)/60),
            'middle': (Q(9, 10), Q(9, 10)),
            'fast': ((55-root_hi)/60, (55-root_lo)/60)}


def interval_ratio(numerator, denominator):
    a, b = numerator
    c, d = denominator
    require(0 < c <= d and 0 <= a <= b, 'positive denominator interval required')
    return a/d, b/c


def rounded_interval(bounds, digits=24):
    unit = 10**digits
    lo, hi = bounds
    return {'lower': str(Q((lo*unit).__floor__(), unit)),
            'upper': str(Q((hi*unit).__ceil__(), unit)),
            'decimal_midpoint_display_only': f'{float((lo+hi)/2):.14g}'}


def rate_ratio_error(gap_a, gap_b, error_a, error_b):
    """Absolute ratio error for independently bounded positive gap errors."""
    require(all(type(x) is Q for x in (gap_a, gap_b, error_a, error_b)),
            'exact rationals required')
    require(gap_a >= 0 and gap_b > error_b >= 0 and error_a >= 0,
            'denominator error must be smaller than its positive gap')
    return (error_a*gap_b+gap_a*error_b)/(gap_b*(gap_b-error_b))


@lru_cache(maxsize=1)
def certificate():
    p, b, laplacian = source_matrices()
    verify_source(p, b, laplacian)
    # The intertwining proves every horizon by induction. These finite powers
    # also check the implementation's order and complete normalization.
    for n in (0, 1, 2, 3):
        require(p**n*b == b*(s.eye(12)-laplacian/60)**n,
                'complete native covariance lost an attempt')
    intervals = mode_intervals()
    energies = {key: neg_log_bounds(*value) for key, value in intervals.items()}
    ratios = {key+'_over_slow': interval_ratio(energies[key], energies['slow'])
              for key in ('middle', 'fast')}
    return {
        'scope': 'complete-attempt native occupation correlations and reconstructed logarithmic generator',
        'complete_input_columns': 66,
        'raw_seam_coin_outcomes_per_input': 60,
        'intertwining': 'P B = B (I-L/60)',
        'all_horizons_identity': 'P^n B = B (I-L/60)^n, n>=0',
        'laplacian_characteristic_polynomial': 'z (z-6)^5 ((z-5)^2-5)^3',
        'uniform_reference': '1/66 on every native configuration',
        'complete_endpoint_covariance': 'C_lambda(n)=(5 rank(Pi_lambda)/33)(1-lambda/60)^n',
        'positive_mode_residues': {'slow': '5/11', 'middle': '25/33', 'fast': '5/11'},
        'mode_rates': {'slow': '-log((55+sqrt(5))/60)', 'middle': '-log(9/10)',
                       'fast': '-log((55-sqrt(5))/60)'},
        'certified_dimensionless_rates': {k: rounded_interval(v) for k, v in energies.items()},
        'certified_same_clock_ratios': {k: rounded_interval(v) for k, v in ratios.items()},
        'same_reference_modular_gap': 'L times the corresponding H=-log(P) gap; additive trace scalar cancels',
        'common_block_clock': '-log(P^b)=b*(-log(P)), b>=1; ratios unchanged',
        'common_positive_clock_and_energy_scale_cancel': True,
        'empirical_mass_or_phase_used': False,
        'native_unitary_clock_implementation_proved': False,
        'reference_loop_length_selected': False,
        'physical_particle_identity_proved': False,
        'physical_GeV_conversion_derived': False,
        'calibrated_charged_action_attached_to_this_source': False,
    }
