"""Certified finite-sample identifiability within the supplied eight-line model.

This is a prospective resolution theorem, not an executed long-time detector
or a claim that arbitrary spectra are identifiable from finitely many samples.
The actual observer circuit has its own execution and error certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

from ..certified_reference import interval as iv
from ..particle_identification import vacuum


SPACING = 39
DETERMINANT_LOWER = F(1400)
# Each 7x7 minor has modulus at most (sqrt(7))**7 < 343*8/3.
INVERSE_INFINITY_UPPER = F(8*343*8, 3*1400)
# p_0 errors at the two quadratures give |delta C| <= 2sqrt(2)e < 3e;
# removing the elastic atom multiplies by 16/15.
WEIGHT_ERROR_MULTIPLIER = INVERSE_INFINITY_UPPER*F(16, 5)


def exact_nonnegative(value):
    if type(value) not in (int, F) or value < 0:
        raise ValueError('an exact nonnegative integer or Fraction is required')
    return F(value)


def interval_pair(low, high):
    if low > high:
        raise ValueError('ordered rational endpoints required')
    return iv.Interval(iv.Interval.point(low).lo, iv.Interval.point(high).hi)


def absolute_sine_interval(low, high):
    """Outward enclosure using exact pi reduction and the Lipschitz bound."""
    source = interval_pair(low, high)
    pi = iv.pi()
    midpoint = F(source.lo+source.hi, 2*iv.SCALE)
    pi_midpoint = F(pi.lo+pi.hi, 2*iv.SCALE)
    ratio = midpoint/pi_midpoint+F(1, 2)
    periods = ratio.numerator//ratio.denominator
    reduced = source-periods*pi  # |sin(x-n*pi)| = |sin(x)| for integer n.
    center = F(reduced.lo+reduced.hi, 2*iv.SCALE)
    radius = F(reduced.hi-reduced.lo, 2*iv.SCALE)
    if abs(center)+radius > 4:
        raise ValueError('reduced angle is outside the certified Taylor domain')
    sine, _ = iv.sine_cosine(center)
    lower, upper = F(sine.lo, iv.SCALE)-radius, F(sine.hi, iv.SCALE)+radius
    if lower <= 0 <= upper:
        return F(), max(-lower, upper)
    return min(abs(lower), abs(upper)), max(abs(lower), abs(upper))


def vandermonde_determinant_bounds(frequencies, spacing):
    """For V[n,j]=exp(-i*n*spacing*omega_j), n=0,...,m-1."""
    spacing = exact_nonnegative(spacing)
    if not spacing or len(frequencies) < 2:
        raise ValueError('positive time spacing and at least two lines required')
    low = high = F(1)
    for left, right in combinations(frequencies, 2):
        a = F(right[0])-F(left[1]); b = F(right[1])-F(left[0])
        s0, s1 = absolute_sine_interval(spacing*a/2, spacing*b/2)
        low *= 2*s0; high *= 2*s1
    return low, high


@lru_cache(maxsize=1)
def spectral_resolution_certificate():
    rows = sorted(vacuum.certificate()['spectral_rows'], key=lambda r: F(r['omega_interval'][0]))
    intervals = [tuple(map(F, row['omega_interval'])) for row in rows]
    lower, upper = vandermonde_determinant_bounds(intervals, SPACING)
    if lower <= DETERMINANT_LOWER:
        raise ArithmeticError('the exact eight-line sampling matrix lost its determinant bound')
    display = vacuum.outward_decimal((lower, upper))
    weight_lower = min(F(row['positive_normalized_probe_weight_interval'][0]) for row in rows)
    return {
        'source_support': 'eight certified greater-correlation frequencies of the centered N6 density probe',
        'sample_model_times': [SPACING*n for n in range(8)],
        'Ramsey_quadratures_per_time': 2,
        'complex_signal': 'C(t)=1/16+(15/16)*sum_j w_j exp(-i omega_j t)',
        'elastic_removal': 'y_n=(16*C(39*n)-1)/15=sum_j V[n,j]*w_j',
        'sampling_matrix': 'V[n,j]=exp(-i*39*n*omega_j), n=0,...,7',
        'determinant_identity': '|det V|=product_(j<k) 2|sin(39*(omega_k-omega_j)/2)|',
        'determinant_certified_interval': [str(x) for x in display],
        'determinant_strict_lower': str(DETERMINANT_LOWER),
        'inverse_infinity_norm_strict_upper': str(INVERSE_INFINITY_UPPER),
        'weight_error_per_uniform_probability_error_strict_upper': str(WEIGHT_ERROR_MULTIPLIER),
        'smallest_source_weight_strict_lower': str(weight_lower),
        'identification_scope': 'unique weights on this source-specified support, with both Ramsey quadratures',
        'arbitrary_unknown_spectrum_identified_by_finitely_many_samples': False,
        'long_time_sampling_protocol_executed_here': False,
        'physical_species_or_rest_mass_assignment': False,
    }


def weight_error_bound(probability_error):
    """Systematic plus statistical error per pointer-zero probability."""
    spectral_resolution_certificate()
    return WEIGHT_ERROR_MULTIPLIER*exact_nonnegative(probability_error)


def prospective_sampling_budget():
    """A rigorous, deliberately explicit sufficient shot/accuracy schedule.

Hoeffding plus a union bound over 16 Bernoulli marginals: failure probability
<=32 exp(-2*N*eps**2). The rational e>8/3 bound replaces a floating log.
All source/repair/export outcomes count; there is no survival postselection.
Independent repeated preparations and controlled systematic bias are premises.
    """
    systematic = F(1, 65536)
    statistical = F(1, 32768)
    shots = 9*(1 << 29)
    exponent = 2*shots*statistical**2
    if exponent != 9:
        raise ArithmeticError('shot schedule exponent changed')
    failure = 32*F(3, 8)**9
    error = weight_error_bound(systematic+statistical)
    minimum = F(spectral_resolution_certificate()['smallest_source_weight_strict_lower'])
    if failure >= F(1, 100) or 2*error >= minimum:
        raise ArithmeticError('the prospective schedule no longer separates every weight from zero')
    return {
        'systematic_probability_error_at_each_setting_at_most': str(systematic),
        'sampling_probability_error_at_each_setting_at_most': str(statistical),
        'independent_shots_per_setting': shots,
        'number_of_settings': 16,
        'total_shots': 16*shots,
        'joint_failure_probability_strict_upper': str(failure),
        'reconstructed_weight_error_strict_upper': str(error),
        'all_eight_positive_weights_separated_from_zero': True,
        'requires_independent_repeated_preparations': True,
        'requires_systematic_error_bound_at_all_16_settings': True,
        'abort_records_are_retained_without_postselection': True,
        'sampling_budget_executed': False,
    }
