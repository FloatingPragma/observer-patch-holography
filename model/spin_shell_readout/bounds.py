"""Exact clock resolution, unnormalized tube weights and sampling errors."""
from fractions import Fraction as F


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational required')
    value = F(value)
    if positive and value <= 0:
        raise ValueError('positive exact rational required')
    return value


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('bounded exact integer required')
    return value


def ceil(value):
    return -(-value.numerator // value.denominator)


def tail_bound(clock_size, turn_margin):
    """Complete circular Fejer tail, summed over every distant clock bin."""
    integer(clock_size, 2)
    margin = rational(turn_margin, positive=True)
    if margin > F(1, 2):
        raise ValueError('circular margin must not exceed half a turn')
    x = clock_size * margin
    return min(F(1), F(1, 2) / x**2 + F(1, 2) / x)


def clock_plan(turn_margin, tail_error):
    margin = rational(turn_margin, positive=True)
    error = rational(tail_error, positive=True)
    if margin > F(1, 2) or error >= 1:
        raise ValueError('half-turn margin and error below one required')
    required = max(2, ceil(1 / (margin * error)))
    size = 1 << (required - 1).bit_length()
    bound = tail_bound(size, margin)
    if bound > error:
        raise ArithmeticError('complete clock tail exceeds its request')
    return {'clock_size': size, 'clock_qubits': size.bit_length() - 1,
            'turn_margin': margin, 'complete_clock_tail_upper': bound,
            'maximum_clock_time_over_tau': size - 1,
            'controlled_binary_duration_sum_over_tau': size - 1}


def bin_ranges(clock_size, center_interval, turn_width, turn_margin):
    """Inner/outer inclusive bin ranges; empty range is represented by (1, 0).

    Every possible true center in the certified interval is covered. The
    no-wrap condition is explicit, and a separate spectrum bound must place
    every excitation within this one phase period.
    """
    integer(clock_size, 2)
    if type(center_interval) is not tuple or len(center_interval) != 2:
        raise ValueError('one exact center interval required')
    lo, hi = map(rational, center_interval)
    width = rational(turn_width, positive=True)
    margin = rational(turn_margin, positive=True)
    if not 0 <= lo <= hi <= 1 or not margin < width:
        raise ValueError('ordered center and margin below tube width required')
    if not (0 < lo-width-margin and hi+width+margin < 1):
        raise ValueError('tube and collar must lie strictly inside one phase period')
    if width + margin >= F(1, 2):
        raise ValueError('tube and collar must occupy less than half a turn')

    def indices(a, b):
        first, last = ceil(clock_size*a), (clock_size*b).numerator // (clock_size*b).denominator
        return (1, 0) if first > last else (first, last)

    return {'inner': indices(hi-width+margin, lo+width-margin),
            'outer': indices(lo-width-margin, hi+width+margin),
            'complete_clock_tail_upper': tail_bound(clock_size, margin),
            'all_other_clock_bins_retained': True}


def tube_bounds(success_interval, inner_interval, outer_interval, tail,
                output_halftrace_error=F()):
    """Bounds on original submeasure weight; failures are never normalized away."""
    def interval(pair):
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError('ordered normalized probability interval required')
        a, b = map(rational, pair)
        if not 0 <= a <= b <= 1:
            raise ValueError('probability intervals must lie in [0,1]')
        return a, b
    slo, shi = interval(success_interval)
    ilo, ihi = interval(inner_interval)
    olo, ohi = interval(outer_interval)
    leakage = rational(tail)
    error = rational(output_halftrace_error)
    if not 0 <= leakage <= 1 or not 0 <= error <= 1:
        raise ValueError('normalized error allowances required')
    # The nested event relations need only have a consistent realization;
    # externally certified intervals can overlap more widely than the truth.
    if ilo > shi or olo > shi or ilo > ohi:
        raise ValueError('success and nested event intervals are inconsistent')
    smin, smax = max(F(), slo-error), min(F(1), shi+error)
    lower = max(F(), ilo-error-leakage*smax)
    upper = min(smax, ohi+error+leakage*smax)
    if lower > upper:
        raise ValueError('tube observations have inconsistent bounds')
    return {'original_tube_weight_lower': lower,
            'original_tube_weight_upper': upper,
            'original_success_weight_lower': smin,
            'original_success_weight_upper': smax,
            'fraction_of_original_success_lower': lower/smax if smax else F(),
            'conditional_probability_inferred_without_positive_success': False,
            'all_failure_and_nonselected_outcomes_retained': True}


def sampled_intervals(trials, success, inner, outer, error, failure, searches=1):
    """Simultaneous elementary Bernoulli intervals for a predeclared finite scan.

    Three event means per candidate; union Chebyshev uses variance <=1/4.
    It assumes independent normalized trials. No simulated trials are implied.
    """
    integer(trials, 1)
    for count in (success, inner, outer):
        integer(count)
    integer(searches, 1)
    eps = rational(error, positive=True)
    beta = rational(failure, positive=True)
    if not 0 <= inner <= outer <= success <= trials or not beta < 1 or eps > 1:
        raise ValueError('complete nested counts and normalized tolerances required')
    required = ceil(F(3*searches, 4)/(beta*eps**2))
    if trials < required:
        raise ValueError('too few trials for the simultaneous requested confidence')
    def around(count):
        p = F(count, trials)
        return max(F(), p-eps), min(F(1), p+eps)
    return {'success_interval': around(success), 'inner_interval': around(inner),
            'outer_interval': around(outer), 'all_other_outcome_count': trials-success,
            'required_trials': required, 'simultaneous_failure_upper': beta,
            'predeclared_candidate_count': searches,
            'independent_normalized_trials_assumed': True}
