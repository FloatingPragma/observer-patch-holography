"""Conditional optical comparison for every distinct regional event clock.

Regions use one native carrier, not independent replicas. A region ticks once
per incident seam event. The finite family includes every nonempty port mask;
two masks exposing exactly the same events define one clock and are excluded
from the distinct-clock rejection. No atomic or servo dynamics is asserted.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from math import comb, isqrt
from pathlib import Path

from . import optical_clock as optical
from ..laboratory_clock import regional_clock as regional
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = optical.ROOT
SOURCE_PATHS = tuple(sorted(set((*optical.SOURCE_PATHS, *regional.SOURCE_PATHS,
    'model/predictions/regional_optical.py'))))


def _fraction(value, *, positive=False):
    # Use the existing canonical input check, including mutable Fraction slots.
    value = regional.e.exact_core._fraction(value)
    if value < 0 or (positive and value == 0):
        raise ValueError('a canonical nonnegative exact fraction is required')
    return value


def bernoulli_variance(a, b, common):
    """Choose an independent Bernoulli factor of the exact renewal count.

    X=Bernoulli(c/a)+Geom0(a/(a+b-c)). A geometric variable with failure q
    factors into an independent Bernoulli(q/(1+q)) and twice Geom0(1-q²).
    The stronger of the two available Bernoulli variances improves the bound.
    """
    law = regional.parameters(a, b, common)
    terminal, p = law['terminal_B_probability'], law['geometric_success_probability']
    parity = (1-p)/(2-p)
    return max(terminal*(1-terminal), parity*(1-parity))


def point_mass_bound(a, b, common, count):
    """Exact outward bound on sup_k Pr[B(count)=k] for distinct event sets.

    Fourier inversion of the independent Binomial(count,r) factor bounds its
    largest mass by sqrt(pi/(8*r*(1-r)*count)). Use pi<22/7, an outward
    rational square root, and floor(sqrt(count)); no normal approximation.
    """
    if type(count) is not int or count < 2:
        raise ValueError('an original integer count of at least two is required')
    variance = bernoulli_variance(a, b, common)
    if variance == 0:
        raise ValueError('identical event sets are one clock, not a distinct pair')
    scale = 10**6
    square = F(22*scale*scale, 56)/variance
    integer = isqrt(square.numerator//square.denominator)
    if integer*integer < square:
        integer += 1
    return min(F(1), F(integer, scale*isqrt(count)))


def decimal_ceiling(probability):
    """Smallest nonpositive integer e with min(1,probability)<=10**e."""
    probability = _fraction(probability, positive=True)
    probability = min(F(1), probability)
    exponent, power = 0, F(1)
    while probability <= power/10:
        exponent, power = exponent-1, power/10
    return exponent


def selection_bound(total, retained, mass, *, choice_budget=1):
    """Anchor-and-subset bound, evaluated exactly before outward rounding.

    Any K retained independent values in an interval require an anchor and K-1 other
    values in that anchor's interval. Their union probability is at most
    T*choose(T-1,K-1)*mass**(K-1). The interval centre and arbitrary removals
    can depend on all values. A finite window/count/model search is charged
    separately. The bound uses the observed width monotonically.
    """
    mass = _fraction(mass, positive=True)
    if (type(total) is not int or type(retained) is not int
            or not 2 <= retained <= total or mass > 1
            or type(choice_budget) is not int or choice_budget < 1):
        raise ValueError('valid integer counts and a positive exact interval mass are required')
    bound = F(total*choice_budget*comb(total-1, retained-1))*mass**(retained-1)
    return {'finite_selection_union_budget': choice_budget,
            'bound_formula': 'min(1,T*choice_budget*choose(T-1,K-1)*p**(K-1))',
            'conservative_probability_upper_power_of_ten': decimal_ceiling(bound)}


def evaluate(row, *, error=F(1, 10**12), gain_lower=F(1)):
    error, gain_lower = _fraction(error), _fraction(gain_lower, positive=True)
    integer_fields = ('block_seconds', 'first_UTC_second', 'last_UTC_second',
                      'possible_calendar_blocks', 'complete_blocks')
    if any(type(row[key]) is not int for key in integer_fields):
        raise ValueError('calendar and block counts must be original integers')
    low = (optical.REFERENCE_HZ-1)*optical.BLOCK_SECONDS
    high = (optical.REFERENCE_HZ+1)*optical.BLOCK_SECONDS
    span = F(row['maximum_block_ratio'])-F(row['minimum_block_ratio'])
    if (span < 0 or row['block_seconds'] != optical.BLOCK_SECONDS
            or not optical.CALENDAR_START <= row['first_UTC_second'] <= row['last_UTC_second'] < optical.CALENDAR_END
            or not 2 <= row['possible_calendar_blocks'] <= optical.CALENDAR_BLOCKS):
        raise ValueError('ordered extrema and a candidate window within the fixed acquisition year are required')
    classes = regional.pair_classes()
    worst = None
    for a, b, common, _ in classes:
        # ratio = gain*a*B(N)/(b*N)+error, so undo the derived mean-rate ratio.
        width = F(b, a)*high*(span+2*error)/gain_lower
        integer_width = -(-width.numerator//width.denominator)
        mass = min(F(1), (integer_width+1)*point_mass_bound(a, b, common, low))
        candidate = (mass, a, b, common, integer_width)
        if worst is None or candidate > worst:
            worst = candidate
    mass, a, b, common, integer_width = worst
    pairs = sum(item[3] for item in classes)
    # Uniformity already covers an unknown pair. This extra family-size charge
    # also covers choosing among distinct pairs after inspecting their records.
    budget = optical.CALENDAR_BLOCKS**3*pairs
    result = selection_bound(row['possible_calendar_blocks'], row['complete_blocks'], mass,
                             choice_budget=budget)
    return {'file': row['file'], 'complete_blocks': row['complete_blocks'],
            'possible_calendar_blocks': row['possible_calendar_blocks'],
            'observed_ratio_range': str(span),
            'worst_incidence_counts': [a, b, common],
            'count_interval_width_upper_at_worst_class': integer_width,
            'single_window_interval_mass_upper_for_every_distinct_pair': str(mass),
            'ordered_distinct_event_set_pair_budget': pairs,
            'calendar_window_and_count_budget': optical.CALENDAR_BLOCKS**3,
            **result,
            'all_distinct_regional_direct_attachments_rejected_at_one_percent':
                result['conservative_probability_upper_power_of_ten'] <= -2}


def servo_interval_mass_bound(a, b, common, count, integer_width):
    """Exact interval-mass bound using binary geometric digits and [-1,1] Fourier split.

    Geom0(p) with failure q=1-p expands into mutually independent binary digits
    sum_{j>=0} 2^j Bernoulli(q^(2^j)/(1+q^(2^j))), independent of Bernoulli(c/a).
    Conditioning on all summands except step s=2^j leaves at most floor(w/s)+1
    values in any interval of width w. On |t|<=1, sin(|t|/2)>=(23/48)|t| and
    8*pi*(23/48)^2>23/4; on 1<=|t|<=pi, sin(t/2)^2>=529/2304 and exp(-x)<=1/(1+x+x²/2+x³/6).
    """
    if (type(count) is not int or count < 2
            or type(integer_width) is not int or integer_width < 0):
        raise ValueError('an integer count of at least two and nonnegative integer width are required')
    law = regional.parameters(a, b, common)
    terminal, p = law['terminal_B_probability'], law['geometric_success_probability']
    q = 1-p
    v_term = terminal*(1-terminal)
    candidates = [min(F(1), (integer_width+1)*point_mass_bound(a, b, common, count))]
    scale = 10**6
    for j in (0, 1, 2):
        step = 1 << j
        qj = q**step
        rj = qj/(1+qj)
        vj = rj*(1-rj)+(v_term if j == 0 else F(0))
        if vj == 0:
            continue
        square = F(4*scale*scale, 23*vj)
        integer = isqrt(square.numerator//square.denominator)
        if integer*integer < square:
            integer += 1
        x = F(529, 1152)*count*vj
        tail = F(1)/(1+x+x*x/2+x*x*x/6)
        pm = min(F(1), F(integer, scale*isqrt(count))+tail)
        candidates.append(min(F(1), (integer_width//step+1)*pm))
    return min(candidates)


def servo_and_allan_analysis(row, *, replicas=8, servo_seconds=30,
                             error=F(1, 10**12), gain_lower=F(1)):
    """Conditional regional range test under R replicas and a causal servo filter,
    plus population-moment targets using the finite-carrier regional floors.
    These sample-derived targets are not empirical confidence bounds.
    """
    error, gain_lower = _fraction(error), _fraction(gain_lower, positive=True)
    if (type(replicas) is not int or replicas < 1
            or type(servo_seconds) is not int or not 0 <= servo_seconds < optical.BLOCK_SECONDS):
        raise ValueError('positive integer replica count and bounded integer servo seconds are required')
    integer_fields = ('block_seconds', 'first_UTC_second', 'last_UTC_second',
                      'possible_calendar_blocks', 'complete_blocks', 'adjacent_complete_block_pairs')
    if (any(type(row.get(key)) is not int for key in integer_fields)
            or row['adjacent_complete_block_pairs'] < 1
            or type(row.get('nonoverlapping_Allan_variance_diagnostic')) is not str):
        raise ValueError('complete integer block counts and a positive Allan diagnostic string are required')
    observed_allan = _fraction(F(row['nonoverlapping_Allan_variance_diagnostic']), positive=True)
    span = F(row['maximum_block_ratio'])-F(row['minimum_block_ratio'])
    if (span < 0 or row['block_seconds'] != optical.BLOCK_SECONDS
            or not optical.CALENDAR_START <= row['first_UTC_second'] <= row['last_UTC_second'] < optical.CALENDAR_END
            or not 2 <= row['possible_calendar_blocks'] <= optical.CALENDAR_BLOCKS
            or not 2 <= row['complete_blocks'] <= row['possible_calendar_blocks']):
        raise ValueError('ordered extrema and a candidate window within the fixed acquisition year are required')
    interior_seconds = optical.BLOCK_SECONDS-servo_seconds
    n_int_per_replica = (optical.REFERENCE_HZ-1)*interior_seconds
    total_interior_ticks = replicas*n_int_per_replica
    high = (optical.REFERENCE_HZ+1)*optical.BLOCK_SECONDS
    classes = regional.pair_classes()
    worst = None
    for a, b, common, _ in classes:
        width = F(b, a)*replicas*high*(span+2*error)/gain_lower
        integer_width = -(-width.numerator//width.denominator)
        mass = servo_interval_mass_bound(a, b, common, total_interior_ticks, integer_width)
        candidate = (mass, a, b, common, integer_width)
        if worst is None or candidate > worst:
            worst = candidate
    mass, a, b, common, integer_width = worst
    pairs = sum(item[3] for item in classes)
    budget = optical.CALENDAR_BLOCKS**3*pairs
    result = selection_bound(row['possible_calendar_blocks'], row['complete_blocks'], mass,
                             choice_budget=budget)
    unequal_floor_min = F(1, 30*high)
    equal_floor_min = F(2, 29*high)
    unequal_ratio = unequal_floor_min/observed_allan
    equal_ratio = equal_floor_min/observed_allan
    r_min_unequal = -(-unequal_ratio.numerator//unequal_ratio.denominator)
    r_min_equal = -(-equal_ratio.numerator//equal_ratio.denominator)
    return {'file': row['file'], 'complete_blocks': row['complete_blocks'],
            'possible_calendar_blocks': row['possible_calendar_blocks'],
            'adjacent_complete_block_pairs': row['adjacent_complete_block_pairs'],
            'replicas_tested': replicas,
            'causal_servo_memory_seconds_upper': servo_seconds,
            'interior_uncontaminated_seconds_per_block': interior_seconds,
            'total_interior_A_ticks_lower': total_interior_ticks,
            'observed_ratio_range': str(span),
            'worst_incidence_counts': [a, b, common],
            'count_interval_width_upper_at_worst_class': integer_width,
            'single_block_conditional_interval_mass_upper_for_every_distinct_pair': str(mass),
            'ordered_distinct_event_set_pair_budget': pairs,
            'calendar_window_and_count_budget': optical.CALENDAR_BLOCKS**3,
            **result,
            'all_distinct_regional_servo_attachments_rejected_at_one_percent':
                result['conservative_probability_upper_power_of_ten'] <= -2,
            'observed_nonoverlapping_Allan_variance': str(observed_allan),
            'single_carrier_regional_Allan_variance_lower': str(unequal_floor_min),
            'single_carrier_equal_rate_regional_Allan_variance_lower': str(equal_floor_min),
            'regional_to_observed_Allan_variance_ratio_floor':
                unequal_ratio.numerator//unequal_ratio.denominator,
            'equal_rate_regional_to_observed_Allan_variance_ratio_floor':
                equal_ratio.numerator//equal_ratio.denominator,
            'moment_matching_independent_regional_replica_count_lower': r_min_unequal,
            'moment_matching_independent_equal_rate_regional_replica_count_lower': r_min_equal,
            'moment_target_requires_unfiltered_zero_error_population_variance_matching': True,
            'sample_Allan_value_is_a_population_upper_confidence_bound': False,
            'thirty_seam_regional_floor_transferred_to_refining_graph': False}


def run():
    summary = optical.strict_json.load(optical.SUMMARY)
    if (set(summary) != {'schema', 'datasets'}
            or summary['schema'] != 'oph.optical_clock.public_summary.v1'
            or [row['file'] for row in summary['datasets']] != list(optical.SOURCES)
            or any(row['sha256'] != optical.SOURCES[row['file']] for row in summary['datasets'])):
        raise ValueError('the retained summaries must bind both named public datasets')
    return {'schema': 'oph.optical_clock.regional_attachment.v1',
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'public_dataset': 'https://zenodo.org/records/6901524',
        'measurement_protocol': 'https://arxiv.org/abs/2205.13053',
        'reference_frequency': 'https://arxiv.org/abs/2010.06620',
        'hypothesis': {
            'source': 'one native twelve-port thirty-seam carrier with IID uniform seam selection',
            'clock': 'one tick per event incident to any port in a declared nonempty region',
            'regional_union_order_and_event_coincidence_access_admitted': True,
            'initial_load_or_quantum_coherence_affects_clock_law': False,
            'attachment': 'reported block ratio = constant_gain*a*B(N)/(b*N) + bounded_error',
            'a_b_common': 'numbers of A events, B events, and common events among thirty seams',
            'constant_gain_lower': '1', 'absolute_ratio_error_per_block_upper': '1/1000000000000',
            'A_ticks_per_block': 'one fixed integer N in the declared frequency interval times 300 seconds',
            'all_calendar_blocks': 'each consecutive 300-second candidate block contains exactly N consecutive A intervals, including unobserved blocks',
            'reference_frequency_Hz_interval': [optical.REFERENCE_HZ-1, optical.REFERENCE_HZ+1],
            'frequency_interval_is_a_declared_hypothesis_not_a_hard_measurement_bound': True,
            'extra_temporal_filtering_or_independent_replica_averaging_included': False,
            'arbitrary_missing_or_discarded_blocks_allowed': True,
            'calendar_window_count_and_regional_pair_selection_charged': True,
            'fixed_calendar_year': 2021, 'search_outside_year_or_other_block_duration_covered': False,
            'actual_atomic_transition_and_servo_dictionary_established': False},
        'public_summary_source_sha256': hashlib.sha256(optical.SUMMARY.read_bytes()).hexdigest(),
        'family': {'nonempty_port_regions': 4095, 'distinct_event_sets': len(regional._catalogue()),
            'incidence_classes': len(regional.pair_classes()),
            'clock_aliases_excluded_from_distinct_pair_test': True,
            'alias_boundary': 'B(N)=N pathwise for identical event sets; the test does not reject those aliases'},
        'conditional_tests': [evaluate(row) for row in summary['datasets']],
        'finite_memory_servo_hypothesis': optical.servo_hypothesis(),
        'servo_and_allan_tests': [servo_and_allan_analysis(row) for row in summary['datasets']],
        'assessment': {'retrospective_conditional_attachment_test': True,
            'variance_alone_used_as_rejection': False,
            'shared_same_carrier_overlap_fully_enumerated': True,
            'atomic_clocks_or_OPH_excluded': False,
            'coherence_or_independent_atom_resource_lower_bound_derived': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        optical.compare_receipt(optical.strict_json.load(args.check), result)
        print('Conditional collective regional clock comparison replayed.')
    else:
        print(optical.canonical(result), end='')


if __name__ == '__main__':
    main()
