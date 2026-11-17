"""A conditional public-data test of a direct native-event clock attachment.

A clock is defined internally by its counted events. The tested hypothesis
identifies two such counts with reported optical phase, up to one constant
gain and a bounded readout error. It is deliberately stronger than defining
a clock. Rejection concerns this attachment, not every OPH clock model.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import re

from ..laboratory_clock import event_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT/'model/predictions/data/optical_clock_summary.json'
SOURCES = {
    '36ar_yb_ratio.txt': '4852780f7d506f811c1bd5eaa52596949d89107d395f6a40c549db97f71e701f',
    '40ar_yb_ratio.txt': '21e297071cacf567fc7ce0f9041fabe21b8a1365c97d93b08d55e93351d5e32d',
}
BLOCK_SECONDS = 300
REFERENCE_HZ = 642121496772645
CALENDAR_START = int(datetime(2021, 1, 1, tzinfo=timezone.utc).timestamp())
CALENDAR_END = int(datetime(2022, 1, 1, tzinfo=timezone.utc).timestamp())
CALENDAR_BLOCKS = (CALENDAR_END-CALENDAR_START)//BLOCK_SECONDS
SOURCE_PATHS = tuple(sorted(set((*event_geometry.SOURCE_PATHS,
    'model/predictions/optical_clock.py', 'model/predictions/__init__.py',
    'model/predictions/data/optical_clock_summary.json',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py'))))


def summarize(content, name):
    """Reduce all published rows; retain only complete fixed UTC-grid blocks.

    No value-based filter is applied. Published removals remain gaps. The
    range test below allows arbitrary removal of whole candidate blocks.
    """
    if name not in SOURCES or type(content) is not bytes:
        raise ValueError('a named public ratio file as original bytes is required')
    lines = content.decode('ascii').splitlines()
    isotope = name[:2]
    if not lines or lines[0] != f'Date Time(UTC) Ratio_{isotope}ar_yb':
        raise ValueError('the published ratio header is required')
    bins, previous, count, first = defaultdict(list), None, 0, None
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 3 or not re.fullmatch(r'[0-9]+\.[0-9]+', fields[2]):
            raise ValueError('each row needs a UTC second and a finite decimal ratio')
        stamp = fields[0]+'T'+fields[1]
        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}', stamp):
            raise ValueError('original whole-second UTC timestamps are required')
        second = int(datetime.fromisoformat(stamp).replace(tzinfo=timezone.utc).timestamp())
        ratio = F(fields[2])
        if not CALENDAR_START <= second < CALENDAR_END:
            raise ValueError('the tested acquisition calendar is the fixed year 2021')
        if previous is not None and second <= previous:
            raise ValueError('timestamps must be unique and strictly increasing')
        if not 1 <= ratio <= 2:
            raise ValueError('the declared optical-ratio range is one through two')
        if first is None:
            first = second
        bins[second//BLOCK_SECONDS].append((second, ratio))
        previous, count = second, count+1
    full = {i: sum((r for _, r in rows), F(0))/BLOCK_SECONDS
            for i, rows in bins.items() if len(rows) == BLOCK_SECONDS
            and rows[0][0] == BLOCK_SECONDS*i
            and rows[-1][0] == BLOCK_SECONDS*i+BLOCK_SECONDS-1}
    if len(full) < 2:
        raise ValueError('at least two complete UTC-grid blocks are required')
    adjacent = [(full[i+1]-value)**2 for i, value in full.items() if i+1 in full]
    mean = sum(full.values(), F(0))/len(full)
    return {'file': name, 'sha256': hashlib.sha256(content).hexdigest(),
            'published_rows': count, 'first_UTC_second': first, 'last_UTC_second': previous,
            'block_seconds': BLOCK_SECONDS, 'complete_blocks': len(full),
            'possible_calendar_blocks': max(bins)-min(bins)+1,
            'minimum_block_ratio': str(min(full.values())),
            'maximum_block_ratio': str(max(full.values())),
            'mean_block_ratio': str(mean),
            'adjacent_complete_block_pairs': len(adjacent),
            'nonoverlapping_Allan_variance_diagnostic': (
                str(sum(adjacent, F(0))/(2*len(adjacent)*mean**2)) if adjacent else None),
            'additional_value_based_filter': False}


def verify_raw(directory):
    rows = []
    for name, expected in SOURCES.items():
        content = (Path(directory)/name).read_bytes()
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError('the public dataset bytes differ from the verified source')
        rows.append(summarize(content, name))
    return {'schema': 'oph.optical_clock.public_summary.v1', 'datasets': rows}


def max_point_mass_bound(count):
    """Uniform bound for either native pair: sup_k P[B(count)=k] <= 2/sqrt(count).

    X=Bernoulli(delta/5)+Geom0(5/(10-delta)). The geometric parity bit is
    independent of its even remainder and has success r=q/(1+q), where
    q is 4/9 or 1/2. Thus B(count) is a convolution with Binomial(count,r).
    Fourier inversion, 1-x<=exp(-x), and sin(|t|/2)>=|t|/pi give its maximum
    mass <= sqrt(pi)/(2*sqrt(2*count*r*(1-r))) < 2/sqrt(count).
    Replacing sqrt(count) by its integer floor preserves an exact upper bound.
    """
    if type(count) is not int or count < 2:
        raise ValueError('an original integer count of at least two is required')
    for delta in (0, 1):
        parameters = event_geometry.renewal_parameters(delta)
        p = parameters['geometric_success_probability']
        r = (1-p)/(2-p)
        if (parameters['terminal_B_probability'] != F(delta, 5)
                or p != F(5, 10-delta) or parameters['mean'] != 1
                or parameters['variance'] != 2-F(2*delta, 5) or 8*r*(1-r) <= 1):
            raise ArithmeticError('the native renewal decomposition lost its uniform bound')
    return min(F(1), F(2, isqrt(count)))


def _fraction(value, *, positive=False):
    value = event_geometry.exact_core._fraction(value)
    if value < 0 or (positive and value == 0):
        raise ValueError('a canonical nonnegative exact fraction is required')
    return value


def selection_bound(total, retained, interval_mass, *, choice_budget=1):
    """Any retained subset of independent draws with a uniform interval-mass bound.

    Anchor the interval at its smallest retained value, union over all total
    anchors and (retained-1)-subsets, then use binomial(n,k)<=(3n/k)^k.
    The constant gain/centre may be chosen after seeing the values. A further
    union budget can cover choosing the calendar window and retained count.
    """
    interval_mass = _fraction(interval_mass, positive=True)
    if (type(total) is not int or type(retained) is not int
            or not 2 <= retained <= total or interval_mass > 1
            or type(choice_budget) is not int or choice_budget < 1):
        raise ValueError('valid integer window counts and exact positive mass are required')
    base = min(F(1), 3*(total-1)*interval_mass/(retained-1))
    power = 0
    while base*2**(power+1) <= 1:
        power += 1
    exponent = power*(retained-1)
    if exponent == 0 and base < 1:
        product = base**(retained-1)
        exponent = max(0, product.denominator.bit_length()-product.numerator.bit_length()-1)
        while product*2**(exponent+1) <= 1:
            exponent += 1
    prefactor = total*choice_budget
    # 2**10>10**3; exact integer digit count bounds the prefactor.
    decimal_exponent = min(0, len(str(prefactor))-3*(exponent//10))
    return {'interval_occupancy_union_base_upper': str(base),
            'calendar_and_count_selection_union_budget': choice_budget,
            'binary_power_per_retained_draw': power,
            'probability_upper_form': f'min(1,{prefactor}*2^(-{exponent}))',
            'conservative_probability_upper_power_of_ten': decimal_exponent}


def evaluate(row, *, error=F(1, 10**12), gain_lower=F(1)):
    error, gain_lower = _fraction(error), _fraction(gain_lower, positive=True)
    integer_fields = ('block_seconds', 'first_UTC_second', 'last_UTC_second',
                      'possible_calendar_blocks', 'complete_blocks')
    if any(type(row[key]) is not int for key in integer_fields):
        raise ValueError('calendar and block counts must be original integers')
    low = (REFERENCE_HZ-1)*BLOCK_SECONDS
    high = (REFERENCE_HZ+1)*BLOCK_SECONDS
    span = F(row['maximum_block_ratio'])-F(row['minimum_block_ratio'])
    if span < 0 or row['block_seconds'] != BLOCK_SECONDS:
        raise ValueError('ordered extrema for the fixed observation window are required')
    width = high*(span+2*error)/gain_lower
    integer_width = -(-width.numerator//width.denominator)
    mass = min(F(1), (integer_width+1)*max_point_mass_bound(low))
    if (not CALENDAR_START <= row['first_UTC_second'] <= row['last_UTC_second'] < CALENDAR_END
            or not 2 <= row['possible_calendar_blocks'] <= CALENDAR_BLOCKS):
        raise ValueError('the candidate calendar window must lie within the fixed acquisition year')
    # At most Y² calendar endpoint pairs and Y retained-count choices. This
    # charges selection of the observed span, start, stop and subset size.
    result = selection_bound(row['possible_calendar_blocks'], row['complete_blocks'], mass,
                             choice_budget=CALENDAR_BLOCKS**3)
    return {'file': row['file'], 'complete_blocks': row['complete_blocks'],
            'possible_calendar_blocks': row['possible_calendar_blocks'],
            'observed_ratio_range': str(span), 'count_interval_width_upper': integer_width,
            'single_window_interval_mass_upper': str(mass), **result,
            'attachment_rejected_at_one_percent_under_stated_hypotheses':
                result['conservative_probability_upper_power_of_ten'] <= -2}


def servo_and_refinement_analysis(row, *, replicas=8, servo_seconds=30,
                                  error=F(1, 10**12), gain_lower=F(1)):
    """Conditional range test under R replicas and a finite-memory causal filter,
    plus a population-moment target calculated from the sample Allan diagnostic.

    The compact impulse-response support, equal tick spacing, independent
    replicas and fixed preselected filter are hypotheses; a servo time constant
    does not establish any of them. The moment targets are not confidence bounds.
    """
    error, gain_lower = _fraction(error), _fraction(gain_lower, positive=True)
    if (type(replicas) is not int or replicas < 1
            or type(servo_seconds) is not int or not 0 <= servo_seconds < BLOCK_SECONDS):
        raise ValueError('positive integer replica count and bounded integer servo seconds are required')
    integer_fields = ('block_seconds', 'first_UTC_second', 'last_UTC_second',
                      'possible_calendar_blocks', 'complete_blocks', 'adjacent_complete_block_pairs')
    if (any(type(row.get(key)) is not int for key in integer_fields)
            or row['adjacent_complete_block_pairs'] < 1
            or type(row.get('nonoverlapping_Allan_variance_diagnostic')) is not str):
        raise ValueError('complete integer block counts and a positive Allan diagnostic string are required')
    observed_allan = _fraction(F(row['nonoverlapping_Allan_variance_diagnostic']), positive=True)
    span = F(row['maximum_block_ratio'])-F(row['minimum_block_ratio'])
    if (span < 0 or row['block_seconds'] != BLOCK_SECONDS
            or not CALENDAR_START <= row['first_UTC_second'] <= row['last_UTC_second'] < CALENDAR_END
            or not 2 <= row['possible_calendar_blocks'] <= CALENDAR_BLOCKS
            or not 2 <= row['complete_blocks'] <= row['possible_calendar_blocks']):
        raise ValueError('ordered extrema and a valid candidate calendar window in 2021 are required')
    interior_seconds = BLOCK_SECONDS-servo_seconds
    n_int_per_replica = (REFERENCE_HZ-1)*interior_seconds
    total_interior_ticks = replicas*n_int_per_replica
    high = (REFERENCE_HZ+1)*BLOCK_SECONDS
    width = replicas*high*(span+2*error)/gain_lower
    integer_width = -(-width.numerator//width.denominator)
    mass = min(F(1), (integer_width+1)*max_point_mass_bound(total_interior_ticks))
    result = selection_bound(row['possible_calendar_blocks'], row['complete_blocks'], mass,
                             choice_budget=CALENDAR_BLOCKS**3)
    adjacent_min = F(8, 5*high)
    nonadjacent_min = F(2, high)
    adj_ratio = adjacent_min/observed_allan
    nonadj_ratio = nonadjacent_min/observed_allan
    r_min_allan = -(-adj_ratio.numerator//adj_ratio.denominator)
    nu_ratio = F(8, 5*BLOCK_SECONDS)/observed_allan
    nu_eff_min = -(-nu_ratio.numerator//nu_ratio.denominator)
    return {'file': row['file'], 'complete_blocks': row['complete_blocks'],
            'possible_calendar_blocks': row['possible_calendar_blocks'],
            'adjacent_complete_block_pairs': row['adjacent_complete_block_pairs'],
            'replicas_tested': replicas,
            'causal_servo_memory_seconds_upper': servo_seconds,
            'interior_uncontaminated_seconds_per_block': interior_seconds,
            'interior_A_ticks_lower_per_replica': n_int_per_replica,
            'total_interior_A_ticks_lower': total_interior_ticks,
            'observed_ratio_range': str(span),
            'interior_count_interval_width_upper': integer_width,
            'single_block_conditional_interval_mass_upper': str(mass), **result,
            'servo_filtered_multi_replica_attachment_rejected_at_one_percent':
                result['conservative_probability_upper_power_of_ten'] <= -2,
            'observed_nonoverlapping_Allan_variance': str(observed_allan),
            'single_carrier_adjacent_Allan_variance_lower': str(adjacent_min),
            'single_carrier_nonadjacent_Allan_variance_lower': str(nonadjacent_min),
            'adjacent_to_observed_Allan_variance_ratio_floor':
                adj_ratio.numerator//adj_ratio.denominator,
            'nonadjacent_to_observed_Allan_variance_ratio_floor':
                nonadj_ratio.numerator//nonadj_ratio.denominator,
            'moment_matching_independent_native_replica_count_lower': r_min_allan,
            'moment_matching_unfiltered_A_tick_frequency_Hz_lower': nu_eff_min,
            'moment_target_requires_unfiltered_zero_error_population_variance_matching': True,
            'sample_Allan_value_is_a_population_upper_confidence_bound': False,
            'continuum_refinement_level_inferred_from_public_data': False}


def clock_field_consistency():
    """Recover the native density transfer from coincident incident-clock moments.

    Every attempt marks its two incident ports, including zero-current attempts.
    With p=E[I], Sigma=Cov(I), fair endpoint repair has mean matrix
    P=Id-diag(p)+(Sigma+p*p^T)/2. All ordered coincident port tapes are admitted;
    the union of those tapes supplies the event index, not host time.
    """
    seams = event_geometry._source().seams
    size, count = event_geometry.PORTS, len(seams)
    means = [F(sum(i in edge for edge in seams), count) for i in range(size)]
    second = [[F(sum(i in edge and j in edge for edge in seams), count)
               for j in range(size)] for i in range(size)]
    covariance = [[second[i][j]-means[i]*means[j] for j in range(size)] for i in range(size)]
    transfer = [[F(i == j)*(1-means[i])+second[i][j]/2
                 for j in range(size)] for i in range(size)]
    branches = 0
    for state in event_geometry.native.states(binary=True):
        actual = [F(0) for _ in range(size)]
        for edge in seams:
            for coin in (0, 1):
                after, _ = event_geometry.native.native(state, edge, coin)
                actual = [value+F(int(after[i]), 2*count) for i, value in enumerate(actual)]
                branches += 1
        predicted = [sum((transfer[i][j]*int(state[j]) for j in range(size)), F(0))
                     for i in range(size)]
        if actual != predicted:
            raise ArithmeticError('clock moments did not predict the complete native density transfer')
    adjacent = set(seams) | {(b, a) for a, b in seams}
    if any(covariance[i][j] != (F(5, 36) if i == j else F(1, 180) if (i, j) in adjacent else -F(1, 36))
           for i in range(size) for j in range(size)):
        raise ArithmeticError('the complete incident-clock covariance changed')
    return {'mean_transfer_from_clock_moments': 'P=Id-diag(p)+(Sigma+p*p^T)/2',
            'clock_tick': 'every incident event including zero current',
            'event_order': 'ordered union of all coincident incident-clock tapes',
            'port_mean_tick_probability': list(map(str, means)),
            'clock_covariance_by_pair': {'same': '5/36', 'seam': '1/180', 'distinct_nonseam': '-1/36'},
            'recovered_density_transfer_by_pair': {'same': '11/12', 'seam': '1/60', 'distinct_nonseam': '0'},
            'native_binary_histories_checked': branches,
            'finite_sample_example': clock_field_sampling_bound(3900000, F(1, 100)),
            'second_field_rate_fitted': False,
            'complete_clock_coincidence_access_and_fair_repair_required': True,
            'public_optical_ratio_files_supply_these_joint_records': False}


def clock_field_sampling_bound(events, tolerance):
    """IID clock-coincidence sampling bound for the inferred mean field matrix.

    Estimate Q_ij=E[I_i I_j] from one fixed window, with Q_ii=E[I_i].
    There are 78 distinct bounded indicators. Chebyshev and a union bound give
    simultaneous entry error <=delta except with probability78/(4*N*delta²).
    P=Id-diag(diag Q)+Q/2 then has entry error <=delta/2 and spectral-norm
    error <=6*delta. No independence between different indicator pairs is used.
    """
    tolerance = _fraction(tolerance, positive=True)
    if type(events) is not int or events < 1 or tolerance > 1:
        raise ValueError('a positive integer sample size and rational tolerance at most one are required')
    return {'complete_union_events': events, 'second_moment_entry_tolerance': str(tolerance),
            'failure_probability_upper': str(min(F(1), F(78, 4*events*tolerance**2))),
            'density_transfer_operator_error_upper': str(6*tolerance),
            'stationary_density_preparation_needed_for_mean_transfer': False,
            'IID_seam_events_and_fixed_unselected_window_required': True}


def refined_incident_clock_scale(level, degree):
    """Expected local ticks in the declared density-continuum coordinate.

    The source has M_l=630*4**l seams and a_l=1260*16**l attempts per unit.
    A fixed port touches degree seams, so its incidence probability shrinks
    with refinement. This mean is not a deterministic count or a lab rate.
    """
    if type(level) is not int or not 0 <= level <= 100 or type(degree) is not int or degree not in (5, 6):
        raise ValueError('a bounded refinement level and native port degree five or six are required')
    seams, attempts = 630*4**level, 1260*16**level
    probability = F(degree, seams)
    return {'level': level, 'degree': degree, 'source_seams': seams,
            'native_attempts_per_density_continuum_unit': attempts,
            'incident_probability': str(probability),
            'mean_incident_ticks_per_density_continuum_unit': str(attempts*probability)}


def servo_hypothesis(replicas=8, seconds=30):
    return {'independent_identically_distributed_native_replicas': replicas,
            'replica_outputs_averaged_with_equal_weights': True,
            'deterministic_filter_fixed_before_data': True,
            'filter': 'linear causal nonnegative unit-integral impulse response',
            'impulse_response_support_seconds': [0, seconds],
            'A_ticks_equally_spaced_in_declared_calendar_coordinate': True,
            'all_outside_interior_tick_increments_conditioned_on': True,
            'conditional_blocks_independent_with_arbitrary_outside_shifts': True,
            'published_servo_time_constant_establishes_compact_support': False,
            'eight_published_servo_outputs_establish_independent_native_replicas': False,
            'actual_atomic_servo_dictionary_established': False}


def run():
    summary = strict_json.load(SUMMARY)
    if (set(summary) != {'schema', 'datasets'} or summary['schema'] != 'oph.optical_clock.public_summary.v1'
            or [row['file'] for row in summary['datasets']] != list(SOURCES)
            or any(row['sha256'] != SOURCES[row['file']] for row in summary['datasets'])):
        raise ValueError('the retained summaries must bind both named public datasets')
    return {'schema': 'oph.optical_clock.direct_attachment.v1',
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
            'public_dataset': 'https://zenodo.org/records/6901524',
            'measurement_protocol': 'https://arxiv.org/abs/2205.13053',
            'reference_frequency': 'https://arxiv.org/abs/2010.06620',
            'data': summary,
            'hypothesis': {'clock_definition': 'each incident repair ticks, including zero current',
                'comparison': 'count B between successive A ticks, including coincident terminal ticks',
                'scheduler': 'independent uniform native seam selection',
                'optical_attachment': 'reported block ratio = constant_gain*B(N)/N + bounded_error',
                'reference_frequency_Hz_interval': [REFERENCE_HZ-1, REFERENCE_HZ+1],
                'A_ticks_per_block': 'one fixed integer N in the frequency interval times 300 seconds',
                'constant_gain_lower': '1', 'absolute_ratio_error_per_block_upper': '1/1000000000000',
                'arbitrary_missing_or_discarded_blocks_allowed': True,
                'window_endpoints_and_retained_count_selection_charged_within_calendar_year': 2021,
                'searches_outside_declared_calendar_year_covered': False,
                'extra_temporal_filtering_or_collective_clock_mechanism_included': False,
                'bounded_error_and_attachment_are_tested_hypotheses_not_established_physics': True},
            'conditional_tests': [evaluate(row) for row in summary['datasets']],
            'joint_clock_field_constraint': clock_field_consistency(),
            'finite_memory_servo_hypothesis': servo_hypothesis(),
            'servo_and_refinement_tests': [servo_and_refinement_analysis(row) for row in summary['datasets']],
            'density_continuum_clock_normalization': {
                'mean_fixed_port_ticks_per_coordinate_unit': '2*d*4**l',
                'physical_seconds_per_coordinate_unit_calibrated': False,
                'primitive_renewal_coefficient_for_equal_degree_adjacent_ports': '2-2/d',
                'examples': [refined_incident_clock_scale(level, degree)
                             for level in (0, 1, 2) for degree in (5, 6)]},
            'predicted_relative_variance_per_N_A_ticks': {'adjacent': '8/(5*N)', 'nonadjacent': '2/N'},
            'assessment': {'retrospective_conditional_attachment_test': True,
                'OPH_or_every_native_clock_excluded': False,
                'native_atomic_transition_or_interacting_QFT_derived': False,
                'new_distinction_from_SM_GR_established': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, help='independently re-extract both original public files')
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.raw_dir:
        compare_receipt(strict_json.load(SUMMARY), verify_raw(args.raw_dir))
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Conditional optical event-clock attachment comparison replayed.')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
