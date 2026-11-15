"""Exact native current waiting/count laws, with complete CQ execution checks.

The count parameter is the number of native attempts. The current counter is
an observer-accessible statistic, not a calibrated physical duration. No
independent clock process, Hamiltonian, fit, or field action is introduced.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from ..observer import native
from ..primitive_observer import runtime, transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/native_clock/clock.py', 'model/native_clock/__init__.py',
    'model/native_clock/__main__.py', 'model/common_source/protocol.py',
    *transfer.SOURCE_PATHS,
    'Lean/QFT/SourceRecurrenceClock.lean',
    'Lean/Time/SourceCountClockEnclosure.lean',
)))
TOLERANCE = 2e-11


def _bound(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= TOLERANCE:
        raise ArithmeticError(f'{name} exceeds the native clock tolerance')
    return value


def _horizon(value, maximum=256):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(f'an integer native-attempt horizon in [0,{maximum}] required')
    return value


def _row_times(row, matrix):
    return tuple(sum((weight*matrix[i][j] for i, weight in enumerate(row) if weight), F(0))
                 for j in range(len(matrix)))


def _add(a, b):
    return tuple(x+y for x, y in zip(a, b, strict=True))


def _sum_matrix(a, b):
    return tuple(_add(x, y) for x, y in zip(a, b, strict=True))


def _probability(row):
    return sum(row, F(0))


def _prior(prior, size=3):
    if (type(prior) is not tuple or len(prior) != size
            or any(type(x) is not F or x < 0 for x in prior) or sum(prior) != 1):
        raise ValueError('a complete exact normalized probability row is required')
    return prior


def lumping(kernel):
    """Prove marked strong lumping by all actual source rows, not a fit."""
    transfer.exact_audit(kernel)
    adjacency = {i: set() for i in range(12)}
    for a, b in kernel.seams:
        adjacency[a].add(b)
        adjacency[b].add(a)
    classes = []
    for state in kernel.states:
        a, b = (i for i, value in enumerate(state) if value)
        classes.append(0 if b in adjacency[a] else 1 if adjacency[a] & adjacency[b] else 2)
    sizes = tuple(classes.count(c) for c in range(3))
    if sizes != (30, 30, 6):
        raise ArithmeticError('the source token-separation partition changed')
    marked = []
    for matrix in kernel.marked_counts:
        reduced = []
        for c in range(3):
            rows = {tuple(sum(matrix[i][j] for j in range(66) if classes[j] == d)
                          for d in range(3)) for i in range(66) if classes[i] == c}
            if len(rows) != 1:
                raise ArithmeticError('marked transition rows are not strongly lumpable')
            reduced.append(tuple(F(x, 60) for x in rows.pop()))
        marked.append(tuple(reduced))
    stationary = tuple(F(size, 66) for size in sizes)
    matrix = _sum_matrix(*marked)
    if any(sum(row) != 1 for row in matrix) or _row_times(stationary, matrix) != stationary:
        raise ArithmeticError('the reduced native stationary law changed')
    move_mass = _probability(_row_times(stationary, marked[1]))
    palm = tuple(x/move_mass for x in _row_times(stationary, marked[1]))
    if move_mass != F(5, 33):
        raise ArithmeticError('the native stationary moving current rate changed')
    return {'classes': tuple(classes), 'sizes': sizes, 'marked': tuple(marked),
            'matrix': matrix, 'stationary': stationary, 'palm': palm, 'move_mass': move_mass}


def word_law(reduction, word, prior=None):
    """Unnormalized terminal class populations of every retained bit word."""
    if (type(word) is not tuple or len(word) > 256
            or any(type(bit) is not int or bit not in (0, 1) for bit in word)):
        raise ValueError('a bounded tuple of native current bits is required')
    row = _prior(reduction['stationary'] if prior is None else prior)
    for bit in word:
        row = _row_times(row, reduction['marked'][bit])
    return row


def count_law(reduction, attempts):
    """All count probabilities, including the zero-current outcome."""
    attempts = _horizon(attempts)
    states = [reduction['stationary']]
    zero = (F(0),)*3
    for _ in range(attempts):
        future = [zero]*(len(states)+1)
        for count, row in enumerate(states):
            future[count] = _add(future[count], _row_times(row, reduction['marked'][0]))
            future[count+1] = _add(future[count+1], _row_times(row, reduction['marked'][1]))
        states = future
    probabilities = tuple(map(_probability, states))
    mean = sum((F(i)*p for i, p in enumerate(probabilities)), F(0))
    variance = sum(((F(i)-mean)**2*p for i, p in enumerate(probabilities)), F(0))
    if sum(probabilities) != 1 or mean != attempts*reduction['move_mass']:
        raise ArithmeticError('the complete native count law lost mass or mean')
    return probabilities, mean, variance


def waiting_law(reduction, attempts, *, start='stationary'):
    """First next current, with the complete finite-horizon no-current tail."""
    attempts = _horizon(attempts)
    if start not in ('stationary', 'palm'):
        raise ValueError('start must be stationary or palm')
    prior = reduction[start]
    probabilities = tuple(_probability(word_law(reduction, (0,)*(n-1)+(1,), prior))
                          for n in range(1, attempts+1))
    tail = _probability(word_law(reduction, (0,)*attempts, prior))
    if sum(probabilities)+tail != 1:
        raise ArithmeticError('the no-current tail was discarded')
    return probabilities, tail


def interval_statistics(reduction):
    """Geometric waiting conditional on source class; not a renewal ansatz."""
    no, move = reduction['marked']
    escape = tuple(sum(row) for row in move)
    if any(no[i][j] != (1-escape[i] if i == j else 0)
           for i in range(3) for j in range(3)):
        raise ArithmeticError('a no-current branch changed configuration class')
    jump = tuple(tuple(value/escape[i] for value in row) for i, row in enumerate(move))
    palm = reduction['palm']
    if _row_times(palm, jump) != palm:
        raise ArithmeticError('the event-conditioned law is not jump-stationary')
    means = tuple(1/p for p in escape)
    mean = sum((p*m for p, m in zip(palm, means, strict=True)), F(0))
    second = sum((p*(2*m*m-m) for p, m in zip(palm, means, strict=True)), F(0))
    adjacent_product = sum((palm[i]*means[i]*jump[i][j]*means[j]
                            for i in range(3) for j in range(3)), F(0))
    residual_mean = sum((p*m for p, m in zip(reduction['stationary'], means, strict=True)), F(0))
    conditional_after_one = _probability(word_law(reduction, (1, 1)))/reduction['move_mass']
    conditional_after_two = (_probability(word_law(reduction, (1, 1, 1)))
                             /_probability(word_law(reduction, (1, 1))))
    if mean != 1/reduction['move_mass'] or conditional_after_one == conditional_after_two:
        raise ArithmeticError('the exact native event conditioning witness changed')
    return {'escape': escape, 'jump': jump, 'mean_interval': mean,
            'interval_variance': second-mean*mean,
            'adjacent_interval_covariance': adjacent_product-mean*mean,
            'stationary_start_wait_mean': residual_mean,
            'P_next_1_given_last_1': conditional_after_one,
            'P_next_1_given_last_two_11': conditional_after_two}


def duration_scale_audit(reduction):
    """Certify that native attempt records contain no seconds scale.

    Candidate seconds-per-attempt values are deliberately hypothetical. Every
    positive scale assigns different durations and rates to the same exact
    dimensionless record law; the source supplies no criterion selecting one.
    """
    probabilities, mean, _ = count_law(reduction, 1)
    expected = (F(28, 33), F(5, 33))
    if probabilities != expected or mean != F(5, 33):
        raise ArithmeticError('the native one-attempt law changed during scale audit')
    candidates = (F(1, 100), F(1), F(100))
    rows = []
    for seconds_per_attempt in candidates:
        rows.append({
            'illustrative_seconds_per_attempt': str(seconds_per_attempt),
            'implied_stationary_event_rate_per_second': str(F(5, 33)/seconds_per_attempt),
            'implied_Palm_mean_interval_seconds': str(F(33, 5)*seconds_per_attempt),
            'implied_stationary_start_wait_seconds': str(F(147, 22)*seconds_per_attempt),
            'exact_native_one_attempt_record_law': list(map(str, probabilities)),
        })
    if len({row['implied_stationary_event_rate_per_second'] for row in rows}) != len(rows):
        raise ArithmeticError('the duration rescaling witness lost distinct physical rates')
    return {
        'native_parameter': 'integer attempt index n',
        'dimensionful_source_quantity_with_seconds_units': False,
        'same_exact_record_law_for_every_positive_seconds_per_attempt': True,
        'candidate_scales_are_estimates': False,
        'illustrative_rescalings': rows,
        'seconds_per_attempt_identified': False,
        'matched_public_lab_dataset_identified': False,
        'reason_no_direct_dataset_test_is_defined': (
            'no laboratory apparatus has been operationally mapped to these native attempts and current bits'),
        'minimal_calibration_if_a_constant_attempt_clock_is_validated': (
            'one nonzero paired native-attempt difference and independently measured lab interval; '
            'a single pair fixes only that interval unless constant-rate t=t0+tau*n is justified'),
    }


def weighted_native_marks(kernel, weights, *, activity=F(1)):
    """Actual native branches for a fixed seam law, with optional idle control."""
    transfer.exact_audit(kernel)
    weights = _prior(weights, 30)
    if type(activity) is not F or not 0 < activity <= 1:
        raise ValueError('an exact activity probability in (0,1] is required')
    index = {state: i for i, state in enumerate(kernel.states)}
    marked = [[[F(0)]*66 for _ in kernel.states] for _ in (0, 1)]
    for i, state in enumerate(kernel.states):
        marked[0][i][i] = 1-activity
        for seam, weight in zip(kernel.seams, weights, strict=True):
            for coin in (0, 1):
                output, record = native.native(state, seam, coin)
                marked[abs(record[1])][i][index[output]] += activity*weight/2
    return tuple(tuple(tuple(row) for row in matrix) for matrix in marked)


def scheduler_controls(kernel):
    """What event normalization removes, and what an observable still detects."""
    uniform = (F(1, 30),)*30
    biased = (F(2, 31),)+(F(1, 31),)*29
    ordinary = weighted_native_marks(kernel, uniform)
    lazy = weighted_native_marks(kernel, uniform, activity=F(1, 2))
    alternate = weighted_native_marks(kernel, biased)
    prior = (F(1, 66),)*66
    def facts(marked):
        no, move = marked
        matrix = _sum_matrix(no, move)
        if _row_times(prior, matrix) != prior or any(sum(row) != 1 for row in matrix):
            raise ArithmeticError('the actual weighted native kernel lost stationarity')
        escape = tuple(sum(row) for row in move)
        if min(escape) <= 0:
            raise ArithmeticError('the proposed event clock may never advance')
        jump = tuple(tuple(value/escape[i] for value in row) for i, row in enumerate(move))
        rate = sum(escape, F(0))/66
        palm = tuple(value/rate for value in _row_times(prior, move))
        if _row_times(palm, jump) != palm:
            raise ArithmeticError('the Palm distribution changed under jump evolution')
        both = _probability(_row_times(_row_times(prior, move), move))
        return rate, both, jump, sum((p/r for p, r in zip(palm, escape, strict=True)), F(0))
    first, second, third = map(facts, (ordinary, lazy, alternate))
    if first[2] != second[2] or first[2] == third[2] or first[0] != third[0]:
        raise ArithmeticError('the scheduler counterexample or lazy invariance changed')
    # The rate identity is checked on every individual seam. It therefore
    # holds for every fixed normalized convex combination of these laws.
    seam_rates = []
    for seam in kernel.seams:
        moving = sum(abs(native.native(state, seam, coin)[1][1])
                     for state in kernel.states for coin in (0, 1))
        seam_rates.append(F(moving, 132))
    if set(seam_rates) != {F(5, 33)}:
        raise ArithmeticError('a native seam violates the scheduler-free stationary mean')
    single = weighted_native_marks(kernel, (F(1),)+(F(0),)*29)
    single_escape = tuple(sum(row) for row in single[1])
    inactive = sum(rate == 0 for rate in single_escape)
    single_rate = sum(single_escape, F(0))/66
    # Event conditioning cannot recover stationary mass in components that
    # never produce an event. This control prevents an overbroad 1/r claim.
    single_palm = tuple(p/single_rate for p in _row_times(prior, single[1]))
    single_interval = sum((p/r for p, r in zip(single_palm, single_escape, strict=True) if r), F(0))
    if (inactive, single_rate, single_interval) != (46, F(5, 33), F(2)):
        raise ArithmeticError('the inactive-component clock counterexample changed')
    return {'all_30_single_seam_stationary_moving_rates': sorted(set(map(str, seam_rates))),
            'fixed_seam_weights_need_not_be_uniform_for_mean_rate': True,
            'stationary_mean_is_not_a_full_process_identification': True,
            'uniform_rate': str(first[0]), 'half_active_rate': str(second[0]),
            'uniform_P_11': str(first[1]), 'half_active_P_11': str(second[1]),
            'double_first_seam_weights': list(map(str, biased)),
            'double_first_seam_rate': str(third[0]), 'double_first_seam_P_11': str(third[1]),
            'double_first_seam_mean_Palm_interval': str(third[3]),
            'single_seam_reducible_negative_control': {
                'never_moving_configurations': inactive,
                'stationary_mean_move_rate': str(single_rate),
                'Palm_mean_interval_attempts': str(single_interval),
                'reciprocal_stationary_rate': str(1/single_rate),
                'reciprocal_rate_identity_fails_without_eventual_events_from_every_state': True},
            'global_lazification_preserves_all_stopped_configuration_transitions': first[2] == second[2],
            'changed_seam_weights_preserve_stopped_transition_law': first[2] == third[2],
            'no_current_records_and_attempt_cost_preserved_by_lazification': False}


def runtime_checks(kernel, reduction):
    """One actual five-attempt CQ history; no-current records are retained."""
    visible, _ = transfer.instruments(kernel)
    events = tuple(runtime.Event(f'native_current_{n+1}', instruments=(visible,)) for n in range(5))
    program = runtime.Program('native_event_counter_five_attempts', 66, np.eye(66)/66, events,
        binding={'parameter': 'native attempt count',
                 'initial_law': 'uniform stationary native mixture',
                 'readout': 'absolute native current, with every zero retained',
                 'Palm_law': 'conditional on first actual retained current bit being one'})
    histories = [runtime.execute(program, stop_after=0)]
    for n in range(1, 6):
        histories.append(runtime.execute(program, stop_after=n, checkpoint=histories[-1]))
    native_matrices = tuple(tuple(tuple(F(value, 60) for value in row) for row in matrix)
                            for matrix in kernel.marked_counts)
    errors = {'probability': 0., 'full_terminal_density': 0., 'reduced_class_law': 0.}
    word_rows = []
    for branch in histories[-1].branches:
        population = (F(1, 66),)*66
        for bit in branch.record:
            population = _row_times(population, native_matrices[bit])
        probability = _probability(population)
        reduced = word_law(reduction, branch.record)
        projected = tuple(sum((population[i] for i, c in enumerate(reduction['classes']) if c == j), F(0))
                          for j in range(3))
        if projected != reduced:
            raise ArithmeticError('the reduced process lost a full native bit word')
        errors['probability'] = max(errors['probability'], abs(branch.probability-float(probability)))
        errors['full_terminal_density'] = max(errors['full_terminal_density'],
            float(np.linalg.norm(branch.density-np.diag(list(map(float, population))), 'fro')))
        errors['reduced_class_law'] = max(errors['reduced_class_law'],
            abs(float(_probability(reduced))-branch.probability))
        if dict(branch.resources) != {'native_attempts': 5, 'retained_current_bits': 5}:
            raise ArithmeticError('complete native clock resources changed')
        word_rows.append({'record': list(branch.record), 'exact_probability': str(probability),
                          'runtime_probability': branch.probability})
    stopping = []
    for start, prefix, normalizer in (('stationary', (), F(1)), ('palm', (1,), reduction['move_mass'])):
        expected, tail = waiting_law(reduction, 4, start=start)
        observed = []
        stopped_density_error = 0.
        for n in range(1, 5):
            record = prefix+(0,)*(n-1)+(1,)
            history = histories[len(record)]
            branch = next(branch for branch in history.branches if branch.record == record)
            observed.append(branch.probability/float(normalizer))
            population = (F(1, 66),)*66
            for bit in record:
                population = _row_times(population, native_matrices[bit])
            stopped_density_error = max(stopped_density_error,
                float(np.linalg.norm(branch.density/float(normalizer)
                    -np.diag([float(p/normalizer) for p in population]), 'fro')))
        no_record = prefix+(0,)*4
        no_branch = next(branch for branch in histories[len(no_record)].branches if branch.record == no_record)
        observed_tail = no_branch.probability/float(normalizer)
        no_population = (F(1, 66),)*66
        for bit in no_record:
            no_population = _row_times(no_population, native_matrices[bit])
        stopped_density_error = max(stopped_density_error,
            float(np.linalg.norm(no_branch.density/float(normalizer)
                -np.diag([float(p/normalizer) for p in no_population]), 'fro')))
        error = max(abs(a-float(b)) for a, b in zip((*observed, observed_tail), (*expected, tail), strict=True))
        stopping.append({'start': start, 'conditioning_prefix': list(prefix),
                         'conditioning_probability': str(normalizer),
                         'first_move_at_attempt_1_to_4': list(map(str, expected)),
                         'no_move_through_four_attempts': str(tail),
                         'all_five_stopping_outcomes_mass': '1',
                         'runtime_first_move_probabilities': observed,
                         'runtime_no_move_tail': observed_tail,
                         'runtime_residual_diagnostic': _bound(error, 'stopped prefix law'),
                         'complete_stopped_density_residual_diagnostic': _bound(stopped_density_error, 'stopped prefix density'),
                         'physical_halt_device_implemented': False,
                         'stopped_law_obtained_from_complete_executed_prefixes': True})
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(histories[2]))
    resumed = runtime.execute(program, checkpoint=restored)
    checkpoint_error = _bound(runtime.complete_distance(resumed, histories[-1]), 'clock checkpoint replay')
    return {'native_attempts': 5, 'complete_bit_records': word_rows,
            'terminal_quantum_dimension_per_record': 66,
            'complete_terminal_states_checked': True,
            'errors_diagnostic': {name: _bound(value, name) for name, value in errors.items()},
            'stopping_laws': stopping,
            'all_no_current_records_retained': True,
            'checkpoint_replay_cq_distance_diagnostic': checkpoint_error,
            'expanded_tree_resources': dict(histories[-1].expanded_tree_resources)}


def run():
    kernel = transfer.native_counts()
    reduction = lumping(kernel)
    intervals = interval_statistics(reduction)
    count_rows = []
    for n in (1, 2, 8, 32):
        probabilities, mean, variance = count_law(reduction, n)
        p = reduction['move_mass']
        independent = tuple(F(math.comb(n, k))*p**k*(1-p)**(n-k) for k in range(n+1))
        tv = sum((abs(a-b) for a, b in zip(probabilities, independent, strict=True)), F(0))/2
        count_rows.append({'native_attempts': n, 'probabilities_by_complete_count': list(map(str, probabilities)),
                           'mean_current_count': str(mean), 'variance_current_count': str(variance),
                           'independent_tick_variance': str(n*p*(1-p)),
                           'TV_from_independent_equal_rate_ticks': str(tv)})
    waited = []
    for start in ('stationary', 'palm'):
        mass, tail = waiting_law(reduction, 16, start=start)
        waited.append({'start': start, 'attempts': 16, 'first_move_probabilities': list(map(str, mass)),
                       'no_move_tail_probability': str(tail)})
    return {'schema': 'oph.native_current_clock.v1', 'all_pass': True,
            'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'exact_marked_lumping': {
                'original_configurations': 66, 'primitive_branch_addresses_checked': 3960,
                'reduced_classes': ['adjacent occupied ports', 'distance-two occupied ports', 'antipodal occupied ports'],
                'class_sizes': list(reduction['sizes']), 'configuration_classes': list(reduction['classes']),
                'no_current_matrix': [list(map(str, row)) for row in reduction['marked'][0]],
                'moving_current_matrix': [list(map(str, row)) for row in reduction['marked'][1]],
                'stationary_class_law': list(map(str, reduction['stationary'])),
                'Palm_class_law_after_actual_moving_current': list(map(str, reduction['palm'])),
                'all_marked_rows_verified_within_each_class': True,
                'observable_scope': 'all finite absolute-current bit words and their count/stopping functions',
                'claim_of_reducing_every_quantum_observable': False},
            'event_clock': {
                'counter': 'N(n)=sum_{k=1}^n abs(j_k)',
                'parameter': 'native attempt index n, not seconds',
                'increments': 'zero or one; adjacent disjoint attempt blocks add exactly',
                'stationary_move_probability': str(reduction['move_mass']),
                'class_escape_probabilities': list(map(str, intervals['escape'])),
                'stationary_start_wait_mean_attempts': str(intervals['stationary_start_wait_mean']),
                'Palm_mean_interval_attempts': str(intervals['mean_interval']),
                'Palm_interval_variance_attempts_squared': str(intervals['interval_variance']),
                'adjacent_Palm_interval_covariance_attempts_squared': str(intervals['adjacent_interval_covariance']),
                'stopped_moving_transition_matrix': [list(map(str, row)) for row in intervals['jump']],
                'P_next_1_given_last_1': str(intervals['P_next_1_given_last_1']),
                'P_next_1_given_last_two_11': str(intervals['P_next_1_given_last_two_11']),
                'independent_renewal_intervals': False, 'current_bits_form_a_first_order_Markov_chain': False},
            'complete_finite_horizon_waiting_laws': waited,
            'complete_current_count_laws': count_rows,
            'duration_scale_identifiability': duration_scale_audit(reduction),
            'scheduler_controls': scheduler_controls(kernel),
            'actual_observer_execution': runtime_checks(kernel, reduction),
            'theorem_attachment': {
                'attached_count_theorems': ['OPH.QFT.SourceRecurrenceClock.returnCount_block_additive',
                                            'OPH.QFT.SourceRecurrenceClock.returnCount_mono',
                                            'OPH.QFT.SourceRecurrenceClock.returnCount_le'],
                'count_instantiation': 'i=1 and omega(k)=abs(native current j_k), viewed in Fin 2; no Markov premise',
                'nearest_waiting_theorem': 'OPH.QFT.SourceRecurrenceClock.meanReturn_hasSum',
                'waiting_theorem_domain': 'a distinct committed two-state chain; its numerical return laws are not imported',
                'new_checked_attachment': 'actual binary native kernel, exact marked strong lumping, geometric waiting sums and joint CQ readout laws',
                'count_clock_enclosure_not_attached': 'Time.SourceCountClockEnclosure requires independent count-volume bounds and fourth-power proper-time scaling',
                'new_Lean_proof_claimed': False},
            'scope': {
                'carrier_and_native_rule_supplied': True,
                'uniform_fixed_seam_scheduler_needed_for_three_state_reduction': True,
                'uniform_scheduler_needed_for_stationary_mean_rate': False,
                'positive_seam_weights_sufficient_for_Palm_mean_interval_identity': True,
                'stationary_native_preparation_supplied': True,
                'state_dependent_or_adaptive_scheduler_covered': False,
                'clock_signal_is_a_current_of_the_same_process': True,
                'unobserved_classes_need_not_be_accessed_to_measure_reported_bits': True,
                'no_event_record_and_finite_horizon_tail_retained': True,
                'independent_stochastic_clock_added': False,
                'physical_duration_or_proper_time_identified': False,
                'new_lab_prediction_or_empirical_validation': False,
                'full_quantum_instrument_selected_by_classical_record_laws': False,
                'roundoff_interval_certificate': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if result.get('all_pass') is not True:
        raise ArithmeticError('native current clock checks failed')
    if args.check:
        expected = strict_json.loads(args.check.read_text(encoding='utf-8'))
        compare_receipt(expected, result)
        print('Native current clock receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')
