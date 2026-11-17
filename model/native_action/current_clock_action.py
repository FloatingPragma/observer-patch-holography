"""Native current action in one port's own movement clock.

The observer reads current records across one carrier and uses movements of
one selected port as its clock. No seam/coin rates or physical units are fitted.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp

from tools import strict_json
from ..event_circulation import construction as event
from ..native_circulation import construction as circulation
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical, compare_receipt
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*event.SOURCE_PATHS, *quantum.SOURCE_PATHS,
    'model/native_action/current_clock_action.py',
))))
MARKS = ((0, 0, 0), (1, 0, -1), (1, 0, 0), (1, 0, 1),
         (1, 1, -1), (1, 1, 0), (1, 1, 1))


def _port(port):
    if type(port) is not int or not 0 <= port < 12:
        raise ValueError('one exact source port from zero through eleven required')
    return port


@lru_cache(maxsize=12)
def stopped_source(port=0):
    """Exact event kernels with derivatives marking the triangular current.

    q[b][r] contains g**r for r=0,1,2 and b indicating a local clock tick.
    Every movement has g in {-1,0,1}; higher derivatives reuse odd/even rows.
    """
    _port(port)
    kernel, triangle, weights, branches = circulation.source()
    degree, palm, _, _ = event.embedded()
    q = [[sp.zeros(66) for _ in range(3)] for _ in range(2)]
    local_degree = [0]*66
    for i, j, _, _, current, score in branches:
        if not current:
            continue
        tick = int(kernel.states[i][port] != kernel.states[j][port])
        local_degree[i] += tick
        for order in range(3):
            q[tick][order][i, j] += sp.Rational(score**order, degree[i])
    q = tuple(tuple(sp.ImmutableMatrix(matrix) for matrix in row) for row in q)
    one = sp.ones(66, 1)
    mu = sp.Matrix(1, 66, [sp.Rational(x) for x in palm])
    nu = sp.Matrix(1, 66, [sp.Rational(d, 100) for d in local_degree])
    if sum(local_degree) != 100 or nu != 6*mu*q[1][0]:
        raise ArithmeticError('the native local-tick Palm law changed')
    if q[0][0]*one+q[1][0]*one != one:
        raise ArithmeticError('movement/tick decomposition lost probability')
    beta = max(q[0][0]**3*one)
    if not 0 <= beta < 1:
        raise ArithmeticError('three-movement no-tick tail is not uniformly contractive')
    inverse = (sp.eye(66)-q[0][0]).inv()
    if nu*inverse != 6*mu:
        raise ArithmeticError('stopped Palm occupation measure differs from six native movements')
    transition = inverse*q[1][0]
    if transition*one != one or nu*transition != nu:
        raise ArithmeticError('complete tick-to-tick transition lost its stationary law')
    support = tuple(i for i, value in enumerate(nu) if value)
    reached = {support[0]}
    while True:
        following = reached | {j for i in reached for j in support if transition[i, j]}
        if following == reached:
            break
        reached = following
    if reached != set(support):
        raise ArithmeticError('tick-to-tick chain is not irreducible on the Palm support')
    occupancy = sp.Matrix([2*state[port]-1 for state in kernel.states])
    if transition*occupancy != -occupancy:
        raise ArithmeticError('a local tick must flip the native binary occupation')
    return {'port': port, 'kernel': kernel, 'triangle': triangle, 'q': q,
            'inverse': inverse, 'transition': transition, 'mu': mu, 'nu': nu,
            'local_degree': tuple(local_degree), 'tail_three_movements': F(beta),
            'Palm_support': support}


def stopped_moments(source, maximum=4, boundary=None):
    """E_i[S_T**r b(X_T)] by exact reward recursions, including unbounded waits."""
    if type(maximum) is not int or not 0 <= maximum <= 4:
        raise ValueError('an exact moment order from zero through four required')
    boundary = sp.ones(66, 1) if boundary is None else sp.Matrix(boundary)
    if boundary.shape != (66, 1):
        raise ValueError('one boundary value for every native configuration required')
    q, inverse = source['q'], source['inverse']
    derivative = lambda order: 0 if order == 0 else 1 if order % 2 else 2
    moments = [inverse*q[1][0]*boundary]
    for order in range(1, maximum+1):
        rhs = q[1][derivative(order)]*boundary
        for previous in range(order):
            rhs += math.comb(order, previous)*q[0][derivative(order-previous)]*moments[previous]
        moments.append(inverse*rhs)
    return tuple(moments)


def local_tick_law(port=0):
    source = stopped_source(port)
    moments = stopped_moments(source)
    nu = source['nu']
    values = tuple(F((nu*moment)[0]) for moment in moments)
    if moments[1] != sp.zeros(66, 1) or values[2] != F(3, 5):
        raise ArithmeticError('stopped cycle martingale or its exact variance changed')
    joint = F((nu*stopped_moments(source, 2, moments[2])[2])[0])
    fourth_connected = values[4]-3*values[2]**2
    squared_dependence = joint-values[2]**2
    if not fourth_connected or not squared_dependence:
        raise ArithmeticError('finite local-clock records no longer distinguish independent Gaussian blocks')
    previous = circulation.exact_geometry()
    projection = sp.Matrix(previous['cycle_projection'])
    # The expected pre-event occupation before the next local tick is six
    # times the original movement Palm law. Each edge has second moment1/30.
    stopped_edge_second = [sp.Rational(0)]*30
    occupation = nu*source['inverse']
    degree = event.embedded()[0]
    for i, _, edge, _, current, _ in circulation.source()[3]:
        if current:
            stopped_edge_second[edge] += occupation[i]*sp.Rational(current*current, degree[i])
    if set(stopped_edge_second) != {sp.Rational(1, 5)}:
        raise ArithmeticError('native local-tick current covariance changed')
    covariance = projection*sp.diag(*stopped_edge_second)*projection
    if covariance != projection/5:
        raise ArithmeticError('projected local-tick covariance is not P/5')
    return {
        'clock_port': port, 'triangle': list(source['triangle']),
        'clock_parameter': 'm recorded movements at the selected port; all carrier currents are retained',
        'local_tick_Palm_configuration_probabilities': [str(F(x)) for x in nu],
        'Palm_support_size': len(source['Palm_support']),
        'one_attempt_herald_success': '5/198', 'one_attempt_herald_failure': '193/198',
        'faithful_full_binary_sector_reference': False,
        'uniform_three_movement_no_local_tick_bound': str(source['tail_three_movements']),
        'waiting_tail_bound': 'Pr(T>3k) <= beta**k for every initial native configuration',
        'exact_Palm_mean_global_movements_per_local_tick': str(F((occupation*sp.ones(66, 1))[0])),
        'exact_cycle_covariance_per_local_tick': 'P/5',
        'exact_triangle_variance_after_m_local_ticks': '3*m/5',
        'one_local_tick_triangle_moments_0_through_4': list(map(str, values)),
        'one_local_tick_triangle_fourth_cumulant': str(fourth_connected),
        'adjacent_local_tick_squared_current_product_mean': str(joint),
        'adjacent_local_tick_squared_current_covariance': str(squared_dependence),
        'finite_block_increments_independent_or_Gaussian': False,
        'Gaussian_limit_current_action_coefficient': '5/2',
        'Gaussian_limit_scalar_triangle_action_coefficient': '5/6',
        'local_clock_occupation_mode_multiplier': '-1',
        'local_clock_occupation_alternation_is_spatial_wave': False,
        'proof_scope': 'exact stopped moments for local-Palm preparation; one-history martingale functional limit with bounded-moment tails',
    }


def joint_instrument(law='C', *, herald=False, feedback_seam=False):
    """The unchanged native Kraus maps, refined by movement/port-zero tick/current.

    feedback_seam selects the actual source seam (0,1) with the same fair coin;
    this is an explicit control, not the uniform schedule used by the theorem.
    """
    if law not in ('C', 'E') or type(herald) is not bool or type(feedback_seam) is not bool:
        raise ValueError('C/E and explicit Boolean protocol choices required')
    kernel, _, _, branches = circulation.source()
    selected = kernel.seams.index((0, 1))
    groups, marks = {}, {}
    for i, j, edge, coin, current, score in branches:
        if feedback_seam and edge != selected:
            continue
        a, b = kernel.seams[edge]
        label = current if law == 'C' else kernel.states[i][a]-kernel.states[i][b]
        key = edge, coin, label
        mark = (int(current != 0), abs(kernel.states[j][0]-kernel.states[i][0]), score)
        if key in marks and marks[key] != mark:
            raise ArithmeticError('joint current/clock readout splits a hidden native Kraus label')
        marks[key] = mark
        groups.setdefault(key, np.zeros((66, 66), complex))[j, i] = 1/math.sqrt(2 if feedback_seam else 60)
    count = 2 if herald else len(MARKS)
    outcome = lambda key: marks[key][1] if herald else MARKS.index(marks[key])
    maps = [tuple(groups[key] for key in sorted(groups) if outcome(key) == k) for k in range(count)]
    maps = tuple(group or (np.zeros((66, 66), complex),) for group in maps)
    return runtime.Instrument(kraus=maps, visible=True,
        resources=(('native_attempts', 1), ('joint_current_clock_records', 1)),
        binding={'candidate': law, 'port_clock': 0, 'joint_marks': MARKS,
            'source_seam': [0, 1] if feedback_seam else 'uniform thirty seams',
            'herald_port_zero_tick': herald})


def exact_joint_columns(*, feedback_seam=False):
    """Independent floor/ceiling branch counts; no source endpoint call."""
    kernel = circulation.source()[0]
    index = {state: i for i, state in enumerate(kernel.states)}
    weights = circulation.source()[2]
    counts = [[[0]*66 for _ in range(66)] for _ in MARKS]
    for i, state in enumerate(kernel.states):
        for edge, (a, b) in enumerate(kernel.seams):
            if feedback_seam and (a, b) != (0, 1):
                continue
            total = state[a]+state[b]
            low, high = total//2, (total+1)//2
            for pair in ((low, high), (high, low)):
                output = list(state)
                output[a], output[b] = pair
                current = state[a]-pair[0]
                mark = (int(current != 0), abs(output[0]-state[0]), weights[edge]*current)
                counts[MARKS.index(mark)][i][index[tuple(output)]] += 1
    return tuple(tuple(tuple(row) for row in matrix) for matrix in counts)


def _row_times(row, matrix, denominator):
    return tuple(sum((row[i]*matrix[i][j] for i in range(66) if row[i]), F(0))/denominator
                 for j in range(66))


def finite_execution():
    """Actual herald and a fixed two-attempt window, retaining failure/no ticks.

    A separate record-controlled seam experiment supplies a finite negative
    control: the uniform-schedule action cannot be reused after that feedback.
    This window does not halt at its first local tick and does not execute the
    infinite stopped quantum channel used only for the diagonal moment law.
    """
    source = stopped_source()
    counts = exact_joint_columns()
    biased_counts = exact_joint_columns(feedback_seam=True)
    palm = tuple(F(v) for v in source['nu'])
    successful_initial = tuple(F(d, 3960) for d in source['local_degree'])
    failed_initial = tuple(F(60-d, 3960) for d in source['local_degree'])
    rows = []
    for law in ('C', 'E'):
        read, herald = joint_instrument(law), joint_instrument(law, herald=True)
        original = quantum.native_instrument(law)
        if sorted(k.tobytes() for group in read.kraus for k in group) != sorted(
                k.tobytes() for group in original.kraus for k in group):
            raise ArithmeticError('joint readout changed the original unread channel')
        inactive = runtime.Instrument(kraus=((np.eye(66),),)+tuple(
            (np.zeros((66, 66)),) for _ in MARKS[1:]), visible=True)
        program = runtime.Program('heralded_local_clock_current', 66, np.eye(66)/66, (
            runtime.Event('local_tick_herald', instruments=(herald,)),
            runtime.Event('first_native_attempt', control=0, instruments=(inactive, read)),
            runtime.Event('second_native_attempt', control=0, instruments=(inactive, read)),
        ))
        history = runtime.execute(program)
        probability_rows, timeout = [], F(0)
        for branch in history.branches:
            flag, first, second = branch.record
            if flag:
                row = _row_times(_row_times(successful_initial, counts[first], 60), counts[second], 60)
                if MARKS[first][1] == MARKS[second][1] == 0:
                    timeout += sum(row)
            else:
                row = failed_initial if first == second == 0 else (F(0),)*66
            if np.linalg.norm(branch.density-np.diag(list(map(float, row))), 'fro') > 2e-12:
                raise ArithmeticError('complete local-clock native output disagrees with its exact law')
            probability_rows.append([list(branch.record), str(sum(row))])
        if sum(F(row[1]) for row in probability_rows) != 1 or not 0 < timeout < F(5, 198):
            raise ArithmeticError('herald failure or stopped timeout was discarded')
        checkpoint = runtime.execute(program, stop_after=2)
        resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
            program, runtime.checkpoint_to_json(checkpoint)))
        if runtime.complete_distance(history, resumed) > 2e-12:
            raise ArithmeticError('recorded local clock continuation changed')
        biased = joint_instrument(law, feedback_seam=True)
        controlled = runtime.Program('current_selected_native_schedule', 66, np.diag(list(map(float, palm))), (
            runtime.Event('first_native_read', instruments=(read,)),
            runtime.Event('record_controlled_seam', control=0,
                instruments=tuple(biased if mark[2] > 0 else read for mark in MARKS)),
        ))
        actual = runtime.execute(controlled)
        distance = F(0)
        tilted = F(0)
        for branch in actual.branches:
            first, second = branch.record
            prefix = _row_times(palm, counts[first], 60)
            chosen = biased_counts if MARKS[first][2] > 0 else counts
            denominator = 2 if MARKS[first][2] > 0 else 60
            row = _row_times(prefix, chosen[second], denominator)
            unchanged = _row_times(prefix, counts[second], 60)
            if np.linalg.norm(branch.density-np.diag(list(map(float, row))), 'fro') > 2e-12:
                raise ArithmeticError('feedback history changed the declared source instrument')
            mass = sum(row)
            distance += abs(mass-sum(unchanged))/2
            # A fitted-free rational evaluation of the complete record
            # generating functional, not a selected or normalized subensemble.
            score = MARKS[first][2]+MARKS[second][2]
            ticks = MARKS[first][1]+MARKS[second][1]
            tilted += mass*F(2)**score*F(3)**ticks
        if not distance:
            raise ArithmeticError('feedback stopped distinguishing the uniform-schedule action')
        rows.append({'candidate': law, 'herald_and_two_attempt_complete_records': len(history.branches),
            'readout_window': 'two native attempts after successful herald, including attempts after any first local tick',
            'preparation_failure_mass': '193/198',
            'conditional_no_local_tick_timeout_after_two_native_attempts': str(timeout/F(5, 198)),
            'complete_native_record_probabilities': probability_rows,
            'feedback_two_attempt_complete_records': len(actual.branches),
            'feedback_vs_unchanged_schedule_record_TV': str(distance),
            'feedback_generating_function_E_2_power_current_3_power_ticks': str(tilted),
            'complete_terminal_densities_and_serialized_continuation_checked': True,
            'quantum_Kraus_family_unchanged_before_explicit_schedule_control': True})
    return rows


def run():
    return {'schema': 'oph.native_action.current_clock_action.v1',
        'source_derived_local_tick_action': local_tick_law(),
        'finite_current_clock_execution': finite_execution(),
        'scope': {'carrier_count': 1, 'binary_total_two_states': 66,
            'uniform_independent_native_seam_and_coin_schedule': True,
            'current_access': 'carrier-wide seam/current records; the clock uses only port-zero movements',
            'extra_independent_clock_or_replicas': False,
            'feedback_preserves_stationary_action_without_rechecking': False,
            'full_quantum_stopped_channel_identified_with_diagonal_jump_chain': False,
            'laboratory_seconds_or_relativistic_wave_action_derived': False},
        'source_bindings': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native current action in a recorded local clock reproduced.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
