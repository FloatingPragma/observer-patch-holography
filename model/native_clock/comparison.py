"""Exact comparison and spatial covariance of twelve clocks on one carrier.

Only actual native currents tick the clocks. The comparison parameter n counts
completed carrier movements, not attempts, modular phase, or laboratory seconds.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp

from ..observer import native
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime, transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*transfer.SOURCE_PATHS, *quantum.SOURCE_PATHS,
    'model/native_clock/comparison.py',
    'Lean/Screen/RegionalContinuity.lean',
    'Lean/Time/ClockReadout.lean',
))))
MARKS = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 0, 1), (1, 1, 1))


def _integer(value, lower, upper, name):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be an integer in [{lower}, {upper}]')
    return value


def jump_source(kernel):
    """Eliminate genuine no-movement holds; no new stochastic law is fitted."""
    transfer.exact_audit(kernel)
    degree = tuple(sum(row) for row in kernel.marked_counts[1])
    transitions = tuple(tuple((j, F(count, degree[i])) for j, count in enumerate(row) if count)
                        for i, row in enumerate(kernel.marked_counts[1]))
    palm = tuple(F(d, sum(degree)) for d in degree)
    for i, row in enumerate(transitions):
        if sum(p for _, p in row) != 1:
            raise ArithmeticError('the movement transition lost mass')
        for j, p in row:
            changed = tuple(v for v in range(12) if kernel.states[i][v] != kernel.states[j][v])
            if len(changed) != 2 or tuple(sorted(changed)) not in kernel.seams:
                raise ArithmeticError('a movement does not tick its two seam endpoints')
            reverse = dict(transitions[j]).get(i, F(0))
            if palm[i]*p != palm[j]*reverse:
                raise ArithmeticError('movement Palm detailed balance failed')
    return transitions, palm


def clock_correctors(kernel):
    """Solve (I-J)h=f with Palm mean zero, then verify every exact equation."""
    transitions, palm = jump_source(kernel)
    size = len(kernel.states)
    J = sp.zeros(size)
    f = sp.zeros(size, 12)
    for i, row in enumerate(transitions):
        for j, p in row:
            J[i, j] = sp.Rational(p.numerator, p.denominator)
        for v in range(12):
            f[i, v] = sum(J[i, j]*abs(kernel.states[i][v]-kernel.states[j][v])
                          for j, _ in row)-sp.Rational(1, 6)
    nu = sp.Matrix(1, size, [sp.Rational(p.numerator, p.denominator) for p in palm])
    if nu*f != sp.zeros(1, 12):
        raise ArithmeticError('local clocks do not have the derived mean rate 1/6')
    H = (sp.eye(size)-J+sp.ones(size, 1)*nu).inv()*f
    if (sp.eye(size)-J)*H != f or nu*H != sp.zeros(1, 12):
        raise ArithmeticError('the exact clock Poisson equations failed')
    h = tuple(tuple(F(H[i, v]) for v in range(12)) for i in range(size))
    covariance = [[F(0)]*12 for _ in range(12)]
    state_variances = [[F(0)]*12 for _ in range(size)]
    increment_bound = F(0)
    for i, row in enumerate(transitions):
        drift = [F(0)]*12
        for j, p in row:
            delta = tuple(F(abs(kernel.states[i][v]-kernel.states[j][v]))-F(1, 6)
                          + h[j][v]-h[i][v] for v in range(12))
            if sum(delta) != 0:
                raise ArithmeticError('the exact sum of twelve clock increments changed')
            for v in range(12):
                drift[v] += p*delta[v]
                state_variances[i][v] += p*delta[v]**2
                for w in range(12):
                    covariance[v][w] += palm[i]*p*delta[v]*delta[w]
            increment_bound = max(increment_bound, *map(abs, delta))
        if any(drift):
            raise ArithmeticError('the corrected native clock is not a martingale')
    ranges = tuple(max(row[v] for row in h)-min(row[v] for row in h) for v in range(12))
    if len(set(ranges)) != 1 or any(sum(row) for row in covariance):
        raise ArithmeticError('carrier symmetry or the total-clock identity failed')
    return {'transitions': transitions, 'palm': palm, 'correctors': h,
            'covariance': tuple(map(tuple, covariance)), 'range': ranges[0],
            'variance_ceiling': max(map(max, state_variances)), 'increment_bound': increment_bound}


def confidence_bound(certificate, movements, error):
    """Uniform Pr(|K_v/n-1/6| >= error), from the bounded martingale correction."""
    n = _integer(movements, 1, 10**12, 'completed movement count')
    if type(error) is not F or error <= 0:
        raise ValueError('a positive exact Fraction error is required')
    margin = n*error-certificate['range']
    return F(1) if margin <= 0 else min(F(1), n*certificate['variance_ceiling']/margin**2)


def distance_shells(kernel, covariance):
    """Check all unordered port pairs; covariance classes recover source hops."""
    adjacency = {v: set() for v in range(12)}
    for a, b in kernel.seams:
        adjacency[a].add(b)
        adjacency[b].add(a)
    distances = {}
    for start in range(12):
        reached, queue = {start: 0}, [start]
        for port in queue:
            for neighbour in sorted(adjacency[port]):
                if neighbour not in reached:
                    reached[neighbour] = reached[port]+1
                    queue.append(neighbour)
        if len(reached) != 12 or max(reached.values()) != 3:
            raise ArithmeticError('the native port graph is not connected with diameter three')
        distances[start] = reached
    rows = []
    for distance in (1, 2, 3):
        pairs = [(v, w) for v in range(12) for w in range(v+1, 12)
                 if distances[v][w] == distance]
        values = {covariance[v][w] for v, w in pairs}
        if len(values) != 1:
            raise ArithmeticError('clock covariance is not constant on a distance shell')
        value = values.pop()
        differences = {covariance[v][v]+covariance[w][w]-2*value for v, w in pairs}
        if len(differences) != 1:
            raise ArithmeticError('clock difference variance changed inside a shell')
        rows.append({'native_hop_distance': distance, 'unordered_pairs': [list(p) for p in pairs],
                     'pair_count': len(pairs), 'asymptotic_covariance_per_movement': str(value),
                     'asymptotic_difference_variance_per_movement': str(differences.pop())})
    if len({row['asymptotic_covariance_per_movement'] for row in rows}) != 3:
        raise ArithmeticError('clock covariance cannot distinguish the three distance shells')
    return rows


def instrument(kernel, law='C', ports=(2, 3)):
    """Refine existing native records into global and two local movement bits."""
    if law not in ('C', 'E') or type(ports) is not tuple or len(ports) != 2:
        raise ValueError('native C/E law and two distinct source ports required')
    for v in ports:
        _integer(v, 0, 11, 'clock port')
    if ports[0] == ports[1]:
        raise ValueError('the two clock ports must be distinct')
    index = {state: i for i, state in enumerate(kernel.states)}
    groups, outcomes = {}, {}
    for edge, seam in enumerate(kernel.seams):
        for coin in (0, 1):
            for i, state in enumerate(kernel.states):
                output, record = native.native(state, seam, coin)
                current = record[1]
                label = current if law == 'C' else state[seam[0]]-state[seam[1]]
                key = (edge, coin, label)
                mark = (abs(current), *(abs(output[v]-state[v]) for v in ports))
                if key in outcomes and outcomes[key] != mark:
                    raise ArithmeticError('one native Kraus branch has conflicting clock marks')
                outcomes[key] = mark
                groups.setdefault(key, {})[i] = index[output]
    matrices = [[] for _ in MARKS]
    coverage = [0]*66
    for key, mapping in sorted(groups.items()):
        if len(set(mapping.values())) != len(mapping):
            raise ArithmeticError('the native branch is not a partial injection')
        matrix = np.zeros((66, 66), complex)
        for i, j in mapping.items():
            matrix[j, i] = 1/math.sqrt(60)
            coverage[i] += 1
        matrices[MARKS.index(outcomes[key])].append(matrix)
    if set(coverage) != {60}:
        raise ArithmeticError('the complete clock instrument lost a native branch')
    # Some port pairs cannot tick together; retain their zero outcome explicitly.
    for group in matrices:
        if not group:
            group.append(np.zeros((66, 66), complex))
    return runtime.Instrument(kraus=tuple(tuple(g) for g in matrices), visible=True,
        resources=(('native_attempts', 1), ('joint_clock_reads', 1)),
        binding={'native_law': law, 'clock_ports': list(ports), 'marks': [list(m) for m in MARKS],
                 'scheduler': 'uniform thirty seams and fair independent placement coin',
                 'no_event_outcome_retained': True})


def exact_attempt_records(kernel, horizon, ports=(2, 3), initial=0):
    """Independent rational source-branch propagation, retaining all zero words."""
    n = _integer(horizon, 0, 4, 'attempt horizon')
    _integer(initial, 0, 65, 'initial configuration')
    index = {state: i for i, state in enumerate(kernel.states)}
    populations = {(): tuple(F(int(i == initial)) for i in range(66))}
    for _ in range(n):
        future = {}
        for word, row in populations.items():
            outputs = [[F(0)]*66 for _ in MARKS]
            for i, mass in enumerate(row):
                if not mass:
                    continue
                state = kernel.states[i]
                for seam in kernel.seams:
                    for coin in (0, 1):
                        output, record = native.native(state, seam, coin)
                        mark = (abs(record[1]), *(abs(output[v]-state[v]) for v in ports))
                        outputs[MARKS.index(mark)][index[output]] += mass/60
            for mark, out in enumerate(outputs):
                future[word+(mark,)] = tuple(out)
        populations = future
    return populations


def runtime_check(kernel):
    exact = exact_attempt_records(kernel, 3)
    rows = []
    for law in ('C', 'E'):
        initial = np.zeros((66, 66), complex)
        initial[0, 0] = 1
        read = instrument(kernel, law)
        original = quantum.native_instrument(law)
        if sorted(k.tobytes() for group in read.kraus for k in group) != sorted(
                k.tobytes() for group in original.kraus for k in group):
            raise ArithmeticError('the clock readout changed the source Kraus family')
        program = runtime.Program('single_carrier_clock_comparison', 66, initial,
            tuple(runtime.Event(str(i), instruments=(read,)) for i in range(3)))
        history = runtime.execute(program)
        residual = max(float(np.linalg.norm(branch.density-np.diag(list(map(float, exact[branch.record]))), 'fro'))
                       for branch in history.branches)
        if not math.isfinite(residual) or residual > 2e-12 or len(history.branches) != 125:
            raise ArithmeticError('the complete native clock execution disagrees with its oracle')
        if sum(sum(row) for row in exact.values()) != 1:
            raise ArithmeticError('a finite clock history was discarded')
        rows.append({'native_law': law, 'native_attempts': 3, 'complete_records': 125,
                     'nonstationary_initial_configuration': 0,
                     'forgotten_record_kraus_family_matches_native_source': True,
                     'maximum_terminal_density_residual_diagnostic': residual,
                     'zero_movement_record_probability': str(sum(exact[(0, 0, 0)]))})
    return rows


def run():
    kernel = transfer.native_counts()
    certificate = clock_correctors(kernel)
    rows = distance_shells(kernel, certificate['covariance'])
    if (certificate['range'], certificate['variance_ceiling'], certificate['increment_bound']) != (
            F(38, 25), F(30341677, 47045000), F(241, 150)):
        raise ArithmeticError('the exact clock correction constants changed')
    return {
        'schema': 'oph.single_carrier_clock_comparison.v1', 'all_pass': True,
        'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source': {'carriers': 1, 'ports': 12, 'seams': 30, 'binary_total_two_states': 66,
                   'movement_degree_sum': 600, 'uniform_independent_seam_coin_scheduler': True,
                   'stationary_initial_state_required_for_rate_and_limit': False,
                   'quantum_clock_record_law_depends_only_on_initial_diagonal': True,
                   'interventions_during_comparison': False},
        'clock_comparison': {
            'parameter': 'n completed native movements, including remote port movements',
            'local_count': 'K_v(n), number of the n movements changing occupation at port v',
            'exact_pathwise_identity': 'sum_v K_v(n)=2*n',
            'asymptotic_local_count_per_movement': '1/6',
            'global_identity_holds_leave_complete_event_count_law_unchanged': True,
            'corrector_equation': '(I-J)h_v=f_v, f_v(i)=E_i[tick_v]-1/6; Palm mean(h_v)=0',
            'correctors_by_state_and_port': [list(map(str, row)) for row in certificate['correctors']],
            'corrector_range': str(certificate['range']),
            'martingale_increment_absolute_bound': str(certificate['increment_bound']),
            'conditional_variance_ceiling': str(certificate['variance_ceiling']),
            'uniform_tail_bound': 'P(abs(K_v(n)/n-1/6)>=eps) <= min(1,n*V/(n*eps-R)^2) when n*eps>R',
            'examples_absolute_rate_error_1_over_60': [
                {'movements': n, 'error': '1/60', 'failure_probability_upper_bound':
                 str(confidence_bound(certificate, n, F(1, 60)))} for n in (10000, 100000, 1000000)]},
        'clock_spatial_covariance': {
            'definition': 'Sigma_vw=lim_(n->infinity) Cov(K_v(n),K_w(n))/n',
            'asymptotic_for_every_initial_native_density': True,
            'derivation': 'bounded Poisson corrector plus finite irreducible chain Cesaro occupation law',
            'matrix': [list(map(str, row)) for row in certificate['covariance']],
            'single_port_variance_per_movement': str(certificate['covariance'][0][0]),
            'distance_shells': rows,
            'all_66_unordered_pairs_recover_native_distance_shell': True,
            'finite_sample_inference_or_lab_distance_calibration_established': False},
        'actual_native_execution': runtime_check(kernel),
        'lab_attachment': {
            'units_are_observer_conventions': True,
            'existing_lab_apparatus_or_dataset_identified': False,
            'minimum_raw_records': 'all global movements and both local port bits at shared retained endpoints; independent lab timestamps for those endpoints',
            'calibration': 'fix one positive lab interval per observed calibration-block count of one clock; retain its uncertainty and do not fit the other clock',
            'held_out_tests': 'local clock count ratios, full event marks, distance-class covariance and stationarity of paired lab/count rate over disjoint validation blocks',
            'required_operational_match': 'same 12-port total-two source, scheduler, native current detector, no unmodeled interventions, complete records or quantified missing-event/censoring errors',
            'finite_bound_clock': 'global movement count n is the declared comparison parameter, not an assumed constant lab-time rate',
            'unproved_bridge': 'a source-identified lab apparatus must realize these currents and support its lab/count-rate validation; one calibration interval alone is insufficient'},
        'scope': {'independent_clock_dynamics_added': False, 'federation_required': False,
                  'native_modular_phase_advance_established': False,
                  'physical_spacetime_or_empirical_agreement_established': False,
                  'proof_type': 'exact rational finite source identities and analytic bounded-martingale argument; no new Lean theorem'}
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
        print('Single-carrier clock comparison receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
