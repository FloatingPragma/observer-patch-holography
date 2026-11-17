"""Exact joint conversion of native incident-hit and movement clocks.

The retained host capture is replayed, not recaptured. Its total-two sector and
software clock are not identified with a laboratory field apparatus or with
the full-binary-sector refinement protocol.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import sympy as sp

from . import comparison, host_clock
from ..observer import native
from ..primitive_observer import transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*host_clock.SOURCE_PATHS,
    'model/native_clock/clock_bridge.py',
    'model/native_clock/host_clock_records.json',
    'model/native_clock/host_clock_results.json'))))
PORT = host_clock.VALIDATION_PORT
RATIO = F(33, 5)
HELD_OUT_FAILURE_BUDGET = F(1, 20)


def _int(value, low, high, name):
    return comparison._integer(value, low, high, name)


def marked_source(kernel, port=PORT):
    """Partition all 60 actual choices by incident hit and moving current."""
    _int(port, 0, 11, 'port')
    transfer.exact_audit(kernel)
    index = {state: i for i, state in enumerate(kernel.states)}
    matrices = [sp.zeros(66) for _ in range(3)]  # remote, incident zero, incident tick
    nonmoving_incident = [F(0)]*66
    for i, state in enumerate(kernel.states):
        for seam in kernel.seams:
            for coin in (0, 1):
                output, (_, current) = native.native(state, seam, coin)
                incident, tick = int(port in seam), int(bool(current) and port in seam)
                matrices[incident+tick][i, index[output]] += sp.Rational(1, 60)
                if not current:
                    if output != state:
                        raise ArithmeticError('a zero-current source branch changed configuration')
                    nonmoving_incident[i] += F(incident, 60)
    A, B0, B1 = matrices
    ones = sp.ones(66, 1)
    if (A*ones != sp.Rational(5, 6)*ones
            or (B0+B1)*ones != sp.Rational(1, 6)*ones
            or 60*(A+B0+B1) != sp.Matrix(kernel.counts)):
        raise ArithmeticError('the complete incident partition changed the native source')
    # Row convention: any number of remote attempts precede one incident hit.
    inverse = (sp.eye(66)-A).inv()
    Q0, Q1 = inverse*B0, inverse*B1
    if (Q0+Q1)*ones != ones or any(x < 0 for Q in (Q0, Q1) for x in Q):
        raise ArithmeticError('the stopped incident law lost positive probability mass')
    stationary = sp.ones(1, 66)/66
    if stationary*(Q0+Q1) != stationary or (stationary*Q1*ones)[0] != F(5, 33):
        raise ArithmeticError('the source-derived incident-to-current rate changed')
    return {'port': port, 'A': A, 'B0': B0, 'B1': B1, 'Q0': Q0, 'Q1': Q1,
            'nonmoving_incident': tuple(nonmoving_incident)}


def count_field_law(kernel, source, hits, initial=0, *, feedback=False):
    """Joint local current count, field reads and actual incident-hit count.

    Feedback after the first scheduled hit reads that port's occupation. A one
    causes one extra existing native seam repair, with its fair coin and tick
    retained. Remote attempts are summed, not suppressed or timed externally.
    """
    _int(hits, 1, 8, 'scheduled incident hits')
    if type(feedback) is not bool:
        raise ValueError('feedback must be Boolean')
    if initial == 'stationary':
        prior = sp.ones(1, 66)/66
    else:
        _int(initial, 0, 65, 'initial configuration')
        prior = sp.zeros(1, 66)
        prior[initial] = 1
    port = source['port']
    seam = next(edge for edge in kernel.seams if port in edge)
    index = {state: i for i, state in enumerate(kernel.states)}
    # key: (current count, first field bit, extra feedback hits)
    rows = {(0, -1, 0): prior}
    for hit in range(hits):
        future = {}
        for (count, first, extra), row in rows.items():
            for tick, Q in enumerate((source['Q0'], source['Q1'])):
                propagated = row*Q
                for j, mass in enumerate(propagated):
                    if not mass:
                        continue
                    bit = kernel.states[j][port] if hit == 0 else first
                    if hit == 0 and feedback and bit:
                        for coin in (0, 1):
                            output, (_, current) = native.native(kernel.states[j], seam, coin)
                            key = (count+tick+abs(current), bit, 1)
                            future.setdefault(key, sp.zeros(1, 66))[index[output]] += mass/2
                    else:
                        key = (count+tick, bit, extra)
                        future.setdefault(key, sp.zeros(1, 66))[j] += mass
        rows = future
    law = {}
    for (count, first, extra), row in rows.items():
        for j, mass in enumerate(row):
            key = (count, first, kernel.states[j][port], hits+extra)
            law[key] = law.get(key, F(0))+F(mass)
    law = {key: value for key, value in law.items() if value}
    if sum(law.values()) != 1 or any(p < 0 for p in law.values()):
        raise ArithmeticError('the field/count/feedback law discarded probability')
    return law


def movement_generating_kernel(kernel, source, z, y):
    """Exact joint PGF at the next global movement, including every wait."""
    if any(type(x) is not F or not 0 <= x <= 1 for x in (z, y)):
        raise ValueError('exact probability-generating arguments in [0,1] required')
    matrix = sp.zeros(66)
    for i in range(66):
        escape = F(sum(kernel.marked_counts[1][i]), 60)
        q1 = source['nonmoving_incident'][i]
        q0 = 1-escape-q1
        denominator = 1-q0-q1*z
        if denominator <= 0:
            raise ArithmeticError('the full waiting series does not converge')
        for j, count in enumerate(kernel.marked_counts[1][i]):
            if count:
                tick = abs(kernel.states[j][source['port']]-kernel.states[i][source['port']])
                matrix[i, j] = F(count, 60)*(z*y)**tick/denominator
    return matrix


def conversion_certificate(kernel, source):
    """Poisson/martingale certificate at any fixed number of global movements.

    Between moves the configuration is fixed. If p is its escape probability
    and q the probability of an incident zero-current attempt, the number G of
    those incident holds is geometric on {0,1,...}, mean q/p and variance
    q(p+q)/p^2, independently of the next moving destination. The next moving
    seam contributes a local tick k and an incident hit k. Therefore one
    increment of I-(33/5)K is G-(28/5)k, not its mean.
    """
    jumps, palm = comparison.jump_source(kernel)
    J = sp.zeros(66)
    f = sp.zeros(66, 1)
    means, variances = [], []
    for i, row in enumerate(jumps):
        escape = F(sum(kernel.marked_counts[1][i]), 60)
        q = source['nonmoving_incident'][i]
        means.append(q/escape)
        variances.append(q*(escape+q)/escape**2)
        for j, probability in row:
            J[i, j] = probability
            tick = abs(kernel.states[j][source['port']]-kernel.states[i][source['port']])
            f[i] += probability*(means[i]-(RATIO-1)*tick)
    nu = sp.Matrix(1, 66, list(palm))
    if nu*f != sp.zeros(1, 1):
        raise ArithmeticError('the exact incident/current conversion has nonzero Palm drift')
    h = (sp.eye(66)-J+sp.ones(66, 1)*nu).inv()*f
    if (sp.eye(66)-J)*h != f or nu*h != sp.zeros(1, 1):
        raise ArithmeticError('the clock-conversion Poisson identity failed')
    state_variance = []
    for i, row in enumerate(jumps):
        drift, variance = F(0), variances[i]
        for j, probability in row:
            tick = abs(kernel.states[j][source['port']]-kernel.states[i][source['port']])
            corrected = means[i]-(RATIO-1)*tick+F(h[j]-h[i])
            drift += probability*corrected
            variance += probability*corrected**2
        if drift:
            raise ArithmeticError('the corrected joint clock increment is not a martingale')
        state_variance.append(variance)
    return {'correctors': tuple(map(F, h)), 'range': F(max(h)-min(h)),
            'variance_ceiling': max(state_variance),
            'asymptotic_variance_per_global_movement': sum(p*v for p, v in zip(palm, state_variance)),
            'geometric_means': tuple(means), 'geometric_variances': tuple(variances),
            'jump': jumps, 'palm': palm}


def error_radius(certificate, movements, failure_probability):
    """Integer d with Pr(|I-(33/5)K|>d) <= failure_probability.

    Uniform in the initial configuration, including starts at a preceding
    stopping boundary. This bounds counts, not host-clock accuracy or entropy.
    """
    _int(movements, 1, 10**12, 'global movements')
    if type(failure_probability) is not F or not 0 < failure_probability < 1:
        raise ValueError('an exact failure probability strictly between zero and one is required')
    variance_budget = movements*certificate['variance_ceiling']/failure_probability
    root = math.isqrt(variance_budget.numerator//variance_budget.denominator)
    if F(root*root) < variance_budget:
        root += 1
    return math.ceil(certificate['range'])+root


def retained_comparison(kernel, certificate, raw_path=host_clock.RAW_PATH):
    """New incident marks of the existing trace; never re-fit host timing."""
    raw = strict_json.loads(Path(raw_path).read_text(encoding='utf-8'))
    old = host_clock.analyse(raw)  # all currents, seeded choices, capture hashes
    compare_receipt(strict_json.loads((ROOT/'model/native_clock/host_clock_results.json').read_text()), old)
    per_tick = F(old['calibration']['nanoseconds_per_local_tick'])
    delta = host_clock.TIMING_TOLERANCE
    radius = error_radius(certificate, host_clock.MOVEMENTS, HELD_OUT_FAILURE_BUDGET/3)
    rows = []
    for run, replay in zip(raw['runs'], old['native_record_replays'], strict=True):
        incident = zero = boundary = 0
        for attempt, code in enumerate(bytes.fromhex(run['records_hex']), 1):
            address, current_code = divmod(code, 3)
            seam = kernel.seams[address//2]
            if PORT in seam:
                incident += 1
                zero += int(current_code == 1)
            if attempt != run['block_endpoints'][boundary+1]['attempts']:
                continue
            block = replay['blocks'][boundary]
            ticks = block['local_counts'][PORT]
            if incident != zero+ticks:
                raise ArithmeticError('incident hits do not split into zero and moving records')
            if boundary > 0:  # calibration block never becomes a held-out test
                low = (1-delta)*per_tick*max(F(0), F(incident-radius, 1)/RATIO)
                high = (1+delta)*per_tick*F(incident+radius, 1)/RATIO
                duration = block['host_duration_ns']
                rows.append({'run': run['name'], 'block': boundary,
                    'global_movements': block['movements'], 'incident_hits': incident,
                    'zero_current_incident_hits': zero, 'local_movement_ticks': ticks,
                    'incident_minus_scaled_current': str(incident-RATIO*ticks),
                    'within_native_count_bound': abs(incident-RATIO*ticks) <= radius,
                    'host_duration_ns': duration,
                    'host_duration_interval_ns': [str(low), str(high)],
                    'host_interval_verdict': 'within_tolerance' if low <= duration <= high else 'reject'})
            boundary += 1
            incident = zero = 0
    if len(rows) != 6:
        raise ArithmeticError('the complete paired held-out blocks were not retained')
    return {'capture_and_previous_analysis_replayed': True,
        'new_timestamps_or_fitted_parameters': False,
        'joint_incident_analysis_is_retrospective': True,
        'calibration': old['calibration'],
        'comparison_port': PORT,
        'count_error_radius': radius,
        'ideal_three_held_out_block_union_bound': str(HELD_OUT_FAILURE_BUDGET),
        'identical_control_records_are_independent_samples': False,
        'pseudorandom_source_has_physical_iid_certificate': False,
        'host_relative_tolerance_retained': str(delta),
        'host_tolerance_is_statistical_confidence_or_accuracy_bound': False,
        'held_out_comparisons': rows}


def run():
    kernel = transfer.native_counts()
    source = marked_source(kernel)
    certificate = conversion_certificate(kernel, source)
    one = sp.ones(66, 1)
    stationary = sp.ones(1, 66)/66
    p11 = F((stationary*source['Q1']*source['Q1']*one)[0])
    M = movement_generating_kernel(kernel, source, F(1, 2), F(2, 3))
    pgf_values = []
    row = sp.zeros(1, 66)
    row[0] = 1
    for n in range(1, 5):
        row = row*M
        pgf_values.append({'global_movements': n, 'value': str(F((row*one)[0]))})
    laws = []
    for initial, feedback in ((0, False), (0, True), ('stationary', False)):
        law = count_field_law(kernel, source, 2, initial, feedback=feedback)
        laws.append({'initial': initial, 'feedback': feedback,
            'scheduled_incident_hits': 2,
            'outcomes': [{'local_movement_ticks': k, 'first_field_bit': c,
                'final_field_bit': d, 'actual_incident_hits': i, 'probability': str(p)}
                for (k, c, d, i), p in sorted(law.items())]})
    return {'schema': 'oph.native_joint_clock_bridge.v1',
        'exact_source_checks_passed': True,
        'sources': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'source_domain': {'ports': 12, 'seams': 30, 'binary_total_occupation': 2,
            'configurations': 66, 'clock_port': PORT,
            'scheduler': 'independent uniform native seam and fair placement coin',
            'nonstationary_initial_law_allowed': True,
            'full_binary_Bernoulli_refinement_preparation_used': False},
        'incident_hit_law': {'row_kernel': '(I-A)^(-1) (B0+z B1)',
            'A': 'remote-attempt subkernel', 'B0': 'incident zero-current subkernel',
            'B1': 'incident moving-current subkernel',
            'stationary_ticks_per_incident_hit': '5/33',
            'two_consecutive_incident_ticks_probability': str(p11),
            'independent_tick_prediction': str(F(5, 33)**2),
            'incident_ticks_are_independent_bernoulli': p11 == F(5, 33)**2,
            'joint_count_and_field_laws': laws,
            'feedback_included_in_actual_clock': True},
        'movement_clock_joint_law': {
            'generating_kernel': 'M(z,y)[i,j]=p_ij*(z*y)^k_ij/(1-q0_i-q1_i*z)',
            'parameters': 'p_ij=moving transition probability per native attempt; q0/q1=nonmoving remote/incident probabilities; k_ij=local movement bit',
            'n_movement_joint_count_generating_function': 'initial_row*M(z,y)^n*ones',
            'checked_generating_values': {'initial_configuration': 0,
                'incident_argument_z': '1/2', 'local_tick_argument_y': '2/3',
                'values': pgf_values},
            'all_nonmoving_incident_hits_and_infinite_waiting_tails_retained': True,
            'incident_hits_per_local_tick_asymptotic_ratio': str(RATIO),
            'ratio_is_exact_pathwise_clock_conversion': False},
        'finite_error_certificate': {
            'parameter': 'n completed global movements',
            'difference': 'Z_n=I_v(n)-(33/5)*K_v(n)',
            'martingale': 'Z_n+h(X_n)-h(X_0)',
            'corrector_range': str(certificate['range']),
            'conditional_variance_ceiling': str(certificate['variance_ceiling']),
            'asymptotic_variance_per_global_movement': str(certificate['asymptotic_variance_per_global_movement']),
            'tail_bound': 'Pr(|Z_n|>=d) <= min(1,n*V/(d-R)^2) for d>R',
            'stationary_preparation_required_for_bound': False,
            'same_bound_conditional_on_preceding_stopping_boundary': True,
            'feedback_or_changed_scheduler_covered_by_this_bound': False},
        'retained_host_attachment': retained_comparison(kernel, certificate),
        'scope': {'new_independent_clock_dynamics_added': False,
            'source_joint_count_distribution_derived': True,
            'full_binary_refinement_protocol_or_H_lab_calibrated': False,
            'native_real_time_unitary_or_modular_phase_derived': False,
            'traceable_lab_clock_or_physical_proper_time_established': False,
            'empirical_validation_of_OPH_physics': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
        print('Joint native clock laws and retained host comparison reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
