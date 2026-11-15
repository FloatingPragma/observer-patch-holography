"""Compile a common modular action from inferred native charge weights.

The inputs are three positive exact weights inferred from complete records,
and a simultaneous confidence box at the same (possibly subnormalized)
scale. No analytical reference weights enter compilation or its certificate.
Dividing by the first weight fixes only the irrelevant scalar phase.

If the true weights lie in the box, log monotonicity bounds the difference
of the two diagonal modular generators modulo identity. Duhamel gives a
reference-uniform channel bound. Two complete cap instruments compose that
bound, retaining feedback, disturbance, all source addresses and export
aborts. The statistical confidence and preparation errors are separate.
Coherent source controls remain declared architecture, not passive repair.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np
from tools import strict_json

from . import cap_local_clock as cap, common_caps as common
from . import local_history as history, local_preparation as local

full = common.full
WORD_GATE = F(1, 10**9)
ACCEPTANCE_CAP = F(9, 100)
TRAINING_TRIALS = 1 << 24
OUTPUT = Path(__file__).with_name('learned_clock_results.json')


def _weights(values):
    if (type(values) is not tuple or len(values) != 3
            or any(type(x) is not F or x <= 0 for x in values)):
        raise ValueError('three original positive Fraction charge weights required')
    return values


def ratios(weights):
    """One shared reference on the native shell, identity on all padding."""
    weights = _weights(weights)
    result = [F(1)]*512
    for code, state in zip(full.support()[0], local.native.states(binary=True), strict=True):
        result[code] = weights[sum(state[p] for p in local.CORE)]/weights[0]
    return tuple(result)


def inference_bound(weights, intervals, s=F(1)):
    """Exact log-spread bound, invariant under a common positive rescaling.

    h_q=log(w_hat_q/w_q) lies in [log(w_hat_q/u_q),log(w_hat_q/l_q)].
    Subtract the midpoint of the combined range from H_hat-H; this changes
    only a scalar phase. Its norm is at most half the range. Channel
    half-diamond distance is at most |s| times that radius, capped at one.
    An unresolved or nonpositive box cannot be promoted to a clock.
    """
    weights = _weights(weights)
    s = full.parameter(s)
    if (type(intervals) is not tuple or len(intervals) != 3
            or any(type(box) is not tuple or len(box) != 2
                   or any(type(x) is not F for x in box)
                   or not 0 < box[0] <= w <= box[1]
                   for w, box in zip(weights, intervals, strict=True))):
        raise ValueError('three faithful exact confidence intervals containing their centers required')
    lower = [common.logarithm_interval(w/hi)[0] for w, (_, hi) in zip(weights, intervals, strict=True)]
    upper = [common.logarithm_interval(w/lo)[1] for w, (lo, _) in zip(weights, intervals, strict=True)]
    return min(F(1), abs(s)*(max(upper)-min(lower))/2)


@lru_cache(maxsize=4)
def _clock(weights, s):
    array = full.compiler.execute(common.diagonal_word(s, ratios(weights)), 9)
    return array.dtype.str, array.tobytes()


def clock(weights, s=F(1)):
    dtype, data = _clock(_weights(weights), full.parameter(s))
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


@lru_cache(maxsize=1)
def _rotation_certificate():
    for index in (0, 1):
        for sign in (-1, 1):
            full.corner.report._gate(cap.certify_rotation(index, sign), cap.ROTATION_GATE,
                                     'complete cap analyzer')
    return True


@lru_cache(maxsize=4)
def _certificate(weights, s):
    full.corner.report._gate(common.certify_diagonal(s, ratios(weights), clock(weights, s)),
                             WORD_GATE, 'record-inferred modular source word')
    _rotation_certificate()
    error, gain = source_bounds()
    word = common.diagonal_word(s, ratios(weights))
    pulses = sum(2*len(cap.rotation_word(i, 1))+len(word) for i in (0, 1))
    return tuple({'source_history_halfdiamond_upper': error,
                  'source_channel_norm_upper': gain,
                  'source_pulses': pulses, 'proper_code_service_events': 2*pulses,
                  'current_reads': 2, 'native_exports': 2,
                  'clock_word_sha256': hashlib.sha256(full.corner.canonical(word).encode()).hexdigest()
                  }.items())


def certificate(weights, s=F(1)):
    return dict(_certificate(_weights(weights), full.parameter(s)))


def source_bounds():
    """Fixed precision caps, independent of the training records."""
    error, gain = F(), F(1)
    for epsilon in (cap.ROTATION_GATE, WORD_GATE, cap.ROTATION_GATE)*2:
        error = (1+epsilon)**2*error+epsilon*(2+epsilon)/2
        gain *= (1+epsilon)**2
    return error, gain


def execute(density, weights, law, *, s=F(1), raw=True, feedback=True):
    density = history._density(density)
    weights, s = _weights(weights), full.parameter(s)
    if (type(law) is not str or law not in ('C', 'E')
            or type(raw) is not bool or type(feedback) is not bool):
        raise ValueError('an original native law and Boolean switches required')
    if raw:
        certificate(weights, s)
        unitary = clock(weights, s)
    else:
        values = ratios(weights)
        phases = {w: np.exp(-1j*float(s)*float(sum(common.logarithm_interval(w))/2))
                  for w in set(values)}
        unitary = np.diag([phases[w] for w in values])
    rotation = cap.rotation if raw else cap.ideal_rotation
    stages = {(index, sign): rotation(index, sign)@unitary@rotation(index, -1)
              for index in (0, 1) for sign in (-1, 1)}
    return cap._execute_stages(density, law, lambda index, sign: stages[index, sign], feedback=feedback)


def prepared(weights, law, *, s=F(1), raw=True, initial=0):
    """Continuing validation history; finite preparation kernel, all flags."""
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('one native C/E law required')
    if type(initial) is not int or not 0 <= initial < 66:
        raise ValueError('an original native input index required')
    maps = local.population_maps(law=law)
    output = {('preparation', flag): np.diag(a[:, initial]).astype(complex)
              for flag, a in maps.items() if flag != 'success'}
    output.update(execute(np.diag(maps['success'][:, initial]).astype(complex), weights, law, s=s, raw=raw))
    return output


def comparison_bound(weights, intervals, s=F(1)):
    """Complete learned history vs true ideal reference, on confidence event.

    The unnormalized inference target is b=q*sigma. Its scalar cancels in
    the action. Each ideal read is contractive, so the two clock errors add
    even with adaptive analyzer selection and arbitrary untouched references.
    Statistical failure probability is not included in this conditional bound.
    Raw execution must first pass certificate(weights, s); the executable
    record-admission protocol enforces that premise or halts explicitly.
    """
    source_error, gain = source_bounds()
    error = inference_bound(weights, intervals, s)
    prep = local.budget()
    return (gain*prep['complete_flagged_reference_halftrace_upper']
            +prep['uniform_ideal_herald_probability']*(source_error+2*error))


def admission(counts, *, sampling_bias=F()):
    """Fixed-budget protocol: learn at s=1 or halt on the training carrier.

    No data-dependent error is integrated as if it were a uniform bound.
    On accepted confidence-event transcripts the error is <=9/100; on the
    complementary event it is <=(G+1)/2. Both comparison protocols retain
    the identical training transcript and halt state when admission fails.
    This joint bound does not apply after conditioning on acceptance.
    """
    from . import reference_tomography as readout
    inferred = readout.infer_counts(counts, sampling_bias=sampling_bias)
    decision = {'status': 'stop', 'reason': 'unresolved_reference',
                'acceptance_cap': ACCEPTANCE_CAP, 'inference': inferred}
    if inferred['trials'] != TRAINING_TRIALS:
        decision['reason'] = 'fixed_trial_budget_required'
    elif inferred['status'] == 'faithful':
        bound = comparison_bound(inferred['central_subnormalized_weights'],
                                 inferred['subnormalized_weight_intervals'])
        decision['observed_confidence_event_bound'] = bound
        if bound <= ACCEPTANCE_CAP:
            decision.update(status='accept', reason='fixed_error_cap_satisfied')
        else:
            decision['reason'] = 'confidence_box_too_wide'
    decision['complete_protocol_halftrace_upper'] = (
        ACCEPTANCE_CAP+inferred['confidence_failure_upper']*(source_bounds()[1]+1)/2)
    return decision


def prepared_from_records(counts, law, *, initial=0, sampling_bias=F()):
    """Executable continuation, including explicit training/compilation stops."""
    if type(initial) is not int or not 0 <= initial < 66:
        raise ValueError('an original native input index required')
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('one native C/E law required')
    decision = admission(counts, sampling_bias=sampling_bias)
    rho = np.diag([float(i == initial) for i in range(66)]).astype(complex)
    if decision['status'] != 'accept':
        return decision, {('training_stop',): rho}
    weights = decision['inference']['central_subnormalized_weights']
    try:
        certificate(weights)
    except ArithmeticError:
        decision.update(status='stop', reason='source_certification_failed')
        return decision, {('compilation_stop',): rho}
    return decision, prepared(weights, law, initial=initial)


def run():
    from . import reference_tomography as readout
    training = readout.run()
    data = training['continuing_execution']
    inference = data['inference']
    if inference['status'] != 'faithful':
        raise ArithmeticError('retained training cannot certify a faithful clock')
    weights = inference['central_subnormalized_weights']
    intervals = inference['subnormalized_weight_intervals']
    error = comparison_bound(weights, intervals)
    phase_error = inference_bound(weights, intervals)
    source = certificate(weights)
    rows = []
    for law in ('C', 'E'):
        # Start future preparation on the last sampled terminal state, not a
        # silently supplied independent copy. Its uniform certificate also
        # covers any correlated quantum input outside this diagonal example.
        decision, actual = prepared_from_records(data['execution']['counts'], law,
            initial=data['execution']['final_state'], sampling_bias=readout.SAMPLER_GATE)
        if decision['status'] != 'accept':
            raise ArithmeticError('retained training does not pass the fixed acceptance protocol')
        residual = history.distance(actual, cap.target_history(law))
        if residual > float(error):
            raise ArithmeticError('learned continued history exceeds the confidence-event bound')
        probabilities = history.probabilities(actual)
        rows.append({'law': law, 'complete_output_halftrace_diagnostic': residual,
                     'complete_output_mass_diagnostic': sum(probabilities.values()),
                     'record_probabilities': {str(key): value for key, value in probabilities.items()}})
    paths = set(training['source_sha256']) | set(common.SOURCE_PATHS) | set(local.SOURCE_PATHS)
    paths.update(('model/native_modular_reference/learned_clock.py',
                  'model/native_modular_reference/cap_local_clock.py',
                  'model/native_modular_reference/local_history.py'))
    return {'schema': 'oph.native_modular_reference.learned_clock.v1',
            'training_trials': data['execution']['trials'],
            'training_counts': data['execution']['counts'],
            'training_record_and_state_sha256': data['execution']['record_and_terminal_state_sha256'],
            'continuation_initial_native_state': data['execution']['final_state'],
            'inferred_subnormalized_weights': weights,
            'inferred_weight_intervals': intervals,
            'parameter': '1', 'source_certificate': source,
            'inference_clock_halfdiamond_upper': phase_error,
            'statistical_failure_upper': inference['confidence_failure_upper'],
            'complete_ideal_reference_history_halftrace_upper_on_confidence_event': error,
            'acceptance_protocol': {
                'status': decision['status'], 'fixed_confidence_event_error_cap': ACCEPTANCE_CAP,
                'complete_joint_protocol_halftrace_upper': decision['complete_protocol_halftrace_upper'],
                'same_training_transcript_and_halt_carrier_in_comparison': True,
                'bound_conditioned_on_acceptance': False},
            'two_read_channel_vs_actual_prepared_modular_upper_on_confidence_event':
                source['source_history_halfdiamond_upper']+2*phase_error
                +2*local.budget()['prepared_modular_generator_operator_error_upper'],
            'training_native_attempts_upper': inference['source_native_attempts_all_trials_upper'],
            'future_preparation_native_attempts_upper': local.budget()['total_native_attempts_upper'],
            'histories': rows,
            'source_sha256': {p: hashlib.sha256((common.ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)},
            'scope': {'source_derived_family_required': True,
                      'known_weights_used_to_center_estimator_or_compile_clock': False,
                      'prepared_reference_used_only_for_source_certificate_and_validation': True,
                      'coherent_source_services_required_for_enactment': True,
                      'confidence_theorem_authenticates_fixed_seed_PRNG_randomness': False,
                      'intermediate_pulse_locality_or_physical_units_established': False,
                      'entire_microscopic_training_tape_executed': False}}


def compare_receipt(supplied, actual):
    full.corner.compare_receipt(supplied, actual)
    for receipt in (supplied, actual):
        for row in receipt['histories']:
            for key in ('complete_output_halftrace_diagnostic', 'complete_output_mass_diagnostic'):
                x = row[key]
                if type(x) is not float or not np.isfinite(x) or x < 0:
                    raise ValueError('nonnegative finite history diagnostics required')
            if abs(row['complete_output_mass_diagnostic']-1) > 2e-10:
                raise ValueError('complete learned history mass changed')
            if any(type(p) is not float or not np.isfinite(p) or p < 0 for p in row['record_probabilities'].values()):
                raise ValueError('nonnegative complete branch probabilities required')
            if (any(p > 1+2e-10 for p in row['record_probabilities'].values())
                    or abs(sum(row['record_probabilities'].values())-row['complete_output_mass_diagnostic']) > 2e-10):
                raise ValueError('branch probabilities disagree with complete history mass')
            if row['complete_output_halftrace_diagnostic'] > float(F(
                    receipt['complete_ideal_reference_history_halftrace_upper_on_confidence_event'])):
                raise ValueError('history diagnostic exceeds the confidence-event certificate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', type=Path)
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    payload = strict_json.loads(full.preparation.canonical(run()))
    if args.check:
        compare_receipt(strict_json.load(args.check), payload)
        print('Record-inferred common source clock and continuing histories verified.')
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(full.corner.canonical(payload).rstrip('\n')+'\n')


if __name__ == '__main__':
    main()
