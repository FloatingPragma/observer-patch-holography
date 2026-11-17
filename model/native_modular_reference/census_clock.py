"""An accurate successful common-clock history from native census records.

Admission uses the normalized successful history, before its herald mass is
applied. The census estimates the reference; it does not enact its modular
action. That action still uses the existing certified coherent source words.
All training stops, preparation failures, current outcomes and exports remain
in the complete protocol. Statistical failure is separate from the error on
the simultaneous confidence event; no admission-conditioned guarantee follows.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from . import learned_clock as learned, reference_tomography as prior
from . import native_census as census

common, full, local = learned.common, learned.full, learned.local
history, cap = learned.history, learned.cap
TRAINING_TRIALS = 1 << 20
SUCCESSFUL_ERROR_CAP = F(1, 10)
OUTPUT = Path(__file__).with_name('census_clock_results.json')
SOURCE_PATHS = tuple(sorted(set(census.SOURCE_PATHS) | set(common.SOURCE_PATHS) | set(local.SOURCE_PATHS) | {
    'model/native_modular_reference/census_clock.py',
    'model/native_modular_reference/learned_clock.py',
    'model/native_modular_reference/cap_local_clock.py',
    'model/native_modular_reference/local_history.py'}))


def _scope():
    return {
        'source_derived_charge_family_required': True,
        'known_weights_used_to_center_estimator_compile_or_admit_clock': False,
        'successful_channel_comparison_covers_arbitrary_work_and_reference_inputs': True,
        'successful_state_history_uses_uniform_conditional_preparation_theorem': True,
        'preparation_success_diagnostics_are_conditioned_but_raw_arrays_are_not_normalized': True,
        'source_numerical_error_is_separate_from_physical_service_accuracy': True,
        'coherent_source_services_required_for_enactment': True,
        'native_event_count_identified_with_modular_parameter_or_physical_seconds': False,
        'strictly_monotone_individual_random_record_clock_established': False,
        'confidence_theorem_authenticates_fixed_seed_PRNG_randomness': False,
        'entire_microscopic_training_tape_executed': False}


def _record_alphabet(law):
    """Complete branch labels from the two native instruments, without arrays."""
    keys = {('preparation', flag) for flag in ('herald_failure', 'warmup_timeout', 'completion_timeout')}
    keys.add(('abort', 0))
    first = cap.interface.fixed_maps(law, cap.SEAMS[0])
    later = cap.interface.fixed_maps(law, cap.SEAMS[1])
    for branch in first:
        record = branch.coin, branch.label, branch.current
        keys.add(('abort', 1, *record))
        for second in later:
            keys.add(('read', *record, second.coin, second.label, second.current))
    return {str(key) for key in keys}


def inference_bound(weights, intervals, s=F(1)):
    """Reference-uniform half-diamond bound modulo one common phase.

    For h_i=log(w_hat_i/w_i), osc(h)=max_{i != j}(h_i-h_j).
    The two extrema cannot be different endpoints of the SAME coordinate.
    This distinction matters when one confidence interval dominates the box.
    Every possible generator difference has norm <=osc(h)/2 after subtracting
    its scalar midpoint. Duhamel supplies |s| osc(h)/2 on arbitrary references.
    The 446 padding addresses are an implementation extension, not a fourth
    weight of the native reference or part of this native channel comparison.
    """
    def source_fraction(value):
        if (type(value) is not F or type(value.numerator) is not int
                or type(value.denominator) is not int or value.denominator <= 0):
            raise ValueError('original Fraction integer slots and a positive denominator required')
        return F(value.numerator, value.denominator)
    if type(weights) is not tuple or len(weights) != 3:
        raise ValueError('three original positive Fraction charge weights required')
    weights = learned._weights(tuple(source_fraction(w) for w in weights))
    s = full.parameter(source_fraction(s) if type(s) is F else s)
    if (type(intervals) is not tuple or len(intervals) != 3
            or any(type(box) is not tuple or len(box) != 2 for box in intervals)):
        raise ValueError('three faithful exact intervals containing their centers required')
    intervals = tuple(tuple(source_fraction(x) for x in box) for box in intervals)
    if any(not 0 < lo <= w <= hi for w, (lo, hi) in zip(weights, intervals, strict=True)):
        raise ValueError('three faithful exact intervals containing their centers required')
    spread = max(common.logarithm_interval(
        weights[i]*intervals[j][1]/(intervals[i][0]*weights[j]))[1]
        for i in range(3) for j in range(3) if i != j)
    return min(F(1), abs(s)*max(F(), spread)/2)


def bounds(weights, intervals, s=F(1)):
    """Successful channel and state bounds without any herald dilution.

    The preparer theorem gives half-trace error B on the complete flagged
    output and ideal herald probability q. Its normalized successful output
    is within 2B/(q-B) of sigma tensor the original reference marginal. Apply
    the raw source history (norm <=G), then its source error and the two
    ideal clock comparisons. The latter compose through adaptive cap reads.
    Only source budget scalars enter this contract, not target weight lookup.
    """
    phase = inference_bound(weights, intervals, s)
    source_error, gain = learned.source_bounds()
    prep_error, q, cost = prior._record_contract()
    if prep_error >= q:
        raise ArithmeticError('the source preparer has no certified success floor')
    conditional_prep = 2*prep_error/(q-prep_error)
    channel = source_error+2*phase
    return {
        'inference_clock_halfdiamond_upper': phase,
        'successful_two_read_channel_halfdiamond_upper': channel,
        'conditional_preparation_halftrace_upper': conditional_prep,
        'successful_reference_history_halftrace_upper': gain*conditional_prep+channel,
        'complete_flagged_history_halftrace_upper': gain*prep_error+q*channel,
        'future_preparation_success_probability_interval': (q-prep_error, q+prep_error),
        'future_preparation_native_attempts_upper': cost,
        'source_channel_norm_upper': gain,
    }


def admission(counts, *, sampling_bias=F()):
    """Fixed successful-output cap; all unsuccessful training explicitly stops.

    On confidence-event transcripts admission bounds the channel error by
    S-G*2B/(q-B), using FIXED source constants B,q,G and successful cap S.
    Thus complete error is <=G*B+q*max(0,S-G*2B/(q-B)). Off the event its
    raw/ideal mass bound is (G+1)/2. Averaging adds alpha*(G+1)/2,
    keeping the identical training transcript and halt carrier in both sides.
    Neither this average nor alpha is conditioned on acceptance.
    """
    inferred = census.infer_counts(counts, sampling_bias=sampling_bias)
    decision = {'status': 'stop', 'reason': 'unresolved_reference',
                'fixed_successful_error_cap': F(SUCCESSFUL_ERROR_CAP.numerator,
                                                SUCCESSFUL_ERROR_CAP.denominator),
                'inference': inferred}
    if inferred['trials'] != TRAINING_TRIALS:
        decision['reason'] = 'fixed_trial_budget_required'
    elif inferred['status'] == 'faithful':
        certificate = bounds(inferred['central_subnormalized_weights'],
                             inferred['subnormalized_weight_intervals'])
        decision['confidence_event_bounds'] = certificate
        if certificate['successful_reference_history_halftrace_upper'] <= SUCCESSFUL_ERROR_CAP:
            decision.update(status='accept', reason='successful_history_error_cap_satisfied')
        else:
            decision['reason'] = 'successful_history_confidence_box_too_wide'
    prep_error, q, _ = prior._record_contract()
    gain = learned.source_bounds()[1]
    channel_cap = max(F(), SUCCESSFUL_ERROR_CAP-gain*2*prep_error/(q-prep_error))
    complete_cap = gain*prep_error+q*channel_cap
    decision['fixed_successful_channel_error_cap'] = channel_cap
    decision['fixed_complete_confidence_event_error_cap'] = complete_cap
    decision['complete_joint_protocol_halftrace_upper'] = (
        complete_cap+inferred['confidence_failure_upper']*(gain+1)/2)
    return decision


def prepared_from_records(counts, law, *, initial=0, sampling_bias=F()):
    """Continue on the terminal training carrier, with source-certification stop."""
    if type(initial) is not int or not 0 <= initial < 66:
        raise ValueError('an original native input index required')
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('one original native C/E law required')
    decision = admission(counts, sampling_bias=sampling_bias)
    rho = np.diag([float(i == initial) for i in range(66)]).astype(complex)
    if decision['status'] != 'accept':
        return decision, {('training_stop',): rho}
    weights = decision['inference']['central_subnormalized_weights']
    try:
        learned.certificate(weights)
    except ArithmeticError:
        decision.update(status='stop', reason='source_certification_failed')
        return decision, {('compilation_stop',): rho}
    return decision, learned.prepared(weights, law, initial=initial)


def run():
    training = census.run()
    data = training['continuing_execution']
    execution = data['execution']
    counts, initial = execution['counts'], execution['final_state']
    decision = admission(counts, sampling_bias=census.SAMPLER_GATE)
    if decision['status'] != 'accept':
        raise ArithmeticError('the retained census does not meet the successful history cap')
    inferred = decision['inference']
    weights = inferred['central_subnormalized_weights']
    budget = decision['confidence_event_bounds']
    # Fraction's private slots remain writable even though normal arithmetic
    # treats it as immutable. Detach the cached provider's scalar objects.
    source = {key: F(value.numerator, value.denominator) if type(value) is F else value
              for key, value in learned.certificate(weights).items()}
    rows = []
    for law in ('C', 'E'):
        checked, actual = prepared_from_records(counts, law, initial=initial,
                                                sampling_bias=census.SAMPLER_GATE)
        if checked['status'] != 'accept':
            raise ArithmeticError('accepted census failed source-word admission')
        # Conditioning here is a reported state diagnostic. No executed source
        # array, successful output mass or export distribution is normalized.
        p = float(local.population_maps(law=law)['success'][:, initial].sum())
        successful = {key: value/p for key, value in actual.items() if key[0] != 'preparation'}
        sigma = np.diag(list(map(float, local.reference()))).astype(complex)
        target_success = cap.execute(sigma, law, raw=False)
        residual = history.distance(actual, cap.target_history(law))
        success_residual = history.distance(successful, target_success)
        if (residual > float(budget['complete_flagged_history_halftrace_upper'])
                or success_residual > float(budget['successful_reference_history_halftrace_upper'])):
            raise ArithmeticError('continued census-clock history exceeds its exact budget')
        probabilities = history.probabilities(actual)
        rows.append({'law': law, 'preparation_success_probability_diagnostic': p,
                     'complete_output_halftrace_diagnostic': residual,
                     'successful_output_halftrace_diagnostic': success_residual,
                     'complete_output_mass_diagnostic': sum(probabilities.values()),
                     'record_probabilities': {str(key): value for key, value in probabilities.items()}})
    return {
        'schema': 'oph.native_modular_reference.census_clock.v1',
        'training_trials': execution['trials'], 'training_counts': counts,
        'training_record_and_state_sha256': execution['record_and_terminal_state_sha256'],
        'continuation_initial_native_state': initial,
        'inferred_subnormalized_weights': weights,
        'inferred_weight_intervals': inferred['subnormalized_weight_intervals'],
        'parameter': '1', 'source_certificate': source,
        'confidence_event_bounds': budget,
        'statistical_failure_upper': inferred['confidence_failure_upper'],
        'acceptance_protocol': {
            'status': decision['status'],
            'fixed_successful_error_cap': F(SUCCESSFUL_ERROR_CAP.numerator,
                                            SUCCESSFUL_ERROR_CAP.denominator),
            'fixed_successful_channel_error_cap': decision['fixed_successful_channel_error_cap'],
            'fixed_complete_confidence_event_error_cap': decision['fixed_complete_confidence_event_error_cap'],
            'complete_joint_protocol_halftrace_upper': decision['complete_joint_protocol_halftrace_upper'],
            'same_training_transcript_and_halt_carrier_in_comparison': True,
            'bound_conditioned_on_acceptance': False,
            'successful_error_multiplied_by_herald_probability_for_admission': False},
        'training_native_attempts_upper': inferred['source_native_attempts_all_trials_upper'],
        'two_read_channel_vs_actual_prepared_modular_upper_on_confidence_event':
            budget['successful_two_read_channel_halfdiamond_upper']
            +2*local.budget()['prepared_modular_generator_operator_error_upper'],
        'histories': rows,
        'source_sha256': {p: hashlib.sha256((common.ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'scope': _scope()}


def compare_receipt(supplied, actual=None):
    """Verify against a fresh source replay, never a caller-supplied oracle."""
    owned = strict_json.loads(full.preparation.canonical(run()))
    _compare_receipt(supplied, owned)
    if actual is not None:
        _compare_receipt(actual, owned)


def _compare_receipt(supplied, actual):
    """Exact proof fields; independent consistency checks on float diagnostics."""
    full.corner.compare_receipt(supplied, actual)
    for receipt in (supplied, actual):
        if receipt['schema'] != 'oph.native_modular_reference.census_clock.v1' or receipt['parameter'] != '1':
            raise ValueError('the declared census-clock schema and parameter one are required')
        initial = receipt['continuation_initial_native_state']
        digest = receipt['training_record_and_state_sha256']
        if type(initial) is not int or not 0 <= initial < 66:
            raise ValueError('an original terminal training state is required')
        if type(digest) is not str or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('an original SHA256 training stream digest is required')
        full.corner.compare_receipt(receipt['scope'], _scope())
        if (type(receipt['histories']) is not list or len(receipt['histories']) != 2
                or [row.get('law') for row in receipt['histories']] != ['C', 'E']):
            raise ValueError('exactly one complete history for each native law is required')
        # Reconstruct the proof-bearing census and admission fields even if a
        # caller supplies the same forged object as both comparison arguments.
        decision = admission(receipt['training_counts'], sampling_bias=census.SAMPLER_GATE)
        if decision['status'] != 'accept' or receipt['training_trials'] != TRAINING_TRIALS:
            raise ValueError('retained records do not admit the declared fixed-budget clock')
        inferred = decision['inference']
        weights = inferred['central_subnormalized_weights']
        word = common.diagonal_word(F(1), learned.ratios(weights))
        pulses = sum(2*len(cap.rotation_word(i, 1))+len(word) for i in (0, 1))
        source_error, gain = learned.source_bounds()
        expected = {
            'inferred_subnormalized_weights': inferred['central_subnormalized_weights'],
            'inferred_weight_intervals': inferred['subnormalized_weight_intervals'],
            'confidence_event_bounds': decision['confidence_event_bounds'],
            'statistical_failure_upper': inferred['confidence_failure_upper'],
            'training_native_attempts_upper': inferred['source_native_attempts_all_trials_upper'],
            'two_read_channel_vs_actual_prepared_modular_upper_on_confidence_event':
                decision['confidence_event_bounds']['successful_two_read_channel_halfdiamond_upper']
                +2*local.budget()['prepared_modular_generator_operator_error_upper'],
            'source_certificate': {
                'source_history_halfdiamond_upper': source_error,
                'source_channel_norm_upper': gain, 'source_pulses': pulses,
                'proper_code_service_events': 2*pulses, 'current_reads': 2, 'native_exports': 2,
                'clock_word_sha256': hashlib.sha256(full.corner.canonical(word).encode()).hexdigest()},
            'acceptance_protocol': {
                'status': 'accept', 'fixed_successful_error_cap': SUCCESSFUL_ERROR_CAP,
                'fixed_successful_channel_error_cap': decision['fixed_successful_channel_error_cap'],
                'fixed_complete_confidence_event_error_cap': decision['fixed_complete_confidence_event_error_cap'],
                'complete_joint_protocol_halftrace_upper': decision['complete_joint_protocol_halftrace_upper'],
                'same_training_transcript_and_halt_carrier_in_comparison': True,
                'bound_conditioned_on_acceptance': False,
                'successful_error_multiplied_by_herald_probability_for_admission': False}}
        expected = strict_json.loads(full.preparation.canonical(expected))
        full.corner.compare_receipt({key: receipt[key] for key in expected}, expected)
        if set(receipt['source_sha256']) != set(SOURCE_PATHS):
            raise ValueError('the complete declared source custody is required')
        for path, digest in receipt['source_sha256'].items():
            if hashlib.sha256((common.ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('retained census-clock source digest changed')
        budget = receipt['confidence_event_bounds']
        for row in receipt['histories']:
            for key in ('preparation_success_probability_diagnostic',
                        'complete_output_halftrace_diagnostic', 'successful_output_halftrace_diagnostic',
                        'complete_output_mass_diagnostic'):
                x = row[key]
                if type(x) is not float or not np.isfinite(x) or x < 0:
                    raise ValueError('nonnegative finite history diagnostics required')
            if abs(row['complete_output_mass_diagnostic']-1) > 2e-10:
                raise ValueError('complete census-clock output mass changed')
            lo, hi = map(F, budget['future_preparation_success_probability_interval'])
            if not float(lo)-2e-10 <= row['preparation_success_probability_diagnostic'] <= float(hi)+2e-10:
                raise ValueError('successful preparation diagnostic escaped its exact interval')
            probabilities = row['record_probabilities']
            if type(probabilities) is not dict or set(probabilities) != _record_alphabet(row['law']):
                raise ValueError('the complete preparation/current/export alphabet is required')
            if any(type(p) is not float or not np.isfinite(p) or p < 0 or p > 1+2e-10
                   for p in probabilities.values()):
                raise ValueError('nonnegative complete branch probabilities required')
            if abs(sum(probabilities.values())-row['complete_output_mass_diagnostic']) > 2e-10:
                raise ValueError('record probabilities disagree with complete output mass')
            preparation_mass = sum(probabilities[str(('preparation', flag))]
                                   for flag in ('herald_failure', 'warmup_timeout', 'completion_timeout'))
            if abs(preparation_mass+row['preparation_success_probability_diagnostic']-1) > 2e-10:
                raise ValueError('preparation success and retained failures disagree')
            for diagnostic, bound in (
                    ('complete_output_halftrace_diagnostic', 'complete_flagged_history_halftrace_upper'),
                    ('successful_output_halftrace_diagnostic', 'successful_reference_history_halftrace_upper')):
                if row[diagnostic] > float(F(budget[bound])):
                    raise ValueError('history diagnostic exceeds its confidence-event certificate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', type=Path)
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    payload = strict_json.loads(full.preparation.canonical(run()))
    if args.check:
        _compare_receipt(strict_json.load(args.check), payload)
        print('Native census common-clock successful and complete histories verified.')
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(full.corner.canonical(payload).rstrip('\n')+'\n')


if __name__ == '__main__':
    main()
