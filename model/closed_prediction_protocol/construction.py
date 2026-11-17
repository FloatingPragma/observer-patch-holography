"""Finite sequential C/E test with no exact reset-state or IID premise.

The exact-source model is distinguished from floating execution and from a
physical device. All parameter counts denote native repair attempts.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..compiled_native_interface import construction as compiled
from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..readout_controls import prospective
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import compare_receipt

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
RECORDS = tuple(product(range(3), range(2)))
SOURCE_PATHS = tuple(sorted(set((*compiled.SOURCE_PATHS,
    'model/native_equilibration/construction.py',
    'model/native_equilibration/README.md',
    'model/readout_controls/prospective.py',
    'model/calibrated_measurement/bounds.py',
    'model/closed_prediction_protocol/__init__.py',
    'model/closed_prediction_protocol/__main__.py',
    'model/closed_prediction_protocol/construction.py',
    'model/closed_prediction_protocol/protocol.json'))))


def preparation_count(error):
    """Smallest n certified by 33*(229/240)**n <= an exact allowance.

    The gap follows from (9/4)^2>5 and the existing primitive gap theorem.
    This is a sharper use of the same bound, not a new mixing theorem.
    """
    if type(error) is not F or not 0 < error < 1:
        raise ValueError('an exact preparation allowance in (0,1) is required')
    if F(9, 4)**2 <= 5:
        raise ArithmeticError('the rational gap comparison failed')
    count, bound = 0, F(33)
    while bound > error:
        count += 1
        bound *= F(229, 240)
        if count > 10000:
            raise ValueError('preparation allocation exceeds this finite protocol')
    return count


@lru_cache(maxsize=1)
def design():
    protocol = strict_json.load(FOLDER/'protocol.json')
    if (protocol['candidate_laws'] != ['C', 'E'] or protocol['work_dimension'] != 66
            or protocol['binary_total'] != 2 or protocol['primitive_gap_lower'] != '11/240'
            or protocol['unread_attempt_designs'] != list(range(1, 9))
            or protocol['selection'] != 'minimum worst-case total native attempts; tie by smaller unread count'
            or protocol['executed_sequential_trials'] != 2):
        raise ValueError('the locked source design changed')
    ideal, prepare, services, target = (F(protocol[key]) for key in (
        'ideal_decision_error', 'whole_experiment_preparation_error',
        'whole_experiment_additional_service_allowance', 'total_decision_error_target'))
    if not all(0 < x < 1 for x in (ideal, prepare, services, target)) or ideal+prepare+services > target:
        raise ValueError('positive allowances must fit the declared error target')
    rows = []
    for attempts in protocol['unread_attempt_designs']:
        first, second = (prospective.exact_law(law, attempts) for law in ('C', 'E'))
        certificate = prospective.trial_design(first, second, ideal)
        trials = certificate['fixed_total_trials']
        epsilon = prepare/trials
        count = preparation_count(epsilon)
        rows.append({'unread_attempts': attempts, 'fixed_trials': trials,
            'discarded_preparation_attempts_per_trial': count,
            'preparation_trace_distance_upper': str(epsilon),
            'maximum_native_attempts': trials*(count+attempts+2),
            'likelihood_certificate': certificate})
    selected = min(rows, key=lambda row: (row['maximum_native_attempts'], row['unread_attempts']))
    return {'locked_design_search': rows, 'selected': selected,
        'error_budget': {'ideal_product_decision': str(ideal),
            'sequential_preparation_allowance': str(prepare),
            'exact_source_model_decision_upper': str(ideal+prepare),
            'additional_conditional_trial_TV_allowance': str(services/selected['fixed_trials']),
            'decision_upper_if_additional_allowance_is_met': str(target),
            'additional_physical_allowance_authenticated': False},
        'selection_uses_no_measurement_outcomes': True}


@lru_cache(maxsize=2)
def native_sparse(law):
    """All source matrix-unit paths of one unread native attempt.

    Count native partial-injection paths exactly, then apply their integer
    counts/60. No spectral exponential or finite-time reset is substituted.
    """
    law = quantum.candidate(law)
    counts = Counter()
    groups = 0
    for seam in quantum.mixing.seams():
        for branch in interface.fixed_maps(law, seam):
            groups += 1
            for a, out_a in branch.mapping:
                for b, out_b in branch.mapping:
                    counts[out_a*66+out_b, a*66+b] += 1
    rows = sorted(counts)
    output = np.array([row[0] for row in rows], dtype=int)
    inputs = np.array([row[1] for row in rows], dtype=int)
    weights = np.array([counts[row] for row in rows], dtype=float)/60
    # Counts, not floating sums, establish preservation of every diagonal
    # input trace. The source maps are partial injections, so off-diagonal
    # inputs have zero output trace separately.
    traces = Counter()
    for (out, source), count in counts.items():
        if out//66 == out % 66:
            traces[source] += count
    if any(traces[i*66+j] != (60 if i == j else 0)
           for i in range(66) for j in range(66)):
        raise ArithmeticError('the actual unread source channel lost its exact trace')
    return output, inputs, weights, groups


def prepare(density, law, count):
    """Execute every discarded source attempt, retaining reference coherences."""
    if type(count) is not int or not 0 <= count <= 10000:
        raise ValueError('bounded nonnegative preparation attempt count required')
    if (type(density) is not np.ndarray or density.ndim != 2 or density.shape[0] != density.shape[1]
            or density.shape[0] not in (66, 132, 198) or not np.isfinite(density).all()):
        raise ValueError('a finite 66D work matrix with up to three reference levels is required')
    reference = density.shape[0]//66
    # Minor reference factors remain untouched. Bincount adds all Kraus
    # alternatives sharing an output matrix unit, with no amplitude cutoff.
    values = np.asarray(density, dtype=complex).reshape(66, reference, 66, reference).transpose(0, 2, 1, 3).reshape(4356, -1).copy()
    output, inputs, weights, _ = native_sparse(law)
    for _ in range(count):
        following = np.empty_like(values)
        for column in range(reference*reference):
            weighted = weights*values[inputs, column]
            following[:, column] = (np.bincount(output, weights=weighted.real, minlength=4356)
                + 1j*np.bincount(output, weights=weighted.imag, minlength=4356))
        values = following
    return values.reshape(66, 66, reference, reference).transpose(0, 2, 1, 3).reshape(density.shape)


def trial_outputs(density, law, attempts):
    """Execute the existing complete instrument on unnormalized CQ branches.

    This calls the shared runtime's native Kraus application at each event.
    It neither resets nor normalizes branches between consecutive trials.
    """
    program = interface.program(law, attempts)
    reference = density.shape[0]//66
    branches = {(): density}
    for event in program.events:
        following = {}
        for record, rho in branches.items():
            instrument = event.instruments[0 if event.control is None else record[event.control]]
            for outcome, result in enumerate(runtime._apply(instrument, rho, reference)):
                following[(*record, outcome) if event.visible else record] = result
        branches = following
    if tuple(branches) != RECORDS:
        raise ArithmeticError('the complete six-record experiment changed')
    return branches


def _distance(first, second):
    return float(np.linalg.norm(first-second, 'nuc')/2)


def classify(counts):
    """The locked future decision counts every trial, including all failures."""
    selected = design()['selected']
    first, second = (prospective.exact_law(law, selected['unread_attempts']) for law in ('C', 'E'))
    return prospective.classify(counts, first, second, selected['fixed_trials'])


def _serialize(density):
    return json.dumps({'shape': list(density.shape), 'real': density.real.tolist(),
                       'imag': density.imag.tolist()}, allow_nan=False, separators=(',', ':'))


def _restore(payload, expected):
    """Replay-validate a local preparation checkpoint against its full prefix."""
    value = strict_json.loads(payload)
    if (type(value) is not dict or set(value) != {'shape', 'real', 'imag'}
            or value['shape'] != list(expected.shape)):
        raise ValueError('checkpoint shape/schema mismatch')
    real, imag = np.asarray(value['real']), np.asarray(value['imag'])
    if (real.shape != expected.shape or imag.shape != expected.shape
            or real.dtype.kind not in 'iuf' or imag.dtype.kind not in 'iuf'
            or not np.isfinite(real).all() or not np.isfinite(imag).all()):
        raise ValueError('finite replay-shaped checkpoint entries required')
    result = real+1j*imag
    if not np.array_equal(result, expected):
        raise ValueError('checkpoint does not replay its preparation prefix')
    return result


def execute(law, selected):
    """Two consecutive prepared trials on one carrier; no fresh-state reset."""
    attempts, count = selected['unread_attempts'], selected['discarded_preparation_attempts_per_trial']
    vector = np.zeros(66, complex)
    vector[0], vector[1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    initial = np.outer(vector, vector.conj())
    prepared = prepare(initial, law, count)
    half = count//2
    prefix = prepare(initial, law, half)
    restored = _restore(_serialize(prefix), prepare(initial, law, half))
    continued = prepare(restored, law, count-half)
    checkpoint_error = _distance(continued, prepared)
    distance = _distance(prepared, np.eye(66)/66)
    if distance > float(F(selected['preparation_trace_distance_upper']))+2e-10:
        raise ArithmeticError('the executed preparation violates its uniform certificate')
    # Independently run the shared validated full instrument and its native
    # herald checkpoint, then compare every unnormalized terminal density.
    program = interface.program(law, attempts)
    actual = runtime.execute(program, initial=prepared)
    checkpoint = runtime.execute(program, initial=prepared, stop_after=2)
    resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
        program, runtime.checkpoint_to_json(checkpoint)))
    first = trial_outputs(prepared, law, attempts)
    runtime_error = max(_distance(branch.density, first[branch.record]) for branch in actual.branches)
    exact = prospective.exact_law(law, attempts)
    first_law = [float(np.trace(first[key]).real) for key in RECORDS]
    rows, joint_tv, conditional_errors = [], 0., []
    for index, key in enumerate(RECORDS):
        # This is the actual previous terminal branch, not the stationary
        # state or a new density. Its trace carries the first-trial weight.
        next_input = prepare(first[key], law, count)
        second = trial_outputs(next_input, law, attempts)
        for j, second_key in enumerate(RECORDS):
            probability = float(np.trace(second[second_key]).real)
            expected = float(exact[index]*exact[j])
            joint_tv += abs(probability-expected)/2
            rows.append({'record': list((*key, *second_key)), 'probability': probability,
                         'ideal_product_probability': str(exact[index]*exact[j])})
        if first_law[index] > 0:
            conditional = [float(np.trace(second[k]).real)/first_law[index] for k in RECORDS]
            conditional_errors.append(sum(abs(p-float(q)) for p, q in zip(conditional, exact))/2)
    first_tv = sum(abs(p-float(q)) for p, q in zip(first_law, exact))/2
    epsilon = float(F(selected['preparation_trace_distance_upper']))
    if max(conditional_errors, default=0) > epsilon+2e-10 or joint_tv > 2*epsilon+2e-10:
        raise ArithmeticError('the sequential complete law violated the conditional bound')
    if max(checkpoint_error, runtime_error, runtime.complete_distance(actual, resumed)) > 2e-10:
        raise ArithmeticError('the full output/checkpoint execution disagrees')
    return {'candidate': law, 'executed_trials_on_same_carrier': 2,
        'original_input': '(|configuration 0> + i |configuration 1>)/sqrt(2)',
        'initial_preparation_trace_distance_diagnostic': distance,
        'first_complete_law': first_law,
        'first_law_TV_from_ideal_diagnostic': first_tv,
        'complete_two_trial_law': rows,
        'two_trial_TV_from_ideal_product_diagnostic': joint_tv,
        'maximum_second_conditional_law_TV_diagnostic': max(conditional_errors, default=0),
        'preparation_checkpoint_distance_diagnostic': checkpoint_error,
        'shared_runtime_all_terminal_densities_distance_diagnostic': runtime_error,
        'shared_runtime_checkpoint_distance_diagnostic': runtime.complete_distance(actual, resumed),
        'normalization_residual_diagnostic': abs(sum(row['probability'] for row in rows)-1),
        'initial_state_replaced_between_trials': False,
        'all_36_records_including_impossible_records_retained': True,
        'full_terminal_matrix_dimension': 66}


def reference_execution(law, selected):
    initial = interface.reference_input()
    reference = 2
    marginal = np.trace(initial.reshape(66, reference, 66, reference), axis1=0, axis2=2)
    count = selected['discarded_preparation_attempts_per_trial']
    # One attempt independently checks the sparse source execution against
    # all native Kraus maps lifted over a genuinely entangled reference.
    oracle, = runtime._apply(quantum.native_instrument(law), initial, reference)
    one_step = _distance(prepare(initial, law, 1), oracle)
    prepared = prepare(initial, law, count)
    target = np.kron(np.eye(66)/66, marginal)
    error = _distance(prepared, target)
    if one_step > 2e-10 or error > float(F(selected['preparation_trace_distance_upper']))+2e-10:
        raise ArithmeticError('source execution/reference-uniform preparation failed')
    return {'candidate': law, 'reference_dimension': 2,
        'one_step_full_native_Kraus_distance_diagnostic': one_step,
        'prepared_distance_from_product_diagnostic': error,
        'untouched_reference_marginal_distance_diagnostic': _distance(
            np.trace(prepared.reshape(66, reference, 66, reference), axis1=0, axis2=2), marginal)}


def run():
    chosen = design()
    selected = chosen['selected']
    attempts, n, count = (selected[k] for k in ('unread_attempts', 'fixed_trials',
        'discarded_preparation_attempts_per_trial'))
    source = prospective.execute(attempts)
    # This fresh source compilation is a raw-array comparison. It is not
    # charged as if it were a normalized physical law or a hardware error.
    full = compiled.full_matrix()
    if not np.array_equal(full, np.kron(np.eye(64), compiled.local_matrix())):
        raise ArithmeticError('the all-padding source lift changed')
    source_bound = compiled.certify()
    if source_bound['two_call_half_diamond_upper'] > F(1, 10**10):
        raise ArithmeticError('the numerical source envelope changed')
    success = 2*interface.SUCCESS
    mixer_resources = compiled.resources()
    return {'schema': 'oph.closed_prediction_protocol.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'design': chosen,
        'exact_selected_complete_laws': {law: list(map(str, prospective.exact_law(law, attempts))) for law in ('C', 'E')},
        'sequential_execution': [execute(law, selected) for law in ('C', 'E')],
        'entangled_reference_execution': [reference_execution(law, selected) for law in ('C', 'E')],
        'fresh_coherent_source': {'resources_per_mixer': mixer_resources,
            'full_source_dimension': 512, 'two_call_raw_array_half_diamond_upper': '1/10000000000',
            'all_padding_retained': True,
            'normalized_physical_instrument_error_certified': False,
            'ideal_runtime_is_target_surrogate_with_separate_source_comparison': True,
            'fresh_selected_trial_source_execution': source['candidates']},
        'resource_budget': {'native_attempts_per_failed_trial': count+1,
            'native_attempts_per_useful_trial': count+attempts+2,
            'maximum_total_native_attempts': selected['maximum_native_attempts'],
            'exact_product_reference_expected_native_attempts': str(n*(count+1+success*(attempts+1))),
            'maximum_total_mixers': 2*n,
            'maximum_source_pulses': 2*n*mixer_resources['native_pulses'],
            'maximum_load_export_events': 2*n*mixer_resources['proper_code_service_events'],
            'useful_herald_reference_probability': str(success),
            'all_failures_count_toward_fixed_trial_budget': True,
            'physical_duration_supplied': False},
        'theorem_attachment': {
            'existing_preparation_bound': 'native_equilibration; quantum mixing theorem, uniform over untouched references',
            'conditional_trial_law_TV': 'at most preparation allowance plus any additional conditional physical trial error',
            'sequential_joint_law_TV': 'at most N times the conditional trial allowance; by chain-rule telescoping',
            'iid_assumed_for_executed_process': False,
            'external_memory_is_untouched_during_native_preparation': True,
            'native_random_schedule_remains_declared_input': True,
            'coherent_proper_code_services_remain_declared_input': True,
            'fixed_total_selected_from_OPH_axioms': False,
            'physical_clock_or_laboratory_mapping_supplied': False,
            'new_empirical_observations': False}}


def compare(expected, actual):
    compare_receipt(expected, actual)
