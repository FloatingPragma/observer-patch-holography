"""Complete native-record phase calibration on a continuing finite carrier.

The decoder sees training and clock records, never the producer's parameter.
The numerical demonstration samples a declared dyadic population process;
raw full-source quantum histories and their error budgets are certified
separately. No microscopic native tape or autonomous laboratory clock is
claimed. Both preparation failures and numerical aborts remain public.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from tools import strict_json
from . import reference_tomography as training, learned_clock as learned
from . import phase_readout as readout, phase_inference as inverse

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('phase_clock_results.json')
GRID = 1 << 53
TRIALS = 1 << 22
SAMPLER_GATE = F(1, 10**9)
NUMERICAL_ABORT = ('numerical_abort',)
SCOPE = {'decoder_receives_programmed_parameter': False,
         'reference_amplitude_fitted_or_supplied_to_decoder': False,
         'raw_source_certificates_and_statistical_sampler_are_separate': True,
         'single_random_record_strict_monotonicity': False,
         'unrestricted_phase_winding_identified': False,
         'autonomous_native_dynamics_or_laboratory_seconds_derived': False,
         'coherent_source_services_still_required': True,
         'public_empirical_dataset_agreement_established': False}


def source_paths():
    paths = set(training.local.SOURCE_PATHS) | set(learned.common.SOURCE_PATHS)
    paths.update(('model/native_modular_reference/reference_tomography.py',
                  'model/native_modular_reference/learned_clock.py',
                  'model/native_modular_reference/cap_local_clock.py',
                  'model/native_modular_reference/local_history.py',
                  'model/native_modular_reference/phase_readout.py',
                  'model/native_modular_reference/phase_inference.py',
                  'model/native_modular_reference/phase_clock.py',
                  'model/native_equilibration/common_variation.py', 'model/common_source/halfline.py'))
    return tuple(sorted(paths))


def sampler_certificate():
    total = training._record_contract()[0]+SAMPLER_GATE
    return {'exact_conditional_columns_checked': 66,
            'conditional_score_mean_bias_upper': 4*total,
            'conditional_score_second_moment_upper': F(9, 20)*training._record_contract()[1]+8*total,
            'entire_multitrial_path_distance_claimed': False,
            'raw_compiled_quantum_carrier_sampled': False}


def resources(law, inference, rows):
    """Conservative fixed protocol costs, including every failed trial."""
    certificates = [row['source_certificates'][q] for row in rows for q in readout.QUADRATURES]
    alphabets = sum(len(readout.record_alphabet(law, i, prepared=True))+1 for i in (0, 1))
    return {'training_native_attempts_upper': inference['source_native_attempts_all_trials_upper'],
            'readout_native_attempts_upper': 4*TRIALS*readout.complete_trial_budget()['native_attempts_per_trial_upper'],
            'readout_source_pulses_upper': TRIALS*sum(int(c['source_pulses']) for c in certificates),
            'readout_proper_code_service_events_upper': TRIALS*sum(int(c['proper_code_service_events']) for c in certificates),
            'all_retained_record_count_bits_upper': 2*alphabets*23+7*25,
            'all_retained_record_state_count_bits_upper': 2*alphabets*66*23+7*66*25,
            'counter_memory_includes_source_program_private_tapes_stream_hashes_or_arithmetic_workspace': False}


def quantized_kernel(weights, s, index, quadrature, law):
    """Declared CDF quantization, retaining every fine record and final state."""
    records, instrument = readout.population_instrument(weights, s, index=index,
                                                       quadrature=quadrature, law=law)
    maps = training.local.population_maps(law=law)
    core, _ = training.population_counts()
    success = (core/2)@maps['success']
    arrays = [a@success for a in instrument.reshape(len(records), 66, 66)]
    labels = list(records)
    for flag in ('herald_failure', 'warmup_timeout', 'completion_timeout'):
        labels.append(('preparation', flag))
        arrays.append(maps[flag])
    raw = np.vstack(arrays)
    if (not np.isfinite(raw).all() or np.any(raw < 0)
            or float(np.max(np.abs(raw.sum(axis=0)-1))) > 5e-10):
        raise ArithmeticError('complete numerical phase kernel failed its mass gate')
    result = np.zeros((len(raw)+66, 66), dtype=np.uint64)
    for x in range(66):
        earlier = 0
        for y, value in enumerate(np.cumsum(raw[:, x])):
            scaled = F(float(value))*GRID
            threshold = min(GRID, -(-scaled.numerator//scaled.denominator))
            if threshold < earlier:
                raise ArithmeticError('nonmonotone complete CDF')
            result[y, x] = threshold-earlier
            earlier = threshold
        result[len(raw)+x, x] = GRID-earlier
    labels.append(NUMERICAL_ABORT)
    _validate_kernel(tuple(labels), result)
    return tuple(labels), result


def _validate_kernel(records, kernel):
    if (type(records) is not tuple or len(set(records)) != len(records)
            or type(kernel) is not np.ndarray or kernel.dtype != np.uint64
            or kernel.shape != (len(records)*66, 66) or np.any(kernel > GRID)
            or any(sum(map(int, kernel[:, x])) != GRID for x in range(66))):
        raise ValueError('a complete original dyadic record/state kernel is required')
    for record in records:
        readout.score(record)


def audit_kernel(records, kernel, weights, s, quadrature):
    """Exact score moments of all 66 columns against analytic phase enclosures.

    This is a conditional-mean/variance certificate, not a claim that the
    entire sampled path is close to an ideal quantum path. The phase target
    and source reference enter validation only, never inverse.infer_counts.
    """
    _validate_kernel(records, kernel)
    source_error, q, _ = training._record_contract()
    sigma = training.local.reference()
    charges = training.variation.charges()
    w = tuple(next(v for v, c in zip(sigma, charges, strict=True) if c == j) for j in range(3))
    amplitude = q*F(9, 2)*(w[1]-w[0])
    re, im, radius = learned.common.phase_enclosure(s, weights[1]/weights[0])
    center = amplitude*(re if quadrature == 'cos' else -im)
    enclosure = abs(amplitude)*radius
    scores = tuple(readout.score(record) for record in records)
    bias, second = F(), F()
    for x in range(66):
        masses = tuple(sum(map(int, kernel[66*j:66*(j+1), x])) for j in range(len(records)))
        mean = sum((F(k*mass, GRID) for k, mass in zip(scores, masses, strict=True)), F())
        moment = sum((F(k*k*mass, GRID) for k, mass in zip(scores, masses, strict=True)), F())
        bias = max(bias, abs(mean-center)+enclosure)
        second = max(second, moment)
    # These portable bounds include finite preparation and all numerical
    # operations; exact checks above are rerun on every receipt replay.
    total = source_error+SAMPLER_GATE
    mean_cap, second_cap = 4*total, F(9, 20)*q+8*total
    if bias > mean_cap or second > second_cap:
        raise ArithmeticError('dyadic phase sampler violates its exact score certificate')
    return sampler_certificate()


def sample_kernel(records, kernel, trials, *, seed, initial, chunk_size=1 << 16):
    """Sample one continuing carrier; batching preserves the exact dyadic law."""
    _validate_kernel(records, kernel)
    training._integer(trials, 1, 1 << 30, 'trials')
    training._integer(seed, 0, (1 << 63)-1, 'seed')
    training._integer(initial, 0, 65, 'initial')
    training._integer(chunk_size, 1, 1 << 20, 'chunk_size')
    common = kernel.min(axis=1)
    common_mass = int(common.sum())
    common_cdf = np.cumsum(common)
    residual_cdf = np.cumsum(kernel-common[:, None], axis=0)
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(kernel), dtype=np.int64)
    digest, state, residual, prefix = hashlib.sha256(), initial, 0, []
    for start in range(0, trials, chunk_size):
        draws = rng.integers(0, GRID, size=min(chunk_size, trials-start), dtype=np.uint64)
        selected = np.searchsorted(common_cdf, draws, side='right')
        for position in np.flatnonzero(draws >= common_mass):
            prior = state if position == 0 else int(selected[position-1]) % 66
            selected[position] = np.searchsorted(residual_cdf[:, prior],
                                                   int(draws[position])-common_mass, side='right')
            residual += 1
        if np.any(selected >= len(kernel)):
            raise ArithmeticError('record/state draw omitted mass')
        state = int(selected[-1]) % 66
        counts += np.bincount(selected, minlength=len(kernel))
        digest.update(selected.astype('<u4').tobytes())
        if len(prefix) < 12:
            prefix.extend((records[int(i)//66], int(i)%66) for i in selected[:12-len(prefix)])
    state_counts = counts.reshape(len(records), 66)
    return {'records': records,
            'counts': tuple(int(row.sum()) for row in state_counts),
            'terminal_state_counts': tuple((int(i), int(j), int(state_counts[i, j]))
                                           for i, j in np.argwhere(state_counts != 0)),
            'trials': trials, 'seed': seed, 'initial_state': initial, 'final_state': state,
            'record_and_terminal_state_sha256': digest.hexdigest(), 'prefix': prefix,
            'residual_draws': residual, 'same_carrier_reused': True,
            'microscopic_native_attempt_tapes_executed': False,
            'seeded_PRNG_demonstration_proves_independent_uniform_randomness': False}


def record_counts(sample):
    return dict(zip(sample['records'], sample['counts'], strict=True))


def _json(value):
    return json.loads(json.dumps(value, default=str, allow_nan=False))


def _same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def inference_summary(result):
    """Retain derived bounds once; their full input records live in sampling."""
    result = deepcopy(result)
    del result['reference_inference']
    for observed in result['observed_scores'].values():
        del observed['complete_counts']
    return result


def _pretty(value, depth=0):
    """Readable compact arrays avoid thousands of one-number receipt lines."""
    pad = '  '*depth
    if type(value) is dict:
        return '{\n'+',\n'.join('  '*(depth+1)+json.dumps(key)+': '+_pretty(value[key], depth+1)
                               for key in sorted(value))+'\n'+pad+'}'
    if type(value) is list:
        compact = json.dumps(value, separators=(',', ':'), ensure_ascii=False, allow_nan=False)
        if len(compact) <= 240:
            return compact
        return '[\n'+',\n'.join('  '*(depth+1)+_pretty(item, depth+1) for item in value)+'\n'+pad+']'
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def _detach_rationals(value):
    """Copy returned containers without sharing Fraction's writable slots."""
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is dict:
        return {key: _detach_rationals(item) for key, item in value.items()}
    if type(value) in (tuple, list):
        return type(value)(_detach_rationals(item) for item in value)
    return value


def run():
    """Fixed source/training/readout protocol; its target is validation only."""
    proofs = readout.exact_effects()
    budget = readout.complete_trial_budget()
    systematic = training._record_contract()[0]+readout.source_bounds()[0]+SAMPLER_GATE
    alternatives = []
    for law, target in (('C', F(3, 4)), ('E', F(-3, 4))):
        trained = training.execute_trials(law=law)
        inferred = trained['inference']
        if inferred['status'] != 'faithful':
            raise ArithmeticError('retained native training did not identify a faithful reference')
        weights = inferred['central_subnormalized_weights']
        state = trained['execution']['final_state']
        rows = []
        for index in (0, 1):
            samples, histories, certificates, audits = {}, {}, {}, {}
            for qi, quadrature in enumerate(readout.QUADRATURES):
                certificate = readout.certificate(weights, target, index=index, quadrature=quadrature)
                records, kernel = quantized_kernel(weights, target, index, quadrature, law)
                audits[quadrature] = audit_kernel(records, kernel, weights, target, quadrature)
                sample = sample_kernel(records, kernel, TRIALS,
                                       seed=20261025+100*index+qi, initial=state)
                # Full raw source execution is checked from the SAME batch
                # input. The population stream itself is explicitly ideal-
                # stage dyadic simulation, not microscopic pulse execution.
                actual = readout.prepared(weights, law, initial=state, index=index,
                                          quadrature=quadrature, s=target, raw=True)
                ideal = readout.prepared(weights, law, initial=state, index=index,
                                         quadrature=quadrature, s=target, raw=False)
                distance = readout.history.distance(actual, ideal)
                probabilities = readout.history.probabilities(actual)
                mass = sum(probabilities.values())
                if distance > float(certificate['source_history_halfdiamond_upper']):
                    raise ArithmeticError('raw complete phase history missed its source bound')
                if abs(mass-1) > 2e-10:
                    raise ArithmeticError('raw complete phase history lost mass')
                histories[quadrature] = {
                    'raw_vs_ideal_complete_halftrace_diagnostic': distance,
                    'complete_mass_diagnostic': mass,
                    'record_probabilities': tuple((record, probability)
                                                 for record, probability in probabilities.items())}
                samples[quadrature], certificates[quadrature] = sample, certificate
                state = sample['final_state']
            result = inverse.infer_counts(
                {key: record_counts(value) for key, value in samples.items()},
                trained['execution']['counts'], law=law, index=index,
                systematic_halftrace_upper=systematic,
                reference_sampling_bias=training.SAMPLER_GATE)
            if result['status'] != 'identified':
                raise ArithmeticError('retained complete records did not identify the phase')
            interval = result['parameter_interval']
            if not interval[0] <= target <= interval[1] or interval[1]-interval[0] > F(1, 5):
                raise ArithmeticError('retained blind interval missed target or declared precision gate')
            rows.append({'index': index, 'sampling': samples, 'inference': inference_summary(result),
                         'sampler_certificates': audits, 'source_certificates': certificates,
                         'raw_source_histories': histories,
                         'validation_only_target_in_interval': True})
        intersection = (max(row['inference']['parameter_interval'][0] for row in rows),
                        min(row['inference']['parameter_interval'][1] for row in rows))
        if not intersection[0] <= target <= intersection[1]:
            raise ArithmeticError('two common-reference caps do not admit their shared parameter')
        alternatives.append({'law': law, 'validation_only_programmed_parameter': target,
            'training': {'execution': trained['execution'],
                         'sampler_certificate': _detach_rationals(trained['sampler_certificate']),
                         'inference': inferred},
            'readouts': rows, 'common_cap_parameter_intersection': intersection,
            'final_carrier_state': state,
            'resources': resources(law, inferred, rows),
            # One common training event plus four quadrature-mean tails.
            'both_caps_joint_confidence_failure_upper': inferred['confidence_failure_upper']+F(8, 2**12)})
    return {'schema': 'oph.native_modular_reference.phase_clock.v1',
            'source_effect_proof': proofs, 'complete_trial_budget': budget,
            'fixed_trials_per_quadrature': TRIALS,
            'fixed_parameter_window': (F(-1), F(1)), 'fixed_interval_width_gate': F(1, 5),
            'alternatives': alternatives,
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths()},
            'scope': SCOPE.copy()}


def _validate_sample(sample, law, index, expected_initial):
    alphabet = readout.record_alphabet(law, index, prepared=True)+(NUMERICAL_ABORT,)
    records = tuple(tuple(row) for row in sample['records'])
    if records != alphabet or sample['trials'] != TRIALS:
        raise ValueError('complete fixed-budget record alphabet required')
    counts = sample['counts']
    if (len(counts) != len(records) or any(type(x) is not int or x < 0 for x in counts)
            or sum(counts) != TRIALS or type(sample['trials']) is not int):
        raise ValueError('complete original integer counts required')
    for key in ('initial_state', 'final_state'):
        training._integer(sample[key], 0, 65, key)
    if sample['initial_state'] != expected_initial:
        raise ValueError('phase batches did not continue the previous carrier')
    training._integer(sample['seed'], 0, (1 << 63)-1, 'seed')
    training._integer(sample['residual_draws'], 0, TRIALS, 'residual draws')
    if (sample['same_carrier_reused'] is not True
            or sample['microscopic_native_attempt_tapes_executed'] is not False
            or sample['seeded_PRNG_demonstration_proves_independent_uniform_randomness'] is not False
            or type(sample['record_and_terminal_state_sha256']) is not str
            or not re.fullmatch('[0-9a-f]{64}', sample['record_and_terminal_state_sha256'])):
        raise ValueError('record custody or sampler scope changed')
    support = {record: {after for _, after in mapping}
               for record, mapping in readout._combined_branches(law, index)}
    support[('abort', index)] = set()  # Ideal-stage sampler has no export.
    reconstructed, seen, cells = [0]*len(records), set(), {}
    for row in sample['terminal_state_counts']:
        if type(row) is not list or len(row) != 3:
            raise ValueError('sparse complete state counts required')
        i, j, count = row
        training._integer(i, 0, len(records)-1, 'record index')
        training._integer(j, 0, 65, 'terminal index')
        training._integer(count, 1, TRIALS, 'terminal count')
        if records[i] in support and j not in support[records[i]]:
            raise ValueError('terminal state is impossible for its native record')
        if (i, j) in seen:
            raise ValueError('duplicate terminal state cell')
        seen.add((i, j)); reconstructed[i] += count
        cells[(records[i], j)] = count
    if reconstructed != counts:
        raise ValueError('terminal state counts disagree with public records')
    if len(sample['prefix']) != 12:
        raise ValueError('retained record prefix missing')
    for record, state in sample['prefix']:
        if tuple(record) not in alphabet:
            raise ValueError('invalid prefix record')
        training._integer(state, 0, 65, 'prefix state')
        if (records.index(tuple(record)), state) not in seen:
            raise ValueError('prefix state missing from retained state counts')
    if not any(j == sample['final_state'] for _, j in seen):
        raise ValueError('final carrier absent from state counts')
    if any(count > cells.get(key, 0) for key, count in
           Counter((tuple(record), state) for record, state in sample['prefix']).items()):
        raise ValueError('prefix multiplicities exceed retained state counts')
    return dict(zip(records, counts, strict=True))


def _validate_training(trained, law):
    execution = trained['execution']
    n = learned.TRAINING_TRIALS
    training._integer(execution['seed'], 0, (1 << 63)-1, 'training seed')
    training._integer(execution['initial_state'], 0, 65, 'training initial')
    if (type(execution['trials']) is not int or execution['trials'] != n
            or execution['seed'] != 20261024 or execution['initial_state'] != 0
            or execution['same_carrier_reused'] is not True
            or execution['microscopic_attempt_tapes_executed'] is not False
            or execution['seeded_PRNG_demonstration_proves_independent_uniform_randomness'] is not False
            or execution['common_subkernel_probability_lower'] != str(1-training.SAMPLER_GATE)
            or not _same(trained['sampler_certificate'], _json(training._portable_sampler_certificate()))
            or type(execution['record_and_terminal_state_sha256']) is not str
            or not re.fullmatch('[0-9a-f]{64}', execution['record_and_terminal_state_sha256'])):
        raise ValueError('reference training custody or scope changed')
    training._integer(execution['final_state'], 0, 65, 'training final state')
    training._integer(execution['residual_draws'], 0, n, 'training residual draws')
    rows = execution['terminal_state_counts']
    if len(rows) != 7 or any(len(row) != 66 for row in rows):
        raise ValueError('complete reference record/state counts required')
    for label, row in zip(training.LABELS, rows, strict=True):
        if any(type(x) is not int or x < 0 for x in row) or sum(row) != execution['counts'][label]:
            raise ValueError('reference state counts disagree with records')
        if label in training.LABELS[:3]:
            current = training.LABELS.index(label)-1
            support = {after for branch in readout.interface.fixed_maps(law, training.SEAM)
                       if branch.current == current for _, after in branch.mapping}
            if any(count and state not in support for state, count in enumerate(row)):
                raise ValueError('reference terminal state is impossible for its native current')
    if len(execution['prefix']) != 32:
        raise ValueError('reference prefix missing')
    for label, state in execution['prefix']:
        training._integer(state, 0, 65, 'training prefix state')
        if label not in training.LABELS or rows[training.LABELS.index(label)][state] <= 0:
            raise ValueError('reference prefix incompatible with retained counts')
    if not any(row[execution['final_state']] > 0 for row in rows):
        raise ValueError('training final carrier absent from state counts')
    if any(count > rows[training.LABELS.index(label)][state] for (label, state), count in
           Counter(tuple(row) for row in execution['prefix']).items()):
        raise ValueError('training prefix multiplicities exceed retained counts')


def _validated_copy(original):
    value = deepcopy(original)
    if (value['schema'] != 'oph.native_modular_reference.phase_clock.v1'
            or not _same(value['scope'], SCOPE) or type(value['fixed_trials_per_quadrature']) is not int
            or value['fixed_trials_per_quadrature'] != TRIALS
            or value['fixed_parameter_window'] != ['-1', '1']
            or value['fixed_interval_width_gate'] != '1/5'
            or not _same(value['source_effect_proof'], _json(readout.exact_effects()))
            or not _same(value['complete_trial_budget'], _json(readout.complete_trial_budget()))):
        raise ValueError('canonical phase proof or scope changed')
    bindings = value['source_sha256']
    if (set(bindings) != set(source_paths()) or any(
                type(p) is not str or Path(p).is_absolute() or '..' in Path(p).parts
            or type(digest) is not str or hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest
            for p, digest in bindings.items())):
        raise ValueError('phase source custody is incomplete or stale')
    if len(value['alternatives']) != 2:
        raise ValueError('both repair-law alternatives required')
    systematic = training._record_contract()[0]+readout.source_bounds()[0]+SAMPLER_GATE
    for alternative, law, target in zip(value['alternatives'], ('C', 'E'), ('3/4', '-3/4'), strict=True):
        if alternative['law'] != law or alternative['validation_only_programmed_parameter'] != target:
            raise ValueError('declared validation alternative changed')
        trained = alternative['training']
        _validate_training(trained, law)
        fresh_training = training.infer_counts(trained['execution']['counts'], sampling_bias=training.SAMPLER_GATE)
        if not _same(trained['inference'], _json(fresh_training)):
            raise ValueError('reference inference disagrees with retained training records')
        if trained['execution']['trials'] != learned.TRAINING_TRIALS:
            raise ValueError('fixed reference training budget changed')
        state = trained['execution']['final_state']
        if len(alternative['readouts']) != 2:
            raise ValueError('both common-reference caps required')
        for row, index in zip(alternative['readouts'], (0, 1), strict=True):
            if type(row['index']) is not int or row['index'] != index or row['validation_only_target_in_interval'] is not True:
                raise ValueError('cap identity or validation changed')
            counts = {}
            for qi, quadrature in enumerate(readout.QUADRATURES):
                sample = row['sampling'][quadrature]
                counts[quadrature] = _validate_sample(sample, law, index, state)
                batch_initial = state
                if sample['seed'] != 20261025+100*index+qi:
                    raise ValueError('fixed phase seed changed')
                state = sample['final_state']
                cert = row['source_certificates'][quadrature]
                if not _same(row['sampler_certificates'][quadrature], _json(sampler_certificate())):
                    raise ValueError('exact conditional sampler certificate changed')
                if not _same(cert, _json(readout.certificate(fresh_training['central_subnormalized_weights'],
                                                    F(target), index=index, quadrature=quadrature))):
                    raise ValueError('raw source certificate changed')
                probe = row['raw_source_histories'][quadrature]
                mass, distance = probe['complete_mass_diagnostic'], probe['raw_vs_ideal_complete_halftrace_diagnostic']
                if (type(mass) is not float or not math.isfinite(mass) or abs(mass-1) > 2e-10
                        or type(distance) is not float or not math.isfinite(distance) or distance < 0
                        or distance > float(F(cert['source_history_halfdiamond_upper']))):
                    raise ValueError('raw phase history failed its diagnostic gates')
                probabilities = probe['record_probabilities']
                expected = readout.record_alphabet(law, index, prepared=True)
                if (len(probabilities) != len(expected) or set(tuple(r) for r, _ in probabilities) != set(expected)
                        or any(type(p) is not float or not math.isfinite(p) or p < 0 or p > 1+2e-10
                               for _, p in probabilities)
                        or abs(sum(p for _, p in probabilities)-mass) > 2e-10):
                    raise ValueError('raw phase probabilities are incomplete')
                ideal_records, ideal_kernel = readout.population_instrument(
                    fresh_training['central_subnormalized_weights'], F(target), index=index,
                    quadrature=quadrature, law=law)
                maps = training.local.population_maps(law=law)
                prepared = (training.population_counts()[0]/2)@maps['success'][:, batch_initial]
                target_probabilities = dict(zip(ideal_records,
                    (ideal_kernel@prepared).reshape(len(ideal_records), 66).sum(axis=1), strict=True))
                target_probabilities.update({('preparation', flag): maps[flag][:, batch_initial].sum()
                                             for flag in readout.PREPARATION_FLAGS})
                if sum(abs(p-target_probabilities[tuple(record)]) for record, p in probabilities)/2 > (
                        float(F(cert['source_history_halfdiamond_upper']))+5e-12):
                    raise ValueError('retained raw records violate the certified ideal-source instrument')
                probe['complete_mass_diagnostic'] = 'checked'
                probe['raw_vs_ideal_complete_halftrace_diagnostic'] = 'checked'
                probe['record_probabilities'] = sorted((r, 'checked') for r, _ in probabilities)
            fresh = inverse.infer_counts(counts, trained['execution']['counts'], law=law, index=index,
                                          systematic_halftrace_upper=systematic,
                                          reference_sampling_bias=training.SAMPLER_GATE)
            if not _same(row['inference'], _json(inference_summary(fresh))) or fresh['status'] != 'identified':
                raise ValueError('phase inverse does not follow complete public records')
            lo, hi = fresh['parameter_interval']
            if not lo <= F(target) <= hi or hi-lo > F(1, 5):
                raise ValueError('declared phase interval gate failed')
        training._integer(alternative['final_carrier_state'], 0, 65, 'alternative final carrier')
        if alternative['final_carrier_state'] != state:
            raise ValueError('final continued carrier mismatch')
        intersection = (max(F(row['inference']['parameter_interval'][0]) for row in alternative['readouts']),
                        min(F(row['inference']['parameter_interval'][1]) for row in alternative['readouts']))
        if (alternative['common_cap_parameter_intersection'] != _json(intersection)
                or not _same(alternative['resources'], _json(resources(law, fresh_training, alternative['readouts'])))
                or alternative['both_caps_joint_confidence_failure_upper'] != str(
                    fresh_training['confidence_failure_upper']+F(8, 2**12))):
            raise ValueError('common-cap confidence or complete resource budget changed')
    return value


def _compare_receipts(supplied, actual):
    """Compare with a producer result owned by this module, after validation."""
    try:
        left, right = _validated_copy(supplied), _validated_copy(actual)
        for a, b in zip(supplied['alternatives'], actual['alternatives'], strict=True):
            for ar, br in zip(a['readouts'], b['readouts'], strict=True):
                for quadrature in readout.QUADRATURES:
                    for field in ('complete_mass_diagnostic', 'raw_vs_ideal_complete_halftrace_diagnostic'):
                        if abs(ar['raw_source_histories'][quadrature][field]
                               -br['raw_source_histories'][quadrature][field]) > 2e-10:
                            raise ValueError('raw source history diagnostic drifted')
                    ap = dict((tuple(record), probability) for record, probability in
                              ar['raw_source_histories'][quadrature]['record_probabilities'])
                    bp = dict((tuple(record), probability) for record, probability in
                              br['raw_source_histories'][quadrature]['record_probabilities'])
                    if any(abs(ap[key]-bp[key]) > 2e-10 for key in ap):
                        raise ValueError('raw complete-record probabilities drifted')
        if json.dumps(left, sort_keys=True, allow_nan=False) != json.dumps(right, sort_keys=True, allow_nan=False):
            raise ValueError('phase receipt drift')
    except (KeyError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('malformed phase receipt') from exc


def compare_receipt(supplied, actual=None):
    """Verify against an owned canonical replay, never a caller's authority.

    The optional second receipt is retained for callers of the older API;
    it must also match the fresh replay and cannot replace that replay.
    """
    fresh = _json(run())
    _compare_receipts(supplied, fresh)
    if actual is not None:
        _compare_receipts(actual, fresh)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', type=Path)
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    payload = _json(run())
    if args.check:
        _compare_receipts(strict_json.load(args.check), payload)
        print('Complete-record native phase inference verified.')
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(_pretty(payload)+'\n')


if __name__ == '__main__':
    main()
