"""Infer an unread passive repair age from current-derived endpoint fields.

The controller discovers a configuration through committed currents, chooses
an antipodal pair containing a witnessed occupied port, lets the same carrier
undergo an unread interval, and discovers its endpoint. No reset, stationary
preparation, coherent operation, supplied weight, or elapsed-count input enters
the decoder. The interval has a fixed unknown integer length across trials.
Its scheduler and stopping gate are outside the endpoint observer's access.
Native auxiliary coarse coins implement bounded uniform seam selection.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import random
import re

from tools import strict_json
from . import configuration as discovery
from ..repair_scheduler import coins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('passive_age_results.json')
TRIALS, WINDOW, FAILURE_BITS, SELECTOR_ROUNDS = 1 << 16, 4, 12, 8
REASONS = ('seed_timeout', 'flip_timeout', 'probe_timeout')
SCHEMA = 'oph.native_clock.passive_age.v1'
SCOPE = {
    'decoder_receives_interval_attempt_count': False,
    'decoder_reads_passive_scheduler_records': False,
    'scheduler_records_and_gate_are_private_producer_inputs': True,
    'record_hashes_authenticate_physical_execution': False,
    'stationary_preparation_reset_or_independent_copies_required': False,
    'coherent_source_services_required': False,
    'auxiliary_reset_or_independent_auxiliary_state_required': False,
    'auxiliary_current_and_state_in_observer_inputs': False,
    'fixed_unknown_integer_window_length_required': True,
    'uniform_fair_native_source_contract_required': True,
    'external_start_end_boundary_and_record_memory_required': True,
    'autonomous_time_generator_or_physical_seconds_claimed': False,
    'seeded_execution_proves_physical_randomness': False,
}


def _integer(value, lower, upper, name):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be an original integer in [{lower},{upper}]')
    return value


@lru_cache(maxsize=1)
def pairs():
    seams = discovery.quantum.seams()
    adjacency = {i: {b if a == i else a for a, b in seams if i in (a, b)} for i in range(12)}
    result = []
    for a in range(12):
        opposite = [b for b in range(12) if b != a and b not in adjacency[a]
                    and not adjacency[a].intersection(adjacency[b])]
        if len(opposite) != 1:
            raise ArithmeticError('source antipodal pairing changed')
        b, = opposite
        if a < b:
            result.append((a, b))
    if len(result) != 6 or set(sum(result, ())) != set(range(12)):
        raise ArithmeticError('source antipodes do not partition the ports')
    return tuple(result)


def source_paths():
    return tuple(sorted(set((*discovery.SOURCE_PATHS,
        'model/repair_scheduler/coins.py', 'model/native_clock/passive_age.py'))))


def exact_identity():
    """Actual native Heisenberg identity, including every matrix unit."""
    states, seams = discovery.native.states(binary=True), discovery.quantum.seams()
    rows = []
    for law in ('C', 'E'):
        for pair in pairs():
            expectation = [F()]*66
            for seam in seams:
                for _, _, mapping in coins.coin_maps(law, seam):
                    if len({j for _, j in mapping}) != len(mapping):
                        raise ArithmeticError('off-diagonal native effect is not zero')
                    for i, j in mapping:
                        expectation[i] += F(sum(states[j][a] for a in pair), 60)
            target = tuple(F(1, 3)+F(9, 10)*(sum(x[a] for a in pair)-F(1, 3)) for x in states)
            if tuple(expectation) != target:
                raise ArithmeticError('native antipodal occupation is not the declared eigenfield')
            rows.append({'law': law, 'pair': pair, 'diagonal_columns': 66,
                         'complete_operator_matrix_units': 4356})
    return {'rows': rows, 'eigenvalue': '9/10', 'equilibrium_offset': '1/3',
            'identity': 'P*(Y-1/3)=(9/10)*(Y-1/3)',
            'untouched_reference_extension': True,
            'fine_conclusive_discovery_words_are_rank_one': True,
            'discovery_source_audit': discovery.source_audit(),
            'native_auxiliary_coin_effect': coins.effect_certificate()}


@lru_cache(maxsize=1)
def alphabet():
    return tuple([('start_failure', why) for why in REASONS]
        +[('scheduler_timeout', p, y) for p in range(6) for y in (1, 2)]
        +[('read_failure', p, y, why) for p in range(6) for y in (1, 2) for why in REASONS]
        +[('read', p, y, z) for p in range(6) for y in (1, 2) for z in (0, 1, 2)])


def score(record):
    if (type(record) is not tuple or any(type(x) not in (str, int) for x in record)
            or record not in alphabet()):
        raise ValueError('an original complete passive-age record is required')
    return F(3*record[3]-1, 3*record[2]-1) if record[0] == 'read' else F()


def choose_pair(outcome):
    if type(outcome) is not discovery.Outcome or not outcome.success:
        raise ValueError('a conclusive current-only discovery is required')
    state = discovery._configuration(outcome.final, binary=True)
    occupied = next(i for i, value in enumerate(state) if value)
    pair_index = next(i for i, pair in enumerate(pairs()) if occupied in pair)
    return pair_index, sum(state[i] for i in pairs()[pair_index])


def decode_current_tape(tape, *, full=False, seed_scans=4, retries=32):
    """Replay only public seam/current replies; no source state or age input."""
    if type(tape) is not tuple:
        raise ValueError('an original current tape is required')
    controller = discovery.policy(seed_scans, retries, full=full)
    request, offset = next(controller), 0
    while True:
        currents = []
        for _ in range(request.limit):
            if offset == len(tape):
                raise ValueError('current tape ended before its complete discovery outcome')
            entry = tape[offset]
            if (type(entry) is not tuple or len(entry) != 2 or entry[0] != request.seam
                    or type(entry[0]) is not tuple or any(type(v) is not int for v in entry[0])
                    or type(entry[1]) is not int or entry[1] not in (-1, 0, 1)):
                raise ValueError('current tape disagrees with the record-only policy')
            currents.append(entry[1]); offset += 1
            if entry[1]:
                break
        try:
            request = controller.send(tuple(currents))
        except StopIteration as end:
            if offset != len(tape):
                raise ValueError('current tape extends beyond its terminal outcome')
            return end.value


def endpoint_record(start, end=None, *, scheduler_timeout=False):
    """Public policy outcomes only; their labels are reconstructed from currents."""
    if type(start) is not discovery.Outcome or type(scheduler_timeout) is not bool:
        raise ValueError('original current-only policy outcomes are required')
    if (not start.success and (end is not None or scheduler_timeout)) or (scheduler_timeout and end is not None):
        raise ValueError('a failed start or selector cannot carry an endpoint read')
    if not start.success:
        record = ('start_failure', start.reason)
    else:
        p, y = choose_pair(start)
        if scheduler_timeout:
            record = ('scheduler_timeout', p, y)
        elif type(end) is not discovery.Outcome:
            raise ValueError('complete endpoint discovery is required')
        elif not end.success:
            record = ('read_failure', p, y, end.reason)
        else:
            z = sum(end.initial[i] for i in pairs()[p])
            record = ('read', p, y, z)
    score(record)
    return record


def budget(trials=TRIALS, window=WINDOW, *, seed_scans=4, retries=32, selector_rounds=SELECTOR_ROUNDS):
    _integer(trials, 1, 1 << 30, 'trials')
    _integer(window, 1, 128, 'window')
    read = discovery.budget(seed_scans, retries)
    selector = coins.selector_certificate(30, selector_rounds)
    failure, timeout = read['failure_upper'], F(selector['timeout_probability'])
    # Start failure loses at most theta<=1. Endpoint failure can delete a
    # score of magnitude <=5/2; each scheduler timeout has the same bound.
    bias = F(7, 2)*failure+F(5, 2)*window*timeout
    main = trials*(2*read['native_attempts_max']+window)+30
    auxiliary = trials*window*selector['maximum_auxiliary_native_attempts']
    return {'trials': trials, 'window': window, 'seed_scans': seed_scans, 'retries': retries,
            'selector_rounds': selector_rounds, 'read_failure_upper': failure,
            'selector_timeout_per_attempt': timeout, 'conditional_score_bias_upper': bias,
            'score_range': (F(-1, 2), F(5, 2)),
            'main_native_attempts_upper': main, 'auxiliary_native_attempts_upper': auxiliary,
            'total_native_attempts_upper': main+auxiliary,
            'discovery_calls_upper': 2*trials, 'passive_main_attempts_upper': trials*window,
            'coherent_source_pulses': 0, 'auxiliary_reset_calls': 0,
            'record_and_state_counter_bits_upper': len(alphabet())*66*trials.bit_length(),
            'counter_storage_is_complete_hardware_memory': False}


def radius(trials, failure_bits=FAILURE_BITS):
    _integer(trials, 1, 1 << 30, 'trials')
    _integer(failure_bits, 1, 64, 'failure bits')
    square = (trials*(failure_bits+1)+1)//2
    t = isqrt(square)
    t += int(t*t < square)
    result = F(3*t, trials)
    if 2*trials*result**2 < 9*(failure_bits+1):
        raise ArithmeticError('rational martingale exponent failed')
    return result


def admission(trials=TRIALS, window=WINDOW, *, failure_bits=FAILURE_BITS, **service):
    limits = budget(trials, window, **service)
    total = radius(trials, failure_bits)+limits['conditional_score_bias_upper']
    separation = F(1, 10)*F(9, 10)**(window-1)
    return {'budget': limits, 'statistical_radius_upper': radius(trials, failure_bits),
            'total_mean_radius_upper': total, 'minimum_mean_separation': separation,
            'all_window_ages_uniquely_decodable_on_coverage': 2*total < separation,
            'coverage_failure_upper': F(1, 2**failure_bits)}


def infer_counts(counts, *, window=WINDOW, failure_bits=FAILURE_BITS, **service):
    if (type(counts) is not dict or set(counts) != set(alphabet())
            or any(type(v) is not int or v < 0 for v in counts.values())):
        raise ValueError('all original integer record counts including zero failures are required')
    for record in counts:
        score(record)
    n = sum(counts.values())
    contract = admission(n, window, failure_bits=failure_bits, **service)
    mean = sum((score(record)*count for record, count in counts.items()), F())/n
    error = contract['total_mean_radius_upper']
    candidates = tuple(age for age in range(window+1) if abs(mean-F(9, 10)**age) <= error)
    return {'status': 'identified' if len(candidates) == 1 else 'unresolved' if candidates else 'incompatible',
            'candidate_attempt_counts': candidates, 'observed_mean': mean,
            'mean_interval': (mean-error, mean+error), 'admission': contract,
            'all_trials_including_failures': n, 'successful_records_conditioned_on': False,
            'age_or_hidden_configuration_input': False,
            'coverage_scope': 'fixed unknown integer interval; conditional fair native source; no physical duration'}


def random_window_obstruction():
    probability = F(9, 19)
    mixture = probability+(1-probability)*F(9, 10)**2
    if mixture != F(9, 10):
        raise ArithmeticError('random-window obstruction changed')
    return {'ages': (0, 2), 'probabilities': (probability, 1-probability),
            'mean_score': mixture, 'same_mean_as_fixed_age': 1,
            'mean_attempt_count': F(20, 19),
            'equality_of_complete_record_distributions_claimed': False,
            'unknown_random_or_adaptively_stopped_windows_identified': False}


def _coin_stream(rng):
    while True:
        yield rng.getrandbits(1)


def _selector(auxiliary, rng, law, rounds, digest):
    """Execute each auxiliary native repair; never replace it by randrange."""
    seam = discovery.quantum.seams()[0]
    calls = 0
    for attempt in range(rounds):
        word = []
        for _ in range(5):
            coin = rng.getrandbits(1)
            auxiliary, _ = discovery.native.native(auxiliary, seam, coin)
            calls += 1
            word.append(coin)
        digest.update(bytes((3, *word)))  # coarse coins only; hidden auxiliary labels are unread
        address = coins.decode_coin_word(tuple(word))
        if address is not None:
            return address, auxiliary, calls, attempt
    return None, auxiliary, calls, rounds


def execute_trials(age, *, trials=TRIALS, law='C', seed=20261026, seed_scans=4,
                   retries=32, selector_rounds=SELECTOR_ROUNDS, initial_state=None, auxiliary_state=None):
    """Actual native source simulation. Age belongs to the opaque producer only."""
    _integer(age, 0, 128, 'producer interval length')
    _integer(trials, 1, 1 << 30, 'trials')
    _integer(seed, 0, (1 << 63)-1, 'seed')
    coins._law(law)
    discovery._parameters(seed_scans, retries, False, 'measure', None)
    coins.selector_certificate(30, selector_rounds)
    state = ((2,)+(0,)*11 if initial_state is None
             else discovery._configuration(initial_state, binary=False))
    auxiliary = (discovery.native.states(binary=True)[-1] if auxiliary_state is None
                 else discovery._configuration(auxiliary_state, binary=False))
    initial, initial_auxiliary = state, auxiliary
    states = discovery.native.states(binary=True)
    state_index = {x: i for i, x in enumerate(states)}
    seams = discovery.quantum.seams()
    seam_index = {seam: i for i, seam in enumerate(seams)}
    records = alphabet(); record_index = {record: i for i, record in enumerate(records)}
    counts, cells = Counter(dict.fromkeys(records, 0)), Counter()
    rng, digest = random.Random(seed), hashlib.sha256()
    main_calls = auxiliary_calls = discovery_calls = passive_calls = rejections = 0
    discovery_attempts = 0
    prefix = []
    def retain(execution, kind):
        digest.update(bytes((kind,)))
        digest.update(bytes(v for seam, coin, label in execution.records
                            for v in (seam_index[seam], coin, label+2)))
        digest.update(b'\xff')
    for trial in range(trials):
        start = discovery.execute(state, _coin_stream(rng), seed_scans, retries,
                                  full=trial == 0, law=law)
        state = start.output
        discovery_calls += 1; main_calls += start.outcome.attempts
        discovery_attempts += start.outcome.attempts
        retain(start, 1)
        end, timed_out = None, False
        if start.outcome.success:
            for _ in range(age):
                address, auxiliary, calls, rejected = _selector(auxiliary, rng, law, selector_rounds, digest)
                auxiliary_calls += calls; rejections += rejected
                if address is None:
                    timed_out = True
                    break
                seam, coin = seams[address], rng.getrandbits(1)
                output, (_, current) = discovery.native.native(state, seam, coin)
                label = current if law == 'C' else state[seam[0]]-state[seam[1]]
                state = output
                digest.update(bytes((4, address, coin, label+2)))
                main_calls += 1; passive_calls += 1
            if not timed_out:
                end = discovery.execute(state, _coin_stream(rng), seed_scans, retries, law=law)
                state = end.output
                discovery_calls += 1; main_calls += end.outcome.attempts
                discovery_attempts += end.outcome.attempts
                retain(end, 2)
        record = endpoint_record(start.outcome, None if end is None else end.outcome,
                                 scheduler_timeout=timed_out)
        counts[record] += 1
        cells[record_index[record], state_index[state]] += 1
        digest.update(bytes((5, record_index[record], state_index[state])))
        if trial < 2:
            prefix.append({'record': record, 'start_full': trial == 0,
                'start_currents': start.outcome.currents,
                'end_currents': () if end is None else end.outcome.currents,
                'scheduler_timeout': timed_out, 'terminal_state': state})
    return {'trials': trials, 'seed': seed, 'law': law,
            'records_and_counts': tuple((r, counts[r]) for r in records),
            'terminal_state_counts': tuple((*cell, count) for cell, count in sorted(cells.items())),
            'initial_state': initial, 'final_state': state,
            'initial_auxiliary_state': initial_auxiliary, 'final_auxiliary_state': auxiliary,
            'main_native_attempts': main_calls, 'auxiliary_native_attempts': auxiliary_calls,
            'passive_main_attempts': passive_calls, 'discovery_calls': discovery_calls,
            'discovery_native_attempts': discovery_attempts,
            'selector_rejections': rejections, 'source_record_sha256': digest.hexdigest(),
            'prefix': prefix, 'microscopic_native_calls_executed': True,
            'same_main_and_auxiliary_carriers_reused': True}


def _json(value):
    return json.loads(json.dumps(value, default=str, allow_nan=False))


def _same(a, b):
    return json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False)


def run():
    proof = exact_identity()
    alternatives = []
    for law, age in (('C', 3), ('E', 1)):
        execution = execute_trials(age, law=law)
        inferred = infer_counts(dict(execution['records_and_counts']))
        if inferred['candidate_attempt_counts'] != (age,):
            raise ArithmeticError('blind passive-age record inference failed its fixed gate')
        alternatives.append({'law': law, 'validation_only_interval_attempts': age,
                             'execution': execution, 'inference': inferred})
    return {'schema': SCHEMA, 'scope': dict(SCOPE), 'source_proof': proof,
            'executed_admission': admission(),
            'larger_window_analytic_admission': admission(1 << 18, 12),
            'random_window_obstruction': random_window_obstruction(),
            'alternatives': alternatives,
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths()})}


def _tape(value):
    if type(value) is not list:
        raise ValueError('current tapes must be original arrays')
    return tuple((tuple(seam), current) for seam, current in value)


@lru_cache(maxsize=256)
def _reachable(initial, age):
    values = {initial}
    for _ in range(age):
        values = {discovery.native.native(state, seam, coin)[0]
                  for state in values for seam in discovery.quantum.seams() for coin in (0, 1)}
    return frozenset(values)


@lru_cache(maxsize=16)
def _terminal_support(age):
    """Source support only: successes include every finite discovery seed path."""
    paths = discovery._discovery_paths()
    final_sets = {initial: frozenset(row[2] for row in rows) for initial, rows in paths.items()}
    starting = set().union(*final_sets.values())
    states = discovery.native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    result = []
    for record in alphabet():
        if record[0] == 'start_failure':
            possible = set(states)
        else:
            _, pair_index, y, *tail = record
            pair = pairs()[pair_index]
            selected = {state for state in starting if sum(state[i] for i in pair) == y
                        and min(i for i, value in enumerate(state) if value) in pair}
            if record[0] == 'scheduler_timeout':
                possible = set().union(*(_reachable(state, age-1) for state in selected)) if age else set()
            elif record[0] == 'read_failure':
                possible = set(states) if selected else set()
            else:
                z, = tail
                endpoints = set().union(*(_reachable(state, age) for state in selected)) if selected else set()
                possible = set().union(*(final_sets[state] for state in endpoints
                                          if sum(state[i] for i in pair) == z))
        result.append(frozenset(index[state] for state in possible))
    return tuple(result)


def _after_currents(initial, tape):
    state = initial
    maps = discovery._current_maps()
    for seam, current in tape:
        try:
            state = maps[seam, current][state]
        except KeyError as exc:
            raise ValueError('current prefix cannot arise from the continued input state') from exc
    return state


def _validate(receipt):
    if type(receipt) is not dict or set(receipt) != {'schema', 'scope', 'source_proof',
        'executed_admission', 'larger_window_analytic_admission', 'random_window_obstruction',
        'alternatives', 'source_sha256'}:
        raise ValueError('canonical passive-age receipt schema required')
    expected = {'schema': SCHEMA, 'scope': SCOPE, 'source_proof': exact_identity(),
        'executed_admission': admission(), 'larger_window_analytic_admission': admission(1 << 18, 12),
        'random_window_obstruction': random_window_obstruction(),
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths()})}
    if any(not _same(receipt[key], _json(value)) for key, value in expected.items()):
        raise ValueError('passive-age source, proof, scope or resource contract changed')
    if type(receipt['alternatives']) is not list or len(receipt['alternatives']) != 2:
        raise ValueError('both source-law alternatives required')
    limits = budget()
    for row, law, target in zip(receipt['alternatives'], ('C', 'E'), (3, 1), strict=True):
        if type(row) is not dict or set(row) != {'law', 'validation_only_interval_attempts', 'execution', 'inference'}:
            raise ValueError('canonical alternative fields required')
        if row['law'] != law or type(row['validation_only_interval_attempts']) is not int or row['validation_only_interval_attempts'] != target:
            raise ValueError('fixed validation alternative changed')
        sample = row['execution']
        if type(sample) is not dict or set(sample) != {
            'trials', 'seed', 'law', 'records_and_counts', 'terminal_state_counts',
            'initial_state', 'final_state', 'initial_auxiliary_state', 'final_auxiliary_state',
            'main_native_attempts', 'auxiliary_native_attempts', 'passive_main_attempts',
            'discovery_calls', 'discovery_native_attempts', 'selector_rejections',
            'source_record_sha256', 'prefix', 'microscopic_native_calls_executed',
            'same_main_and_auxiliary_carriers_reused'}:
            raise ValueError('canonical continuing execution fields required')
        if (sample['trials'] != TRIALS or type(sample['trials']) is not int
                or sample['seed'] != 20261026 or type(sample['seed']) is not int or sample['law'] != law
                or sample['microscopic_native_calls_executed'] is not True
                or sample['same_main_and_auxiliary_carriers_reused'] is not True
                or type(sample['source_record_sha256']) is not str
                or not re.fullmatch('[0-9a-f]{64}', sample['source_record_sha256'])):
            raise ValueError('passive-age execution custody changed')
        entries = sample['records_and_counts']
        records = tuple(tuple(r) for r, _ in entries)
        if records != alphabet():
            raise ValueError('complete ordered record counts required')
        counts = dict((tuple(r), count) for r, count in entries)
        fresh = infer_counts(counts)
        if sum(counts.values()) != TRIALS or not _same(row['inference'], _json(fresh)) or fresh['candidate_attempt_counts'] != (target,):
            raise ValueError('blind inverse does not follow complete records')
        for name, full in (('initial_state', True), ('final_state', False),
                           ('initial_auxiliary_state', False), ('final_auxiliary_state', False)):
            discovery._configuration(tuple(sample[name]), binary=not full)
        if sample['initial_state'] != [2]+[0]*11 or sample['initial_auxiliary_state'] != list(discovery.native.states(binary=True)[-1]):
            raise ValueError('initial continuing carriers changed')
        reconstructed, seen = [0]*len(records), set()
        supports = _terminal_support(target)
        for i, j, count in sample['terminal_state_counts']:
            _integer(i, 0, len(records)-1, 'record index'); _integer(j, 0, 65, 'state index')
            _integer(count, 1, TRIALS, 'state count')
            if (i, j) in seen:
                raise ValueError('duplicate terminal state cell')
            if j not in supports[i]:
                raise ValueError('terminal state is impossible for the native readout record')
            seen.add((i, j)); reconstructed[i] += count
        if reconstructed != [count for _, count in entries]:
            raise ValueError('terminal states disagree with public records')
        final_index = discovery.native.states(binary=True).index(tuple(sample['final_state']))
        if not any(j == final_index for _, j in seen):
            raise ValueError('final continuing state absent from retained cells')
        maxima = {'main_native_attempts': limits['main_native_attempts_upper'],
                  'auxiliary_native_attempts': limits['auxiliary_native_attempts_upper'],
                  'passive_main_attempts': target*TRIALS, 'discovery_calls': 2*TRIALS,
                  'selector_rejections': target*TRIALS*SELECTOR_ROUNDS,
                  'discovery_native_attempts': 2*TRIALS*discovery.budget()['native_attempts_max']+30}
        for key, maximum in maxima.items():
            _integer(sample[key], 0, maximum, key)
        start_failures = sum(n for r, n in counts.items() if r[0] == 'start_failure')
        timeouts = sum(n for r, n in counts.items() if r[0] == 'scheduler_timeout')
        complete = TRIALS-start_failures-timeouts
        passive, rejected = sample['passive_main_attempts'], sample['selector_rejections']
        if (sample['discovery_calls'] != TRIALS+complete
                or sample['main_native_attempts'] != sample['discovery_native_attempts']+passive
                or sample['auxiliary_native_attempts'] != 5*(passive+rejected)
                or not target*complete <= passive <= target*complete+(target-1)*timeouts
                or not SELECTOR_ROUNDS*timeouts <= rejected <= (SELECTOR_ROUNDS-1)*passive+SELECTOR_ROUNDS*timeouts
                or sample['discovery_native_attempts'] < sample['discovery_calls']+30):
            raise ValueError('complete native source call identities failed')
        if len(sample['prefix']) != 2:
            raise ValueError('current-only prefix missing')
        prefix_initial = tuple(sample['initial_state'])
        for index, prefix in enumerate(sample['prefix']):
            if type(prefix) is not dict or set(prefix) != {'record', 'start_full', 'start_currents',
                    'end_currents', 'scheduler_timeout', 'terminal_state'}:
                raise ValueError('canonical prefix fields required')
            if type(prefix['start_full']) is not bool or prefix['start_full'] != (index == 0):
                raise ValueError('initial full-sector read changed')
            start_tape, end_tape = _tape(prefix['start_currents']), _tape(prefix['end_currents'])
            start = decode_current_tape(start_tape, full=prefix['start_full'])
            after_start = _after_currents(prefix_initial, start_tape)
            if start.success and (start.initial != prefix_initial or start.final != after_start):
                raise ValueError('prefix discovery lost its actual input or final state')
            end = decode_current_tape(end_tape) if end_tape else None
            record = endpoint_record(start, end, scheduler_timeout=prefix['scheduler_timeout'])
            if list(record) != prefix['record'] or counts[record] == 0:
                raise ValueError('prefix is not current-derived')
            terminal = discovery._configuration(tuple(prefix['terminal_state']), binary=True)
            if not start.success:
                possible_terminal = {after_start}
            elif prefix['scheduler_timeout']:
                possible_terminal = _reachable(after_start, target-1) if target else set()
            else:
                endpoints = _reachable(after_start, target)
                if end is not None and end.success and end.initial not in endpoints:
                    raise ValueError('endpoint cannot follow the declared unread native interval')
                possible_terminal = set()
                for endpoint in endpoints:
                    try:
                        possible_terminal.add(_after_currents(endpoint, end_tape))
                    except ValueError:
                        pass
            if terminal not in possible_terminal or (end is not None and end.success and terminal != end.final):
                raise ValueError('prefix lost its disturbed final state')
            prefix_initial = terminal
            if (records.index(record), discovery.native.states(binary=True).index(terminal)) not in seen:
                raise ValueError('prefix terminal cell absent')


def _compare_replayed_receipt(supplied, actual):
    """Internal comparison against a replay owned by this module's caller."""
    try:
        _validate(supplied); _validate(actual)
        if not _same(supplied, actual):
            raise ValueError('fresh passive-age source execution drifted')
    except (KeyError, TypeError, IndexError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('malformed passive-age receipt') from exc


def compare_receipt(supplied, actual=None):
    """Replay the source; caller-supplied matching receipts are not authority."""
    canonical = _json(run())
    _compare_replayed_receipt(supplied, canonical)
    if actual is not None:
        _compare_replayed_receipt(actual, canonical)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    value = _json(run())
    if args.check:
        _compare_replayed_receipt(strict_json.load(args.check), value)
        print('Current-only passive native age inference verified.')
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')


if __name__ == '__main__':
    main()
