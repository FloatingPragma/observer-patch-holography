"""Measure a native execution against a host clock, then replay its raw records.

This is a software apparatus attachment, not a derivation of physical proper
or modular time. Timing acceptance is a declared engineering tolerance.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import random
import re
import time

from . import comparison
from ..observer import native
from ..common_source import slow_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
RAW_PATH = Path(__file__).with_name('host_clock_records.json')
BLOCKS = 4
MOVEMENTS = 10000
ATTEMPT_CAP = 500000
SEED = 1729
CALIBRATION_PORT = 2
VALIDATION_PORT = 3
TIMING_TOLERANCE = F(1, 4)
RATE_ERROR = F(1, 12)
SOURCE_PATHS = tuple(sorted(set((*comparison.SOURCE_PATHS,
    'model/native_clock/host_clock.py', 'model/native_clock/comparison.py'))))
POLICY = {
    'blocks': BLOCKS, 'movements_per_block': MOVEMENTS, 'attempt_cap_per_run': ATTEMPT_CAP,
    'seed': SEED, 'calibration_port': CALIBRATION_PORT, 'validation_port': VALIDATION_PORT,
    'relative_host_duration_tolerance': str(TIMING_TOLERANCE),
    'native_rate_error': str(RATE_ERROR), 'control_host_work_iterations': 64,
    'calibration_block': 0, 'validation_blocks': [1, 2, 3],
}


def _strict_keys(value, keys, name):
    if type(value) is not dict or set(value) != set(keys):
        raise ValueError(f'{name} has missing or unrecognized fields')


def _int(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f'{name} must be an integer in [{low},{high}]')
    return value


def _hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}


def _timed_run(name, work, clock=time.perf_counter_ns):
    seams = tuple(slow_geometry.source().seams())
    state = native.states(binary=True)[0]
    rng = random.Random(SEED)
    records = bytearray()
    endpoints = [{'attempts': 0, 'elapsed_ns': 0}]
    moves = sink = 0
    started = clock()
    for attempt in range(1, ATTEMPT_CAP+1):
        address = rng.randrange(60)
        edge, coin = divmod(address, 2)
        state, (_, current) = native.native(state, seams[edge], coin)
        records.append(address*3+current+1)
        # Declared host-only workload: never changes native state or choices.
        for _ in range(work):
            sink = (1664525*sink+1013904223) & 0xffffffff
        moves += abs(current)
        if current and moves % MOVEMENTS == 0:
            endpoints.append({'attempts': attempt, 'elapsed_ns': clock()-started})
            if moves == BLOCKS*MOVEMENTS:
                break
    complete = moves == BLOCKS*MOVEMENTS
    if not complete and endpoints[-1]['attempts'] != len(records):
        endpoints.append({'attempts': len(records), 'elapsed_ns': clock()-started})
    return {'name': name, 'host_work_iterations': work, 'records_hex': records.hex(),
            'block_endpoints': endpoints, 'complete': complete, 'work_sink': sink}


def capture():
    info = time.get_clock_info('perf_counter')
    if not info.monotonic:
        raise RuntimeError('this attachment requires a monotonic host performance counter')
    # Import/source initialization occurs before timed intervals. There is no
    # selection or rerun based on the measured performance.
    slow_geometry.source()
    return {'schema': 'oph.native_host_clock_records.v1', 'policy': POLICY,
            'capture_sources': _hashes(),
            'clock': {'api': 'time.perf_counter_ns', 'monotonic': info.monotonic,
                      'adjustable': info.adjustable, 'reported_resolution_seconds': info.resolution,
                      'implementation': info.implementation,
                      'traceable_accuracy_bound_supplied': False},
            'scheduler': 'Python random.Random(1729).randrange(60); pseudorandom software service',
            'initial_configuration': 0,
            'record_encoding': 'one byte per attempt: 3*(2*seam_index+coin)+current+1; lowercase hex',
            'runs': [_timed_run('baseline', 0), _timed_run('host_work_control', 64)]}


def _replay_run(run, seams):
    _strict_keys(run, ('name', 'host_work_iterations', 'records_hex', 'block_endpoints', 'complete', 'work_sink'), 'run')
    if type(run['complete']) is not bool:
        raise ValueError('completion must be Boolean')
    work = _int(run['host_work_iterations'], 0, 64, 'host workload')
    _int(run['work_sink'], 0, 2**32-1, 'host workload sink')
    encoded = run['records_hex']
    if type(encoded) is not str or len(encoded) > 2*ATTEMPT_CAP or re.fullmatch(r'(?:[0-9a-f]{2})*', encoded) is None:
        raise ValueError('bounded canonical hexadecimal source records required')
    records = bytes.fromhex(encoded)
    endpoints = run['block_endpoints']
    if type(endpoints) is not list or not 2 <= len(endpoints) <= BLOCKS+1:
        raise ValueError('complete bounded endpoint list required')
    previous = (-1, -1)
    for point in endpoints:
        _strict_keys(point, ('attempts', 'elapsed_ns'), 'clock endpoint')
        attempts = _int(point['attempts'], 0, ATTEMPT_CAP, 'endpoint attempt')
        elapsed = _int(point['elapsed_ns'], 0, 10**18, 'host elapsed nanoseconds')
        if attempts <= previous[0] or elapsed <= previous[1]:
            raise ValueError('clock endpoint order must increase strictly')
        previous = attempts, elapsed
    if endpoints[0] != {'attempts': 0, 'elapsed_ns': 0} or endpoints[-1]['attempts'] != len(records):
        raise ValueError('clock endpoints must retain exactly the complete source prefix')
    if run['complete'] and len(endpoints) != BLOCKS+1:
        raise ValueError('completed trace lacks a declared clock block')
    state = native.states(binary=True)[0]
    rng = random.Random(SEED)
    boundary = 1
    total_moves = 0
    counts = [0]*12
    block_counts = [0]*12
    blocks = []
    for attempt, code in enumerate(records, 1):
        address, observed = divmod(code, 3)
        _int(address, 0, 59, 'source address')
        if address != rng.randrange(60):
            raise ValueError('the recorded source address differs from its declared seeded service')
        seam_number, coin = divmod(address, 2)
        seam = seams[seam_number]
        state, (_, current) = native.native(state, seam, coin)
        if current != observed-1:
            raise ValueError(f'native current mismatch at attempt {attempt}')
        if current:
            total_moves += 1
            for port in seam:
                counts[port] += 1
                block_counts[port] += 1
        if attempt == endpoints[boundary]['attempts']:
            expected_complete = run['complete'] or boundary < len(endpoints)-1
            movements = sum(block_counts)//2
            if sum(block_counts) % 2 or (expected_complete and movements != MOVEMENTS):
                raise ValueError('a declared movement boundary has the wrong event count')
            if not expected_complete and movements >= MOVEMENTS:
                raise ValueError('a timeout block contains a discarded completed boundary')
            if expected_complete and current == 0:
                raise ValueError('a complete block must stop on its actual moving event')
            blocks.append({'block': boundary-1, 'native_attempts': attempt-endpoints[boundary-1]['attempts'],
                           'movements': movements, 'local_counts': block_counts,
                           'host_duration_ns': endpoints[boundary]['elapsed_ns']-endpoints[boundary-1]['elapsed_ns'],
                           'complete': expected_complete})
            block_counts = [0]*12
            boundary += 1
    sink, a, b, steps = 0, 1664525, 1013904223, work*len(records)
    while steps:
        if steps & 1:
            sink = (a*sink+b) & 0xffffffff
        a, b = (a*a) & 0xffffffff, (b*(a+1)) & 0xffffffff
        steps //= 2
    if sink != run['work_sink']:
        raise ValueError('the host workload sink differs from its declared computation')
    if not run['complete'] and len(records) != ATTEMPT_CAP:
        raise ValueError('an incomplete capture must retain its declared attempt cap')
    return {'name': run['name'], 'native_attempts': len(records), 'movements': total_moves,
            'complete': run['complete'], 'all_native_currents_replayed': True,
            'host_work_iterations': work, 'blocks': blocks,
            'seeded_source_addresses_and_work_sink_replayed': True,
            'records_sha256': hashlib.sha256(records).hexdigest(),
            'all_twelve_counts': counts}


def analyse(raw):
    _strict_keys(raw, ('schema', 'policy', 'capture_sources', 'clock', 'scheduler',
                      'initial_configuration', 'record_encoding', 'runs'), 'raw capture')
    if (raw['schema'] != 'oph.native_host_clock_records.v1'
            or canonical(raw['policy']) != canonical(POLICY)
            or type(raw['initial_configuration']) is not int or raw['initial_configuration'] != 0):
        raise ValueError('the declared host-clock protocol changed')
    if raw['capture_sources'] != _hashes():
        raise ValueError('capture source bindings changed')
    if raw['scheduler'] != 'Python random.Random(1729).randrange(60); pseudorandom software service':
        raise ValueError('the declared pseudorandom scheduler changed')
    if raw['record_encoding'] != 'one byte per attempt: 3*(2*seam_index+coin)+current+1; lowercase hex':
        raise ValueError('the raw source-record encoding changed')
    clock = raw['clock']
    _strict_keys(clock, ('api', 'monotonic', 'adjustable', 'reported_resolution_seconds', 'implementation',
                        'traceable_accuracy_bound_supplied'), 'clock metadata')
    if (clock['api'] != 'time.perf_counter_ns' or clock['monotonic'] is not True
            or type(clock['adjustable']) is not bool or clock['traceable_accuracy_bound_supplied'] is not False
            or type(clock['implementation']) is not str
            or re.fullmatch(r'[A-Za-z0-9_(). +\-]{1,100}', clock['implementation']) is None
            or type(clock['reported_resolution_seconds']) not in (int, float)
            or not math.isfinite(clock['reported_resolution_seconds'])
            or not 0 < clock['reported_resolution_seconds'] <= 1):
        raise ValueError('valid bounded clock metadata without an invented accuracy guarantee required')
    if type(raw['runs']) is not list or len(raw['runs']) != 2:
        raise ValueError('both unselected baseline and workload-control captures are required')
    seams = tuple(slow_geometry.source().seams())
    runs = [_replay_run(row, seams) for row in raw['runs']]
    if [(row['name'], row['host_work_iterations']) for row in runs] != [('baseline', 0), ('host_work_control', 64)]:
        raise ValueError('baseline and host-only workload protocols changed')
    if runs[0]['records_sha256'] != runs[1]['records_sha256']:
        raise ValueError('the timing control did not preserve the same native record stream')
    if not runs[0]['blocks'][0]['complete']:
        raise ValueError('the calibration block did not complete')
    cal = runs[0]['blocks'][0]
    count = cal['local_counts'][CALIBRATION_PORT]
    if count == 0:
        raise ValueError('the chosen calibration clock never ticked')
    duration_per_tick = F(cal['host_duration_ns'], count)
    tests = []
    for run in runs:
        for block in run['blocks']:
            if run['name'] == 'baseline' and block['block'] == 0:
                continue
            prediction = duration_per_tick*block['local_counts'][VALIDATION_PORT]
            if prediction == 0:
                relative_error = None
                verdict = 'reject'
            else:
                relative_error = abs(F(block['host_duration_ns'])-prediction)/prediction
                verdict = 'within_tolerance' if block['complete'] and relative_error <= TIMING_TOLERANCE else 'reject'
            tests.append({'run': run['name'], 'block': block['block'], 'clock_port': VALIDATION_PORT,
                          'observed_host_duration_ns': block['host_duration_ns'],
                          'predicted_host_duration_ns': str(prediction),
                          'relative_error': None if relative_error is None else str(relative_error),
                          'engineering_verdict': verdict})
    # The exact certificate constants are independently checked by comparison's
    # source receipt. This test checks the retained finite rate envelope only.
    rate_tests = []
    for block in runs[0]['blocks']:
        port = CALIBRATION_PORT if block['block'] == 0 else VALIDATION_PORT
        n, k = block['movements'], block['local_counts'][port]
        deviation = abs(F(k, n)-F(1, 6)) if n else None
        rate_tests.append({'block': block['block'], 'clock_port': port,
                           'absolute_rate_deviation': None if deviation is None else str(deviation),
                           'within_declared_native_rate_envelope': bool(block['complete'] and deviation <= RATE_ERROR)})
    # Chebyshev: these are model bounds for ideal independent uniform choices,
    # not confidence guarantees for a single deterministic pseudorandom stream.
    R, V = F(38, 25), F(30341677, 47045000)
    tail = min(F(1), MOVEMENTS*V/(MOVEMENTS*RATE_ERROR-R)**2)
    control_ratios = [str(F(b['host_duration_ns'], a['host_duration_ns']))
                      for a, b in zip(runs[0]['blocks'], runs[1]['blocks'], strict=True)]
    return {'schema': 'oph.native_host_clock_attachment.v1', 'replay_integrity_passed': True,
            'sources': _hashes(), 'clock_metadata': clock,
            'raw_records_canonical_sha256': hashlib.sha256(canonical(raw).encode()).hexdigest(),
            'calibration': {'run': 'baseline', 'block': 0, 'clock_port': CALIBRATION_PORT,
                            'local_tick_count': count, 'host_duration_ns': cal['host_duration_ns'],
                            'nanoseconds_per_local_tick': str(duration_per_tick),
                            'other_clock_and_later_blocks_used_in_fit': False},
            'native_record_replays': runs,
            'held_out_clock_tests': tests,
            'native_rate_tests': rate_tests,
            'declared_native_rate_envelope': {'absolute_error': str(RATE_ERROR),
                'per_block_failure_probability_upper_bound_for_ideal_scheduler': str(tail),
                'four_block_union_bound_for_ideal_scheduler': str(min(F(1), 4*tail)),
                'pseudorandom_sample_has_independent_randomness_certificate': False},
            'timing_control': {'same_complete_native_records': True,
                'control_to_baseline_duration_ratios_by_block': control_ratios,
                'control_rejected_in_every_measured_block': all(t['engineering_verdict'] == 'reject'
                    for t in tests if t['run'] == 'host_work_control'),
                'timing_tolerance_is_a_confidence_interval': False},
            'scope': {'actual_native_execution_measured_by_host_clock': True,
                      'replay_reproduces_timestamps': False,
                      'new_capture_requires_new_measurement': True,
                      'host_clock_traceable_lab_calibration_established': False,
                      'clock_resolution_is_an_accuracy_bound': False,
                      'host_service_latency_or_rate_error_bounded': False,
                      'physical_proper_time_or_modular_phase_derived': False,
                      'empirical_confirmation_of_oph_as_physics': False}}


def run(path=RAW_PATH):
    return analyse(strict_json.loads(Path(path).read_text(encoding='utf-8')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--capture', type=Path)
    modes.add_argument('--write', type=Path)
    modes.add_argument('--check', type=Path)
    parser.add_argument('--records', type=Path, default=RAW_PATH)
    args = parser.parse_args()
    if args.capture:
        with args.capture.open('x', encoding='utf-8') as stream:
            stream.write(canonical(capture()))
        print(f'Captured native host-clock records at {args.capture}')
        return
    result = run(args.records)
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
        print('Native host-clock records and retained analysis reproduced; timestamps were not remeasured')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
