"""Prospectively registered exact census and complete native record histories."""
import argparse
from collections import deque
import gzip
import hashlib
import json
from pathlib import Path
import random
from source import (CANDIDATE, N, SEAMS, native, body_event, decode_probe,
                    decode_body_path, run_history)

HERE = Path(__file__).resolve().parent


def compositions(total, n=N):
    if n == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, n - 1):
                yield (first,) + rest


def exact_census():
    sectors, counts = [], {'states': 0, 'native_branches': 0,
                          'binary_states': 0, 'binary_exit_branches': 0}
    for total in range(5):
        states = list(compositions(total))
        index = {x: i for i, x in enumerate(states)}
        reverse = [[] for _ in states]
        branches = []
        for i, x in enumerate(states):
            for e, edge in enumerate(SEAMS):
                for coin in (0, 1):
                    y, (_, current) = native(x, edge, coin)
                    target = index[y]
                    branches.append([i, e, coin, current, target])
                    reverse[target].append(i)
                    if max(x) <= 1 and max(y) > 1:
                        counts['binary_exit_branches'] += 1
        distance = [-1] * len(states)
        queue = deque(i for i, x in enumerate(states) if max(x) <= 1)
        for i in queue:
            distance[i] = 0
        while queue:
            target = queue.popleft()
            for i in reverse[target]:
                if distance[i] < 0:
                    distance[i] = distance[target] + 1
                    queue.append(i)
        sectors.append({'total': total, 'states': states,
                        'branches': branches, 'distance_to_binary': distance})
        counts['states'] += len(states)
        counts['native_branches'] += len(branches)
        counts['binary_states'] += sum(max(x) <= 1 for x in states)
    counts['unreachable_states'] = sum(d < 0 for s in sectors
                                      for d in s['distance_to_binary'])
    counts['max_shortest_entry_events'] = max(d for s in sectors
                                             for d in s['distance_to_binary'])
    return sectors, counts


def normalize(x):
    return [v - x[0] for v in x]


def audit_history(record):
    x = tuple(record['initial_audit_only'])
    path = [x]
    for a, b, c, j in record['body']:
        x, (_, actual) = body_event(x, (a, b), c, record['law'])
        if j != actual:
            raise ValueError('Body current mismatch')
        path.append(x)
    for a, b, c, j in record['probe']:
        x, (_, actual) = native(x, (a, b), c)
        if j != actual:
            raise ValueError('Probe current mismatch')
    decoded = decode_probe(record['probe'], record['probe_cap'])
    reconstructed = decode_body_path(record['body'], decoded)
    if reconstructed is not None:
        if [normalize(y) for y in path] != [normalize(y) for y in reconstructed]:
            raise ValueError('Current-only reconstruction mismatch')
    return reconstructed


def preparations():
    result = []
    for p in range(N):
        for family in ('double', 'quadruple', 'signed', 'shifted'):
            x = [2 if family == 'shifted' else 0] * N
            if family == 'signed':
                x[p], x[(p + 1) % N] = 8, -6
            else:
                x[p] += 4 if family == 'quadruple' else 2
            result.append((family, x))
    return result


def calculate():
    sectors, census = exact_census()
    histories, audited_conclusive = [], 0
    for i, (family, initial) in enumerate(preparations()):
        for law in ('native', 'swap_hold'):
            record = run_history(initial, 72301000 + i, law)
            record['family'] = family
            audited_conclusive += audit_history(record) is not None
            histories.append(record)
    probes = []
    for i, initial in enumerate(([2] + [0] * 11,
                                 [1, 1] + [0] * 10,
                                 [0] * 12,
                                 [2] * 12)):
        r = run_history(initial, 72301900 + i, horizon=0)
        audit_history(r)
        probes.append(r)
    timeout = run_history([2] + [0] * 11, 72301999, horizon=0,
                          cap=2, forced_probe_coin=0)
    audit_history(timeout)
    mutated = json.loads(json.dumps(probes[0]))
    mutated['probe'][0][3] += 1
    try:
        audit_history(mutated)
        mutation_rejected = False
    except ValueError:
        mutation_rejected = True
    native_records = [r for r in histories if r['law'] == 'native']
    controls = [r for r in histories if r['law'] == 'swap_hold']
    binary_exits = sum(any(v > 1 for v in r['body_spans_from_records'][
                       r['first_balanced_body_event']:]) for r in native_records
                       if r['first_balanced_body_event'] is not None)
    degree = sum(0 in e for e in SEAMS)
    receipt = {'schema': 1, 'candidate': CANDIDATE, 'class': 'N+',
        'official_row_status': 'untested', 'full_row_solved': False,
        'coherent_memory_used': False, 'nature_values_used': [],
        'exact_census': census,
        'double_initial_entry_probability': f'{degree}/{len(SEAMS)}',
        'double_initial_mean_entry_events': len(SEAMS) // degree,
        'histories': {'native': len(native_records), 'changed_law': len(controls),
            'native_balanced_within_horizon': sum(r['first_balanced_body_event']
               is not None for r in native_records),
            'changed_law_balanced_within_horizon': sum(r['first_balanced_body_event']
               is not None for r in controls),
            'inconclusive_terminal_readings': sum(r['probe_outcome'] != 'conclusive'
               for r in histories), 'postentry_binary_exits': binary_exits,
            'body_events': sum(len(r['body']) for r in histories),
            'native_terminal_probe_events': sum(len(r['probe']) for r in histories)},
        'detector_controls': {'double_span': probes[0]['body_spans_from_records'],
            'binary_span': probes[1]['body_spans_from_records'],
            'zero_span': probes[2]['body_spans_from_records'],
            'uniform_two_span': probes[3]['body_spans_from_records'],
            'forced_timeout_outcome': timeout['probe_outcome'],
            'forced_timeout_decoded': timeout['decoded_preprobe_relative'],
            'changed_current_rejected_by_independent_event_audit': mutation_rejected},
        'limitations': ['Integer alphabet supplied; no fractional native load',
            'Selection is adjacent occupancy levels, not electric unit derivation',
            'Absolute background invisible to current-only observer',
            '0/1 conclusion requires total in [0,12]',
            'Uniform all-two state is fixed and is retained as counterexample',
            'Proof for unbounded integer inputs exceeds finite enumerated census',
            'No physical charge dictionary, topology derivation, G6 or G9 award']}
    gates = {
        'all_1820_states_enumerated': census['states'] == 1820,
        'all_109200_native_branches_executed': census['native_branches'] == 109200,
        'binary_class_closed': census['binary_exit_branches'] == 0,
        'every_census_state_reaches_binary': census['unreachable_states'] == 0,
        'double_detector_reads_two': probes[0]['body_spans_from_records'] == [2],
        'binary_detector_reads_one': probes[1]['body_spans_from_records'] == [1],
        'timeout_retained': timeout['probe_outcome'] == 'inconclusive' and
                            timeout['decoded_preprobe_relative'] is None,
        'current_mutation_rejected': mutation_rejected,
        'controls_do_not_enter_balanced_class': all(r['first_balanced_body_event']
            is None for r in controls),
        'no_native_postentry_exit': binary_exits == 0,
        'all_conclusive_reconstructions_audited': audited_conclusive ==
            sum(r['probe_outcome'] == 'conclusive' for r in histories)}
    receipt['checks'] = gates
    receipt['passed'] = all(gates.values())
    evidence = {'sectors': sectors, 'histories': histories,
                'detector_trials': probes, 'forced_timeout': timeout,
                'mutation': mutated}
    return receipt, evidence


def encode(obj):
    return (json.dumps(obj, sort_keys=True, separators=(',', ':')) + '\n').encode()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    receipt, evidence = calculate()
    evidence_bytes = gzip.compress(encode(evidence), mtime=0)
    receipt['evidence_sha256'] = hashlib.sha256(evidence_bytes).hexdigest()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / 'RECEIPT.json').write_bytes(encode(receipt))
        (args.output / 'EVIDENCE.json.gz').write_bytes(evidence_bytes)
    print(json.dumps(receipt, indent=2))
    raise SystemExit(0 if receipt['passed'] else 1)


if __name__ == '__main__':
    main()
