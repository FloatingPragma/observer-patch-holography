"""Standard-library independent checker. Imports none of the producer code."""
from collections import deque
from fractions import Fraction
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def endpoint(values, edge, coin, mode):
    a, b = edge
    result = list(values)
    if mode == 'native':
        total = values[a] + values[b]
        low, high = total // 2, -((-total) // 2)
        result[a], result[b] = (low, high) if coin == 0 else (high, low)
    elif mode == 'swap_hold':
        if coin:
            result[a], result[b] = values[b], values[a]
    else:
        raise AssertionError(mode)
    return result, values[a] - result[a]


def check_history(record):
    original = record['initial_audit_only']
    actual = original[:]
    actual_path = [actual[:]]
    for a, b, coin, recorded in record['body']:
        actual, emitted = endpoint(actual, (a, b), coin, record['law'])
        assert emitted == recorded
        actual_path.append(actual[:])
    assert len(record['body']) == record['horizon']
    before_probe = actual[:]
    difference_records, delta, offset = [], [0] * len(actual), 0
    tape = record['probe']
    while offset < len(tape):
        a, b, first_coin, first_current = tape[offset]
        original_delta = delta[a] - delta[b]
        completed = False
        for attempt in range(record['probe_cap']):
            assert offset < len(tape)
            aa, bb, coin, recorded = tape[offset]
            assert (aa, bb) == (a, b)
            actual, emitted = endpoint(actual, (a, b), coin, 'native')
            assert emitted == recorded
            delta[a] -= recorded
            delta[b] += recorded
            offset += 1
            if coin != first_coin:
                difference_records.append((a, b,
                    2 * first_current + recorded - original_delta))
                completed = True
                break
        if not completed:
            assert offset == len(tape)
            assert record['probe_outcome'] == 'inconclusive'
            assert record['decoded_preprobe_relative'] is None
            return 0
    known = {0: 0}
    while len(known) < len(actual):
        before = len(known)
        for a, b, difference in difference_records:
            if a in known:
                known[b] = known[a] - difference
            elif b in known:
                known[a] = known[b] + difference
        assert len(known) > before
    inferred = [known[i] for i in range(len(actual))]
    assert inferred == [v - before_probe[0] for v in before_probe]
    assert inferred == record['decoded_preprobe_relative']
    reverse, last = [inferred[:]], inferred[:]
    for a, b, coin, current in record['body'][::-1]:
        last[a] += current
        last[b] -= current
        reverse.append(last[:])
    observer_path = reverse[::-1]
    for truth, read in zip(actual_path, observer_path):
        assert [v - truth[0] for v in truth] == [v - read[0] for v in read]
    spans = [max(v) - min(v) for v in observer_path]
    assert spans == record['body_spans_from_records']
    assert next((i for i, v in enumerate(spans) if v <= 1), None) == record[
        'first_balanced_body_event']
    return len(tape) + len(record['body'])


def verify(directory):
    edges = json.loads((HERE / 'seams.json').read_text())
    receipt = json.loads((directory / 'RECEIPT.json').read_text())
    raw = (directory / 'EVIDENCE.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == receipt['evidence_sha256']
    data = json.loads(gzip.decompress(raw))
    checked, states_checked = 0, 0
    for sector in data['sectors']:
        states, total = sector['states'], sector['total']
        assert len(states) == comb(total + 11, 11)
        assert len({tuple(x) for x in states}) == len(states)
        assert all(len(x) == 12 and sum(x) == total and min(x) >= 0 for x in states)
        reverse = [[] for _ in states]
        visits = set()
        double_probabilities = []
        for i, e, coin, current, target in sector['branches']:
            assert (i, e, coin) not in visits
            visits.add((i, e, coin))
            expected, emitted = endpoint(states[i], edges[e], coin, 'native')
            assert expected == states[target] and current == emitted
            assert sum(v * v for v in expected) <= sum(v * v for v in states[i])
            if max(states[i]) <= 1:
                assert max(expected) <= 1
            reverse[target].append(i)
            checked += 1
        assert len(visits) == len(states) * len(edges) * 2
        distance = [-1] * len(states)
        queue = deque(i for i, x in enumerate(states) if max(x) <= 1)
        for i in queue:
            distance[i] = 0
        while queue:
            j = queue.popleft()
            for i in reverse[j]:
                if distance[i] == -1:
                    distance[i] = distance[j] + 1
                    queue.append(i)
        assert distance == sector['distance_to_binary'] and min(distance) >= 0
        if total == 2:
            for i, x in enumerate(states):
                if max(x) == 2:
                    exits = sum(max(states[t]) <= 1 for j, _, _, _, t in
                                sector['branches'] if j == i)
                    assert Fraction(exits, 2 * len(edges)) == Fraction(1, 6)
                    double_probabilities.append(Fraction(exits, 2 * len(edges)))
            assert len(double_probabilities) == 12
        states_checked += len(states)
    calls = sum(check_history(r) for r in data['histories'] + data['detector_trials'])
    assert check_history(data['forced_timeout']) == 0
    try:
        check_history(data['mutation'])
    except AssertionError:
        mutation_rejected = True
    else:
        mutation_rejected = False
    assert mutation_rejected
    for r in data['histories']:
        if r['law'] == 'swap_hold':
            assert len(set(r['body_spans_from_records'])) == 1
            assert r['first_balanced_body_event'] is None
    assert checked == 109200 and states_checked == 1820
    assert receipt['passed'] and all(receipt['checks'].values())
    return {'independent_branch_checks': checked, 'states_checked': states_checked,
            'history_records_checked': calls, 'same_group': True,
            'mutated_record_rejected': mutation_rejected, 'passed': True}


if __name__ == '__main__':
    print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE), indent=2))
