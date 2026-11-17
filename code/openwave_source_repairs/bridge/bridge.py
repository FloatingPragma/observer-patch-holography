"""Actual OLD240 multi-load selection and continuing discovered-bond program."""
from collections import Counter
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import random
from native_event import native

HERE = Path(__file__).resolve().parent
EDGES = [tuple(e) for e in json.loads((HERE / 'seams.json').read_text())]
REF = (12, 13)
CALLS = 0


def tree():
    seen, queue, result = {0}, [0], []
    for a in queue:
        for u, v in EDGES:
            b = v if a == u else u if a == v else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                result.append((a, b))
    assert len(result) == 11
    return result


TREE = tree()


def call(x, edge, coin, law='native'):
    global CALLS
    if law == 'native':
        CALLS += 1
        return native(x, edge, coin)
    a, b = edge
    y = list(x)
    j = coin * (x[a] - x[b])
    y[a] -= j
    y[b] += j
    return tuple(y), (coin, j)


def selector(x, rng, forced=False):
    tape = []
    for _ in range(16):
        label = 0
        for _ in range(5):
            c = 1 if forced else rng.randrange(2)
            x, (_, j) = call(x, REF, c)
            tape.append([*REF, c, j])
            label = 2 * label + c
        if label < 30:
            return x, EDGES[label], tape
    return x, None, tape


def selection(x, seed, law):
    rng, blocks = random.Random(seed), []
    for _ in range(1024):
        x, edge, refs = selector(x, rng)
        if edge is None:
            blocks.append({'reference': refs, 'body': None})
            return x, blocks, 'selector_timeout'
        c = rng.randrange(2)
        x, (_, j) = call(x, edge, c, law)
        blocks.append({'reference': refs, 'body': [*edge, c, j]})
    return x, blocks, 'selection_horizon'


def census(x, seed, cap=32):
    rng, tape = random.Random(seed), []
    for a, b in TREE:
        first = None
        for _ in range(cap):
            c = rng.randrange(2)
            x, (_, j) = call(x, (a, b), c)
            tape.append([a, b, c, j])
            if first is None:
                first = c
            elif c != first:
                break
        else:
            return x, tape
    return x, tape


def read_census(tape, cap=32):
    """Current-only reader with explicit total1/first-cell domain promise."""
    relative, delta, cursor = [None] * 12, [0] * 12, 0
    relative[0] = 0
    for a, b in TREE:
        correction = delta[a] - delta[b]
        first = None
        for _ in range(cap):
            if cursor == len(tape):
                return None
            aa, bb, c, j = tape[cursor]
            cursor += 1
            if (aa, bb) != (a, b) or c not in (0, 1):
                raise ValueError('Malformed census')
            delta[a] -= j
            delta[b] += j
            if first is None:
                first, j0 = c, j
            elif c != first:
                relative[b] = relative[a] - (2 * j0 + j - correction)
                break
        else:
            if cursor != len(tape):
                raise ValueError('Records after timeout')
            return None
    if cursor != len(tape):
        raise ValueError('Extra census records')
    offset = F(1 - sum(relative), 12)
    pre = [F(v) + offset for v in relative]
    if any(v.denominator != 1 for v in pre):
        raise ValueError('Total1 integer promise inconsistent')
    pre = [int(v) for v in pre] + [0] * 228
    post = pre[:]
    for i, v in enumerate(delta):
        post[i] += v
    return pre, post


def apply_record(x, record):
    a, b, _, j = record
    y = list(x)
    y[a] -= j
    y[b] += j
    return y


def position(x):
    indices = [i for i, v in enumerate(x) if v]
    if len(indices) != 1 or x[indices[0]] != 1:
        raise ValueError('Not an actual positive source unit')
    return indices[0]


def encoded(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':')) + '\n').encode()


def handoff(x, observed, prefix):
    if list(x) != observed:
        raise ValueError('Actual/observed handoff mismatch')
    return {'actual_source_audit_only': x, 'observed_source': observed,
            'sign': 1, 'occupied': position(observed),
            'prefix_sha256': hashlib.sha256(encoded(prefix)).hexdigest(),
            'encoding': 'declared basis isometry |delta_s> -> |s>; no quantum quotient claim'}


def one_trial(initial, index, law):
    record = {'initial_audit_only': initial, 'index': index, 'selection_law': law,
              'discovery': [], 'holding': [], 'basis_handoff': None,
              'final_handoff': None, 'bond': None}
    x, blocks, status = selection(tuple(initial), 72321000 + index, law)
    record['selection'] = blocks
    x, probe = census(x, 72323000 + index)
    record['census'] = probe
    read = read_census(probe)
    if read is None:
        record['outcome'] = 'census_timeout'
        return record
    pre, observed = read
    z, path = pre[:], [pre[:]]
    for block in blocks[::-1]:
        if block['body'] is not None:
            a, b, _, j = block['body']
            z[a] += j
            z[b] -= j
        path.append(z[:])
    spans = [max(v[:12]) - min(v[:12]) for v in path[::-1]]
    record['selection_spans'] = spans
    record['first_selected'] = next((i for i, s in enumerate(spans) if s <= 1), None)
    if status != 'selection_horizon':
        record['outcome'] = status
        return record
    if max(pre[:12]) - min(pre[:12]) > 1:
        record['outcome'] = 'unsettled_at_horizon'
        return record
    record['basis_handoff'] = handoff(x, observed, {'selection': blocks, 'census': probe})
    rng = random.Random(72322000 + index)
    for _ in range(256):
        x, edge, refs = selector(x, rng)
        for row in refs:
            observed = apply_record(observed, row)
        if edge is None:
            record['discovery'].append({'reference': refs, 'body': None})
            record['outcome'] = 'discovery_selector_timeout'
            return record
        before = position(observed)
        c = rng.randrange(2)
        x, (_, j) = call(x, edge, c)
        row = [*edge, c, j]
        observed = apply_record(observed, row)
        record['discovery'].append({'reference': refs, 'body': row})
        if j:
            after = position(observed)
            record['bond'] = [before, after]
            break
    else:
        record['outcome'] = 'discovery_horizon_timeout'
        return record
    rng = random.Random(72324000 + index)
    edge = tuple(sorted(record['bond']))
    for _ in range(64):
        c = rng.randrange(2)
        x, (_, j) = call(x, edge, c)
        row = [*edge, c, j]
        record['holding'].append(row)
        observed = apply_record(observed, row)
    record['final_handoff'] = handoff(x, observed, {
        'basis': record['basis_handoff']['prefix_sha256'],
        'discovery': record['discovery'], 'holding': record['holding']})
    record['outcome'] = 'selected_and_continued'
    return record


def branch_certificate(edge):
    distributions, rows = {}, []
    a, b = edge
    for law in ('native', 'hold', 'swap'):
        distributions[law] = []
        for start in edge:
            x = tuple(int(i == start) for i in range(240))
            outcomes = []
            for c in (0, 1):
                if law == 'native':
                    y, (_, j) = call(x, edge, c)
                else:
                    j = 0 if law == 'hold' else x[a] - x[b]
                    y = apply_record(x, [a, b, c, j])
                target = position(y)
                outcomes.append(target)
                rows.append([law, start, c, j, target])
            distributions[law].append(Counter(outcomes))
    distances = {law: str(sum(F(abs(ds[0][k] - ds[1][k]), 4) for k in edge))
                 for law, ds in distributions.items()}
    assert distances == {'native': '0', 'hold': '1', 'swap': '1'}
    return {'edge': edge, 'branches': rows, 'total_variation': distances}


def calculate():
    global CALLS
    CALLS = 0
    graph = json.loads((HERE / 'graph.json').read_text())
    assert graph['ports'] == 240 and len(graph['seams']) == 630
    assert [list(e) for e in EDGES] == [e for e in graph['seams'] if max(e) < 12]
    assert list(REF) in graph['seams']
    trials = []
    for p in range(12):
        for large in (2, 8, 17):
            initial = [0] * 240
            initial[p], initial[(p + 1) % 12] = large, 1 - large
            index = len(trials) // 2
            for law in ('native', 'swap_hold'):
                trials.append(one_trial(initial, index, law))
    main = [r for r in trials if r['selection_law'] == 'native']
    controls = [r for r in trials if r['selection_law'] == 'swap_hold']
    certificates = [branch_certificate(tuple(sorted(r['bond']))) for r in main]
    initial = [2, -1] + [0] * 238
    _, detector_tape = census(initial, 72325000)
    detector = read_census(detector_tape)
    _, forced_census = census(initial, 72325001, cap=1)
    _, edge, forced_selector = selector(initial, random.Random(72325002), forced=True)
    gates = {
        'native_selected_36': all(r['first_selected'] is not None for r in main),
        'no_native_postentry_exit': all(max(r['selection_spans'][r['first_selected']:]) <= 1 for r in main),
        'control_selected_0': all(r['first_selected'] is None for r in controls),
        'actual_record_handoffs_36': all(r['basis_handoff'] is not None for r in main),
        'continuations_36': all(r['outcome'] == 'selected_and_continued' for r in main),
        'complete_holding_64': all(len(r['holding']) == 64 for r in main),
        'native_opposite_start_laws_equal': all(r['total_variation']['native'] == '0' for r in certificates),
        'oriented_controls_distinct': all(r['total_variation']['hold'] == r['total_variation']['swap'] == '1' for r in certificates),
        'multiload_detector': detector[0] == initial,
        'census_timeout_not_success': read_census(forced_census, cap=1) is None,
        'selector_timeout_retained': edge is None and len(forced_selector) == 80,
    }
    evidence = {'trials': trials, 'bond_certificates': certificates,
                'detector': {'initial_audit_only': initial, 'census': detector_tape},
                'forced_census': forced_census, 'forced_selector': forced_selector}
    receipt = {'candidate': 'A-old240-selected-unit-handoff-v1', 'gates': gates,
               'all_gates_pass': all(gates.values()), 'literal_native_calls': CALLS,
               'control_body_events': sum(len(r['selection']) for r in controls),
               'outcomes': dict(Counter(r['outcome'] for r in trials)),
               'selected_handoff_addresses': [r['basis_handoff']['occupied'] for r in main],
               'source240_actual_not_difference': True, 'coherent_state_selection_claimed': False,
               'physical_spin_or_energy_claimed': False, 'full_row_solved': False,
               'selector_rejected_attempts': sum((len(b['reference']) // 5 - (b['body'] is not None))
                     for r in trials for b in r['selection'] + r['discovery'])}
    return receipt, evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    prereg = json.loads((HERE / 'PREREG.json').read_text())
    for name, expected in prereg['sha256'].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == expected
    receipt, evidence = calculate()
    raw = gzip.compress(encoded(evidence), mtime=0)
    receipt['prereg_sha256'] = hashlib.sha256((HERE / 'PREREG.json').read_bytes()).hexdigest()
    receipt['evidence_sha256'] = hashlib.sha256(raw).hexdigest()
    if args.record:
        for name, data in [('RECEIPT.json', encoded(receipt)), ('RECORDS.json.gz', raw)]:
            with (HERE / name).open('xb') as f:
                f.write(data)
    else:
        assert json.loads((HERE / 'RECEIPT.json').read_text()) == receipt
        assert (HERE / 'RECORDS.json.gz').read_bytes() == raw
    print(json.dumps(receipt, indent=2))
    if not receipt['all_gates_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
