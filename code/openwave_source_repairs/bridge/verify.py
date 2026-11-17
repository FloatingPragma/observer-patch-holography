"""Separate standard-library audit, without producer imports (same author)."""
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import random
import sys


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def endpoint(x, row, law='native'):
    a, b, c, j = row
    y = list(x)
    if law == 'native':
        total = x[a] + x[b]
        lo, hi = total // 2, -((-total) // 2)
        y[a], y[b] = (lo, hi) if c == 0 else (hi, lo)
    elif law == 'swap_hold':
        if c:
            y[a], y[b] = x[b], x[a]
    else:
        raise AssertionError(law)
    assert c in (0, 1) and j == x[a] - y[a]
    return y


def read_census(tape, tree, cap=32):
    delta, values, cursor = [0] * 12, {0: 0}, 0
    for a, b in tree:
        before = delta[a] - delta[b]
        assert cursor < len(tape)
        first = tape[cursor]
        complete = False
        for _ in range(cap):
            assert cursor < len(tape)
            aa, bb, c, j = tape[cursor]
            assert (aa, bb) == (a, b)
            cursor += 1
            delta[a] -= j
            delta[b] += j
            if c != first[2]:
                d = 2 * first[3] + j - before
                values[b] = values[a] - d
                complete = True
                break
        if not complete:
            assert cursor == len(tape)
            return None
    assert cursor == len(tape)
    base = F(1 - sum(values.values()), 12)
    pre = [F(values[i]) + base for i in range(12)]
    assert all(v.denominator == 1 for v in pre)
    pre = list(map(int, pre)) + [0] * 228
    post = pre[:]
    for i in range(12):
        post[i] += delta[i]
    return pre, post


def verify(path):
    p = Path(path)
    edges = [tuple(e) for e in json.loads((p / 'seams.json').read_text())]
    graph = json.loads((p / 'graph.json').read_text())
    assert graph['ports'] == 240 and len(graph['seams']) == 630
    assert [list(e) for e in edges] == [e for e in graph['seams'] if max(e) < 12]
    assert [12, 13] in graph['seams']
    queue, seen, tree = [0], {0}, []
    for a in queue:
        for u, v in edges:
            b = v if a == u else u if a == v else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                tree.append((a, b))
    pre = json.loads((p / 'PREREG.json').read_text())
    for name, value in pre['sha256'].items():
        assert hashlib.sha256((p / name).read_bytes()).hexdigest() == value
    receipt = json.loads((p / 'RECEIPT.json').read_text())
    raw = (p / 'RECORDS.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == receipt['evidence_sha256']
    assert hashlib.sha256((p / 'PREREG.json').read_bytes()).hexdigest() == receipt['prereg_sha256']
    data = json.loads(gzip.decompress(raw))
    calls, rejected = 0, 0

    def blocks(x, rows, seed, law):
        nonlocal calls, rejected
        rng, states = random.Random(seed), [x[:]]
        for block in rows:
            refs = block['reference']
            assert len(refs) in range(5, 81, 5)
            for start in range(0, len(refs), 5):
                label = 0
                for row in refs[start:start + 5]:
                    assert row[:2] == [12, 13]
                    assert row[2] == rng.randrange(2)
                    x = endpoint(x, row)
                    calls += 1
                    assert row[3] == 0
                    label = 2 * label + row[2]
                if start + 5 != len(refs):
                    assert label >= 30
                    rejected += 1
            body = block['body']
            if body is None:
                assert label >= 30 and len(refs) == 80
                rejected += 1
                return x, states
            assert label < 30 and body[:2] == list(edges[label])
            assert body[2] == rng.randrange(2)
            x = endpoint(x, body, law)
            calls += law == 'native'
            states.append(x[:])
        return x, states

    def check_probe(x, tape, seed, cap=32):
        nonlocal calls
        decoded = read_census(tape, tree, cap)
        before = x[:]
        rng = random.Random(seed)
        for row in tape:
            assert row[2] == rng.randrange(2)
            x = endpoint(x, row)
            calls += 1
        if decoded is not None:
            assert decoded == (before, x)
        return x, decoded

    def check_handoff(actual, observed, handoff, prefix):
        assert actual == observed == handoff['actual_source_audit_only'] == handoff['observed_source']
        occupied = [i for i, v in enumerate(observed) if v]
        assert len(occupied) == 1 and observed[occupied[0]] == 1 and occupied[0] < 12
        assert handoff['occupied'] == occupied[0] and handoff['sign'] == 1
        assert handoff['prefix_sha256'] == hashlib.sha256(canonical(prefix)).hexdigest()

    assert len(data['trials']) == 72
    keys = {(r['index'], r['selection_law']) for r in data['trials']}
    assert keys == {(i, law) for i in range(36) for law in ('native', 'swap_hold')}
    for r in data['trials']:
        i, law = r['index'], r['selection_law']
        paxis, family = divmod(i, 3)
        large = (2, 8, 17)[family]
        initial = [0] * 240
        initial[paxis], initial[(paxis + 1) % 12] = large, 1 - large
        assert initial == r['initial_audit_only'] and sum(initial) == 1
        assert len(r['selection']) == 1024
        x, states = blocks(initial, r['selection'], 72321000 + i, law)
        spans = [max(z[:12]) - min(z[:12]) for z in states]
        assert spans == r['selection_spans']
        first = next((j for j, s in enumerate(spans) if s <= 1), None)
        assert first == r['first_selected']
        x, decoded = check_probe(x, r['census'], 72323000 + i)
        assert decoded is not None
        if law == 'swap_hold':
            assert first is None and r['outcome'] == 'unsettled_at_horizon'
            assert r['discovery'] == r['holding'] == []
            assert r['basis_handoff'] is r['final_handoff'] is None
            continue
        assert first is not None and max(spans[first:]) <= 1
        observed = decoded[1]
        check_handoff(x, observed, r['basis_handoff'], {'selection': r['selection'], 'census': r['census']})
        x, _ = blocks(x, r['discovery'], 72322000 + i, 'native')
        assert len(r['discovery']) in range(1, 257)
        assert all(b['body'][3] == 0 for b in r['discovery'][:-1])
        last = r['discovery'][-1]['body']
        assert abs(last[3]) == 1
        departed = observed.index(1)
        for block in r['discovery']:
            for a, b, _, j in block['reference'] + [block['body']]:
                observed[a] -= j
                observed[b] += j
        arrived = observed.index(1)
        assert r['bond'] == [departed, arrived]
        assert sorted((departed, arrived)) == last[:2]
        assert x == observed
        rng = random.Random(72324000 + i)
        assert len(r['holding']) == 64
        for row in r['holding']:
            assert row[:2] == sorted(r['bond']) and row[2] == rng.randrange(2)
            x = endpoint(x, row)
            calls += 1
            a, b, _, j = row
            observed[a] -= j
            observed[b] += j
            assert x == observed and sum(x) == 1 and not any(x[12:])
        check_handoff(x, observed, r['final_handoff'], {
            'basis': r['basis_handoff']['prefix_sha256'],
            'discovery': r['discovery'], 'holding': r['holding']})
        assert r['outcome'] == 'selected_and_continued'
    main = [r for r in data['trials'] if r['selection_law'] == 'native']
    assert len(data['bond_certificates']) == 36
    for cert, r in zip(data['bond_certificates'], main):
        assert cert['edge'] == sorted(r['bond'])
        a, b = cert['edge']
        distributions = {law: {start: Counter() for start in (a, b)} for law in ('native', 'hold', 'swap')}
        assert len(cert['branches']) == 12
        for law, start, c, j, target in cert['branches']:
            x = [int(i == start) for i in range(240)]
            if law == 'native':
                y = endpoint(x, [a, b, c, j])
                calls += 1
            else:
                expected = 0 if law == 'hold' else x[a] - x[b]
                assert j == expected
                y = x[:]
                y[a] -= j
                y[b] += j
            assert y.index(1) == target and sum(y) == 1
            distributions[law][start][target] += 1
        distances = {law: str(sum(F(abs(d[a][k] - d[b][k]), 4) for k in (a, b)))
                     for law, d in distributions.items()}
        assert distances == cert['total_variation'] == {'native': '0', 'hold': '1', 'swap': '1'}
    d = data['detector']
    assert d['initial_audit_only'] == [2, -1] + [0] * 238
    _, decoded = check_probe(d['initial_audit_only'], d['census'], 72325000)
    assert decoded[0] == d['initial_audit_only']
    _, decoded = check_probe(d['initial_audit_only'], data['forced_census'], 72325001, cap=1)
    assert decoded is None and len(data['forced_census']) == 1
    assert len(data['forced_selector']) == 80
    x = d['initial_audit_only'][:]
    for row in data['forced_selector']:
        assert row == [12, 13, 1, 0]
        x = endpoint(x, row)
        calls += 1
    assert calls == receipt['literal_native_calls']
    assert rejected == receipt['selector_rejected_attempts']
    assert receipt['control_body_events'] == 36 * 1024
    assert receipt['outcomes'] == {'selected_and_continued': 36, 'unsettled_at_horizon': 36}
    assert receipt['selected_handoff_addresses'] == [r['basis_handoff']['occupied'] for r in main]
    assert all(receipt['gates'].values()) and receipt['all_gates_pass']
    assert receipt['source240_actual_not_difference'] is True
    assert receipt['coherent_state_selection_claimed'] is receipt['physical_spin_or_energy_claimed'] is receipt['full_row_solved'] is False
    # Integrity mutation reaches the event checker, not merely a checksum test.
    initial = data['trials'][0]['initial_audit_only']
    wrong = data['trials'][0]['selection'][0]['body'][:]
    wrong[3] += 1
    try:
        endpoint(initial, wrong)
    except AssertionError:
        mutation_rejected = True
    else:
        raise AssertionError('Corrupt current was not rejected')
    return {'verified': True, 'producer_imports': [], 'native_calls_accounted': calls,
            'complete_trials': 72, 'successful_continuations': 36,
            'selector_rejections_retained': rejected, 'current_mutation_rejected': mutation_rejected,
            'record_sha256': hashlib.sha256(raw).hexdigest()}


if __name__ == '__main__':
    print(json.dumps(verify(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent), indent=2))
