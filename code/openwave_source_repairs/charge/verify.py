"""Retrospective independent implementation; imports no producer modules.

Written by the same Codex worker: implementation independence, not other-group
review. Full native events are checked by floor/ceiling endpoint allocation.
"""
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import sys


def step(x, edge, coin, law):
    a, b = edge
    y = list(x)
    if law == 'native':
        total = x[a] + x[b]
        lo, hi = total // 2, -((-total) // 2)
        y[a], y[b] = (lo, hi) if coin == 0 else (hi, lo)
    else:
        assert law == 'swap_hold'
        if coin:
            y[a], y[b] = x[b], x[a]
    return y, x[a] - y[a]


def read_probe(tape, tree, cap=32):
    n = 12
    cumulative, known = [[0] * n, [0] * n], [{0: 0}, {0: 0}]
    cursor = 0
    for a, b in tree:
        initial_delta = [v[a] - v[b] for v in cumulative]
        if cursor == len(tape):
            return None
        first = tape[cursor]
        complete = False
        for _ in range(cap):
            assert cursor < len(tape)
            aa, bb, c, jx, jy = tape[cursor]
            assert (aa, bb) == (a, b) and c in (0, 1)
            cursor += 1
            for side, j in enumerate((jx, jy)):
                assert type(j) is int
                cumulative[side][a] -= j
                cumulative[side][b] += j
            if c != first[2]:
                for side, j in enumerate((jx, jy)):
                    d = 2 * first[3 + side] + j - initial_delta[side]
                    known[side][b] = known[side][a] - d
                complete = True
                break
        if not complete:
            assert cursor == len(tape)
            return None
    assert cursor == len(tape)
    return [[v[i] for i in range(n)] for v in known], cumulative


def relative(x):
    return [v - x[0] for v in x]


def reverse(tape, last, side):
    x, path = last[:], [last[:]]
    for row in tape[::-1]:
        a, b, _, jx, jy = row
        j = (jx, jy)[side]
        x[a] += j
        x[b] -= j
        path.append(x[:])
    return path[::-1]


def verify(directory):
    p = Path(directory)
    edges = [tuple(e) for e in json.loads((p / 'seams.json').read_text())]
    n, m = 12, len(edges)
    queue, seen, tree = [0], {0}, []
    for a in queue:
        for u, v in edges:
            b = v if a == u else u if a == v else None
            if b is not None and b not in seen:
                queue.append(b)
                seen.add(b)
                tree.append((a, b))
    assert len(tree) == 11 and len(edges) == 30
    receipt = json.loads((p / 'RECEIPT.json').read_text())
    data_bytes = (p / 'EVIDENCE.json.gz').read_bytes()
    assert hashlib.sha256(data_bytes).hexdigest() == receipt['evidence_sha256']
    data = json.loads(gzip.decompress(data_bytes))
    pre = json.loads((p / 'PREREG.json').read_text())
    assert hashlib.sha256((p / 'PREREG.json').read_bytes()).hexdigest() == receipt['prereg_sha256']
    for name, expected in pre['sha256'].items():
        assert hashlib.sha256((p / name).read_bytes()).hexdigest() == expected
    native_calls = 0
    assert len(data['local']) == 9**4
    for row, values in zip(data['local'], product(range(-4, 5), repeat=4)):
        assert row[:4] == list(values)
        xa, xb, ya, yb = values
        paired, direct = [], []
        for c in (0, 1):
            xp, jx = step([xa, xb], (0, 1), c, 'native')
            yp, jy = step([ya, yb], (0, 1), c, 'native')
            qp, jq = step([xa - ya, xb - yb], (0, 1), c, 'native')
            paired.append([xp[0] - yp[0], xp[1] - yp[1], jx - jy])
            direct.append(qp + [jq])
            native_calls += 3
        assert row[4] == paired and row[5] == direct and sorted(paired) == sorted(direct)
    rational = data['rational_response']
    lap = [[0] * n for _ in range(n)]
    for a, b in edges:
        lap[a][a] += 1
        lap[b][b] += 1
        lap[a][b] -= 1
        lap[b][a] -= 1
    assert lap == rational['laplacian']
    inv = [[F(v) for v in row] for row in rational['grounded_inverse']]
    for i in range(n - 1):
        for j in range(n - 1):
            assert sum(lap[i][k] * inv[k][j] for k in range(n - 1)) == (i == j)
    responses = {tuple(q): [F(v) for v in v] for q, v in rational['mean_responses']}
    assert len(responses) == 132
    responses[(0,) * n] = [F(0)] * m
    for q, vals in responses.items():
        phi = [sum(row[i] * q[i] for i in range(n - 1)) for row in inv] + [F(0)]
        assert [phi[a] - phi[b] for a, b in edges] == vals
        assert [sum(lap[i][j] * phi[j] for j in range(n)) for i in range(n)] == list(q)
    reward = {q: [F(0)] * m for q in responses if any(q)}
    visited = set()
    for q, e, c, j, qp in rational['native_branches']:
        key = (tuple(q), e, c)
        assert key not in visited
        visited.add(key)
        actual, current = step(q, edges[e], c, 'native')
        native_calls += 1
        assert current == j and actual == qp
        for f in range(m):
            reward[tuple(q)][f] += (responses[tuple(qp)][f] + (j if e == f else 0)) / (2 * m)
    assert len(visited) == 7920 and all(reward[q] == responses[q] for q in reward)

    def check_history(r):
        nonlocal native_calls
        x = list(r['initial_audit_only'])
        y = x[:]
        selection_path = [x[:]]
        assert len(r['selection_body']) == 512
        for a, b, c, jx, jy in r['selection_body']:
            x, j = step(x, (a, b), c, r['selection_law'])
            y, k = step(y, (a, b), c, r['selection_law'])
            native_calls += 2 * (r['selection_law'] == 'native')
            assert (j, k) == (jx, jy)
            selection_path.append(x[:])
        pre = [x[:], y[:]]
        read = read_probe(r['selection_probe'], tree)
        assert read is not None
        for side in range(2):
            assert read[0][side] == relative(pre[side])
        decoded_path = reverse(r['selection_body'], read[0][0], 0)
        assert [relative(z) for z in decoded_path] == [relative(z) for z in selection_path]
        spans = [max(z) - min(z) for z in decoded_path]
        assert spans == r['selection_spans']
        first = next((i for i, span in enumerate(spans) if span <= 1), None)
        assert first == r['first_selected']
        for a, b, c, jx, jy in r['selection_probe']:
            x, j = step(x, (a, b), c, 'native')
            y, k = step(y, (a, b), c, 'native')
            native_calls += 2
            assert (j, k) == (jx, jy)
        if r['selection_law'] == 'swap_hold':
            assert first is None and r['outcome'] == 'unsettled_at_horizon'
            assert r['intervention'] == r['response_body'] == r['response_probe'] == []
            return
        assert first is not None and max(spans[first:]) <= 1
        endpoint = [a + b for a, b in zip(read[0][0], read[1][0])]
        chosen = next(e for e in edges if abs(endpoint[e[0]] - endpoint[e[1]]) == 1)
        assert len(r['intervention']) in range(1, 33)
        for i, (a, b, c, jx) in enumerate(r['intervention']):
            assert (a, b) == chosen
            x, j = step(x, chosen, c, 'native')
            native_calls += 1
            assert j == jx
            assert (j != 0) == (i == len(r['intervention']) - 1)
        assert abs(j) == 1
        qs = [[a - b for a, b in zip(x, y)]]
        response_path = [[x[:]], [y[:]]]
        assert len(r['response_body']) == 4096
        for a, b, c, jx, jy in r['response_body']:
            x, j = step(x, (a, b), c, r['response_law'])
            y, k = step(y, (a, b), c, r['response_law'])
            native_calls += 2 * (r['response_law'] == 'native')
            assert (j, k) == (jx, jy)
            qs.append([a - b for a, b in zip(x, y)])
            response_path[0].append(x[:])
            response_path[1].append(y[:])
        terminal = read_probe(r['response_probe'], tree)
        assert terminal is not None
        reconstructed = [reverse(r['response_body'], terminal[0][s], s) for s in range(2)]
        for side in range(2):
            assert [relative(z) for z in reconstructed[side]] == [relative(z) for z in response_path[side]]
        for xread, yread, q in zip(*reconstructed, qs):
            raw = [a - b for a, b in zip(xread, yread)]
            assert [F(a) - F(sum(raw), n) for a in raw] == q
        for a, b, c, jx, jy in r['response_probe']:
            x, j = step(x, (a, b), c, 'native')
            y, k = step(y, (a, b), c, 'native')
            native_calls += 2
            assert (j, k) == (jx, jy)
        tau = next((i for i, q in enumerate(qs) if not any(q)), None)
        assert r['response']['first_synchronization'] == tau
        assert r['response']['initial_q'] == qs[0]
        assert r['response']['terminal_q'] == qs[-1]
        assert sorted(qs[0]) == [-1] + [0] * 10 + [1]
        if r['response_law'] == 'native':
            assert tau is not None and r['outcome'] == 'synchronized'
            assert not any(any(q) for q in qs[tau:])
            flux = [0] * m
            for a, b, _, jx, jy in r['response_body'][:tau]:
                flux[edges.index((a, b))] += jx - jy
            assert flux == r['response']['stopped_flux']
            divergence = [0] * n
            for (a, b), j in zip(edges, flux):
                divergence[a] += j
                divergence[b] -= j
            assert divergence == qs[0]
        else:
            assert tau is None and r['outcome'] == 'response_horizon_timeout'

    assert len(data['histories']) == 144
    keys = {(r['index'], r['selection_law'], r['response_law']) for r in data['histories']}
    assert keys == {(i, a, b) for i in range(48) for a, b in
                    [('native', 'native'), ('swap_hold', 'native'), ('native', 'swap_hold')]}
    for r in data['histories']:
        paxis, kind = divmod(r['index'], 4)
        family = ('double', 'signed', 'shifted', 'dense')[kind]
        initial = [2 if kind == 2 else 5 if kind == 3 else 0] * n
        if kind == 1:
            initial[paxis], initial[(paxis + 1) % n] = 8, -6
        else:
            initial[paxis] += -2 if kind == 3 else 2
        assert r['family'] == family and r['initial_audit_only'] == initial
        check_history(r)
    for d in data['detectors']:
        x = y = d['initial_audit_only'][:]
        decoded = read_probe(d['probe'], tree)
        assert decoded[0][0] == relative(x)
        assert d['span'] == max(x) - min(x)
        for a, b, c, jx, jy in d['probe']:
            x, j = step(x, (a, b), c, 'native')
            y, k = step(y, (a, b), c, 'native')
            native_calls += 2
            assert (j, k) == (jx, jy)
    assert [d['span'] for d in data['detectors']] == [2, 1]
    assert read_probe(data['forced_timeout'], tree, cap=1) is None
    assert len(data['forced_timeout']) == 1
    a, b, c, j, k = data['forced_timeout'][0]
    assert step([2] + [0] * 11, (a, b), c, 'native')[1] == j == k
    native_calls += 2
    before = native_calls
    try:
        check_history(data['integrity_mutation']['record'])
    except AssertionError:
        pass
    else:
        raise AssertionError('Corruption was not detected')
    native_calls = before + 1  # producer aborts after the first incorrect X current
    assert native_calls == receipt['literal_native_calls'] == 555001
    assert all(receipt['gates'].values()) and receipt['all_gates_pass']
    assert receipt['full_row_solved'] is receipt['g6_confirmed'] is False
    assert receipt['physical_electric_identification'] is receipt['source_pin_changed'] is False
    outcomes = dict(Counter(r['outcome'] for r in data['histories']))
    assert outcomes == receipt['all_outcomes'] == {
        'synchronized': 48, 'unsettled_at_horizon': 48, 'response_horizon_timeout': 48}
    return {'verified': True, 'producer_imports': [], 'local_quadruples': 6561,
            'primitive_branches': 7920, 'complete_histories': 144,
            'literal_native_calls_accounted': native_calls, 'outcomes': outcomes,
            'evidence_sha256': hashlib.sha256(data_bytes).hexdigest()}


if __name__ == '__main__':
    print(json.dumps(verify(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent), indent=2))
