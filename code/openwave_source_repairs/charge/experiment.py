"""Exact native paired-law extension and record-only selection/response join."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import random
from native_event import native

HERE = Path(__file__).resolve().parent
EDGES = tuple(map(tuple, json.loads((HERE / 'seams.json').read_text())))
N, M = 12, len(EDGES)
CALLS = 0


def event(x, edge, coin, law='native'):
    global CALLS
    if law == 'native':
        CALLS += 1
        return native(x, edge, coin)
    if law != 'swap_hold':
        raise ValueError(law)
    a, b = edge
    y = list(x)
    j = coin * (x[a] - x[b])
    y[a] -= j
    y[b] += j
    return tuple(y), (coin, j)


def tree():
    seen, queue, edges = {0}, [0], []
    for a in queue:
        for u, v in EDGES:
            b = v if a == u else u if a == v else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                edges.append((a, b))
    assert len(seen) == N
    return tuple(edges)


TREE = tree()


def pair_body(x, y, seed, law, horizon=512):
    rng, tape = random.Random(seed), []
    for _ in range(horizon):
        a, b = EDGES[rng.randrange(M)]
        c = rng.randrange(2)
        x, (_, jx) = event(x, (a, b), c, law)
        y, (_, jy) = event(y, (a, b), c, law)
        tape.append([a, b, c, jx, jy])
    return x, y, tape


def pair_probe(x, y, seed, cap=32):
    rng, tape = random.Random(seed), []
    for a, b in TREE:
        first = None
        for _ in range(cap):
            c = rng.randrange(2)
            x, (_, jx) = event(x, (a, b), c)
            y, (_, jy) = event(y, (a, b), c)
            tape.append([a, b, c, jx, jy])
            if first is None:
                first = c
            elif c != first:
                break
        else:
            return x, y, tape
    return x, y, tape


def read_probe(tape, cap=32):
    """No source state, initial total, or preparation identity is an input."""
    delta = [[0] * N, [0] * N]
    relative = [[None] * N, [None] * N]
    relative[0][0] = relative[1][0] = 0
    cursor = 0
    for a, b in TREE:
        correction = [d[a] - d[b] for d in delta]
        first = None
        for _ in range(cap):
            if cursor == len(tape):
                return None
            aa, bb, c, jx, jy = tape[cursor]
            cursor += 1
            if (aa, bb) != (a, b) or c not in (0, 1):
                raise ValueError('Malformed probe address')
            js = (jx, jy)
            for r in range(2):
                if type(js[r]) is not int:
                    raise ValueError('Noninteger source current')
                delta[r][a] -= js[r]
                delta[r][b] += js[r]
            if first is None:
                first, first_js = c, js
            elif c != first:
                for r in range(2):
                    contrast = 2 * first_js[r] + js[r] - correction[r]
                    relative[r][b] = relative[r][a] - contrast
                break
        else:
            if cursor != len(tape):
                raise ValueError('Data after timeout')
            return None
    if cursor != len(tape):
        raise ValueError('Extra probe records')
    return relative, delta


def backwards(tape, relative, side):
    x = relative[:]
    result = [x[:]]
    for a, b, _, jx, jy in reversed(tape):
        j = (jx, jy)[side]
        x[a] += j
        x[b] -= j
        result.append(x[:])
    return result[::-1]


def neutral_difference(rx, ry):
    d = [a - b for a, b in zip(rx, ry)]
    offset = F(sum(d), N)
    q = [F(v) - offset for v in d]
    if any(v.denominator != 1 for v in q):
        raise ValueError('Equal-total integer-source promise inconsistent')
    return [int(v) for v in q]


def read_response(tape, probe, cap=32):
    read = read_probe(probe, cap)
    if read is None:
        return None
    relative, _ = read
    paths = [backwards(tape, relative[r], r) for r in range(2)]
    qs = [neutral_difference(x, y) for x, y in zip(*paths)]
    tau = next((i for i, q in enumerate(qs) if not any(q)), None)
    flux = [0] * M
    if tau is not None:
        for a, b, _, jx, jy in tape[:tau]:
            flux[EDGES.index((a, b))] += jx - jy
        divergence = [0] * N
        for (a, b), j in zip(EDGES, flux):
            divergence[a] += j
            divergence[b] -= j
        if divergence != qs[0]:
            raise ValueError('Record-only Gauss failure')
    return {'initial_q': qs[0], 'terminal_q': qs[-1],
            'first_synchronization': tau, 'stopped_flux': flux if tau is not None else None,
            'post_synchronization_exits': 0 if tau is None else
                sum(any(q) for q in qs[tau:])}


def preparations():
    for p in range(N):
        for family in ('double', 'signed', 'shifted', 'dense'):
            x = [2 if family == 'shifted' else 5 if family == 'dense' else 0] * N
            if family == 'signed':
                x[p], x[(p + 1) % N] = 8, -6
            else:
                x[p] += -2 if family == 'dense' else 2
            yield family, tuple(x)


def pipeline(initial, index, selection_law='native', response_law='native'):
    x, y, body = pair_body(initial, initial, 72311000 + index, selection_law)
    x, y, selection_probe = pair_probe(x, y, 72311500 + index)
    r = {'initial_audit_only': initial, 'index': index,
         'selection_law': selection_law, 'response_law': response_law,
         'selection_body': body, 'selection_probe': selection_probe,
         'intervention': [], 'response_body': [], 'response_probe': [],
         'selection_spans': None, 'first_selected': None, 'response': None}
    read = read_probe(selection_probe)
    if read is None:
        r['outcome'] = 'selection_probe_timeout'
        return r
    relative, delta = read
    path = backwards(body, relative[0], 0)
    spans = [max(z) - min(z) for z in path]
    r['selection_spans'] = spans
    r['first_selected'] = next((i for i, span in enumerate(spans) if span <= 1), None)
    if spans[-1] > 1:
        r['outcome'] = 'unsettled_at_horizon'
        return r
    endpoint = [a + b for a, b in zip(relative[0], delta[0])]
    edge = next((e for e in EDGES if abs(endpoint[e[0]] - endpoint[e[1]]) == 1), None)
    if edge is None:
        r['outcome'] = 'uniform_no_relative_excitation'
        return r
    rng = random.Random(72312000 + index)
    for _ in range(32):
        c = rng.randrange(2)
        x, (_, j) = event(x, edge, c)
        r['intervention'].append([*edge, c, j])
        if j:
            break
    else:
        r['outcome'] = 'intervention_timeout'
        return r
    x, y, r['response_body'] = pair_body(x, y, 72313000 + index, response_law, horizon=4096)
    x, y, r['response_probe'] = pair_probe(x, y, 72314000 + index)
    r['response'] = read_response(r['response_body'], r['response_probe'])
    if r['response'] is None:
        r['outcome'] = 'response_probe_timeout'
    else:
        r['outcome'] = ('synchronized' if r['response']['first_synchronization'] is not None
                        else 'response_horizon_timeout')
    return r


def local_closure():
    evidence = []
    for xa, xb, ya, yb in product(range(-4, 5), repeat=4):
        x, y, q = (xa, xb), (ya, yb), (xa - ya, xb - yb)
        paired, direct = [], []
        for c in (0, 1):
            xp, (_, jx) = event(x, (0, 1), c)
            yp, (_, jy) = event(y, (0, 1), c)
            qp, (_, jq) = event(q, (0, 1), c)
            paired.append([xp[0] - yp[0], xp[1] - yp[1], jx - jy])
            direct.append([*qp, jq])
        if sorted(paired) != sorted(direct):
            raise ValueError('All-integer paired closure counterexample')
        evidence.append([xa, xb, ya, yb, paired, direct])
    return evidence


def invert(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [v / scale for v in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [v - scale * w for v, w in zip(a[i], a[j])]
    return [row[n:] for row in a]


def rational_response():
    lap = [[0] * N for _ in range(N)]
    for a, b in EDGES:
        lap[a][a] += 1
        lap[b][b] += 1
        lap[a][b] -= 1
        lap[b][a] -= 1
    inverse = invert([row[:-1] for row in lap[:-1]])

    def h(q):
        phi = [sum(a * b for a, b in zip(row, q[:-1])) for row in inverse] + [F(0)]
        assert [sum(a * b for a, b in zip(row, phi)) for row in lap] == list(q)
        return [phi[a] - phi[b] for a, b in EDGES]

    zero = (0,) * N
    hs, branch_rows, responses = {zero: [F(0)] * M}, [], []
    for p in range(N):
        for r in range(N):
            if p != r:
                q = tuple((i == p) - (i == r) for i in range(N))
                hs[q] = h(q)
    for q, value in hs.items():
        if q == zero:
            continue
        mean = [F(0)] * M
        for e, edge in enumerate(EDGES):
            for c in (0, 1):
                qp, (_, j) = event(q, edge, c)
                reward = hs[qp][:]
                reward[e] += j
                mean = [a + b / (2 * M) for a, b in zip(mean, reward)]
                branch_rows.append([q, e, c, j, qp])
        if mean != value:
            raise ValueError('Stopped reward equation residual')
        responses.append([q, [str(v) for v in value]])
    return {'laplacian': lap, 'grounded_inverse': [[str(v) for v in row] for row in inverse],
            'native_branches': branch_rows, 'mean_responses': responses}


def audit_events(r):
    """Producer audit only: separate from the observer's record decoder."""
    x = y = tuple(r['initial_audit_only'])
    for key, law in [('selection_body', r['selection_law']),
                     ('selection_probe', 'native'), ('intervention', 'native'),
                     ('response_body', r['response_law']), ('response_probe', 'native')]:
        for row in r[key]:
            a, b, c, jx = row[:4]
            x, (_, expected) = event(x, (a, b), c, law)
            if jx != expected:
                raise ValueError('Current record disagrees with source')
            if key != 'intervention':
                y, (_, expected) = event(y, (a, b), c, law)
                if row[4] != expected:
                    raise ValueError('Reference current disagrees with source')


def calculate():
    global CALLS
    CALLS = 0
    local = local_closure()
    response = rational_response()
    histories = []
    for i, (family, initial) in enumerate(preparations()):
        for slaw, rlaw in [('native', 'native'), ('swap_hold', 'native'),
                           ('native', 'swap_hold')]:
            r = pipeline(initial, i, slaw, rlaw)
            r['family'] = family
            histories.append(r)
    main = [r for r in histories if r['selection_law'] == r['response_law'] == 'native']
    selection_control = [r for r in histories if r['selection_law'] == 'swap_hold']
    response_control = [r for r in histories if r['response_law'] == 'swap_hold']
    detectors = []
    for x in ([2] + [0] * 11, [1, 1] + [0] * 10):
        _, _, tape = pair_probe(x, x, 72315900 + len(detectors))
        decoded = read_probe(tape)
        detectors.append({'initial_audit_only': x, 'probe': tape,
                          'span': max(decoded[0][0]) - min(decoded[0][0])})
    _, _, timeout = pair_probe([2] + [0] * 11, [2] + [0] * 11, 72315999, cap=1)
    corrupted = json.loads(json.dumps(main[0]))
    corrupted['selection_body'][0][3] += 1
    calls_before_audit = CALLS
    try:
        audit_events(corrupted)
    except ValueError:
        rejected = True
    else:
        rejected = False
    audit_calls = CALLS - calls_before_audit
    gates = {
        'local_all_pairs_equal': len(local) == 6561,
        'exact_rational_response': len(response['native_branches']) == 7920,
        'native_selection_48': all(r['first_selected'] is not None for r in main),
        'no_native_selection_escape': all(max(r['selection_spans'][r['first_selected']:]) <= 1 for r in main),
        'control_selection_0': all(r['first_selected'] is None for r in selection_control),
        'native_interventions_unit': all(abs(r['intervention'][-1][3]) == 1 for r in main),
        'intervention_sources_primitive': all(sorted(r['response']['initial_q']) == [-1] + [0] * 10 + [1] for r in main),
        'native_responses_48': all(r['outcome'] == 'synchronized' for r in main),
        'control_responses_0': all(r['outcome'] == 'response_horizon_timeout' for r in response_control),
        'no_response_exit': all(r['response']['post_synchronization_exits'] == 0 for r in main),
        'double_detector': [r['span'] for r in detectors] == [2, 1],
        'timeout_not_success': read_probe(timeout, cap=1) is None,
        'corrupt_current_rejected': rejected,
    }
    evidence = {'local': local, 'rational_response': response, 'histories': histories,
                'detectors': detectors, 'forced_timeout': timeout,
                'integrity_mutation': {'record': corrupted, 'rejected': rejected}}
    receipt = {'candidate': 'A-native-multiload-gauss-selection-v2', 'gates': gates,
               'all_gates_pass': all(gates.values()), 'literal_native_calls': CALLS,
               'integrity_audit_native_calls': audit_calls,
               'all_outcomes': dict(Counter(r['outcome'] for r in histories)),
               'local_endpoint_quadruples': len(local), 'primitive_source_branches': 7920,
               'native_synchronization_times': [r['response']['first_synchronization'] for r in main],
               'native_selection_times': [r['first_selected'] for r in main],
               'full_row_solved': False, 'g6_confirmed': False,
               'physical_electric_identification': False, 'source_pin_changed': False}
    return receipt, evidence


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def verify_freeze():
    prereg = json.loads((HERE / 'PREREG.json').read_text())
    for name, digest in prereg['sha256'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen input changed: ' + name)
    return hashlib.sha256((HERE / 'PREREG.json').read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    prereg = verify_freeze()
    receipt, evidence = calculate()
    raw = gzip.compress(encoded(evidence), mtime=0)
    receipt['evidence_sha256'] = hashlib.sha256(raw).hexdigest()
    receipt['prereg_sha256'] = prereg
    if args.record:
        for name, data in [('RECEIPT.json', encoded(receipt)), ('EVIDENCE.json.gz', raw)]:
            with (HERE / name).open('xb') as handle:
                handle.write(data)
    else:
        if json.loads((HERE / 'RECEIPT.json').read_text()) != receipt:
            raise ValueError('Receipt does not replay')
        if (HERE / 'EVIDENCE.json.gz').read_bytes() != raw:
            raise ValueError('Full evidence does not replay byte for byte')
    print(json.dumps(receipt, indent=2))
    if not receipt['all_gates_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
