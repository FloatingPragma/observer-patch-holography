"""Independent full-criterion row1 review; writes only this review directory.

The candidate's topology routines are not used for the independent integer
calculation. Its public admission facade is the implementation under test.
Frozen complete-law evidence is consumed and checked, not called a fresh run.
"""
from collections import deque
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ROW = ROOT / 'openwave-tuning/row_completion/ROW_01'
PROVIDER = ROOT / 'reverse-engineering-reality/code/openwave_row1_integrity'
AUDIT = ROOT / 'openwave-tuning/next/CHAIN3/row1_integrity_audit'
sys.path[:0] = [str(ROW), str(PROVIDER), str(AUDIT)]
import reviewed
import integrity


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def pin_inputs():
    pins = {}
    for name, expected in (
        ('MANIFEST.json', 'c4f6db87b1d2c7ddd4caa4c8ac5a95f09c33679f060d50845aa0bf42f5e54d12'),
        ('REVIEW_FIX_MANIFEST.json', 'f695ee571e25cfe59e808b9606a2b788aa232623f12dc33f0f4cac9db66c47d3'),
    ):
        p = ROW / name
        assert sha(p) == expected, name
        pins[str(p.relative_to(ROOT))] = expected
        data = json.loads(p.read_text())
        for key in ('files', 'source_sha256'):
            for path, h in data.get(key, {}).items():
                assert sha(ROOT / path) == h, path
                pins[path] = h
    # The full dependency audit includes all normalized instrument laws.
    m = json.loads((AUDIT / 'MANIFEST.json').read_text())
    for path, h in m['files_sha256'].items():
        p = AUDIT / path
        assert sha(p) == h, path
        pins[str(p.relative_to(ROOT))] = h
    return pins


def bareiss(matrix):
    """Fraction-free integer determinant, independent of candidate elimination."""
    a = [list(row) for row in matrix]
    previous = sign = 1
    for k in range(len(a) - 1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        p = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = p * a[i][j] - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = p
    return sign * a[-1][-1]


def independent_geometry(edges):
    index = {edge: i for i, edge in enumerate(edges)}
    assert len(index) == 30
    faces = [f for f in combinations(range(12), 3)
             if all(e in index for e in combinations(f, 2))]
    assert len(faces) == 20
    unsigned = []
    for a, b, c in faces:
        col = [0] * 30
        for x, y, s in ((a, b, 1), (b, c, 1), (a, c, -1)):
            col[index[x, y]] = s
        unsigned.append(col)
    signs = {0: 1}
    queue = deque([0])
    while queue:
        f = queue.popleft()
        for e, value in enumerate(unsigned[f]):
            if not value:
                continue
            neighbours = [g for g in range(20) if unsigned[g][e] and g != f]
            assert len(neighbours) == 1
            g = neighbours[0]
            want = -signs[f] * value * unsigned[g][e]
            if g in signs:
                assert signs[g] == want
            else:
                signs[g] = want
                queue.append(g)
    assert len(signs) == 20
    A = [[signs[f] * unsigned[f][e] for f in range(20)] for e in range(30)]
    B = [[int(v == a) - int(v == b) for a, b in edges] for v in range(12)]
    assert all(sum(B[v][e] * A[e][f] for e in range(30)) == 0
               for v in range(12) for f in range(20))
    # A deliberately different spanning tree, rooted at port11 and reverse
    # edge order, tests the lattice rather than the producer's chosen minor.
    tree = []
    reached = {11}
    while len(reached) < 12:
        for e in reversed(range(30)):
            a, b = edges[e]
            if (a in reached) != (b in reached):
                tree.append(e)
                reached.update((a, b))
                break
        else:
            raise AssertionError('disconnected carrier')
    chords = [e for e in range(30) if e not in tree]
    d1 = bareiss([[B[v][e] for e in tree] for v in range(11)])
    d2 = bareiss([[A[e][f] for f in range(1, 20)] for e in chords])
    assert abs(d1) == abs(d2) == 1
    # Every port link is a single closed five-edge curve on the dual graph.
    dual = [(row.index(-1), row.index(1)) for row in A]
    for v in range(12):
        arcs = [(a, b) if B[v][e] == 1 else (b, a)
                for e, (a, b) in enumerate(dual) if B[v][e]]
        step = dict(arcs)
        assert len(arcs) == len(step) == 5 and set(step) == set(step.values())
        p = start = min(step)
        visited = set()
        for _ in range(5):
            visited.add(p)
            p = step[p]
        assert p == start and len(visited) == 5
    return B, A, {'tree_minor': d1, 'face_minor': d2, 'root': 11,
                  'ports': 12, 'seams': 30, 'faces': 20, 'dual_pentagons': 12}


def states(total):
    for occupied in combinations_with_replacement(range(12), total):
        yield tuple(occupied.count(v) for v in range(12))


def native(n, edge, c):
    a, b = edge
    lower, odd = divmod(n[a] + n[b], 2)
    result = list(n)
    result[a], result[b] = lower + c * odd, lower + (1 - c) * odd
    return tuple(result), n[a] - result[a]


def replay_source(G, claim, edges):
    n = tuple(G)
    F = [0] * 30
    for records in (claim['source_records'], claim['sweep_records']):
        for record in records:
            e, c = record['seam'], record['coin']
            a, b = edges[e]
            if record['kind'] == 'E':
                assert record['contrast'] == n[a] - n[b]
            after, j = native(n, edges[e], c)
            assert record['current'] == j
            F[e] += j
            n = after
    return n, F


def main():
    pins = pin_inputs()
    edges = tuple(integrity.EDGES)
    B, A, geometry = independent_geometry(edges)
    branches = pairs = 0
    for total in (1, 2):
        configurations = list(states(total))
        for n in configurations:
            for e in edges:
                for c in (0, 1):
                    after, j = native(n, e, c)
                    actual, (_, actual_j) = integrity.native(n, e, c)
                    assert (after, j) == (actual, actual_j)
                    assert min(after) >= 0 and sum(after) == total
                    if total == 1:
                        assert after.count(1) == 1
                    branches += 1
            for G in configurations:
                delta = [g - v for g, v in zip(G, n)]
                assert sum(delta) == 0
                assert sum(x for x in delta if x > 0) <= total
                assert sum(x for x in delta if x < 0) >= -total
                assert len({w % 6 for w in range(-total, total + 1)}) == 2 * total + 1
                pairs += 1
    print('Independent lattice and native sectors passed', branches, pairs, flush=True)
    examples = json.loads((ROW / 'EXAMPLES.json').read_text())
    vacuum = json.loads((ROW / 'VACUUM.json').read_text())
    cases = examples + [{'G': vacuum['G'], 'region': [0], 'claim': vacuum['claim'],
                         'topological_degree': 0, 'raw_pointer_record_degree': 0}]
    signed = set()
    outcomes = {}
    mutations = 0
    readouts = []
    cycle = integrity.meter.CYCLE
    cycle_sign = {tuple(sorted((a, b))): (1 if a < b else -1)
                  for a, b in zip(cycle, cycle[1:] + cycle[:1])}
    for row in cases:
        claim, G, region = row['claim'], row['G'], row['region']
        n, F = replay_source(G, claim, edges)
        degree = sum(B[v][e] * F[e] for v in region for e in range(30))
        assert degree == sum(G[v] - n[v] for v in region) == row['topological_degree']
        verdict = reviewed.admit_and_decode(claim, G, region, 2)
        assert verdict['admitted'], verdict
        signed.add(degree)
        outcomes[verdict['outcome']] = outcomes.get(verdict['outcome'], 0) + 1
        steps = sum(cycle_sign[tuple(r['seam'])] * r['current']
                    for r in claim['meter']['records'] if r['stage'] == 'measure')
        if claim['meter']['success']:
            assert steps == 6 * degree and verdict['readout'] == -degree
        else:
            assert verdict['readout'] is None
        assert all(Fraction(x) in (0, Fraction(1, 2**verdict['record_attempts']))
                   for x in verdict['pointer_initial_classical_effect'])
        bad = deepcopy(claim)
        bad['meter']['final_loads'] = {'arbitrary': 'not a state'}
        assert reviewed.admit_and_decode(bad, G, region, 2) == verdict
        mutations += 1
        readouts.append({'degree': degree, 'readout': verdict['readout'],
                         'outcome': verdict['outcome'], 'native_attempts': verdict['record_attempts']})
    assert signed == {-2, -1, 0, 1, 2}
    for change in ('empty', 'wrong_charge', 'fractional', 'false_contrast', 'wrong_latch'):
        good = next(r for r in cases if r['claim']['meter']['success'] and r['topological_degree'] == 1)
        bad = deepcopy(good['claim'])
        if change == 'empty': bad['meter']['records'] = []
        if change == 'wrong_charge': bad['decoded_relative_charge'] = Fraction(1, 2)
        if change == 'fractional': bad['source_records'][0]['current'] = Fraction(1, 2)
        if change == 'false_contrast': bad['sweep_records'][0]['contrast'] += 1
        if change == 'wrong_latch': bad['charge_read_at'] = 'L0'
        assert not reviewed.admit_and_decode(bad, good['G'], good['region'], 2)['admitted'], change
        mutations += 1
    # The half-current violates the shared endpoint of the specified local
    # circle map: exp(2pi i/2)=-1, whereas both dual vertex values are +1.
    # Equal daughter integer degrees obey 2d=+/-1, which is impossible.
    assert all(2 * d not in (-1, 1) for d in range(-2, 3))
    laws = json.loads((PROVIDER / 'RESULTS.json').read_text())['complete_bounded_laws']
    opposing = json.loads((AUDIT / 'RESULT.json').read_text())['complete_joint_laws']
    assert len(laws['rows']) == 768 and digest(laws['rows']) == opposing['row_list_sha256']
    for law in laws['rows']:
        assert sum(Fraction(x) for x in law['outcome_mass'].values()) == 1
    assert opposing['every_frozen_row_exactly_equal'] and opposing['all_failures_retained']
    assert pin_inputs() == pins
    save('READOUTS.json', readouts)
    result = {'schema': 'oph.row1.full_criterion.independent_review.v1',
              'reviewer': 'Codex root full-criterion reviewer',
              'producer_independent': True,
              'independence_scope': 'Independent of ROW_01 producer and auditor-root review. Earlier instrument suppliers are reused with explicit pins; no claim of independent authorship of all dependencies.',
              'candidate_id': 'literal-native-relative-spatial-charge-AU5-record-only-v1',
              'criterion_sha256': '3af50af0955c64a0aa2019fa84af3cf0dd36e6c63a5abc76eb2c60c9b413cb6f',
              'base_manifest_sha256': sha(ROW / 'MANIFEST.json'),
              'correction_manifest_sha256': sha(ROW / 'REVIEW_FIX_MANIFEST.json'),
              'full_criterion_met': True,
              'verdict': 'PASS on the declared bounded signed-relative-charge domain',
              'geometry': geometry, 'native_branches': branches, 'preparation_state_pairs': pairs,
              'full_received_instruments_replayed': len(cases), 'signed_degrees': sorted(signed),
              'outcomes': outcomes, 'metadata_and_rejection_controls': mutations,
              'complete_probability_laws_consumed': 768, 'fresh_complete_law_rerun': False,
              'conditions': ['Known actual preparation G and total one or two',
                  'Faithful received closing E records; complete emitted histories for fine likelihoods',
                  'Explicit local circle-path interpretation on the actual dual surface',
                  'Charge refers to L1 after the disturbing closing sweep',
                  'All native pointer failures remain unconditional outcomes'],
              'not_established': ['Official OpenWave award', 'An electron mass or spin identification',
                  'Unrestricted mod6 charge decoding', 'Authenticity against replacement of the entire channel',
                  'Equivalence to the separate coherent field source'],
              'source_sha256': pins, 'provider_bytes_unchanged': True}
    save('REVIEW.json', result)
    print(json.dumps({k: result[k] for k in ('full_criterion_met', 'native_branches',
        'preparation_state_pairs', 'full_received_instruments_replayed', 'metadata_and_rejection_controls')}), flush=True)


if __name__ == '__main__':
    main()
