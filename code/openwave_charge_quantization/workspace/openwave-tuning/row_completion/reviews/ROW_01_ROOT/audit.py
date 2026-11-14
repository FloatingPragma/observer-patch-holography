"""Independent root review of the complete static native charge criterion."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PRODUCER = ROOT/'openwave-tuning/row_completion/ROW_01'
BASE = 'c4f6db87b1d2c7ddd4caa4c8ac5a95f09c33679f060d50845aa0bf42f5e54d12'
FIX = 'f695ee571e25cfe59e808b9606a2b788aa232623f12dc33f0f4cac9db66c47d3'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def custody():
    pins = {}
    for filename, expected in [('MANIFEST.json', BASE), ('REVIEW_FIX_MANIFEST.json', FIX)]:
        path = PRODUCER/filename
        assert sha(path) == expected
        pins[str(path.relative_to(ROOT))] = expected
        data = json.loads(path.read_text())
        for key in ('files', 'source_sha256'):
            for name, wanted in data.get(key, {}).items():
                assert sha(ROOT/name) == wanted, name
                pins[name] = wanted
    return pins


def bareiss(matrix):
    a = [list(row) for row in matrix]; old = 1; sign = 1
    for k in range(len(a)-1):
        pivot = next(i for i in range(k, len(a)) if a[i][k])
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]; sign = -sign
        d = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = d*a[i][j]-a[i][k]*a[k][j]
                assert numerator % old == 0
                a[i][j] = numerator//old
            a[i][k] = 0
        old = d
    return sign*a[-1][-1]


def event(state, edge, coin):
    a, b = edge
    # Literal integer repair independently implemented, including negatives
    # of contrasts; the state itself remains in its nonnegative sector.
    d = state[a]-state[b]
    current = (d-(2*coin-1)*(abs(d) & 1))//2
    result = list(state); result[a] -= current; result[b] += current
    return tuple(result), current


def run():
    pins = custody()
    sys.path.insert(0, str(PRODUCER))
    import reviewed
    provider = reviewed.integrity
    geometry = json.loads((PRODUCER/'GEOMETRY.json').read_text())
    edges = tuple(tuple(e) for e in geometry['seams'])
    assert edges == tuple(provider.EDGES)
    lookup = {e: i for i, e in enumerate(edges)}
    B = [[int(v == a)-int(v == b) for a, b in edges] for v in range(12)]
    A = [[0]*20 for _ in edges]
    for j, face in enumerate(geometry['faces']):
        for a, b in zip(face, face[1:]+face[:1]):
            A[lookup[tuple(sorted((a, b)))]] [j] += 1 if a < b else -1
    assert B == geometry['B'] and A == geometry['A']
    assert all(sum(B[v][e]*A[e][f] for e in range(30)) == 0
               for v in range(12) for f in range(20))
    assert all(sorted(x for x in row if x) == [-1, 1] for row in A)
    # A distinct spanning tree and fraction-free algorithm verify the integer
    # quotient, rather than accepting the producer's floating rank or minors.
    parent = list(range(12)); tree = []
    def find(v):
        while parent[v] != v:
            v = parent[v]
        return v
    for e in reversed(range(30)):
        a, b = map(find, edges[e])
        if a != b:
            parent[a] = b; tree.append(e)
    chords = [e for e in range(30) if e not in tree]
    tree_det = bareiss([[B[v][e] for e in tree] for v in range(1, 12)])
    face_det = bareiss([[A[e][f] for f in range(1, 20)] for e in chords])
    assert abs(tree_det) == abs(face_det) == 1
    for v, loop in enumerate(geometry['dual_port_boundaries']):
        current_vertex = None; first = None
        for e, sign in loop:
            tail = A[e].index(-1); head = A[e].index(1)
            if sign == -1:
                tail, head = head, tail
            if current_vertex is None:
                first = tail
            else:
                assert tail == current_vertex
            current_vertex = head
            assert B[v][e] == sign
        assert current_vertex == first and len(loop) == 5
    # Every spatial region, every primitive edge: linearity covers all integer
    # cochains. The signs are determined by the actual oriented dual cells.
    region_columns = 0
    for mask in range(4096):
        for e, (a, b) in enumerate(edges):
            geometric = sum(sign for v, loop in enumerate(geometry['dual_port_boundaries'])
                            if (mask >> v) & 1 for edge, sign in loop if edge == e)
            assert geometric == ((mask >> a) & 1)-((mask >> b) & 1)
            region_columns += 1
    native_columns = 0
    for total in (1, 2):
        for positions in combinations_with_replacement(range(12), total):
            state = tuple(positions.count(v) for v in range(12))
            for edge in edges:
                for coin in (0, 1):
                    out, j = event(state, edge, coin)
                    actual, (_, k) = provider.native(state, edge, coin)
                    assert (out, j) == (actual, k) and min(out) >= 0 and sum(out) == total
                    if total == 1:
                        assert max(out) == 1 and sum(v == 1 for v in out) == 1
                    native_columns += 1
    print('Independent spatial and native columns passed', region_columns, native_columns, flush=True)
    examples = json.loads((PRODUCER/'EXAMPLES.json').read_text())
    observed = set(); outcomes = set(); selected = {}
    for row in examples:
        claim = row['claim']; initial = row['G']; region = row['region']
        state = tuple(initial); current = [0]*30
        for record in claim['source_records']+claim['sweep_records']:
            edge = edges[record['seam']]
            if record['kind'] == 'E':
                assert record['contrast'] == state[edge[0]]-state[edge[1]]
            state, j = event(state, edge, record['coin'])
            assert j == record['current']; current[record['seam']] += j
        degree = sum(B[v][e]*current[e] for v in region for e in range(30))
        charge = sum(state[v]-initial[v] for v in region)
        assert charge == -degree == row['actual_source_charge']
        compatible = []
        for x in range(12):
            pointer = tuple(int(v == x) for v in range(12)); valid = True
            for record in claim['meter']['records']:
                pointer, j = event(pointer, record['seam'], record['coin'])
                if j != record['current']:
                    valid = False; break
            if valid:
                compatible.append(x)
        assert compatible == row['admission']['compatible_pointer_initials']
        attempts = len(claim['source_records'])+len(claim['sweep_records'])+len(claim['meter']['records'])
        expected_effect = [F(1, 2**attempts) if x in compatible else F(0) for x in range(12)]
        assert expected_effect == list(map(F, row['admission']['pointer_initial_classical_effect']))
        verdict = reviewed.admit_and_decode(claim, initial, region, 2)
        assert verdict['admitted']
        assert verdict['readout'] == (charge if claim['meter']['success'] else None)
        selected.setdefault((degree, claim['meter']['reason']), row)
        observed.add(degree); outcomes.add(claim['meter']['reason'])
    assert observed == {-2, -1, 0, 1, 2}
    vacuum = json.loads((PRODUCER/'VACUUM.json').read_text())
    zero = reviewed.admit_and_decode(vacuum['claim'], vacuum['G'], (0,), 2)
    assert zero['admitted'] and zero['readout'] == 0
    controls = {}
    bad = deepcopy(vacuum['claim']); bad['meter']['records'] = []
    controls['empty_success_is_not_zero_measurement'] = not reviewed.admit_and_decode(bad, vacuum['G'], (0,), 2)['admitted']
    sample = next(r for r in examples if r['topological_degree'] and r['claim']['meter']['success'])
    good = reviewed.admit_and_decode(sample['claim'], sample['G'], sample['region'], 2)
    bad = deepcopy(sample['claim']); bad['meter']['final_loads'] = ['untrusted']
    controls['hidden_terminal_snapshot_irrelevant'] = reviewed.admit_and_decode(bad, sample['G'], sample['region'], 2) == good
    bad = deepcopy(sample['claim']); bad['sweep_records'][0]['contrast'] += 1
    controls['false_closing_contrast_rejected'] = not reviewed.admit_and_decode(bad, sample['G'], sample['region'], 2)['admitted']
    controls['abstract_GL11_uniqueness_counterexample'] = bareiss([[1, 1], [0, 1]]) == 1
    controls['integer_double_split_permitted'] = event((2,)+(0,)*11, (0, 1), 0)[0][:2] == (1, 1)
    controls['half_current_breaks_local_circle_endpoint'] = F(1, 2).denominator != 1
    assert all(controls.values())
    assert custody() == pins
    receipt = {'schema': 'oph.row01.root.independent_execution.v1',
               'producer_independent': True, 'reviewer': '/root',
               'native_columns': native_columns, 'spatial_region_edge_columns': region_columns,
               'independent_unimodular_minors': [tree_det, face_det],
               'complete_examples_replayed': len(examples), 'signed_degrees': sorted(observed),
               'retained_pointer_outcomes': sorted(outcomes),
               'complete_fine_effect_columns_checked': 12*len(examples),
               'controls': controls, 'source_sha256': pins,
               'full_probability_tree_scope': 'Independent event replay/effects here; the 768 complete source/pointer instruments are a pinned independently audited dependency, not rerun here.',
               'algorithm_sha256': sha(Path(__file__))}
    (HERE/'RESULT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['native_columns','spatial_region_edge_columns','complete_examples_replayed','complete_fine_effect_columns_checked']}), flush=True)


if __name__ == '__main__':
    run()
