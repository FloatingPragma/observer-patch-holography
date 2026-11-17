"""Independent endpoint/count-clock audit; no producer imports (same author)."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import sys


def step_sparse(state, edge, c):
    a, b = edge
    left, right = state.get(a, 0), state.get(b, 0)
    total = left + right
    lo, hi = total // 2, -((-total) // 2)
    out = dict(state)
    out[a], out[b] = (lo, hi) if c == 0 else (hi, lo)
    out = {i: v for i, v in out.items() if v}
    return out, left - out.get(a, 0)


def verify(directory):
    p = Path(directory)
    graph = json.loads((p / 'graph.json').read_text())
    edges = [tuple(e) for e in graph['seams']]
    receipt = json.loads((p / 'RECEIPT.json').read_text())
    freeze = json.loads((p / 'PREREG.json').read_text())
    for name, digest in freeze['sha256'].items():
        assert hashlib.sha256((p / name).read_bytes()).hexdigest() == digest
    raw = (p / 'RECORDS.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == receipt['evidence_sha256']
    assert hashlib.sha256((p / 'PREREG.json').read_bytes()).hexdigest() == receipt['prereg_sha256']
    data = json.loads(gzip.decompress(raw))
    calls = Counter()
    assert len(data['source_modes']) == 8
    basis = json.loads((p / 'DEVELOPMENT_BASIS.json').read_text())['basis_by_cell']
    expected_inputs = {}
    for c in range(2):
        for j in range(2):
            expected_inputs[f'S{c}{j}'] = [[12*c+i, v] for i, v in enumerate(basis[c][j]) if v]
    for a in range(2):
        for b in range(2):
            expected_inputs[f'P{a}{b}'] = [[p0*256+n0, v*w] for p0, v in expected_inputs[f'S0{a}']
                                          for n0, w in expected_inputs[f'S1{b}']]
    for mode in data['source_modes']:
        assert mode['inputs'] == expected_inputs[mode['name']]
        coeff = dict(mode['inputs'])
        actual, visited = defaultdict(int), set()
        for code, e, c, current, dest in mode['branches']:
            assert (code, e, c) not in visited
            visited.add((code, e, c))
            if mode['kind'] == 'single':
                state = {code: 1}
            else:
                a, b = divmod(code, 256)
                state = {a: 1, b: -1}
            out, j = step_sparse(state, edges[e], c)
            calls['source_mode_columns'] += 1
            assert current == j
            if mode['kind'] == 'single':
                assert len(out) == 1 and list(out.values()) == [1]
                wanted = next(iter(out))
            elif not out:
                wanted = -1
            else:
                assert sorted(out.values()) == [-1, 1]
                wanted = 256 * next(i for i, v in out.items() if v == 1) + next(i for i, v in out.items() if v == -1)
            assert wanted == dest
            actual[dest] += coeff[code]
        assert len(visited) == len(coeff) * 1260
        actual = {q: v for q, v in actual.items() if v}
        eigen = 1254 if mode['kind'] == 'single' else 1248
        assert actual == {q: eigen * v for q, v in coeff.items()}
        assert mode['mean_numerator'] == sorted([q, v] for q, v in actual.items())
        assert mode['eigen_numerator'] == eigen and mode['hold_generator_eigenvalue'] == 0

    # Independent reciprocal positive-series enclosure of exp(-1).
    ex = sum((F(1, math.factorial(k)) for k in range(201)), F())
    ex_upper = ex + F(202, 201 * math.factorial(201))
    lo, hi = map(F, data['exp_minus_one_interval'])
    assert lo <= 1/ex_upper <= 1/ex <= hi
    cdfs, total = [], F()
    for n in range(33):
        total += F(1, math.factorial(n))
        cdfs.append((lo * total, min(F(1), hi * total)))

    @lru_cache(maxsize=None)
    def limits(bits):
        result, previous = [], F()
        for n, (low, high) in enumerate(cdfs):
            minimum = -((-previous.numerator * (1 << bits)) // previous.denominator)
            maximum = (low.numerator * (1 << bits)) // low.denominator - 1
            result.append((minimum, maximum, n))
            previous = high
        return result

    def classify(prefix, bits):
        for a, b, n in limits(bits):
            if a <= prefix <= b:
                return n
        return -2 if F(prefix, 1 << bits) > cdfs[-1][1] else -1

    tree = data['clock_tree']
    states = {(0, 0): 0}
    visited = set()
    for depth, prefix, pos, c, j, child, target in tree['native_prefix_branches']:
        assert (depth, prefix, c) not in visited
        visited.add((depth, prefix, c))
        assert states[(depth, prefix)] == pos
        out, actual = step_sparse({pos: 1}, (0, 1), c)
        calls['clock_prefix_columns'] += 1
        assert actual == j and list(out) == [target]
        assert child == 2 * prefix + abs(j)
        states[(depth + 1, child)] = target
    assert len(visited) == 131070
    assert all((16, i) in states for i in range(65536))
    leaves = [classify(i, 16) for i in range(65536)]
    assert leaves == tree['leaf_outcomes']
    counts = Counter(leaves)
    assert {str(k): v for k, v in sorted(counts.items())} == tree['leaf_counts']
    fail = F(sum(v for n, v in counts.items() if n < 0), 65536)
    assert fail == F(tree['failure_mass']) == F(1, 8192)
    tail = F(34, 33 * math.factorial(33))
    assert tail == F(data['poisson_count_cap_tail_upper'])
    for n in range(33):
        mass = F(counts[n], 65536)
        assert mass <= hi/math.factorial(n) and mass + fail + tail >= lo/math.factorial(n)
    for f, a, exponent in zip(tree['transfer_factors'], (F(209,210),F(104,105)), (F(1,210),F(1,105))):
        lower = sum(F(v,65536)*a**n for n,v in counts.items() if n>=0)
        assert lower == F(f['retained_factor'])
        exp_lo, exp_hi = map(F,f['ideal_factor_interval'])
        s = sum((exponent**k/math.factorial(k) for k in range(201)),F())
        upper = s + exponent**201/math.factorial(201)/(1-exponent/202)
        assert exp_lo <= 1/upper <= 1/s <= exp_hi
        assert lower <= exp_lo <= exp_hi <= lower+fail

    def aux_record(aux, row, rng):
        c,j,target = row
        assert c == rng.randrange(2)
        out, actual = step_sparse({aux:1},(0,1),c)
        calls['clock_auxiliary'] += 1
        assert actual==j and list(out)==[target]
        return target,abs(j)

    seen_counts=set()
    assert len(data['clock_histories'])==8
    for index,r in enumerate(data['clock_histories']):
        kind='single' if index<4 else 'pair';port=4+index%4
        initial={port:1}
        if kind=='pair':initial[16]=-1
        assert r['kind']==kind and r['port']==port and r['seed']==72331000+index
        assert r['initial_audit_only']==[initial.get(i,0) for i in range(240)]
        x,aux,rng=initial,0,random.Random(r['seed'])
        assert len(r['ticks'])==64 and r['complete']
        for tick in r['ticks']:
            prefix=0
            for k,row in enumerate(tick['clock_records'],1):
                aux,bit=aux_record(aux,row,rng);prefix=2*prefix+bit
                value=classify(prefix,k)
                if k<len(tick['clock_records']):assert value==-1
            assert value==tick['count'] and value>=0
            seen_counts.add(value)
            assert len(tick['events'])==value
            for event in tick['events']:
                for trial_index,trial in enumerate(event['address_records']):
                    assert len(trial)==10
                    label=0
                    for row in trial:
                        aux,bit=aux_record(aux,row,rng);label=2*label+bit
                    assert (label<630)==(trial_index==len(event['address_records'])-1)
                e,c,j=event['source'];assert e==label and c==rng.randrange(2)
                x,actual=step_sparse(x,edges[e],c)
                calls['clocked_source_events']+=1;assert j==actual
        assert r['final_source_audit_only']==[x.get(i,0) for i in range(240)]
        assert r['final_auxiliary']==[int(aux==0),int(aux==1)]
    assert 0 in seen_counts and any(n>=2 for n in seen_counts)
    rng=random.Random(72331999)
    assert len(data['hold_auxiliary'])==16
    for t in data['hold_auxiliary']:
        assert t['count']==0
        prefix=0
        for k,(c,j,pos) in enumerate(t['records'],1):
            assert c==rng.randrange(2) and j==pos==0
            value=classify(prefix,k)
            assert (value==0)==(k==len(t['records']))
    for row,denominator in zip(data['rare_event'],(1,16,256,4096)):
        eps=F(1,denominator);assert F(row['epsilon'])==eps
        certified=[]
        for rate,label in ((6,'single'),(12,'pair')):
            z=rate*eps/1260
            partial=sum((z**k/k for k in range(1,65)),F())
            bounds=(1260/eps*partial,1260/eps*(partial+z**65/(65*(1-z))))
            a,b=map(F,row[label]);assert a<=bounds[0]<=bounds[1]<=b
            certified.append(bounds)
        lower,upper=map(F,row['excess'])
        assert lower<=certified[1][0]-2*certified[0][1]
        assert upper>=certified[1][1]-2*certified[0][0]
        assert 0<lower<=upper<F(row['excess_upper'])==36*eps/(1260-12*eps)
    failure=F(68,2**128)+68*(hi-lo)+tail
    assert failure==F(receipt['per128bit_clock_sampler_failure_upper'])<F(1,10**30)
    assert F(receipt['per_source_address_failure_upper'])==F(197,512)**64
    assert receipt['literal_native_calls_by_kind']==dict(calls)
    assert sum(calls.values())==receipt['literal_native_calls']==242175
    assert receipt['native_single_generator_values']==[6]*4
    assert receipt['native_pair_generator_values']==[12]*4
    assert F(receipt['fixed_one_attempt_factor_gap'])==F(209,210)**2-F(104,105)==F(1,44100)
    assert all(receipt['gates'].values()) and receipt['all_gates_pass']
    for field in ('rare_event_finite_excess_is_zero','physical_rest_energy_identified',
                  'coherent_UL_execution_in_this_packet','full_row_solved','source_pin_changed'):
        assert receipt[field] is False
    return {'verified':True,'producer_imports':[], 'native_calls_accounted':sum(calls.values()),
            'source_columns':calls['source_mode_columns'],'prefix_columns':calls['clock_prefix_columns'],
            'actual_clocked_events':calls['clocked_source_events'],'complete_histories':8,
            'exact_ideal_source_generator_additivity':'12=6+6',
            'retained16bit_failure_mass':str(fail),'per128bit_failure_upper':float(failure),
            'evidence_sha256':hashlib.sha256(raw).hexdigest()}


if __name__=='__main__':
    print(json.dumps(verify(sys.argv[1] if len(sys.argv)>1 else Path(__file__).parent),indent=2))
