"""Independent reader and floor/ceiling audit; imports no producer module."""
from collections import deque
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent


def repair(x, edge, coin, law='native'):
    a, b = edge
    y = list(x)
    if law == 'native':
        low = (x[a]+x[b])//2
        high = x[a]+x[b]-low
        y[a], y[b] = (low, high) if coin == 0 else (high, low)
    elif law == 'swap_hold':
        if coin:
            y[a], y[b] = x[b], x[a]
    else:
        raise AssertionError('Unknown law')
    return tuple(y), x[a]-y[a]


def record_reading(n, edges, body, probe, cap, burn, calibration):
    adjacency = [[] for _ in range(n)]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    order, seen, spanning = [0], {0}, []
    for a in order:
        for b in adjacency[a]:
            if b not in seen:
                seen.add(b)
                order.append(b)
                spanning.append((a, b))
    changes, relative, cursor = [0]*n, [None]*n, 0
    relative[0] = 0
    for a, b in spanning:
        prior = changes[a]-changes[b]
        first = None
        for _ in range(cap):
            if cursor == len(probe):
                return {'outcome': 'inconclusive'}
            aa, bb, coin, current = probe[cursor]
            cursor += 1
            assert (aa, bb) == (a, b) and coin in (0, 1) and type(current) is int
            changes[a] -= current
            changes[b] += current
            if first is None:
                first = (coin, current)
            elif coin != first[0]:
                relative[b] = relative[a]-(2*first[1]+current-prior)
                break
        else:
            assert cursor == len(probe)
            return {'outcome': 'inconclusive'}
    assert cursor == len(probe) and len(seen) == n
    # Independent absolute anchor: integrate ONLY the reference-port currents.
    reference_final = sum(j*((b == 12)-(a == 12)) for a, b, c, j in body)
    terminal = [v + reference_final-relative[12] for v in relative]
    initial = list(terminal)
    for a, b, c, j in reversed(body):
        initial[a] += j
        initial[b] -= j
    state = list(initial)
    first, bad, fitted, cut_sum, nonzero, magnitudes = None, 0, None, 0, 0, set()
    for t in range(len(body)+1):
        units = set(state) <= {0, 1} or set(state) <= {-1, 0}
        if first is None and units:
            first = t
        bad += t >= burn and not units
        if t == len(body):
            break
        a, b, c, j = body[t]
        if burn <= t < burn+calibration and fitted is None and j:
            fitted = abs(j)
        if t >= burn+calibration:
            cut_sum += sum(state[u] != state[v] for u, v in edges)
            nonzero += j != 0
            magnitudes.add(abs(j))
        state[a] -= j
        state[b] += j
    assert state == terminal
    return {'outcome':'conclusive', 'initial_from_records':initial,
            'final_from_records':terminal, 'total_from_records':sum(initial),
            'first_unit_state':first, 'postburn_nonunit_states':bad,
            'calibration_unit':fitted, 'heldout_magnitudes':sorted(magnitudes),
            'heldout_nonzero_events':nonzero, 'heldout_cut_sum':cut_sum,
            'conditional_probability_denominator':2*len(edges)}


def check(directory):
    task = json.loads((HERE/'TASK.json').read_text())
    cell = json.loads((HERE/'seams.json').read_text())
    receipt = json.loads((directory/'RECEIPT.json').read_text())
    raw = (directory/'EVIDENCE.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == receipt['evidence_sha256']
    data = json.loads(gzip.decompress(raw))
    graph_by_bound, census_calls = {}, 0
    assert len(data['censuses']) == len(task['bounds'])
    for census, bound in zip(data['censuses'], task['bounds']):
        n = 12*bound
        edges = [[12*c+a,12*c+b] for c in range(bound) for a,b in cell]
        edges += [[12*c+v,12*(c+1)+v] for c in range(bound-1) for v in range(12)]
        assert census['bound'] == bound and census['edges'] == edges
        graph_by_bound[bound] = edges
        assert len(census['sectors']) == 3
        for total, sector in enumerate(census['sectors']):
            expected_states = []
            for locations in combinations_with_replacement(range(n), total):
                expected_states.append([locations.count(v) for v in range(n)])
            states = sector['states']
            assert sector['total'] == total and states == expected_states
            reverse, incoming = [[] for _ in states], [0]*len(states)
            seen, nonzero = set(), 0
            for i, e, c, j, target in sector['branches']:
                assert 0 <= i < len(states) and 0 <= target < len(states)
                assert 0 <= e < len(edges) and c in (0,1)
                assert (i,e,c) not in seen
                seen.add((i,e,c))
                out, current = repair(states[i], edges[e], c)
                assert list(out) == states[target] and j == current
                reverse[target].append(i)
                if max(states[i]) <= 1:
                    assert max(out) <= 1
                    incoming[target] += 1
                    nonzero += j != 0
            assert len(seen) == len(states)*len(edges)*2
            census_calls += len(seen)
            binary = [i for i,s in enumerate(states) if max(s) <= 1]
            assert all(incoming[i] == 2*len(edges) for i in binary)
            distance = [-1]*len(states)
            queue = deque(binary)
            for i in binary:
                distance[i] = 0
            while queue:
                target = queue.popleft()
                for previous in reverse[target]:
                    if distance[previous] == -1:
                        distance[previous] = distance[target]+1
                        queue.append(previous)
            assert min(distance) >= 0 and distance == sector['distance_to_binary']
            rate = Fraction(nonzero,len(binary)*2*len(edges))
            assert rate == Fraction(total*(n-total),n*(n-1))
            assert str(rate) == sector['stationary_current_rate']
    expected_local = {(B,a,b,c) for B in task['bounds'] for a in range(-B,B+1)
                      for b in range(-B,B+1) for c in (0,1)}
    assert {(B,a,b,c) for B,a,b,c,j,u,v in data['local_branches']} == expected_local
    for B,a,b,c,j,u,v in data['local_branches']:
        out, actual = repair([a,b],(0,1),c)
        assert out == (u,v) and actual == j
    assert len(data['records']) == 72
    native_body = control_body = probes = 0
    readings = []
    for index, record in enumerate(data['records']):
        bound = task['bounds'][(index//2)//18]
        assert record['bound'] == bound
        law = 'native' if index%2 == 0 else 'swap_hold'
        assert record['law'] == law
        assert record['seed'] == task['seed_start']+index//2
        edges, n = graph_by_bound[bound], 12*bound
        state = tuple(record['initial_audit_only'])
        family = ('stack','saturated','signed','negative','neutral','mixed')[(index//2)%18//3]
        rotation = (0,3,7)[(index//2)%3]
        starts = {'stack':[bound]+[0]*11, 'saturated':[bound]*12,
                  'signed':[bound,0]+[bound,-bound]*5,
                  'negative':[-bound,0]+[-bound,bound]*5,
                  'neutral':[bound,-bound]*6,
                  'mixed':[bound,bound,-bound,1]+[0]*8}
        first = starts[family]
        declared = first[rotation:]+first[:rotation]+[0]*(n-12)
        assert record['family'] == family and record['rotation'] == rotation
        assert record['initial_audit_only'] == declared
        assert len(state) == n and all(abs(v) <= bound for v in state[:12])
        assert state[12:] == (0,)*(n-12)
        if index%2:
            assert record['initial_audit_only'] == data['records'][index-1]['initial_audit_only']
        assert len(record['body']) == task['horizon']
        rng = random.Random(record['seed'])
        for a,b,c,j in record['body']:
            assert [a,b] == edges[rng.randrange(len(edges))] and c == rng.randrange(2)
            state, actual = repair(state,(a,b),c,law)
            assert actual == j
        body_final = state
        for a,b,c,j in record['probe']:
            assert [a,b] in edges or [b,a] in edges
            assert c == rng.randrange(2)
            state, actual = repair(state,(a,b),c)
            assert actual == j
        reading = record_reading(n,edges,record['body'],record['probe'],task['probe_cap'],task['burn'],task['calibration'])
        assert reading == record['reading']
        assert reading['initial_from_records'] == record['initial_audit_only']
        assert reading['final_from_records'] == list(body_final)
        assert reading['outcome'] == record['probe_outcome'] == 'conclusive'
        readings.append(reading)
        native_body += len(record['body']) if law == 'native' else 0
        control_body += len(record['body']) if law == 'swap_hold' else 0
        probes += len(record['probe'])
    u = readings[0]['calibration_unit']
    expected = sum((Fraction(r['heldout_cut_sum'],r['conditional_probability_denominator']) for r in readings[::2]),Fraction(0))
    observed = sum(r['heldout_nonzero_events'] for r in readings[::2])
    T = 36*task['heldout']
    threshold = math.sqrt(T*math.log(2/task['martingale_alpha'])/2)
    assert u > 0 and all(r['postburn_nonunit_states'] == 0 for r in readings[::2])
    assert all(set(r['heldout_magnitudes']) <= {0,u} for r in readings[::2])
    assert all(r['first_unit_state'] is None and max(r['heldout_magnitudes'])>u for r in readings[1::2])
    assert abs(observed-float(expected)) <= threshold
    values = {'native_body_calls':native_body,'control_body_calls':control_body,
              'native_probe_calls':probes,'native_census_calls':census_calls,
              'native_local_calls':len(data['local_branches']),'calibrated_unit':u,
              'heldout_trials':T,'heldout_nonzero_observed':observed,
              'heldout_conditional_expected':str(expected),'martingale_threshold':threshold}
    assert receipt['candidate'] == task['candidate']
    assert receipt['native_histories'] == receipt['control_histories'] == 36
    assert receipt['native_first_unit_max'] == max(r['first_unit_state'] for r in readings[::2])
    assert receipt['martingale_alpha'] == task['martingale_alpha']
    assert len(receipt['checks']) == 8
    assert all(receipt[k] == value for k,value in values.items())
    assert receipt['passed'] and all(receipt['checks'].values())
    return dict(values, independent_reader=True, producer_imports=False,
                same_group=True, passed=True)


if __name__ == '__main__':
    print(json.dumps(check(Path(sys.argv[1]) if len(sys.argv)>1 else HERE),indent=2,sort_keys=True))
