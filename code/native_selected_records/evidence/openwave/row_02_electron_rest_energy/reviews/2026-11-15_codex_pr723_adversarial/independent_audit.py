"""Independent review arithmetic; imports no producer source or verifier."""
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'work/2026-11-15_codex_pr723_rest_selection'


def endpoint(position, sign, edge, coin):
    assert type(coin) is int and coin in (0, 1)
    a, b = edge
    values = {position: sign} if position is not None else {}
    old = values.get(a, 0)
    total = old + values.get(b, 0)
    lo = total // 2
    hi = total - lo
    values[a], values[b] = (lo, hi) if coin == 0 else (hi, lo)
    occupied = [v for v, q in values.items() if q]
    assert len(occupied) in (0, 1)
    return occupied[0] if occupied else None, old - values[a]


def distances(edges, targets):
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    d = {v: 0 for v in targets}
    queue = deque(targets)
    while queue:
        a = queue.popleft()
        for b in adj[a]:
            if b not in d:
                d[b] = d[a] + 1
                queue.append(b)
    return d, adj


def audit_static(graph, receipt, records):
    edges = graph['seams']
    local = [e for e in edges if max(e) < 12]
    f = {4: F(1, 2), 5: F(-1, 2), 6: F(-1, 2), 7: F(1, 2)}
    rho = defaultdict(F)
    for edge in edges:
        for c in (0, 1):
            groups = defaultdict(lambda: defaultdict(F))
            for p, amplitude in f.items():
                y, j = endpoint(p, 1, edge, c)
                groups[j][y] += amplitude
            for state in groups.values():
                for a, aa in state.items():
                    for b, bb in state.items():
                        rho[a, b] += aa * bb / (2 * len(edges))
    fidelity = sum(f.get(a, 0) * f.get(b, 0) * weight
                   for (a, b), weight in rho.items())
    outside = sum(v for (a, b), v in rho.items() if a == b and a not in f)
    assert fidelity == F(2503, 2520) and outside == F(1, 315)
    assert receipt['baseline']['compact_fidelity_after_one_unread_event'] == str(fidelity)
    assert receipt['baseline']['compact_probability_outside_after_one_event'] == str(outside)
    K = defaultdict(lambda: [[0, 0], [0, 0]])
    work = defaultdict(lambda: defaultdict(F))
    phase_counts = Counter()
    for phase, n, p, sign, edge, coin, j, dest, meta in records['columns']:
        expected, current = endpoint(p, sign, edge, coin)
        assert (dest, j) == (expected, current)
        phase_counts[phase] += 1
        if phase == 'hold_active':
            eid, h, t, selector = meta
            assert p == t or selector == 0
            assert edge == sorted((h, t))
            K[eid, h, t, sign][int(p == t)][int(dest == t)] += 1
        elif phase == 'hold_idle':
            eid, h, t, s, observed_p, selector = meta
            assert observed_p == h and selector == 1
            assert p is None and sign == 0 and n == 2 and edge == [0, 1]
            K[eid, h, t, s][0][0] += 1
        elif phase == 'return':
            eid, h, t = meta
            d, adj = distances(local, (h, t))
            target = min(v for v in adj[p] if d[v] == d[p] - 1)
            assert edge == sorted((p, target))
        elif phase == 'work':
            eid, h, t = meta
            work[eid, h, t, sign][int(dest == t) - int(p == t)] += (
                F(1, 3) if p == t else F(2, 3)) / 2
    assert len(K) == len(work) == 120
    assert all(k == [[3, 1], [2, 2]] for k in K.values())
    assert all(w == {-1: F(1, 6), 0: F(1, 2), 1: F(1, 3)} for w in work.values())
    assert len(records['columns']) == receipt['record_columns'] == 14520
    assert phase_counts['return'] == 2400
    assert F(receipt['work']['reset_exponential_work_mean']) == F(171, 256)
    assert receipt['physical_rest_energy_established'] is False
    assert not any(receipt['physical_gates'].values())
    return {'static_records_independently_checked': len(records['columns']),
            'compact_fidelity': str(fidelity), 'compact_outside_probability': str(outside),
            'local_policy_contexts': len(K), 'equilibrium_work_contexts': len(work)}


def audit_trajectory(trial, local):
    position, sign = trial['preparation']
    h = tail = gate = None
    positions = [position]
    work = []
    mode = 'SEARCH'
    for i, (edge, coin, j, tag, observed) in enumerate(trial['events']):
        before = position
        assert type(j) is int
        if tag == 'search':
            assert h is None and edge in local
        elif tag == 'gate':
            assert edge == [12, 13] and position in (h, tail)
        elif tag in ('hold_active', 'hold_idle'):
            active = position == tail or gate == 0
            assert tag == ('hold_active' if active else 'hold_idle')
            assert edge == (sorted((h, tail)) if active else [12, 13])
        elif tag in ('kick', 'return'):
            target = [trial['kick_target']] if tag == 'kick' else [h, tail]
            d, adj = distances(local, target)
            nxt = min(v for v in adj[position] if d[v] == d[position] - 1)
            assert edge == sorted((position, nxt))
        elif tag == 'work':
            assert edge == sorted((h, tail)) and position in (h, tail)
        else:
            raise AssertionError('Unknown event purpose')
        if trial['held_control']:
            assert j == 0
        else:
            position, expected_j = endpoint(position, sign, edge, coin)
            assert j == expected_j
        if tag == 'search' and j:
            h, tail, mode = position, before, 'HOLD'
        elif tag == 'gate':
            gate = coin
        elif tag == 'kick':
            mode = 'PERTURBED'
        elif tag == 'return':
            mode = 'HOLD' if position in (h, tail) else 'RETURN'
        assert observed == [position, h, tail, mode, gate]
        if tag == 'work':
            work.append([i, int(before == tail), int(position == tail),
                         int(position == tail) - int(before == tail)])
        positions.append(position)
    assert work == trial['work']
    if 'return_start' in trial:
        assert positions[trial['return_start']] == trial['kick_target']
        assert positions[trial['return_end']] in (h, tail)
        assert 0 < trial['return_end'] - trial['return_start'] <= 32
    return positions


def audit_traces(graph, receipt, traces):
    local = [e for e in graph['seams'] if max(e) < 12]
    assert len(traces['trials']) == len(traces['held_controls']) == 24
    for trial in traces['trials'] + traces['held_controls']:
        audit_trajectory(trial, local)
    assert all(not t['failures'] for t in traces['trials'])
    assert all(t['failures'] == ['discovery_timeout'] for t in traces['held_controls'])
    for pair in traces['pairs']:
        aa, bb = audit_trajectory(pair['a'], local), audit_trajectory(pair['b'], local)
        joint = pair['joint']
        a_count = next(i + 1 for i, e in enumerate(pair['a']['events'])
                       if e[4][1] is not None) + 16
        b_count = next(i + 1 for i, e in enumerate(pair['b']['events'])
                       if e[4][1] is not None) + 16
        assert len(joint) == 2048
        assert len(pair['selector_records']) == (2048 if pair['exclusive'] else 0)
        for i, (before, after, end_a, end_b) in enumerate(joint):
            tail_a = pair['a']['events'][end_a - 1][4][2]
            tail_b = pair['b']['events'][end_b - 1][4][2]
            assert before == [int(aa[a_count] == tail_a), int(bb[b_count] == tail_b)]
            assert after == [int(aa[end_a] == tail_a), int(bb[end_b] == tail_b)]
            if pair['exclusive']:
                c, j = pair['selector_records'][i]
                assert type(c) is int and c in (0, 1) and j == 0
                assert (end_a - a_count, end_b - b_count) == ((2, 0) if c == 0 else (0, 2))
            else:
                assert (end_a - a_count, end_b - b_count) == (2, 2)
            a_count, b_count = end_a, end_b
    work = Counter(w[3] for t in traces['trials'] for w in t['work'])
    assert sum(work.values()) == receipt['work_draws'] == 768
    assert receipt['work_counts'] == {str(k): work[k] for k in (-1, 0, 1)}
    actual_events = sum(len(t['events']) for t in traces['trials']) + sum(
        len(p['a']['events']) + len(p['b']['events']) + len(p['selector_records'])
        for p in traces['pairs'])
    assert receipt['native_event_records'] == actual_events
    assert receipt['physical_rest_energy_established'] is False
    return {'opaque_source_record_events_checked': actual_events,
            'work_counts': dict(work), 'work_draws': sum(work.values()),
            'single_source_burnin_tv_bound': str(F(2, 3) / 4 ** 8),
            'work_law_burnin_tv_bound': str(F(1, 3) / 4 ** 8),
            'independent_macroepoch_native_events': 4,
            'exclusive_macroepoch_native_events_including_selector': 3}


def audit_quarter(receipt, records):
    kernels = defaultdict(lambda: [[F(0), F(0)], [F(0), F(0)]])
    work = defaultdict(lambda: defaultdict(F))
    groups = defaultdict(list)
    for phase, h, t, sign, p, bits, c, j, y in records:
        key = h, t, sign
        assert p in (h, t) and sign in (-1, 1)
        if phase == 'work':
            expected, current = endpoint(p, sign, sorted((h, t)), c)
            assert (y, j) == (expected, current)
            work[key][int(y == t) - int(p == t)] += (F(1, 5) if p == t else F(4, 5)) / 2
        else:
            assert len(bits) == 3 and all(type(v) is int and v in (0, 1) for v in bits)
            groups[h, t, sign, p, tuple(bits)].append((phase, c, j, y))
    for (h, t, sign, p, bits), events in groups.items():
        assert len(events) == 3
        assert events[0] == ('gate', bits[0], 0, None)
        assert events[1] == ('gate', bits[1], 0, None)
        active = p == t or bits[:2] == (0, 0)
        phase, coin, current, y = events[2]
        assert phase == ('active' if active else 'idle') and coin == bits[2]
        expected, j = endpoint(p, sign, sorted((h, t)), coin) if active else (p, 0)
        assert (y, current) == (expected, j)
        kernels[h, t, sign][int(p == t)][int(y == t)] += F(1, 8)
    assert len(records) == receipt['native_calls'] == 6240
    assert len(groups) == 120 * 16
    assert len(kernels) == len(work) == 120
    assert all(k == [[F(7, 8), F(1, 8)], [F(1, 2), F(1, 2)]] for k in kernels.values())
    assert all(w == {-1: F(1, 10), 0: F(1, 2), 1: F(2, 5)} for w in work.values())
    assert receipt['stationary_reference'] == ['4/5', '1/5']
    assert receipt['physical_rest_energy_established'] is False
    return {'quarter_gate_native_records_checked': len(records),
            'quarter_gate_contexts': len(kernels),
            'quarter_gate_stationary_reference': ['4/5', '1/5'],
            'quarter_gate_modular_gap': 'log(4)',
            'quarter_gate_native_events_per_epoch': 3}


def main():
    graph = json.loads((SOURCE / 'inputs/graph.json').read_text())
    receipt = json.loads((SOURCE / 'RECEIPT.json').read_text())
    records = json.loads((SOURCE / 'RECORDS.json').read_text())
    result = audit_static(graph, receipt, records)
    if (SOURCE / 'TUNABILITY_RECORDS.json').exists():
        result.update(audit_quarter(
            json.loads((SOURCE / 'TUNABILITY_RECEIPT.json').read_text()),
            json.loads((SOURCE / 'TUNABILITY_RECORDS.json').read_text())))
    if (SOURCE / 'TRAJECTORIES.json').exists():
        tr = json.loads((SOURCE / 'TRAJECTORY_RECEIPT.json').read_text())
        traces = json.loads((SOURCE / 'TRAJECTORIES.json').read_text())
        result.update(audit_traces(graph, tr, traces))
        mutants = {}
        for target in ('work_summary', 'joint_first_before', 'hidden_position', 'native_current'):
            rr, tt = copy.deepcopy(tr), copy.deepcopy(traces)
            if target == 'work_summary':
                rr['work_counts']['0'] += 1
            elif target == 'joint_first_before':
                tt['pairs'][0]['joint'][0][0][0] ^= 1
            elif target == 'hidden_position':
                tt['trials'][0]['events'][0][4][0] = 13
            else:
                tt['trials'][0]['events'][0][2] += 1
            try:
                audit_traces(graph, rr, tt)
            except AssertionError:
                mutants[target] = 'rejected'
            else:
                mutants[target] = 'MISSED'
        assert all(v == 'rejected' for v in mutants.values())
        result['mutations'] = mutants
    result['same_group'] = True
    result['physical_rest_energy_established'] = False
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
