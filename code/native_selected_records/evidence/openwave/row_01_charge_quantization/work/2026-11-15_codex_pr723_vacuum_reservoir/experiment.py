"""Frozen reservoir-family execution, complete traces and exact low-total census."""
import argparse
from collections import deque
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
from model import graph, history
from native_function import native
from reader import observe

HERE = Path(__file__).resolve().parent


def preparations(bound):
    families = {'stack': [bound]+[0]*11, 'saturated': [bound]*12,
                'signed': [bound, 0]+[bound, -bound]*5,
                'negative': [-bound, 0]+[-bound, bound]*5,
                'neutral': [bound, -bound]*6,
                'mixed': [bound, bound, -bound, 1]+[0]*8}
    for name, first in families.items():
        for rotation in (0, 3, 7):
            rotated = first[rotation:]+first[:rotation]
            yield name, rotation, rotated + [0]*(12*(bound-1))


def exact_census(n, edges):
    sectors = []
    for total in (0, 1, 2):
        states = []
        for occupied in combinations_with_replacement(range(n), total):
            state = [0]*n
            for v in occupied:
                state[v] += 1
            states.append(tuple(state))
        index = {state: i for i, state in enumerate(states)}
        reverse, branches = [[] for _ in states], []
        incoming = [0]*len(states)
        binary = [max(state) <= 1 for state in states]
        nonzero = 0
        for i, state in enumerate(states):
            for e, edge in enumerate(edges):
                for coin in (0, 1):
                    out, (_, current) = native(state, edge, coin)
                    target = index[out]
                    branches.append([i, e, coin, current, target])
                    reverse[target].append(i)
                    if binary[i]:
                        assert binary[target]
                        incoming[target] += 1
                        nonzero += current != 0
        distance = [0 if yes else -1 for yes in binary]
        queue = deque(i for i, yes in enumerate(binary) if yes)
        while queue:
            target = queue.popleft()
            for previous in reverse[target]:
                if distance[previous] < 0:
                    distance[previous] = distance[target]+1
                    queue.append(previous)
        assert min(distance) >= 0
        assert all(incoming[i] == 2*len(edges) for i, yes in enumerate(binary) if yes)
        rate = Fraction(nonzero, sum(binary)*2*len(edges))
        assert rate == Fraction(total*(n-total), n*(n-1))
        sectors.append({'total': total, 'states': states, 'branches': branches,
                        'distance_to_binary': distance, 'stationary_current_rate': str(rate)})
    return sectors


def calculate():
    settings = json.loads((HERE/'TASK.json').read_text())
    cell = json.loads((HERE/'seams.json').read_text())
    records, censuses, local = [], [], []
    for bound in settings['bounds']:
        edges, n = graph(cell, bound), 12*bound
        censuses.append({'bound': bound, 'edges': edges, 'sectors': exact_census(n, edges)})
        for a in range(-bound, bound+1):
            for b in range(-bound, bound+1):
                for coin in (0, 1):
                    state = [a, b]+[0]*(n-2)
                    out, (_, current) = native(state, (0, 1), coin)
                    local.append([bound, a, b, coin, current, out[0], out[1]])
        for family, rotation, initial in preparations(bound):
            index = len(records)//2
            for law in ('native', 'swap_hold'):
                seed = settings['seed_start'] + index
                body, probe, outcome = history(initial, law, seed, edges,
                    settings['horizon'], settings['probe_cap'])
                reading = observe(n, edges, body, probe, settings['probe_cap'],
                                  settings['burn'], settings['calibration'])
                records.append({'bound': bound, 'family': family, 'rotation': rotation,
                    'law': law, 'seed': seed, 'initial_audit_only': initial,
                    'body': body, 'probe': probe, 'probe_outcome': outcome,
                    'reading': reading})
    natives, controls = records[::2], records[1::2]
    unit = natives[0]['reading'].get('calibration_unit')
    expected, observed = Fraction(0), 0
    for record in natives:
        reading = record['reading']
        expected += Fraction(reading.get('heldout_cut_sum', 0),
                             reading.get('conditional_probability_denominator', 1))
        observed += reading.get('heldout_nonzero_events', 0)
    trials = len(natives)*settings['heldout']
    threshold = math.sqrt(trials*math.log(2/settings['martingale_alpha'])/2)
    gates = {
        'all_readings_conclusive': all(r['reading']['outcome']=='conclusive' for r in records),
        'all_initial_states_recovered': all(r['reading'].get('initial_from_records') == r['initial_audit_only'] for r in records),
        'all_native_bodies_select_units_by_burn': all(r['reading'].get('postburn_nonunit_states') == 0 for r in natives),
        'calibration_nonzero': unit is not None and unit > 0,
        'native_heldout_currents_are_zero_or_calibrated_unit': unit is not None and all(set(r['reading'].get('heldout_magnitudes', [])) <= {0, unit} for r in natives),
        'control_retains_multiloads': all(r['reading'].get('first_unit_state') is None for r in controls),
        'control_transports_larger_than_calibrated_unit': unit is not None and all(max(r['reading'].get('heldout_magnitudes', [0])) > unit for r in controls),
        'conditional_transport_prediction': abs(observed-float(expected)) <= threshold}
    receipt = {'candidate': settings['candidate'], 'checks': gates, 'passed': all(gates.values()),
        'native_histories': len(natives), 'control_histories': len(controls),
        'native_body_calls': len(natives)*settings['horizon'],
        'control_body_calls': len(controls)*settings['horizon'],
        'native_probe_calls': sum(len(r['probe']) for r in records),
        'native_census_calls': sum(len(s['branches']) for c in censuses for s in c['sectors']),
        'native_local_calls': len(local), 'calibrated_unit': unit,
        'heldout_trials': trials, 'heldout_nonzero_observed': observed,
        'heldout_conditional_expected': str(expected), 'martingale_threshold': threshold,
        'martingale_alpha': settings['martingale_alpha'],
        'native_first_unit_max': max((r['reading'].get('first_unit_state') or 0) for r in natives),
        'same_group': True, 'physical_row_award': False}
    return receipt, {'censuses': censuses, 'local_branches': local, 'records': records}


def encode(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt, evidence = calculate()
    payload = gzip.compress(encode(evidence), mtime=0)
    receipt['evidence_sha256'] = hashlib.sha256(payload).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'EVIDENCE.json.gz').write_bytes(payload)
    (args.output/'RECEIPT.json').write_bytes(encode(receipt))
    print(json.dumps(receipt, indent=2))
