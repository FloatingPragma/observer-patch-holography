"""Actual controller under an independent floor/ceiling source implementation."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import time

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]/'work/2026-11-16_codex_record_only_admission_v2'


def require(x, why):
    if not x:
        raise ValueError(why)


def load_controller():
    spec = importlib.util.spec_from_file_location('tested_record_observer', PACKET/'observer.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.admit


def execute(admit, graph, first, seed, forced_coin=None):
    n = graph['core_N']
    initial = list(first)+[0]*(graph['N']-12)
    for ring in graph['rings']:
        initial[ring[2*(len(ring)//3)]] = 1
    live = initial.copy()
    rng = random.Random(seed)
    tape = []
    edges = {tuple(sorted(x)) for x in graph['edges']}
    def request(edge, stage):
        a, b = edge
        require(tuple(sorted(edge)) in edges, 'non-native requested seam')
        coin = rng.randrange(2) if forced_coin is None else forced_coin
        old_a, old_b = live[a], live[b]
        total = old_a+old_b
        live[a] = (total+coin)//2
        live[b] = total-live[a]
        j = old_a-live[a]
        d = old_a-old_b
        require(j == (d-(2*coin-1)*(abs(d) % 2))//2, 'independent endpoint disagreement')
        record = [a, b, coin, j]
        tape.append([stage, record])
        return record
    result = admit(request, graph)
    if forced_coin is not None:
        require(result['outcome'] == 'probe_timeout', 'lost timeout')
        return {'outcome': result['outcome'], 'events': len(tape)}
    require(result['initial_core'] == initial[:n], 'hidden initial-state reconstruction')
    require(result['observed_final'] == live, 'actual source/controller disagreement')
    q = sum(first)
    require(sum(live) == sum(initial), 'global conservation')
    if q == 0:
        require(result['outcome'] == 'no_unit' and all(x == 0 for x in live[:n]), 'zero-source output')
    else:
        sign = 1 if q > 0 else -1
        require(result['outcome'] == 'signed_target_unit', 'nonzero source not admitted')
        require([live[v] for v in graph['target']] == [0, sign, 0], 'target not isolated')
        require(sum(live[:n]) == q-sign and set(live[:n]) <= {0, sign}, 'excess reservoir missing')
    expected_outside = initial[n:]
    if q:
        expected_outside[graph['target'][1]-n] = 1 if q > 0 else -1
    require(live[n:] == expected_outside, 'marker/copy/blank region changed')
    require(all(y <= x-2 for x, y in result['strict_decreases']), 'unjustified potential decrease')
    return {'input': first, 'inferred_total': sum(result['initial_core']), 'outcome': result['outcome'],
            'target': [live[v] for v in graph['target']], 'reservoir_total': sum(live[:n]),
            'events': len(tape), 'strict_decreases': result['strict_decreases'],
            'tape_sha256': hashlib.sha256(json.dumps(tape, separators=(',', ':')).encode()).hexdigest()}


def run():
    started = time.perf_counter()
    menu = json.loads((HERE/'INPUTS.json').read_text())
    graphs = json.loads((PACKET/'PUBLIC_GRAPHS.json').read_text())
    admit = load_controller()
    rows = []
    for block in menu['groups']:
        g = graphs[f"{block['B']}:{block['M']}"]
        for first in block['inputs']:
            for seed in menu['public_source_seeds']:
                result = execute(admit, g, first, seed)
                rows.append({'B': block['B'], 'M': block['M'], 'public_stream': seed, **result})
    timeout = execute(admit, graphs['2:8'], [2, -1]+[0]*10, 0, forced_coin=0)
    require(len(rows) == 480, 'held-out case count')
    return {'passed': True, 'cases': rows, 'case_count': len(rows),
            'signed_extractions': sum(x['outcome'] == 'signed_target_unit' for x in rows),
            'zero_outcomes': sum(x['outcome'] == 'no_unit' for x in rows),
            'literal_independent_source_events': sum(x['events'] for x in rows)+timeout['events'],
            'timeout': timeout, 'seconds': time.perf_counter()-started,
            'independence': 'Imports actual observer under test; independently implements source and input generation; no producer source import',
            'external_G6': False}


def main():
    p = argparse.ArgumentParser()
    group = p.add_mutually_exclusive_group()
    group.add_argument('--record-to', type=Path)
    group.add_argument('--verify', type=Path)
    a = p.parse_args()
    if a.record_to:
        require(not a.record_to.exists(), 'refusing to overwrite evidence')
    result = run()
    if a.record_to:
        with a.record_to.open('x') as stream:
            json.dump(result, stream, indent=2); stream.write('\n')
    else:
        saved = json.loads((a.verify or HERE/'AUDIT.json').read_text())
        observed = dict(result)
        saved.pop('seconds'); observed.pop('seconds')
        require(saved == observed, 'saved audit differs')
    print(json.dumps({k: v for k, v in result.items() if k != 'cases'}))


if __name__ == '__main__':
    main()
