"""Actual two-source admission/calibration and complete global joint columns.

Column enumeration is a linearity audit on actual admitted backgrounds, not
a sequence of physical resets or a sampled quantum history.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
CALIBRATION = HERE.parent/'2026-11-16_codex_record_calibration_v1'
sys.path.insert(0, str(CALIBRATION))
from integration import graph, prepare, start

EDGES = ((0, 1), (0, 1), (0, 2), (0, 2), (1, 2))
STREAMS = ((72365001, 72366001, 72365002, 72366002),
           (72365011, 72366011, 72365012, 72366012))
PATTERNS = (([2, -1, 2, -1]+[0]*8, [2, 2, -1]+[0]*9),
            ([1, -2, 2, 2, -1]+[0]*7, [2, -2, 1, -1, 2]+[0]*7))


def require(value, why):
    if not value:
        raise ValueError(why)


def native(q, edge, coin):
    a, b = edge
    out = list(q)
    total = q[a]+q[b]
    out[a] = (total+coin)//2
    out[b] = total-out[a]
    return out, [a, b, coin, q[a]-out[a]]


def index_table(record):
    columns = record['calibration']['columns']
    require(len(columns) == 18, 'incomplete actual calibration')
    table = {(r['input_position'], tuple(r['edge']), r['coin']): r for r in columns}
    require(len(table) == 18, 'duplicate actual calibration')
    return table


def arc(q, g, offset, sign, column, table):
    position, choice = divmod(column, 11)
    if choice == 10:
        out, event = native(q, tuple(offset+x for x in g['blank']), 0)
        require(out == q and event[3] == 0, 'inactive blank changed source')
        return out, column, event
    ticket, coin = divmod(choice, 2)
    edge = EDGES[ticket]
    out, event = native(q, tuple(offset+g['target'][x] for x in edge), coin)
    measured = table[position, edge, coin]
    require(event[3] == measured['current'], 'global current differs from actual calibration')
    y = measured['decoded_output']
    require([out[offset+v] for v in g['target']] == [sign*int(x == y) for x in range(3)],
            'global target does not match measured output')
    cc = coin ^ bool(event[3])
    return out, 11*y+2*ticket+cc, event


def joint_columns(g, left, right, signs, records):
    offset = g['N']
    initial = left+right
    tables = [index_table(r) for r in records]
    ledger = []
    for x in range(33):
        for y in range(33):
            embedded = initial.copy()
            for side, column in enumerate((x, y)):
                off = side*offset
                embedded[off+g['target'][1]] -= signs[side]
                embedded[off+g['target'][column//11]] += signs[side]
            finals = []
            for order in ((0, 1), (1, 0)):
                q = embedded.copy()
                outputs = [x, y]
                events = []
                for side in order:
                    q, outputs[side], event = arc(q, g, side*offset, signs[side], outputs[side], tables[side])
                    events.append(event)
                require(sum(q) == sum(initial), 'joint charge conservation')
                spectator = q.copy()
                for side, column in enumerate(outputs):
                    spectator[side*offset+g['target'][column//11]] -= signs[side]
                    spectator[side*offset+g['target'][1]] += signs[side]
                require(spectator == initial, 'joint spectator changed')
                finals.append((q, outputs))
                ledger.append({'input': [x, y], 'order': list(order), 'output': outputs, 'events': events})
            require(finals[0] == finals[1], 'joint operator depends on address order')
    require(len(ledger) == 2178, 'joint catalog incomplete')
    require(len({tuple(row['output']) for row in ledger if row['order'] == [0, 1]}) == 1089,
            'joint map not a complete permutation')
    # Its full column identity proves action on correlated inputs, not only
    # product densities. Local reflections and complete flags then tensorize.
    return ledger


def run():
    started = time.perf_counter()
    g = graph(2, 32)
    rows = []
    source_calls = 0
    for pattern, pair in enumerate(PATTERNS):
        for sa in (-1, 1):
            for sb in (-1, 1):
                for stream_index, seeds in enumerate(STREAMS):
                    records, audits, initial_states = [], [], []
                    for side, sign in enumerate((sa, sb)):
                        initial = prepare([sign*x for x in pair[side]], g)
                        request, inspect, record = start(initial, g, seeds[2*side], seeds[2*side+1])
                        require(record['outcome'] == 'ready', 'actual common source did not reach continuation')
                        audit = inspect()
                        require(record['admission']['sign'] == sign, 'wrong sign from actual admission')
                        require(audit['final'] == record['admission']['observed_final'], 'calibration reset or changed background')
                        require([audit['final'][v] for v in g['target']] == [0, sign, 0], 'continuation entry is not actual leaf1')
                        require(sum(audit['final'][:g['core_N']]) == sum(pair[side])*sign-sign,
                                'excess reservoir lost')
                        records.append(record)
                        audits.append(audit)
                        initial_states.append(initial)
                        source_calls += sum(audit['literal_endpoint_calls'].values())
                    joint = joint_columns(g, audits[0]['final'], audits[1]['final'], (sa, sb), records)
                    rows.append({'pattern': pattern, 'signs': [sa, sb], 'public_stream_pair': stream_index,
                        'initial_states_audit_only': initial_states, 'observer_records': records,
                        'source_audits': audits, 'joint_columns': joint})
    require(len(rows) == 16, 'missing common input/stream pair')
    failure = F(2*23+6*4*23+24+597, 2**64)
    return {'candidate': 'A-native-amplified-common-pair-v1', 'passed': True,
        'graph': g, 'joint_graph_vertices': 2*g['N'],
        'idle_connecting_edge': [g['blank'][1], g['N']+g['blank'][0]],
        'idle_connecting_edge_events': 0, 'rows': rows, 'pair_count': len(rows),
        'literal_continuing_source_calls': source_calls,
        'literal_global_column_evaluations': len(rows)*2178*2,
        'classical_timeout_union_bound_each': str(failure),
        'classical_success_lower_each': str(1-failure),
        'classical_joint_success_lower': str((1-failure)**2),
        'scope': 'Actual classical continuation and exact joint operator identification; finite quantum instrument outputs are separately bound.',
        'seconds': time.perf_counter()-started, 'external_G6': False}


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--record-to', type=Path)
    group.add_argument('--verify', type=Path)
    args = parser.parse_args()
    if args.record_to:
        require(not args.record_to.exists(), 'refusing to overwrite common evidence')
    result = run()
    if args.record_to:
        with args.record_to.open('xb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as stream:
                stream.write((json.dumps(result, separators=(',', ':'))+'\n').encode())
    else:
        with gzip.open(args.verify or HERE/'SOURCE_JOIN.json.gz', 'rt') as stream:
            saved = json.load(stream)
        left, right = dict(result), dict(saved)
        left.pop('seconds'); right.pop('seconds')
        require(left == right, 'saved common evidence differs')
    print(json.dumps({key: value for key, value in result.items() if key not in ('rows', 'graph')}))


if __name__ == '__main__':
    main()
