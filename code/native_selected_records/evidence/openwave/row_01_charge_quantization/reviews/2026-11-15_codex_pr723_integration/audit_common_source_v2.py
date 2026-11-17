"""Independent finite embeddings into the proposed common instrument.

This is a compatibility audit of already known row results, not a new physical
prediction or a full common-controller implementation.
"""
from collections import defaultdict
from hashlib import sha256
import gzip
import json
from math import sqrt
from pathlib import Path
import common_source_v2 as common

HERE = Path(__file__).resolve().parent
OW = HERE.parents[2]
R1 = OW/'row_01_charge_quantization/work/2026-11-15_codex_pr723_unit_selection'
R2 = OW/'row_02_electron_rest_energy/work/2026-11-15_codex_pr723_rest_selection'
R7 = OW/'row_07_double_cover/work/2026-11-15_codex_pr723_rotation'


def embed(q):
    return tuple(q) + (0,)*(14-len(q))


def verify_input_pins():
    pins = json.loads((HERE/'COMMON_V2_PREREG.json').read_text())
    for relative, expected in pins['inputs'].items():
        assert sha256((OW/relative).read_bytes()).hexdigest() == expected, relative


def main():
    verify_input_pins()
    checked = defaultdict(int)
    for coin in (0, 1):
        seen = set()
        for a in range(-8, 9):
            for b in range(-8, 9):
                q = embed((a, b))
                out, j = common.event(q, (0, 1), coin)
                low, high = (a+b)//2, (a+b)-(a+b)//2
                assert out[:2] == ((low, high) if coin == 0 else (high, low))
                assert max(map(abs, out)) <= 8 and sum(out) == sum(q)
                assert (out, j) not in seen
                seen.add((out, j))
                restored = list(out)
                restored[0] += j
                restored[1] -= j
                assert tuple(restored) == q
                checked['local_bound_inputs'] += 1

    e1 = json.loads(gzip.decompress((R1/'EVIDENCE.json.gz').read_bytes()))
    edges = json.loads((R1/'seams.json').read_text())
    for sector in e1['sectors']:
        states = sector['states']
        for initial, edge, coin, current, target in sector['branches']:
            out, j = common.event(embed(states[initial]), edges[edge], coin)
            assert out == embed(states[target]) and j == current
            checked['row01_exhaustive_branches'] += 1
    for record in e1['histories']:
        if record['law'] != 'native':
            continue
        q = embed(record['initial_audit_only'])
        for a, b, c, current in record['body'] + record['probe']:
            q, j = common.event(q, (a, b), c)
            assert j == current
            checked['row01_native_history_events'] += 1

    for phase, n, p, sign, edge, coin, current, dest, meta in json.loads((R2/'RECORDS.json').read_text())['columns']:
        if n == 240:
            checked['row02_old_240_port_comparators_excluded'] += 1
            continue
        assert n in (2, 12)
        offset = 12 if n == 2 else 0
        q = [0]*14
        if p is not None:
            q[p+offset] = sign
        out, j = common.event(tuple(q), tuple(x+offset for x in edge), coin)
        expected = [0]*14
        if dest is not None:
            expected[dest+offset] = sign
        assert out == tuple(expected) and j == current
        checked['row02_local_native_columns'] += 1

    signed = [(0,)*12] + [tuple(int(i == p)-int(i == r) for i in range(12))
                         for p in range(12) for r in range(12) if p != r]
    for initial, edge, coin, current, target in json.loads((R7/'RECORDS.json').read_text())['native']:
        out, j = common.event(embed(signed[initial]), edges[edge], coin)
        assert out == embed(signed[target]) and j == current
        checked['row07_signed_native_columns'] += 1

    chart = (0, 1, 5)
    triangle = ((0, 1), (1, 2), (2, 0))
    berry = R7/'berry'
    br = json.loads((berry/'RECORDS.json').read_text())['source']['native']
    wr = json.loads((berry/'RECEIPT.json').read_text())['weight_loop']
    rotor = json.loads((berry/'ROTOR_RECORDS.json').read_text())
    rotor_receipt = json.loads((berry/'ROTOR_RECEIPT.json').read_text())
    all_diags = json.loads((berry/'RECEIPT.json').read_text())['instruments']['native']['loop'] + [json.loads((berry/'RECEIPT.json').read_text())['instruments']['native']['reference']]
    all_diags += rotor_receipt['instruments']['polar'] + rotor_receipt['instruments']['apolar']
    datasets = [(w, table) for w, table in zip(wr, br['loop'])]
    datasets.append(((32, 16, 16), br['reference']))
    for coupling in ('polar', 'apolar'):
        weights = [g['native_ticket_weights'] for g in rotor_receipt['geometry_records']
                   if g['coupling'] == coupling]
        datasets.extend(zip(weights, rotor[coupling]))
    for weights, records in datasets:
        tickets = [tuple(chart[i] for i in triangle[e])
                   for e, count in enumerate(weights) for _ in range(count)]
        counts = [[0]*3 for _ in range(3)]
        for initial, ticket, coin, current, target in records:
            q = [0]*14
            q[chart[initial]] = 1
            out, j = common.event(tuple(q), tickets[ticket], coin)
            expected = [0]*14
            expected[chart[target]] = 1
            assert out == tuple(expected) and j == current
            counts[target][initial] += 1
            checked['row07_coherent_native_columns'] += 1
        # Independently compare the full positive instrument and eraser map
        # on a complex superposition, not only the frozen real optical ray.
        amplitudes = {}
        values = (1/sqrt(3), 1j/sqrt(3), -1/sqrt(3))
        for p, value in zip(chart, values):
            q = [0]*14
            q[p] = 1
            amplitudes[tuple(q)] = value
        observed = common.observe(amplitudes, tickets, 'optical_erase')
        assert abs(sum(common.norm2(v) for v in observed['fine'].values())-1) < 1e-12
        for target, p in enumerate(chart):
            q = [0]*14
            q[p] = 1
            expected = sum(counts[target][i]*values[i] for i in range(3))/(128*sqrt(3))
            assert abs(observed['success'].get(tuple(q), 0)-expected) < 1e-12
        assert observed['all_failure_probability'] >= 0
        assert abs(observed['success_probability']+observed['all_failure_probability']-1) < 1e-12
        assert len(observed['failure_states']) >= 383
        # Compare each original Householder failure outcome, not only success.
        original = {}
        for p, value in zip(chart, (1/sqrt(6), 1/sqrt(6), -2/sqrt(6))):
            q = [0]*14
            q[p] = 1
            original[tuple(q)] = value
        full = common.observe(original, tickets, 'optical_erase')
        probabilities = []
        for t in range(64):
            for c in (0,1):
                for j in (-1,0,1):
                    flag = (t,c,j)
                    vector = full['success'] if flag == (0,0,-1) else full['failure_states'][flag]
                    probabilities.append(common.norm2(vector))
        assert abs(sum(probabilities)-1) < 1e-12
        diag = next(d for d in all_diags if tuple(d['weights']) == tuple(weights))
        assert max(abs(x-y) for x,y in zip(probabilities, diag['all_384_analyzer_probabilities_on_initial_ray'])) < 1e-10
        checked['quantum_eraser_vector_checks'] += 1

    q = embed((8, -8))
    wider = common.observe({q: 1}, [(0, 1)], 'optical_erase')
    assert wider['success_probability'] == 0
    assert abs(wider['all_failure_probability']-1) < 1e-12
    try:
        common.event(embed((1,0,0)), (0,2), 0)
    except ValueError:
        nonedge_rejected = True
    else:
        raise AssertionError('Nonedge accepted')
    result = {'nonedge_rejected': nonedge_rejected, 'complete_failure_source_states_returned': True,
              'candidate': 'A-pr723-common-positive-instrument-v2',
              'checks': dict(checked), 'passed': True, 'same_group': True,
              'wider_current_flags_retained_as_failure': True,
              'complete_common_controller_implemented': False,
              'native_coherent_apparatus_compiled': False,
              'physical_rows_awarded': False,
              'scope': 'Common source instrument, observation modes and exact event embeddings; supplied controllers/eraser/analyzer remain explicit.'}
    (HERE/'COMMON_SOURCE_V2_AUDIT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
