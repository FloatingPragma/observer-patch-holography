#!/usr/bin/env python3
"""Native positive transfer products around a scheduler degeneracy.

No phase pulse or instantaneous eigenprojector is applied. The coherent
source preparation and branch eraser are declared apparatus.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from native_endpoint import native

HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (1, 2), (2, 0))
INITIAL = np.array([1., 1., -2.]) / np.sqrt(6)
TICKETS = 64
REPETITIONS = (1, 2, 4)

def canonical(x):
    return json.dumps(x, sort_keys=True, indent=2, allow_nan=False) + '\n'

def num(x):
    return float(format(float(x), '.12g'))

def weight_loop():
    # Start is frozen at (32,16,16); each emitted point advances one of48
    # steps. Last point returns to start. No sign occurs in this schedule.
    vertices = [(32, 16, 16), (16, 32, 16), (16, 16, 32), (32, 16, 16)]
    return [tuple(a[k] + step * (b[k] - a[k]) // 16 for k in range(3))
            for a, b in zip(vertices, vertices[1:]) for step in range(1, 17)]

def sampler_records():
    records = []
    for ticket in range(64):
        bits = []
        for bit in range(6):
            coin = (ticket >> bit) & 1
            out, record = native((0, 0), (0, 1), coin)
            assert out == (0, 0) and record == (coin, 0)
            bits.append(record[0])
            records.append([ticket, bit, record[0], record[1]])
        assert sum(c << k for k, c in enumerate(bits)) == ticket
    return records

def instrument(weights, law):
    tickets = [e for e, count in enumerate(weights) for _ in range(count)]
    assert len(tickets) == 64
    # flag=(ticket,coin,current in{-1,0,1}); output source on dimension3.
    fine = np.zeros((384, 3, 3))
    counts = np.zeros((3, 3), dtype=int)
    records = []
    for initial in range(3):
        q = tuple(int(i == initial) for i in range(3))
        for ticket, e in enumerate(tickets):
            a, b = EDGES[e]
            for coin in (0, 1):
                if law == 'native':
                    out, (_, current) = native(q, (a, b), coin)
                else:
                    out = list(q)
                    current = q[a] - q[b] if law == 'move_always' else 0
                    out[a] -= current; out[b] += current; out = tuple(out)
                # Source position is reconstructed from public prep and current.
                read = list(q); read[a] -= current; read[b] += current
                assert tuple(read) == out
                final = read.index(1)
                flag = (ticket * 2 + coin) * 3 + current + 1
                fine[flag, final, initial] = 1 / np.sqrt(128)
                counts[final, initial] += 1
                records.append([initial, ticket, coin, current, final])
    effect = np.einsum('rai,raj->ij', fine, fine)
    assert np.max(abs(effect - np.eye(3))) < 1e-12
    # A fixed uniform coherent analyzer is a complete Householder basis.
    w = np.ones(384) / np.sqrt(384)
    u = w.copy(); u[0] -= 1; u /= np.linalg.norm(u)
    transformed = fine - 2 * u[:, None, None] * np.einsum('r,rij->ij', u, fine)
    a = transformed[0]
    p = counts / 128
    assert np.max(abs(a - p / np.sqrt(3))) < 1e-12
    complete_effect = np.einsum('rai,raj->ij', transformed, transformed)
    probabilities = np.einsum('rij,j->ri', transformed, INITIAL)
    probabilities = np.sum(probabilities ** 2, axis=1)
    assert abs(sum(probabilities) - 1) < 1e-12
    assert np.max(abs(complete_effect - np.eye(3))) < 1e-12
    # Reading the flag destroys exactly those cross terms the eraser uses.
    dephased = np.einsum('rij,j->ri', fine, INITIAL)
    dephased = np.einsum('ri,rj->ij', dephased, dephased) / 384
    coherent = a @ INITIAL
    return p, records, {'weights': weights, 'all_384_analyzer_probabilities_on_initial_ray': list(map(num, probabilities)),
                        'complete_effect_error': num(np.max(abs(complete_effect - np.eye(3)))),
                        'success_map_error': num(np.max(abs(a - p / np.sqrt(3)))),
                        'coherent_success_probability': num(coherent @ coherent),
                        'flag_committed_success_probability': num(np.trace(dephased)),
                        'flag_committed_normalized_initial_ray_weight': num(INITIAL @ dephased @ INITIAL / np.trace(dephased))}

def transport(matrices, initial=INITIAL):
    # Normalization is numerical bookkeeping only. All removed norm factors
    # remain in log_norm and source-erasure probability3^-number_of_events.
    v = initial.copy()
    log_norm = 0.
    for p in matrices:
        v = p @ v
        # This mathematical zero-sum restriction is invariant exactly. We DO
        # NOT project numerically: roundoff leakage into constant mode is
        # reported below and independently checked in exact 2D coordinates.
        length = np.linalg.norm(v)
        if length == 0:
            return np.zeros(3), -math.inf
        v /= length
        log_norm += math.log(length)
    return v, log_norm

def transport_optical(matrices):
    # Orthonormal chart of an EXACT invariant space, not an applied physical
    # projector. Avoid amplified floating constant-mode roundoff for long runs.
    b = np.array([[1., 1.], [-1., 1.], [0., -2.]]) / np.array([np.sqrt(2), np.sqrt(6)])
    v = np.array([0., 1.]); log_norm = 0.
    for p in matrices:
        reduced = b.T @ p @ b
        v = reduced @ v
        length = np.linalg.norm(v)
        assert length > 0
        v /= length; log_norm += math.log(length)
    return b @ v, log_norm

def interference(matrices, reference):
    a, la = transport_optical(matrices)
    b, lb = transport_optical([reference] * len(matrices))
    ratio = math.exp(-abs(la - lb))
    visibility = 2 * ratio * float(a @ b) / (1 + ratio * ratio)
    logsuccess = 2 * max(la, lb) + math.log1p(ratio * ratio) - math.log(2) - len(matrices) * math.log(3)
    return {'source_events': len(matrices), 'normalized_source_overlap_with_initial': num(a @ INITIAL),
            'optical_constant_component': num(abs(sum(a))),
            'relative_interference_visibility': num(visibility),
            'conditional_minus_probability': num((1 - visibility) / 2),
            'log10_complete_success_probability': num(logsuccess / math.log(10)),
            'log10_original_attempts_per_common_success': num(-logsuccess / math.log(10)),
            'log_norm_loop_without_eraser': num(la), 'log_norm_reference_without_eraser': num(lb)}

def execute():
    scheduler = sampler_records()
    weights = weight_loop()
    all_records = {}
    instruments = {}
    products = {}
    for law in ('native', 'move_always', 'hold'):
        matrices = []
        branch_records = []
        diagnostics = []
        for w in weights:
            p, rec, diag = instrument(w, law)
            matrices.append(p); branch_records.append(rec); diagnostics.append(diag)
        reference, refrecords, refdiag = instrument((32, 16, 16), law)
        # This reference setting is also last point of loop; records are retained
        # separately so the stated source call count includes the actual rerun.
        all_records[law] = {'loop': branch_records, 'reference': refrecords}
        instruments[law] = {'loop': diagnostics, 'reference': refdiag}
        tests = {}
        for dwell in REPETITIONS:
            sequence = [p for p in matrices for _ in range(dwell)]
            tests[str(dwell)] = {'one_loop': interference(sequence, reference),
                                 'two_loops': interference(sequence + sequence, reference),
                                 'retraced_loop': interference(sequence + sequence[::-1], reference),
                                 'fixed_scheduler': interference([reference] * len(sequence), reference)}
        products[law] = tests
    main = products['native']['4']
    checks = {'same_source_native_one_loop_negative': main['one_loop']['normalized_source_overlap_with_initial'] < -.99,
              'same_source_native_two_loops_positive': main['two_loops']['normalized_source_overlap_with_initial'] > .99,
              'retracing_positive': main['retraced_loop']['normalized_source_overlap_with_initial'] > .99,
              'relative_interference_detects_negative': main['one_loop']['relative_interference_visibility'] < 0,
              'hold_control_no_negative_sign': products['hold']['4']['one_loop']['normalized_source_overlap_with_initial'] > .99,
              'fixed_scheduler_no_negative_sign': main['fixed_scheduler']['normalized_source_overlap_with_initial'] > .99,
              'phase_eraser_changes_instrument': abs(instruments['native']['reference']['coherent_success_probability'] - instruments['native']['reference']['flag_committed_success_probability']) > .01}
    return {'candidate': 'A-native-triangle-optical-holonomy-v1', 'class': 'N+',
            'native_calls': 384 + 49 * 384, 'changed_law_calls': 2 * 49 * 384,
            'all_diagnostic_checks': all(checks.values()), 'checks': checks,
            'weight_loop': weights, 'products': products, 'instruments': instruments,
            'physical_spatial_rotation_derived': False, 'native_eraser_compiled': False,
            'row_criterion_passed': False, 'nature_values_used': [],
            'scope': 'Conditional coherent native transfer holonomy around externally programmed scheduler loop; full costs retained.'}, {'native_scheduler': scheduler, 'source': all_records}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--record', action='store_true'); args = parser.parse_args()
    result, records = execute()
    if args.record:
        for name, obj in [('RECEIPT.json', result), ('RECORDS.json', records)]:
            with (HERE / name).open('x') as f: f.write(canonical(obj))
    else:
        assert canonical(records) == (HERE / 'RECORDS.json').read_text()
        assert canonical(result) == (HERE / 'RECEIPT.json').read_text()
    print(canonical({'checks': result['checks'], 'native_calls': result['native_calls'],
                     'native_dwell4': result['products']['native']['4']}), end='')
