#!/usr/bin/env python3
"""Exact rational products, separate endpoint averaging, no producer imports."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (1, 2), (2, 0))
BASIS = ((1, -1, 0), (1, 1, -2))

def rational_matrix(weights, law):
    counts = [[0] * 3 for _ in range(3)]
    records = []
    tickets = [i for i, n in enumerate(weights) for _ in range(n)]
    for z in range(3):
        x = [int(k == z) for k in range(3)]
        for ticket, edge in enumerate(tickets):
            a, b = EDGES[edge]
            for coin in (0, 1):
                y = x.copy()
                if law == 'native':
                    s = x[a] + x[b]
                    y[a], y[b] = (s // 2, (s + 1) // 2) if coin == 0 else ((s + 1) // 2, s // 2)
                elif law == 'move_always':
                    y[a], y[b] = x[b], x[a]
                out = y.index(1); j = x[a] - y[a]
                counts[out][z] += 1
                records.append([z, ticket, coin, j, out])
    reduced = [[F(sum(BASIS[i][a] * counts[a][b] * BASIS[k][b]
                      for a in range(3) for b in range(3)), 128 * (2 if i == 0 else 6))
                for k in range(2)] for i in range(2)]
    return reduced, records

def product(sequence):
    a, b = F(0), F(1)
    for p in sequence:
        a, b = p[0][0] * a + p[0][1] * b, p[1][0] * a + p[1][1] * b
    return a, b

def diagnostics(sequence, reference):
    a, b = product(sequence)
    ar, br = product([reference] * len(sequence))
    assert ar == 0
    norm = a * a / 3 + b * b
    total = norm + br * br
    visibility = 2 * b * br / total
    success = total / (2 * 3 ** len(sequence))
    # Logs of big integers remain finite; no float conversion of tiny success.
    log10success = (math.log(success.numerator) - math.log(success.denominator)) / math.log(10)
    overlap = float(b / (abs(b) if b else 1)) / math.sqrt(1 + float(a * a / (3 * b * b))) if b else 0.
    return {'source_overlap': overlap, 'visibility': float(visibility), 'log10_success': log10success,
            'amplitude_sign_exact': (b > 0) - (b < 0),
            'exact_product_sha256': hashlib.sha256((str(a) + '\n' + str(b) + '\n').encode()).hexdigest(),
            'numerator_bits': [a.numerator.bit_length(), b.numerator.bit_length()],
            'denominator_bits': [a.denominator.bit_length(), b.denominator.bit_length()]}

def check():
    saved = json.loads((HERE / 'RECEIPT.json').read_text())
    records = json.loads((HERE / 'RECORDS.json').read_text())
    checked = 0; exact = {}
    for law in ('native', 'move_always', 'hold'):
        matrices = []
        for i, weights in enumerate(saved['weight_loop']):
            p, table = rational_matrix(weights, law)
            assert table == records['source'][law]['loop'][i]
            checked += len(table); matrices.append(p)
        reference, table = rational_matrix((32, 16, 16), law)
        assert table == records['source'][law]['reference']; checked += len(table)
        exact[law] = {}
        for dwell in (1, 2, 4):
            sequence = [p for p in matrices for _ in range(dwell)]
            exact[law][str(dwell)] = {}
            for name, seq in [('one_loop', sequence), ('two_loops', sequence + sequence),
                              ('retraced_loop', sequence + sequence[::-1]),
                              ('fixed_scheduler', [reference] * len(sequence))]:
                result = diagnostics(seq, reference)
                expected = saved['products'][law][str(dwell)][name]
                assert abs(result['source_overlap'] - expected['normalized_source_overlap_with_initial']) < 2e-10
                assert abs(result['visibility'] - expected['relative_interference_visibility']) < 2e-10
                assert abs(result['log10_success'] - expected['log10_complete_success_probability']) < 2e-8
                exact[law][str(dwell)][name] = result
    for ticket, bit, coin, current in records['native_scheduler']:
        assert coin == ((ticket >> bit) & 1) and current == 0
    return {'checked_body_branches': checked, 'scheduler_records_checked': 384,
            'method': 'exact rational invariant-chart products from separately averaged native endpoints',
            'same_group': True, 'all_checks': True, 'exact_products': exact}

if __name__ == '__main__':
    print(json.dumps(check(), indent=2, sort_keys=True))
