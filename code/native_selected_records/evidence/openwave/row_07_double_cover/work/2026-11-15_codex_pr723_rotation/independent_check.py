#!/usr/bin/env python3
"""Separately implemented record checker. Imports no producer/source functions."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def run():
    edges = json.loads((HERE / 'seams.json').read_text())
    saved = json.loads((HERE / 'RECORDS.json').read_text())
    receipt = json.loads((HERE / 'RECEIPT.json').read_text())
    states = [(0,) * 12]
    for p in range(12):
        for m in range(12):
            if p == m:
                continue
            x = [0] * 12; x[p] = 1; x[m] = -1; states.append(tuple(x))
    index = {s: i for i, s in enumerate(states)}
    checked = 0
    diagonal_margins = {}
    for law in ('native', 'move_always', 'hold'):
        expected = []
        counts = np.zeros((133, 133), dtype=int)
        for z, x in enumerate(states):
            for e, (a, b) in enumerate(edges):
                for c in (0, 1):
                    y = list(x)
                    if law == 'native':
                        total = x[a] + x[b]
                        low, high = total // 2, -((-total) // 2)
                        y[a], y[b] = (low, high) if c == 0 else (high, low)
                    elif law == 'move_always':
                        y[a], y[b] = x[b], x[a]
                    j = x[a] - y[a]
                    expected.append([z, e, c, j, index[tuple(y)]])
                    counts[z, index[tuple(y)]] += 1
                    checked += 1
        assert expected == saved[law]
        assert np.all(counts.sum(axis=1) == 60)
        margin = int(np.min(2 * np.diag(counts) - counts.sum(axis=1)))
        assert margin > 0
        diagonal_margins[law] = margin
        p = counts / 60
        positions = []
        for plus, minus in ((0, 3), (3, 0), (1, 2), (2, 1)):
            x = [0] * 12; x[plus] = 1; x[minus] = -1; positions.append(index[tuple(x)])
        for result in receipt['dynamics'][law]['checkpoints']:
            # Entire independent matrix power per checkpoint, instead of the
            # producer's incremental selected-row propagation.
            pt = np.linalg.matrix_power(p, result['attempts'])
            a, b, c, d = (pt[i] for i in positions)
            sign = np.linalg.norm(a - b, ord=1) / 2
            axis = np.linalg.norm(a + b - c - d, ord=1) / 4
            assert abs(sign - result['antipodal_total_variation']) < 2e-10
            assert abs(axis - result['axis_total_variation_after_explicit_even_projection']) < 2e-10
    turn = receipt['fine_turn']['port_permutation']
    v = list(range(12))
    for _ in range(5):
        v = [turn[x] for x in v]
    assert v == list(range(12))
    assert receipt['fine_turn']['full_turn_minus_probability'] == 0
    assert receipt['fine_turn']['negative_phase_plant_minus_probability'] == 1
    return {'checked_source_branches': checked, 'strict_diagonal_dominance_integer_margin': diagonal_margins,
            'ordinary_rotation_fifth_power_identity': True,
            'comparison_implementation': 'separate endpoint averaging and full matrix powers; same Codex group',
            'no_package_imports': True, 'all_checks': True}

if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
