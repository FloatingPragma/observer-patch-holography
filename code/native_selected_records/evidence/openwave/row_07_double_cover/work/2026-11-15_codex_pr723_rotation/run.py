#!/usr/bin/env python3
"""Execute frozen finite study. Exact records, numerical spectral diagnostics."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import source

HERE = Path(__file__).resolve().parent
TIMES = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)

def canonical(x):
    return json.dumps(x, indent=2, sort_keys=True, allow_nan=False) + '\n'

def number(x):
    return float(format(float(x), '.12g'))

def geometry():
    f = (1 + 5 ** .5) / 2
    return np.array([(-1, f, 0), (1, f, 0), (-1, -f, 0), (1, -f, 0),
                     (0, -1, f), (0, 1, f), (0, -1, -f), (0, 1, -f),
                     (f, 0, -1), (f, 0, 1), (-f, 0, -1), (-f, 0, 1)])

def label(p, r):
    return source.INDEX[tuple(int(i == p) - int(i == r) for i in range(12))]

def moments(records):
    points = geometry()
    anti = [int(np.argmin(np.linalg.norm(points + v, axis=1))) for v in points]
    stats = []
    for p, r in enumerate(anti):
        vectors = [j * (points[source.SEAMS[e][1]] - points[source.SEAMS[e][0]])
                   for z, e, c, j, out in records if z == label(p, r)]
        x = np.array(vectors)
        mean = np.mean(x, axis=0)
        cov = x.T @ x / len(x) - np.outer(mean, mean)
        v = points[p] / np.linalg.norm(points[p])
        par = v @ cov @ v
        per = (np.trace(cov) - par) / 2
        fit = per * np.eye(3) + (par - per) * np.outer(v, v)
        stats.append({'source': [p, r], 'mean': list(map(number, mean)),
                      'covariance': [[number(a) for a in row] for row in cov],
                      'parallel': number(par), 'transverse': number(per),
                      'shape_residual': number(np.max(abs(cov - fit)))})
    return stats

def transition(records):
    counts = np.zeros((133, 133), dtype=np.int64)
    for z, e, c, j, out in records:
        counts[z, out] += 1
    assert np.all(counts.sum(axis=1) == 60)
    return counts

def evolution(counts):
    p = counts / 60
    seeds = [label(0, 3), label(3, 0), label(1, 2), label(2, 1)]
    rows = np.eye(133)[seeds]
    results = []
    previous = 0
    for t in TIMES:
        rows = rows @ np.linalg.matrix_power(p, t - previous)
        previous = t
        tv_sign = np.sum(abs(rows[0] - rows[1])) / 2
        tv_axis = np.sum(abs((rows[0] + rows[1] - rows[2] - rows[3]) / 2)) / 2
        survival = 1 - rows[0, 0]
        # The average here diagnoses axis information after sign erasure. It
        # is explicitly NOT credited as a naturally apolar field decoder.
        results.append({'attempts': t, 'antipodal_total_variation': number(tv_sign),
                        'axis_total_variation_after_explicit_even_projection': number(tv_axis),
                        'survival': number(survival),
                        'numerical_candidate_window': bool(tv_sign <= 1e-6 and tv_axis >= 1e-3
                                                           and survival >= 1e-3)})
    q = p[1:, 1:]
    assert np.max(abs(q - q.T)) < 1e-15
    vals, vecs = np.linalg.eigh(q)
    odd = np.zeros(132)
    even = np.zeros(132)
    odd[seeds[0] - 1], odd[seeds[1] - 1] = 1, -1
    for i, sign in zip(seeds, [1, 1, -1, -1]):
        even[i - 1] += sign / 2
    def active_radius(v):
        weights = (v @ vecs) ** 2
        return number(max(vals[weights > 1e-18]))
    return {'checkpoints': results, 'smallest_transient_eigenvalue': number(vals[0]),
            'largest_transient_eigenvalue': number(vals[-1]),
            'antipodal_difference_active_radius': active_radius(odd),
            'axis_difference_active_radius': active_radius(even),
            'any_candidate_window': any(r['numerical_candidate_window'] for r in results),
            'numerical_zero_policy': 'numbers below resolution are not exact identifications; use PROOF.md'}

def fine_turn_comparator():
    """Conditional diagnostic. Rotation and analyzer are declared apparatus."""
    v = geometry()
    axis = v[0] / np.linalg.norm(v[0])
    theta = 2 * np.pi / 5
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
    r = np.cos(theta) * np.eye(3) + (1 - np.cos(theta)) * np.outer(axis, axis) + np.sin(theta) * cross
    permutation = [int(np.argmin(np.linalg.norm(v - r @ x, axis=1))) for x in v]
    assert len(set(permutation)) == 12
    assert max(np.min(np.linalg.norm(v - r @ x, axis=1)) for x in v) < 1e-12
    assert all(tuple(sorted((permutation[a], permutation[b]))) in source.SEAMS for a, b in source.SEAMS)
    state_perm = []
    for q in source.STATES:
        out = [0] * 12
        for k, value in enumerate(q):
            out[permutation[k]] = value
        state_perm.append(source.INDEX[tuple(out)])
    history = list(range(133))
    outcome = []
    for turn in range(1, 11):
        history = [state_perm[x] for x in history]
        # For computational-basis preparations, overlap is 0 or +1. Analyzer
        # probability is (1-overlap)/2. Negative overlap is absent in this lift.
        minus = [(1 - int(x == start)) / 2 for start, x in enumerate(history)]
        outcome.append({'steps_72_degrees': turn, 'mean_minus_probability': number(np.mean(minus)),
                        'all_columns_return': history == list(range(133))})
    return {'port_permutation': permutation, 'all_133_columns': outcome,
            'full_turn_minus_probability': outcome[4]['mean_minus_probability'],
            'two_turn_minus_probability': outcome[9]['mean_minus_probability'],
            'negative_phase_plant_minus_probability': (1 - (-1)) / 2,
            'dephased_plant_minus_probability': .5,
            'status': 'conditional comparator; no native analyzer or phase-generating intervention claimed'}

def run():
    records = {law: source.all_records(law) for law in ('native', 'move_always', 'hold')}
    counts = {law: transition(table) for law, table in records.items()}
    stats = {law: moments(table) for law, table in records.items()}
    phi = (1 + 5 ** .5) / 2
    formula_error = max(abs(stats['native'][0]['parallel'] - (3 - phi) / 9),
                        abs(stats['native'][0]['transverse'] - (phi + 2) / 15))
    return_error = {law: max(abs(np.array(stats[law][0]['mean']) + stats[law][3]['mean']).max(),
                            abs(np.array(stats[law][0]['covariance']) - stats[law][3]['covariance']).max())
                    for law in ('native', 'move_always')}
    dynamics = {law: evolution(c) for law, c in counts.items()}
    compare = fine_turn_comparator()
    checks = {'complete_branch_rows': all(len(r) == 7980 for r in records.values()),
              'native_formulas_reproduced': formula_error < 1e-10,
              'alternate_law_keeps_old_return_claim': all(e < 1e-10 for e in return_error.values()),
              'alternate_law_changes_magnitudes': abs(stats['native'][0]['parallel'] - stats['move_always'][0]['parallel']) > .05,
              'hold_retains_antipodal_distinction': all(r['antipodal_total_variation'] == 1 for r in dynamics['hold']['checkpoints']),
              'phase_detector_sees_negative_plant': compare['negative_phase_plant_minus_probability'] == 1,
              'positive_port_lift_has_no_spinor_sign': compare['full_turn_minus_probability'] == 0,
              'native_finite_time_mean_transfer_invertible': float(np.linalg.eigvalsh(np.eye(12) - graph_laplacian() / 60)[0]) > 0}
    result = {'schema': 'oph.row07.pr723.rotation.v1', 'candidate': 'A-signed-record-rotation-audit-v1',
              'native_endpoint_calls': len(records['native']), 'changed_law_calls': len(records['move_always']) + len(records['hold']),
              'record_count': sum(map(len, records.values())), 'checks': checks,
              'all_diagnostic_checks': all(checks.values()), 'moment_results': stats,
              'dynamics': dynamics, 'fine_turn': compare, 'row_criterion_passed': False,
              'memory_phase_used_by_classical_observables': False,
              'official_status': 'untested', 'G6_awarded': False, 'G9_awarded': False}
    return result, records

def graph_laplacian():
    lap = np.zeros((12, 12))
    for a, b in source.SEAMS:
        lap[a, a] += 1; lap[b, b] += 1; lap[a, b] -= 1; lap[b, a] -= 1
    return lap

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    result, records = run()
    if args.record:
        for name, obj in [('RECEIPT.json', result), ('RECORDS.json', records)]:
            with (HERE / name).open('x') as stream:
                stream.write(canonical(obj))
    else:
        assert canonical(records) == (HERE / 'RECORDS.json').read_text()
        assert canonical(result) == (HERE / 'RECEIPT.json').read_text()
    print(canonical({k: result[k] for k in ('native_endpoint_calls', 'record_count', 'all_diagnostic_checks', 'row_criterion_passed')}), end='')
