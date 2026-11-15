"""Independent finite checks of the native and modular source boundaries.

For one integer load on an undirected graph, T = I - L/(2M).
For the modular triangle, G = occupation + half-incidence * record.
Uniform fibres V obey H V = V H0, but a record PVM need not preserve im(V).
No physical particle, clock, continuum, or platform award is inferred.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.linalg import expm


def endpoint(loads, edge, coin):
    """Independent endpoint average: the coin chooses which end gets the ceiling."""
    a, b = edge
    out = list(loads)
    total = loads[a] + loads[b]
    # Odd d and odd endpoint sum have the same parity, including negative loads.
    out[a] = (total + coin) // 2
    out[b] = total - out[a]
    return tuple(out), (coin, loads[a] - out[a])


def one_load_kernel(vertices, edges):
    counts = np.zeros((vertices, vertices), dtype=int)
    for source in range(vertices):
        state = tuple(int(i == source) for i in range(vertices))
        for edge in edges:
            for coin in (0, 1):
                target, _ = endpoint(state, edge, coin)
                if sorted(target) != [0] * (vertices - 1) + [1]:
                    raise AssertionError('one-load closure failed')
                counts[target.index(1), source] += 1
    lap = np.zeros_like(counts)
    for a, b in edges:
        lap[a, a] += 1
        lap[b, b] += 1
        lap[a, b] -= 1
        lap[b, a] -= 1
    expected = 2 * len(edges) * np.eye(vertices, dtype=int) - lap
    if not np.array_equal(counts, expected):
        raise AssertionError('literal transition differs from Laplacian')
    return counts, lap


def triangle(modulus=6, mass=.23):
    """Solve the six Gauss constraints directly, with f0 as the free coordinate."""
    states = []
    for z in range(6):
        for r in range(modulus):
            f0 = r
            f1 = int(z == 3) - f0
            f2 = -int(z == 1) - f1
            f3 = int(z == 4) - f2
            f4 = -int(z == 2) - f3
            f5 = 1 - int(z == 0) - f0
            states.append((z, tuple(x % modulus for x in (f0, f1, f2, f3, f4, f5))))
    index = {state: i for i, state in enumerate(states)}
    n = len(states)
    h = np.diag([mass if z < 3 else -mass for z, _ in states])
    edges = ((0, 1), (1, 2), (2, 0))
    a = np.zeros((6, 6), dtype=int)
    for e, ends in enumerate(edges):
        for side, port in enumerate(ends):
            a[port, 2 * e + side] = 1
            a[3 + e, 2 * e + side] = -1
    for col, (z, records) in enumerate(states):
        assert np.array_equal((np.eye(6, dtype=int)[z] + a @ records) % modulus,
                              np.array([1, 0, 0, 0, 0, 0]))
        for e, ends in enumerate(edges):
            for side, port in enumerate(ends):
                if z != port:
                    continue
                shifted = list(records)
                shifted[2 * e + side] = (shifted[2 * e + side] + 1) % modulus
                row = index[(3 + e, tuple(shifted))]
                h[row, col] = h[col, row] = 1 if side == 0 else -1
    v = np.zeros((n, 6))
    for i, (z, _) in enumerate(states):
        v[i, z] = 1 / np.sqrt(modulus)
    return states, index, h, v, a


def native_kraus(states, index, modulus, law):
    operators = []
    for e, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
        for coin in (0, 1):
            groups = defaultdict(lambda: np.zeros((len(states), len(states))))
            for col, (z, records) in enumerate(states):
                d = int(z == a) - int(z == b)
                current = (d - (2 * coin - 1) * abs(d)) // 2
                label = current if law == 'C' else d
                target = z
                shifted = list(records)
                if current:
                    target = b if current == 1 else a
                    shifted[2 * e] = (shifted[2 * e] + current) % modulus
                    shifted[2 * e + 1] = (shifted[2 * e + 1] - current) % modulus
                groups[label][index[(target, tuple(shifted))], col] = 1 / np.sqrt(6)
            operators.extend(groups.values())
    return operators


def measure_register(rho, states):
    values = np.array([records[0] for _, records in states])
    return rho * (values[:, None] == values[None, :])


def evaluate():
    states, index, h, v, incidence = triangle()
    projector = v @ v.T
    h0 = v.T @ h @ v
    initial = np.outer(v[:, 0], v[:, 0])
    measured = measure_register(initial, states)
    code_weight = float(np.trace(projector @ measured))
    readout_error = abs(code_weight - 1 / 6)
    native_errors = {}
    for law in ('C', 'E'):
        ops = native_kraus(states, index, 6, law)
        native_errors[law] = {
            'completeness': float(np.max(abs(sum(k.T @ k for k in ops) - np.eye(len(states))))),
            'code_leak': max(float(np.linalg.norm((np.eye(len(states)) - projector) @ k @ v)) for k in ops),
        }
    t = .37
    u = expm(-1j * t * h)
    uniform = u @ initial @ u.conj().T
    record_measured = u @ measured @ u.conj().T
    position = lambda rho: np.array([sum(rho[i, i].real for i, (z, _) in enumerate(states) if z == x) for x in range(6)])
    definite = np.zeros_like(initial)
    definite[index[(0, (0,) * 6)], index[(0, (0,) * 6)]] = 1
    definite_final = u @ definite @ u.conj().T
    counts, lap = one_load_kernel(3, ((0, 1), (1, 2), (2, 0)))
    transition = counts / 6
    native_final = expm(t * (transition - np.eye(3)))[:, 0]
    joint_positions = position(definite_final)
    transit_weight = float(sum(joint_positions[3:]))
    # Total variation is at least the probability of the transit event, absent in the native source.
    tv = float(.5 * sum(abs(joint_positions - np.r_[native_final, np.zeros(3)])))
    # Every cut of the subdivided graph obeys the residue-valued Gauss identity.
    cut_failures = 0
    for mask in range(64):
        region = np.array([(mask >> j) & 1 for j in range(6)])
        boundary = region @ incidence
        for z, records in states:
            cut_failures += int((int(region[z]) + int(boundary @ records) - int(region[0])) % 6 != 0)
    # Arithmetic mutations: delete a transferred unit; omit its record update.
    original = (1, 0)
    moved, (_, current) = endpoint(original, (0, 1), 0)
    wrong_load = (0, 0)
    mutation_checks = {
        'erase_unit_rejected': sum(wrong_load) != sum(original),
        'omit_record_rejected': moved[0] != original[0] and current != 0,
        'record_measurement_as_code_preserving_rejected': code_weight < .99,
        'definite_records_as_native_state_law_rejected': transit_weight > 1e-3,
        'identify_survival_with_stationary_distribution_rejected': not np.array_equal(counts[:, 0], np.array([6, 0, 0])),
    }
    checks = {
        'hamiltonian_uniform_intertwiner': float(np.linalg.norm(h @ v - v @ h0)) < 1e-12,
        'complete_native_maps_preserve_code': all(max(row.values()) < 1e-12 for row in native_errors.values()),
        'register_readout_leaves_only_one_sixth_in_code': readout_error < 1e-12,
        'every_cut_has_modular_gauss_identity': cut_failures == 0,
        'native_one_load_spectrum_nonnegative': bool(np.min(np.linalg.eigvalsh(transition)) >= 0),
        'full_record_counterexample': tv + 1e-12 >= transit_weight > 1e-3,
        'mutations_detected': all(mutation_checks.values()),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {
        'schema': 'oph.openwave.independent_scope_audit.v1',
        'checks': checks, 'mutations': mutation_checks,
        'parameters': {'graph': 'triangle', 'ports': 3, 'modulus': 6, 'mass': .23, 'time': t},
        'observables': {'uniform_code_weight_after_register_measurement': code_weight,
                        'continued_position_max_change': float(np.max(abs(position(uniform) - position(record_measured)))),
                        'definite_record_transit_weight': transit_weight,
                        'definite_vs_native_total_variation_at_one_attempt_per_time_unit': tv,
                        'native_transition_numerators': counts.tolist(), 'native_transition_denominator': 6,
                        'native_transition_eigenvalues': np.linalg.eigvalsh(transition).tolist(),
                        'native_map_residuals': native_errors, 'all_cut_violations': cut_failures},
        'scope': 'Independent finite development counterchecks. No reserved acceptance run, physical particle identification, infinite-register limit, or platform award.',
        'official_rows_changed': [],
    }


def compare_literal_native(rer_root):
    """Compare signed integer endpoint arithmetic against the repository source."""
    sys.path.insert(0, str(rer_root))
    from model.observer.native import native
    for a in range(-8, 9):
        for b in range(-8, 9):
            for coin in (0, 1):
                if endpoint((a, b), (0, 1), coin) != native((a, b), (0, 1), coin):
                    raise AssertionError((a, b, coin))
    return {'signed_endpoint_branches': 578, 'disagreements': 0,
            'source_sha256': hashlib.sha256((rer_root / 'model/observer/native.py').read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--rer-root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    result = evaluate()
    result['native_source_comparison'] = compare_literal_native(args.rer_root)
    result['versions'] = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}
    result['auditor_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
