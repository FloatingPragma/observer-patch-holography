"""Independent finite-field certificates for the C12 persistence classification.

Only the auditor's prior direct Gauss/CAR reconstruction is imported. Producer
Python is neither imported nor executed. Full column rank modulo an odd prime
certifies full rank over Q for the exact binary64 rational matrix.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from joint_source_audit import EDGES, enumerate_basis, operators, sign

PRIME = 65521


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def modulo(matrix):
    out = np.zeros(matrix.shape, dtype=np.int64)
    for index, value in np.ndenumerate(matrix):
        numerator, denominator = float(value).as_integer_ratio()
        out[index] = (numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME
    return out


def rank_mod(matrix):
    a = matrix.copy() % PRIME
    rank = 0
    pivots = []
    for col in range(a.shape[1]):
        rows = np.flatnonzero(a[rank:, col])
        if not len(rows):
            continue
        pivot = rank + int(rows[0])
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, PRIME) % PRIME
        for row in range(rank + 1, a.shape[0]):
            a[row] = (a[row] - a[row, col] * a[rank]) % PRIME
        pivots.append(col)
        rank += 1
        if rank == a.shape[0]:
            break
    return rank, pivots


def inertia_below(matrix, boundary):
    """Exact LDL congruence; chosen non-root boundaries have nonzero pivots."""
    n = len(matrix)
    a = [[Fraction(float(matrix[i, j])) - (boundary if i == j else 0)
          for j in range(n)] for i in range(n)]
    negative = 0
    for k in range(n):
        pivot = a[k][k]
        if not pivot:
            raise ValueError('zero LDL pivot: select a different rational boundary')
        negative += pivot < 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] -= a[i][k] * a[k][j] / pivot
    return negative


def native_maps(states):
    index = {state: i for i, state in enumerate(states)}
    maps = {}
    moving_columns = set()
    for e, (a, b) in enumerate(EDGES):
        for coin in (0, 1):
            for law in ('C', 'E'):
                for value in (-1, 0, 1):
                    maps[e, coin, law, value] = np.zeros((384, 384), dtype=np.int64)
            for col, (mask, f) in enumerate(states):
                ua, ub = (mask >> a) & 1, (mask >> b) & 1
                current = ua - (ua + ub + coin) // 2
                dest_mask = mask ^ (1 << a) ^ (1 << b) if current else mask
                shifted = list(f)
                shifted[2*e] = (shifted[2*e] + current) % 6
                shifted[2*e+1] = (shifted[2*e+1] - current) % 6
                dest = index[dest_mask, tuple(shifted)]
                phase = sign(mask, a, b) if current else 1
                for law, value in (('C', current), ('E', ua-ub)):
                    maps[e, coin, law, value][dest, col] = phase
                if current:
                    moving_columns.add(col)
    return maps, moving_columns


def run(c10, c12):
    paths = [c10/'JOINT_SOURCE.json', c10/'joint_source_states.npz',
             c12/'pure_persistence.py', c12/'PURE_PERSISTENCE.json', c12/'PURE_PERSISTENCE.md']
    root = c10.resolve().parents[1]
    before = {str(p.resolve().relative_to(root)): digest(p) for p in paths}
    producer = json.loads(paths[3].read_text())
    states = enumerate_basis()
    h, field = operators(states)
    with np.load(paths[1], allow_pickle=False) as archived:
        actual_states = [(int(m), tuple(map(int, f))) for m, f in
                         zip(archived['occupation_masks'], archived['half_records'])]
        assert len(actual_states) == len(set(actual_states)) == len(states)
        assert set(states) == set(actual_states)
        index = {state: i for i, state in enumerate(states)}
        permutation = [index[state] for state in actual_states]
        h = h[np.ix_(permutation, permutation)]
        field = field[np.ix_(permutation, permutation)]
        states = actual_states
        assert np.array_equal(h, archived['hamiltonian'])
    numbers = np.array([m.bit_count() for m, _ in states])
    dark = np.array([(m & 7) in (0, 7) for m, _ in states])
    maps, moving_columns = native_maps(states)
    nilpotent = []
    for (e, coin, law, value), matrix in maps.items():
        current = value if law == 'C' else (value-(2*coin-1)*abs(value))//2
        if current:
            nilpotent.append(not np.any(matrix @ matrix))
            assert np.max(np.count_nonzero(matrix, axis=0)) <= 1
            assert np.max(np.count_nonzero(matrix, axis=1)) <= 1
        assert np.array_equal(matrix[:, dark], np.eye(384, dtype=np.int64)[:, dark]
                              * (value == 0))
    rows = []
    endpoints = []
    for n in range(7):
        ids = np.flatnonzero((numbers == n) & dark)
        outside = np.flatnonzero((numbers == n) & ~dark)
        coupling = modulo(h[np.ix_(outside, ids)])
        internal = modulo(h[np.ix_(ids, ids)])
        rank1, _ = rank_mod(coupling)
        rank2, pivots = rank_mod(np.vstack((coupling, coupling @ internal % PRIME)))
        invariant_dimension = len(ids) if not len(outside) else 0
        assert not len(outside) or rank2 == len(ids)
        rows.append({'number': n, 'native_support_dimension': len(ids),
                     'coupling_rank_mod_prime': rank1,
                     'two_constraint_rank_mod_prime': rank2, 'pivot_columns': pivots,
                     'maximal_H_invariant_dimension': invariant_dimension})
        if n in (0, 6):
            block = h[np.ix_(ids, ids)]
            for value in np.linalg.eigvalsh(block):
                center = round(float(value)*10**7)
                lower, upper = Fraction(center-2, 10**7), Fraction(center+2, 10**7)
                count = inertia_below(block, upper)-inertia_below(block, lower)
                assert count == 1
                endpoints.append({'number': n, 'lower': str(lower), 'upper': str(upper),
                                  'exact_eigenvalue_count': count})
    ordered = sorted(endpoints, key=lambda item: Fraction(item['lower']))
    assert all(Fraction(a['upper']) < Fraction(b['lower']) for a, b in zip(ordered, ordered[1:]))
    psi = np.array([1/np.sqrt(18) if m in (8, 16, 32) else 0 for m, _ in states])
    energy = float(psi @ h @ psi)
    variance = float(np.linalg.norm(h @ psi - energy*psi)**2)
    evolved = expm(-.37j*h) @ psi
    leakage = float(np.sum(abs(evolved[~dark])**2))
    control = producer['former_cycle_control']
    checks = {
        'independently_reconstructed_actual_H_exact': True,
        'moving_maps_nilpotent_and_injective': all(nilpotent),
        'moving_domains_exhaust_native_dark_complement': moving_columns == set(np.flatnonzero(~dark)),
        'all_native_maps_scalar_on_dark_support': True,
        'full_rank_certifies_no_middle_number_H_invariant_subspace': all(
            row['two_constraint_rank_mod_prime'] == row['native_support_dimension'] for row in rows[1:6]),
        'producer_dimensions_match': all(row['maximal_H_invariant_dimension'] == target['exact_H_invariant_pure_fixed_dimension']
                                         for row, target in zip(rows, producer['numbers'])),
        'twelve_simple_distinct_endpoint_eigenvalues_exactly_enclosed': len(endpoints) == 12,
        'former_cycle_variance_and_escape_match': max(abs(variance-control['full_H_energy_variance']),
                                                     abs(leakage-control['port_occupation_leak_probability'])) < 1e-12,
        'audited_artifacts_unchanged': before == {str(p.resolve().relative_to(root)): digest(p) for p in paths},
    }
    assert all(checks.values()), checks
    return {'schema': 'oph.openwave.persistence_independent_audit.v1', 'checks': checks,
            'finite_field_prime': PRIME, 'number_sector_certificates': rows,
            'endpoint_exact_spectral_intervals': endpoints,
            'former_cycle': {'variance': variance, 'port_leakage_t_037': leakage},
            'source_sha256': before, 'auditor_sha256': digest(Path(__file__)),
            'reference_implementation_sha256': digest(Path(__file__).with_name('joint_source_audit.py')),
            'scope': 'Frozen C10 all-number triangle. Native branchwise pure persistence and H-invariant support only. Algebraic-time extension requires the hypotheses in PERSISTENCE_PROOF.md. No approximate lifetime, mixed-state, particle or platform acceptance claim.',
            'official_rows_changed': []}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=Path, required=True)
    parser.add_argument('--c12', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.c10, args.c12)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'number_sectors': result['number_sector_certificates']}))
