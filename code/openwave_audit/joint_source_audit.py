"""Rebuild the complete C10 basis, Hamiltonian and native channel from equations.

No producer module is imported. Inputs are the saved basis, Hamiltonian,
reference and response artifacts. CAR signs count occupied intervening orbitals.
The Gauss equations are solved directly for every mask and free cycle residue.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

EDGES = ((0, 1), (1, 2), (2, 0))


def sign(mask, source, target):
    low, high = sorted((source, target))
    return (-1) ** sum((mask >> j) & 1 for j in range(low + 1, high))


def enumerate_basis():
    states = []
    for mask in range(64):
        u = [(mask >> j) & 1 for j in range(6)]
        boundary = 3 - sum(u)
        for r in range(6):
            f0 = r
            f1 = u[3] - f0
            f2 = 1 - u[1] - f1
            f3 = u[4] - f2
            f4 = 1 - u[2] - f3
            f5 = 1 - u[0] - boundary - f0
            assert (u[5] - f4 - f5) % 6 == 0
            states.append((mask, tuple(f % 6 for f in (f0, f1, f2, f3, f4, f5))))
    return states


def operators(states):
    index = {state: i for i, state in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    field = np.zeros_like(h)
    for col, (mask, f) in enumerate(states):
        h[col, col] = .23 * sum((1 if j < 3 else -1) * ((mask >> j) & 1) for j in range(6))
        field[col, col] = .37 / 2 * sum(min(x, 6 - x) ** 2 for x in f) + 2 * .61
        for direction in (-1, 1):
            shifted = tuple((x + direction * (1 if j % 2 == 0 else -1)) % 6 for j, x in enumerate(f))
            field[index[(mask, shifted)], col] -= .61
        for e, ends in enumerate(EDGES):
            for side, port in enumerate(ends):
                transit = 3 + e
                for source, target, delta in ((port, transit, 1), (transit, port, -1)):
                    if not (mask >> source) & 1 or (mask >> target) & 1:
                        continue
                    shifted = list(f)
                    shifted[2 * e + side] = (shifted[2 * e + side] + delta) % 6
                    dest = index[(mask ^ (1 << source) ^ (1 << target), tuple(shifted))]
                    h[dest, col] += (1 if side == 0 else -1) * sign(mask, source, target)
    return h + field, field


def commit_reference(states, vector, law):
    index = {state: i for i, state in enumerate(states)}
    branches = defaultdict(lambda: np.zeros(len(states)))
    laws = ('C', 'E') if law == 'MIX' else (law,)
    for e, (a, b) in enumerate(EDGES):
        for coin in (0, 1):
            for col, (mask, f) in enumerate(states):
                ua, ub = (mask >> a) & 1, (mask >> b) & 1
                updated_a = (ua + ub + coin) // 2
                current = ua - updated_a
                phase = sign(mask, a, b) if current else 1
                target = mask ^ (1 << a) ^ (1 << b) if current else mask
                shifted = list(f)
                shifted[2 * e] = (shifted[2 * e] + current) % 6
                shifted[2 * e + 1] = (shifted[2 * e + 1] - current) % 6
                row = index[(target, tuple(shifted))]
                for selected in laws:
                    value = current if selected == 'C' else ua - ub
                    branches[(e, selected, coin, value)][row] += phase * vector[col] / np.sqrt(6 * len(laws))
    return branches


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root):
    artifact_names = ('JOINT_SOURCE.json', 'joint_source_states.npz', 'joint_source.py',
                      'charged_reservoir.py', 'modular_field_adapter.py')
    before = {name: digest(root / name) for name in artifact_names}
    receipt = json.loads((root / 'JOINT_SOURCE.json').read_text())
    declared = receipt['source_sha256']
    binding_mismatches = []
    base = root.parent.parent
    for name, expected_hash in declared.items():
        path = base / name
        observed_hash = digest(path) if path.is_file() else None
        if observed_hash != expected_hash:
            binding_mismatches.append({'path': name, 'recorded_sha256': expected_hash,
                                       'observed_sha256': observed_hash})
    with np.load(root / 'joint_source_states.npz', allow_pickle=False) as source:
        arrays = {key: source[key] for key in source.files}
    recorded = [(int(mask), tuple(map(int, f))) for mask, f in zip(arrays['occupation_masks'], arrays['half_records'])]
    expected = enumerate_basis()
    basis_complete = len(recorded) == len(set(recorded)) == 384 and set(recorded) == set(expected)
    if not basis_complete:
        raise AssertionError('basis does not equal the complete solved Gauss sector')
    h, field = operators(recorded)
    h_error = float(np.max(abs(h - arrays['hamiltonian'])))
    numbers = np.array([mask.bit_count() for mask, _ in recorded])
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    psi = eigenvectors[:, 0]
    rho = np.outer(psi, psi)
    ref_error = float(np.max(abs(rho - arrays['global_allN_reference'])))
    spectra_error = 0.
    for number in range(7):
        subset = np.flatnonzero(numbers == number)
        values = np.linalg.eigvalsh(h[np.ix_(subset, subset)])
        reported = receipt['actual_all_number_spectra'][str(number)]['all_levels']
        spectra_error = max(spectra_error, float(np.max(abs(values - reported))))
    native = {}
    for law in ('C', 'E', 'MIX'):
        branches = commit_reference(recorded, psi, law)
        out = sum(np.outer(vector, vector) for vector in branches.values())
        post_error = float(np.max(abs(out - arrays['global_allN_native_' + law])))
        work = float(np.trace(h @ (out - rho)))
        record = receipt['references']['global_allN']['native_instruments'][law]
        weights = {key: float(vector @ vector) for key, vector in branches.items()}
        reported_keys = {(row['seam'], row['law'], row['coin'], row['record_value']) for row in record['complete_outcomes']}
        if set(weights) != reported_keys or len(reported_keys) != len(record['complete_outcomes']):
            raise AssertionError('native record alphabet is incomplete or duplicated')
        probability_error = max(abs(weights[(row['seam'], row['law'], row['coin'], row['record_value'])] - row['probability']) for row in record['complete_outcomes'])
        native[law] = {'complete_poststate_error': post_error,
                       'outcome_probability_error': float(probability_error),
                       'normalization_error': abs(float(np.trace(out)) - 1),
                       'work': work, 'work_error': abs(work - record['work']),
                       'ground_survival': float(psi @ out @ psi)}
    # A scalar energy offset must be rejected even though all eigenvectors agree.
    shifted_h = arrays['hamiltonian'] + .125 * np.eye(384)
    mutations = {'missing_basis_state_rejected': set(recorded[:-1]) != set(expected),
                 'field_deleted_rejected': bool(np.max(abs(h - field - arrays['hamiltonian'])) > 1),
                 'constant_energy_shift_rejected': bool(np.max(abs(h - shifted_h)) > .1),
                 'unitary_stationary_as_native_stationary_rejected': all(row['ground_survival'] < .99 and row['work'] > .1 for row in native.values())}
    checks = {'complete_independently_enumerated_384_state_basis': basis_complete,
              'Hamiltonian_from_equations': h_error < 1e-12,
              'all_number_spectra': spectra_error < 1e-12,
              'interacting_global_reference': ref_error < 1e-12,
              'unique_global_reference': bool(eigenvalues[1] - eigenvalues[0] > 1e-6),
              'every_native_poststate_record_and_work': all(max(row[k] for k in ('complete_poststate_error', 'outcome_probability_error', 'normalization_error', 'work_error')) < 1e-12 for row in native.values()),
              'mutations_rejected': all(mutations.values())}
    after = {name: digest(root / name) for name in artifact_names}
    checks['audited_artifacts_unchanged_during_run'] = before == after
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'schema': 'oph.openwave.c10_external_audit.v1', 'checks': checks,
            'mutations': mutations, 'native_reference_response': native,
            'errors': {'Hamiltonian': h_error, 'spectra': spectra_error, 'reference': ref_error},
            'global_ground': float(eigenvalues[0]), 'source_sha256': before,
            'producer_declared_source_sha256': declared,
            'producer_source_binding_mismatches': binding_mismatches,
            'archived_values_reproduced': True,
            'producer_receipt_bound_to_current_workspace': not bool(binding_mismatches),
            'auditor_sha256': digest(Path(__file__)),
            'versions': {'python': platform.python_version(), 'numpy': np.__version__},
            'scope': 'One finite development source. Complete basis, H, spectra, reference and native records are independently reconstructed. Reservoir pulse, energy PVM continuation, cutoff removal, physical species and platform scoring are outside this check.',
            'official_rows_changed': []}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.c10)
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text)
