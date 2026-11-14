"""Audit the producer's actual compact Luders states and continuing native map.

Projectors are reconstructed by independent number-block eigendecomposition.
All serialized coefficients are compared with the full source input state.
Native continuation is built from index maps, without producer imports.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from joint_source_audit import EDGES, digest, operators, sign


def native_density(states, rho):
    index = {state: i for i, state in enumerate(states)}
    groups = defaultdict(list)
    for e, (a, b) in enumerate(EDGES):
        for coin in (0, 1):
            for col, (mask, records) in enumerate(states):
                ua, ub = (mask >> a) & 1, (mask >> b) & 1
                current = ua - (ua + ub + coin) // 2
                target = mask ^ (1 << a) ^ (1 << b) if current else mask
                f = list(records)
                f[2 * e] = (f[2 * e] + current) % 6
                f[2 * e + 1] = (f[2 * e + 1] - current) % 6
                row = index[(target, tuple(f))]
                phase = sign(mask, a, b) if current else 1
                for law, value in (('C', current), ('E', ua - ub)):
                    groups[(e, law, coin, value)].append((col, row, phase))
    out = np.zeros_like(rho, dtype=complex)
    probabilities = {}
    for label, entries in groups.items():
        src, dest, phases = np.array(entries).T
        if len(set(dest)) != len(dest):
            raise AssertionError('outcome index map is not injective')
        branch = rho[np.ix_(src, src)] * phases[:, None] * phases[None, :] / 12
        out[np.ix_(dest, dest)] += branch
        probabilities[label] = float(np.trace(branch).real)
    return out, probabilities


def load_arrays(path):
    with np.load(path, allow_pickle=False) as source:
        return {key: source[key] for key in source.files}


def run(root):
    names = ('CANONICAL_CONSTRUCTION.json', 'canonical_states.npz', 'joint_source_states.npz', 'spectral_readout.py')
    before = {name: digest(root / name) for name in names}
    result = json.loads((root / names[0]).read_text())
    spectral = result['spectral_readout_instrument']
    saved = load_arrays(root / 'canonical_states.npz')
    source = load_arrays(root / 'joint_source_states.npz')
    states = [(int(mask), tuple(map(int, records))) for mask, records in zip(source['occupation_masks'], source['half_records'])]
    h, _ = operators(states)
    if not np.array_equal(h, source['hamiltonian']):
        raise AssertionError('source H differs from independent equations')
    numbers = np.array([mask.bit_count() for mask, _ in states])
    u = expm(-.37j * h)
    rho = source['global_allN_reference']
    coherent = u @ rho @ u.conj().T
    committed, _ = native_density(states, coherent)
    selected = .71 * coherent + .29 * committed
    selected_error = float(np.max(abs(selected - saved['selected_protocol_poststate'])))
    inputs = {key: value for key, value in source.items() if value.shape == h.shape and key != 'hamiltonian'}
    inputs['selected_protocol_poststate'] = selected
    independent = {}
    for number in range(7):
        ids = np.flatnonzero(numbers == number)
        values, vectors = np.linalg.eigh(h[np.ix_(ids, ids)])
        groups = []
        for i, value in enumerate(values):
            if groups and abs(value - groups[-1][0]) < 1e-9:
                groups[-1][1].append(i)
            else:
                groups.append([value, [i]])
        independent[number] = (ids, values, vectors, groups)
    errors = {'projector': 0., 'coefficient': 0., 'probability': 0., 'normalization': 0., 'eigen_residual': 0.}
    checked = 0
    for label, rows in spectral['complete_branch_records'].items():
        if label not in inputs:
            continue  # The declared doublet control is checked separately below.
        seen = set()
        for row in rows:
            number = row['source_number']
            ids, values, vectors, groups = independent[number]
            match = [j for j, (energy, _) in enumerate(groups) if abs(energy - row['energy']) < 1e-8]
            if len(match) != 1 or (number, match[0]) in seen:
                raise AssertionError('missing, ambiguous or repeated spectral outcome')
            seen.add((number, match[0]))
            expected_columns = groups[match[0]][1]
            archived_ids = saved[row['position_indices_archive_key']]
            basis = saved[row['basis_archive_key']]
            if not np.array_equal(ids, archived_ids) or len(expected_columns) != row['degeneracy']:
                raise AssertionError('incorrect number-sector basis or degeneracy')
            expected = vectors[:, expected_columns]
            errors['projector'] = max(errors['projector'], float(np.max(abs(basis @ basis.conj().T - expected @ expected.conj().T))))
            coefficient = saved[row['coefficient_archive_key']]
            actual = basis.conj().T @ inputs[label][np.ix_(ids, ids)] @ basis
            errors['coefficient'] = max(errors['coefficient'], float(np.max(abs(actual - coefficient))))
            errors['probability'] = max(errors['probability'], abs(float(np.trace(actual).real) - row['probability']))
            errors['eigen_residual'] = max(errors['eigen_residual'], float(np.max(abs(h[np.ix_(ids, ids)] @ basis - row['energy'] * basis))))
            checked += 1
        if len(seen) != sum(len(data[3]) for data in independent.values()):
            raise AssertionError('complete zero outcomes missing')
        errors['normalization'] = max(errors['normalization'], abs(sum(row['probability'] for row in rows) - float(np.trace(inputs[label]).real)))
    def reconstruct(row):
        ids = saved[row['position_indices_archive_key']]
        basis = saved[row['basis_archive_key']]
        coefficient = saved[row['coefficient_archive_key']]
        state = np.zeros_like(h, dtype=complex)
        state[np.ix_(ids, ids)] = basis @ coefficient @ basis.conj().T
        return state
    actual = spectral['actual_populated_degenerate_continuation']
    branch = reconstruct(actual['branch'])
    evolved = u @ branch @ u.conj().T
    native, probabilities = native_density(states, evolved)
    post = .71 * evolved + .29 * native
    continuing_error = float(np.max(abs(post - saved[actual['poststate_archive_key']])))
    continuing_probability_error = 0.
    for row in actual['complete_joint_native_outcomes']:
        if row.get('status') == 'no_commit':
            expected = .71 * float(np.trace(branch).real)
        else:
            expected = .29 * probabilities[(row['seam'], row['law'], row['coin'], row['record_value'])]
        continuing_probability_error = max(continuing_probability_error, abs(expected - row['joint_probability']))
    control = next(row for row in spectral['complete_branch_records']['declared_doublet_control'] if row['probability'] > .99)
    luders = reconstruct(control)
    ids = saved[control['position_indices_archive_key']]
    basis = saved[control['basis_archive_key']]
    fine = np.zeros_like(h, dtype=complex)
    fine[np.ix_(ids, ids)] = basis @ np.diag([.5, .5]) @ basis.conj().T
    distance = float(np.linalg.norm(luders - fine, ord='nuc') / 2)
    controls_error = 0.
    for state, key in ((luders, 'Luders_continuation_archive_key'), (fine, 'fine_forget_continuation_archive_key')):
        evolved = u @ state @ u.conj().T
        native, _ = native_density(states, evolved)
        controls_error = max(controls_error, float(np.max(abs(.71 * evolved + .29 * native - saved[spectral['fine_basis_rejecting_control'][key]]))))
    checks = {'every_archived_input_and_zero_outcome_checked': checked == 5916,
              'projectors_coefficients_and_probabilities': max(errors.values()) < 1e-12,
              'selected_protocol_from_independent_native_map': selected_error < 1e-12,
              'actual_degenerate_branch_continuation': continuing_error < 1e-12 and continuing_probability_error < 1e-12,
              'fine_forgetting_mutation_distinguished': abs(distance - .5) < 1e-12,
              'both_producer_control_poststates': controls_error < 1e-12,
              'artifact_bytes_unchanged': before == {name: digest(root / name) for name in names}}
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'schema': 'oph.openwave.producer_spectral_archive_audit.v1', 'checks': checks,
            'producer_compared': True, 'checked_compact_branches': checked, 'errors': errors,
            'selected_protocol_error': selected_error, 'continuing_poststate_error': continuing_error,
            'continuing_record_probability_error': continuing_probability_error,
            'fine_forget_trace_distance': distance, 'control_poststates_error': controls_error,
            'source_sha256': before, 'auditor_sha256': digest(Path(__file__)),
            'auditor_dependency_sha256': {'joint_source_audit.py': digest(Path(__file__).with_name('joint_source_audit.py'))},
            'scope': 'Actual serialized producer projectors, full complex coefficients, 17 source inputs, zero outcomes and one subsequent native step. No synthesis from native C/E, physical energy apparatus or criterion award follows.',
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
