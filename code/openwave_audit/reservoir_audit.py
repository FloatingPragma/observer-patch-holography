"""Audit the actual C10 reservoir apparatus using a direct signed permutation.

The target module is imported only in a subprocess that exports its operators.
The reference algorithm solves the Gauss basis directly and never calls the
producer's CAR, swap, sector enumeration, reference or readout functions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from joint_source_audit import digest, enumerate_basis, operators


EXPORT = r'''
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, sys.argv[1])
from charged_reservoir import BoundarySource
s = BoundarySource(6)
a = {'masks': np.array([m for m, f in s.states]),
     'records': np.array([f for m, f in s.states])}
for meta, k in s.reservoir_instrument():
    a['k_%d_%d' % (meta['reservoir_initial'], meta['reservoir_final'])] = k
np.savez_compressed(sys.argv[2], **a)
'''


def signed_exchange(states):
    index = {state: j for j, state in enumerate(states)}
    n = len(states)
    blocks = {(i, o): np.zeros((n, n)) for i in (0, 1) for o in (0, 1)}
    boundary_errors = []
    for col, (mask, records) in enumerate(states):
        occupied = mask & 1
        for initial in (0, 1):
            final = occupied
            target = mask if occupied == initial else mask ^ 1
            # Seven ordered CAR modes: source port0, five spectators, reservoir.
            phase = (-1 if initial else 1) if occupied == initial else (-1) ** ((mask >> 1).bit_count())
            row = index[(target, records)]
            blocks[initial, final][row, col] = phase / np.sqrt(2)
            bi, bf = 3 - mask.bit_count(), 3 - target.bit_count()
            boundary_errors.append((initial - bi) - (final - bf))
    return blocks, max(map(abs, boundary_errors))


def energy_groups(h, numbers):
    groups = []
    references = {}
    for number in range(7):
        ids = np.flatnonzero(numbers == number)
        values, vectors = np.linalg.eigh(h[np.ix_(ids, ids)])
        selected = abs(values - values[0]) < 1e-10
        rho = np.zeros_like(h)
        rho[np.ix_(ids, ids)] = vectors[:, selected] @ vectors[:, selected].T / sum(selected)
        references[number] = (float(values[0]), rho)
        start = 0
        while start < len(values):
            end = start + 1
            while end < len(values) and abs(values[end] - values[start]) < 1e-9:
                end += 1
            groups.append((number, float(values[start]), end - start, ids, vectors[:, start:end]))
            start = end
    return groups, references


def audit(c10):
    bundle = c10.parents[1]
    # Pin all code in the target's local dependency directories before import.
    files = set(c10.glob('*.py')) | set((bundle / 'next/C9').glob('*.py'))
    files |= set((bundle / 'candidate/oph_openwave').glob('*.py'))
    files |= set((bundle / 'validation/bridges/providers').glob('*.py'))
    files |= {c10 / 'JOINT_SOURCE.json', c10 / 'joint_source_states.npz'}
    before = {os.path.relpath(p, c10): digest(p) for p in sorted(files)}
    receipt = json.loads((c10 / 'JOINT_SOURCE.json').read_text())
    with np.load(c10 / 'joint_source_states.npz', allow_pickle=False) as source:
        archive = {key: source[key] for key in source.files}
    with tempfile.TemporaryDirectory(prefix='oph-reservoir-audit-') as temp:
        export = Path(temp) / 'operators.npz'
        subprocess.run([sys.executable, '-c', EXPORT, str(c10), str(export)], check=True,
                       env={**os.environ, 'OPENBLAS_NUM_THREADS': '1'})
        with np.load(export, allow_pickle=False) as source:
            actual = {key: source[key] for key in source.files}
    states = [(int(m), tuple(map(int, f))) for m, f in zip(actual['masks'], actual['records'])]
    assert len(states) == len(set(states)) == 384 and set(states) == set(enumerate_basis())
    archive_states = [(int(m), tuple(map(int, f))) for m, f in zip(archive['occupation_masks'], archive['half_records'])]
    assert states == archive_states
    blocks, boundary_error = signed_exchange(states)
    h, field = operators(states)
    numbers = np.array([m.bit_count() for m, f in states])
    groups, refs = energy_groups(h, numbers)
    ground = min(e for e, r in refs.values())
    # Reconstruct the full global ground projector, including possible sectors.
    candidates = [rho for e, rho in refs.values() if abs(e - ground) < 1e-10]
    global_ref = sum(candidates) / len(candidates)
    # This recorded source has a unique ground; refuse an unweighted degeneracy transfer.
    assert len(candidates) == 1 and abs(np.trace(global_ref @ global_ref) - 1) < 1e-10
    references = {'global_allN': (ground, global_ref), 'fixed_N3': refs[3]}
    errors = {'operator': max(float(np.max(abs(actual[f'k_{i}_{o}'] - k))) for (i, o), k in blocks.items()),
              'completeness': float(np.max(abs(sum(k.T @ k for k in blocks.values()) - np.eye(len(states))))),
              'boundary_reservoir_charge': boundary_error,
              'field_commutator': max(float(np.max(abs(field @ k - k @ field))) for k in blocks.values()),
              'hamiltonian': float(np.max(abs(h - archive['hamiltonian']))),
              'reference': 0., 'poststate': 0., 'numbers_energy_records': 0., 'ledger': 0.}
    summaries = {}
    outcome_count = 0
    for label, (energy, rho) in references.items():
        report = receipt['references'][label]
        errors['reference'] = max(errors['reference'], float(np.max(abs(rho - archive[report['reference_archive_key']]))))
        reports = {(r['reservoir_initial'], r['reservoir_final']): r for r in report['complete_reservoir_response']}
        assert set(reports) == set(blocks) and len(report['complete_reservoir_response']) == 4
        work_sum = 0.
        for key, k in blocks.items():
            i, o = key
            recorded = reports[key]
            post = k @ rho @ k.T
            probability = float(np.trace(post))
            final_energy = float(np.trace(h @ post))
            work = final_energy - energy * probability + .23 * (o - i) * probability
            work_sum += work
            values = {'probability': probability, 'reservoir_preparation_probability': .5,
                      'initial_source_energy': energy, 'unnormalized_final_source_energy': final_energy,
                      'unnormalized_field_energy': float(np.trace(field @ post)),
                      'unnormalized_Hm_energy': float(np.trace((h - field) @ post)),
                      'unnormalized_apparatus_work': work, 'source_charge_change': i - o,
                      'boundary_record_change': o - i}
            errors['ledger'] = max(errors['ledger'], *(abs(values[name] - recorded[name]) for name in values))
            errors['poststate'] = max(errors['poststate'], float(np.max(abs(post - archive[recorded['poststate_archive_key']]))))
            records = recorded['complete_source_number_energy_records']
            assert len(records) == len(groups)
            for r, (q, e, multiplicity, ids, vectors) in zip(records, groups):
                assert r['source_number'] == q and r['degeneracy'] == multiplicity
                p = float(np.trace(vectors.T @ post[np.ix_(ids, ids)] @ vectors))
                errors['numbers_energy_records'] = max(errors['numbers_energy_records'], abs(p - r['probability']), abs(e - r['energy']))
            outcome_count += len(groups)
        summaries[label] = {'mean_apparatus_work': work_sum, 'recorded_work_error': abs(work_sum - report['mean_reservoir_apparatus_work']),
                            'probability_total': sum(float(np.trace(k @ rho @ k.T)) for k in blocks.values())}
    # Ground-state branches have fixed number. They do not detect spectator parity.
    fixed_number_parity_error = max(float(np.max(abs(k @ global_ref @ k.T - abs(k) @ global_ref @ abs(k).T))) for k in blocks.values())
    # One source state with N=0 and one with N=1, both port0 empty, witnesses the phase.
    probe = np.zeros(len(states))
    probe[next(j for j, (m, f) in enumerate(states) if m == 0)] = 1 / np.sqrt(2)
    probe[next(j for j, (m, f) in enumerate(states) if m == 2)] = 1 / np.sqrt(2)
    good = blocks[1, 0] @ probe
    bad = abs(blocks[1, 0]) @ probe
    normalized_distance = float(np.linalg.norm(np.outer(good, good) / (good @ good) - np.outer(bad, bad) / (bad @ bad), ord='nuc') / 2)
    after = {name: digest(c10 / name) for name in before}
    checks = {'all_operator_and_archival_errors_below_1e_10': max(errors.values()) < 1e-10,
              'all_work_totals_match': all(r['recorded_work_error'] < 1e-10 for r in summaries.values()),
              'all_four_outcomes_complete': all(abs(r['probability_total'] - 1) < 1e-10 for r in summaries.values()),
              'cross_number_spectator_mutation_rejected': normalized_distance > .99,
              'omitted_preparation_weight_rejected': np.max(abs(sum(2 * k.T @ k for k in blocks.values()) - np.eye(len(states)))) > .99,
              'unchanged_inspected_source': before == after}
    return {'schema': 'oph.audit.reservoir.v1', 'checks': {k: bool(v) for k, v in checks.items()},
            'errors': errors, 'source_dimension': 384, 'joint_dimension': 768,
            'archived_reservoir_branches': 8, 'complete_number_energy_outcomes': outcome_count,
            'references': summaries, 'fixed_number_parity_mutation_poststate_error': fixed_number_parity_error,
            'cross_number_mutation_conditional_trace_distance': normalized_distance,
            'source_sha256': before, 'reviewer_sha256': digest(Path(__file__)),
            'operator_bytes_sha256': {f'{i}_{o}': hashlib.sha256(actual[f'k_{i}_{o}'].tobytes()).hexdigest() for i, o in blocks},
            'scope': 'Actual producer operators extracted in a separate process; independent signed-permutation reference, complete basis, archive and work checks. Fair preparation, supplied pulse, zero extra boundary energy and mode ordering are declared premises. Cross-number control is a finite mathematical probe, not a derived physical preparation. No pulse duration, full apparatus cost, physical species or platform award follows.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.c10.resolve())
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'errors', 'references', 'complete_number_energy_outcomes')}, indent=2))
    raise SystemExit(0 if all(result['checks'].values()) else 1)
