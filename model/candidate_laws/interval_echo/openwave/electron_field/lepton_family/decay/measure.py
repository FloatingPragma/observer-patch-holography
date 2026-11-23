"""Resolve daughter and neutral energy on the native written-field two-cell bridge."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import networkx as nx
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csgraph

HERE = Path(__file__).resolve().parent
TOL = 1e-8
TIMES = (np.pi/8, np.pi/4, np.pi/2)


def norm2(x):
    return float(np.vdot(x, x).real)


def face_pairs(patch, words):
    word_ids = {w: i for i, w in enumerate(words)}
    flips = [[] for _ in patch.faces]
    for i, w in enumerate(words):
        for face, other in patch.ring_moves(w):
            j = word_ids[other]
            if i < j:
                flips[face].append((i, j))
    return [np.array(f, dtype=int).reshape(-1, 2).T for f in flips]


def ring_projection(wave, face):
    projected = np.zeros_like(wave)
    i, j = face
    projected[i] = (wave[i] - wave[j]) / 2
    projected[j] = -projected[i]
    return projected


def ring_read(wave, modes, energies, faces, near, far):
    coeffs = modes.T @ wave
    expected = float(np.dot(energies, abs(coeffs)**2))
    rows = []
    for wait in (0., np.pi/4, np.pi/2):
        psi = modes @ (coeffs * np.exp(-1j * energies * wait))
        local = [ring_projection(psi, face) for face in faces]
        clicks = [norm2(x) for x in local]
        field_energy = sum(clicks)
        assert abs(field_energy - expected) < 1e-10
        both = ring_projection(local[near], faces[far])
        assert np.max(abs(both - ring_projection(local[far], faces[near]))) < 1e-12
        rows.append({'wait': wait, 'near_click': clicks[near], 'far_click': clicks[far],
                     'both_click': norm2(both), 'sum_native_face_energies': field_energy, 'all_face_clicks': clicks})
    return rows


def run_sector(sign, patch, gateway, receiving, perm, charged_sector, lowest, free_flight, action):
    h, states, ports = charged_sector(patch, sign)
    levels, vectors, charged_groups, residual = lowest(h)
    assert all(len(group) == 1 for group in charged_groups[:3])
    energies = levels[:3]
    vec = np.array(vectors[:, :3], order='C')
    idx = np.flatnonzero(ports == gateway)
    adj = -h[idx][:, idx]
    adj.setdiag(0)
    adj.eliminate_zeros()
    assert np.all(adj.data == 1)
    nh = csgraph.laplacian(adj, normed=False).tocsr()
    field_energies, modes = eigh(nh.toarray()/2)
    neutral = np.full(len(idx), 1/np.sqrt(len(idx)))
    np.testing.assert_allclose(abs(modes[:, 0] @ neutral), 1., atol=1e-12)
    assert abs(field_energies[0]) < TOL and field_energies[1] > TOL
    field_groups = np.split(np.arange(len(idx)), np.flatnonzero(np.diff(field_energies) > TOL) + 1)
    grouped_energies = np.array([field_energies[g].mean() for g in field_groups])
    grouped_energies[0] = 0.
    face_flips = face_pairs(patch, [states[i][1] for i in idx])
    near = next(f for f, face in enumerate(patch.faces) if gateway in face)
    distance = nx.single_source_shortest_path_length(patch.dual, near)
    far = min(f for f, d in distance.items() if d == max(distance.values()))
    assert set(patch.faces[near]).isdisjoint(patch.faces[far])
    outputs = []
    ground_transfer = np.zeros((len(TIMES), 3, 3), dtype=complex)
    overlap = neutral @ vec[idx]
    for i, parent in enumerate(vec.T):
        source = np.outer(parent, neutral)
        work = 4 * norm2(parent[idx])
        initial_energy = energies[i] + work
        variance = 8 * work - work**2
        src_h, dest_h = action(h, nh, idx, source, np.zeros_like(source))
        np.testing.assert_allclose(norm2(src_h-initial_energy*source) + norm2(dest_h), variance, atol=1e-10)
        bare = -4 * overlap[:, None] * (modes.T @ parent[idx])[None, :]
        np.testing.assert_allclose(vec.T @ dest_h @ modes, bare, atol=1e-11)
        del src_h, dest_h
        shells = []
        for j in range(i):
            release = float(energies[i] - energies[j])
            k = int(np.argmin(abs(grouped_energies[1:] - release))) + 1
            equal = [g for n, g in enumerate(field_groups) if n and abs(release-grouped_energies[n]) <= TOL]
            shells.append({'daughter': j, 'isolated_energy_release': release, 'nearest_positive_field_energy': grouped_energies[k],
                           'signed_release_minus_neutral_level': release - grouped_energies[k],
                           'bare_squared_bridge_matrix_element_nearest_group': norm2(bare[j, field_groups[k]]),
                           'on_shell_positive_field_eigenstates': sum(len(g) for g in equal),
                           'summed_bare_squared_bridge_matrix_element_on_shell': sum(norm2(bare[j, g]) for g in equal)})
        evolved, precision = free_flight(h, nh, idx, source, TIMES)
        rows = []
        for ti, (left, right) in enumerate(evolved):
            arrival = norm2(right)
            norm_error = norm2(left) + arrival - 1
            hl, hr = action(h, nh, idx, left, right)
            energy_error = float((np.vdot(left, hl) + np.vdot(right, hr)).real - initial_energy)
            del hl, hr
            assert abs(norm_error) < 2e-12 and abs(energy_error) < 2e-10
            heralds = vec.T @ right
            ground_transfer[ti, i] = heralds @ neutral
            channels = []
            for j, field in enumerate(heralds):
                total = norm2(field)
                coeffs = modes.T @ field
                weights = abs(coeffs)**2
                ground = float(weights[0])
                positive = float(sum(weights[1:]))
                assert abs(total - ground - positive) < 1e-12
                channel = {'daughter': j, 'charged_daughter_joint_probability': total, 'neutral_ground_joint_probability': ground,
                           'neutral_positive_joint_probability': positive, 'conditional_mean_neutral_energy': float(weights @ field_energies/total)}
                if ti == 1:
                    channel['neutral_energy_group_joint_probabilities'] = [float(weights[g].sum()) for g in field_groups]
                    channel['isolated_field_face_reads_joint_with_charged_daughter'] = ring_read(field, modes, field_energies, face_flips, near, far)
                channels.append(channel)
            outside = arrival - sum(c['charged_daughter_joint_probability'] for c in channels)
            assert outside > -1e-12
            rows.append({'time': TIMES[ti], 'charge_arrival_probability': arrival, 'arrival_outside_first_three_charged_modes': outside,
                         'norm_error': norm_error, 'energy_error': energy_error, 'channels': channels})
        outputs.append({'parent': i, 'bridge_opening_work': work, 'conserved_mean_open_energy': initial_energy,
                        'open_energy_standard_deviation': float(np.sqrt(variance)), 'separated_cell_shell_checks': shells, 'times': rows})
        del evolved, left, right, source
    reciprocity_error = float(np.max(abs(ground_transfer - ground_transfer.transpose(0, 2, 1))))
    assert reciprocity_error < 1e-10
    # The existing two-load weak service conditions on an absolute cell occupation of two.
    # Each charged basis state has one excess or vacancy and the other cell has unit background.
    background = len(patch.cell.neighbours)
    counts = [(background + sign, background), (background, background + sign)]
    eligible = len(states) * len(idx) * sum(2 in pair for pair in counts)
    assert eligible == 0
    return {'sign': sign, 'gateways': [gateway, receiving], 'receiver_rotation': list(perm),
            'charged_dimension': len(states), 'neutral_dimension': len(idx), 'open_dimension': 2*len(states)*len(idx),
            'charged_energies_native': energies.tolist(), 'max_charged_eigenpair_residual': residual,
            'neutral_energy_groups_native': grouped_energies.tolist(), 'neutral_group_multiplicities': [len(g) for g in field_groups],
            'near_far_face_ids': [near, far], 'near_far_triangle_ports': [patch.faces[near], patch.faces[far]],
            'neutral_flip_pairs_per_face': [pairs.shape[1] for pairs in face_flips],
            'dual_distance_from_near_face': [distance[f] for f in range(len(patch.faces))],
            'dual_face_separation': distance[far], 'precision': precision, 'reciprocal_ground_transfer_max_amplitude_error': reciprocity_error,
            'existing_two_load_weak_service': {'absolute_cell_occupations_left_right': counts,
                                               'eligible_basis_states': eligible, 'checked_basis_dimension': 2*len(states)*len(idx)},
            'cut_bridge': {'derivation': 'two separated Hamiltonians conserve the cell of the charge and the initially neutral ground energy',
                           'receiver_charge_probability_all_t': 0, 'neutral_field_excitation_probability_all_t': 0}, 'parents': outputs}


def measure():
    archive = HERE.parents[1]/'records.zip'
    bundle_dir = 'experiments/interval_echo_electron_flux_binding'
    with TemporaryDirectory() as tmp, ZipFile(archive) as bundle:
        for filename in ('search.py', 'bridge.py', 'free_bridge.py'):
            bundle.extract(f'{bundle_dir}/{filename}', tmp)
        sys.path.insert(0, str(HERE.parents[6]))
        sys.path.insert(0, str(Path(tmp)/bundle_dir))
        from bridge import matched_bridge
        from free_bridge import free_flight, hamiltonian_action
        from search import Patch, charged_sector, lowest
        gateway, receiving, perm = matched_bridge()
        patch = Patch()
        return {'absolute_native_shell_tolerance': TOL, 'ideal_basis_preparations_and_reads': True,
                'field_ground_and_excited_weights_unconditional': True,
                'sectors': [run_sector(s, patch, gateway, receiving, perm, charged_sector, lowest, free_flight, hamiltonian_action) for s in (-1, 1)]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = measure()
    with args.output.open('x') as destination:
        json.dump(report, destination, indent=2, allow_nan=False)
        destination.write('\n')
