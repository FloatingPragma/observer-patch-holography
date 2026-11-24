"""Pair-fixing reflections and the three-triplet singlet on the intact ground state."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier, rotations
from .routings import partitions


def reflection_witnesses(car, sector, anti):
    eye = np.arange(car.n)
    adj = car.adjacency().toarray().astype(int)
    mirrors = []
    for perm in carrier.automorphisms(car.seams):
        g = np.array(perm)
        if not np.array_equal(g[g], eye):
            continue
        # 4 tr(gP) = a + b/sqrt(5) for the calibrated rank-three projector.
        a = int(np.sum(g == eye) - np.sum(g == anti))
        b = int(sum(adj[v, g[v]] - adj[anti[v], g[v]] for v in eye))
        if (a, b) == (4, 0):
            mapping = [sector.index[tuple(sorted(g[list(t)]))] for t in sector.states]
            np.testing.assert_array_equal(sector.agreed[mapping], sector.agreed)
            np.testing.assert_array_equal(sector.K[np.ix_(mapping, mapping)], sector.K)
            mirrors.append(perm)
    choices = []
    for a, b in combinations(eye, 2):
        choices.append(next(i for i, g in enumerate(mirrors) if g[a] == a and g[b] == b))
    return mirrors, choices


def preparations(car, sector):
    states = {}
    for kappa in (0, candidate.KAPPA):
        eigenvalues, vectors = np.linalg.eigh(sector.hamiltonian(kappa))
        ground = vectors[:, 0]
        assert eigenvalues[1] > eigenvalues[0] + 1e-9
        states[f'ground_k{int(kappa)}'] = ground
        if kappa:
            for i, face in enumerate(car.faces()[:2]):
                column = sector.index[face]
                states[f'flight_{i}'] = vectors @ (np.exp(-1j * eigenvalues) * vectors[column, :])
    return states


def run():
    car = carrier.cell()
    sector = BinarySector(car, 3)
    components, _ = connected_components(sector.K < 0, directed=False)
    assert components == 1
    anti = np.array(currents._geometry(car)[0])
    mirrors, witnesses = reflection_witnesses(car, sector, anti)
    energies, vectors = np.linalg.eigh(car.adjacency().toarray())
    triplet = vectors[:, np.isclose(energies, np.sqrt(5))]
    assert triplet.shape == (12, 3)
    grid = np.indices((12, 12, 12)).transpose(1, 2, 3, 0)
    epsilon = np.linalg.det(triplet[grid]) / np.sqrt(6)
    # Rotate only the first destination port as a physical relative-frame control.
    rotation = min(tuple(g) for g in rotations() if tuple(g) != tuple(range(12)))
    frames = {'common': epsilon, 'rotate_first': epsilon[list(rotation)]}
    states = preparations(car, sector)
    names = list(states)
    inputs = np.column_stack(list(states.values()))
    triples = np.array(sector.states)
    peaks = {frame: np.zeros(len(names)) for frame in frames}
    count = 0
    for groups in partitions():
        owner = np.empty(12, int)
        for i, group in enumerate(groups):
            owner[list(group)] = i
        labels = owner[triples]
        ok = (labels[:, 0] != labels[:, 1]) & (labels[:, 1] != labels[:, 2]) & (labels[:, 0] != labels[:, 2])
        ordered = np.take_along_axis(triples[ok], np.argsort(labels[ok], axis=1), axis=1)
        for frame, eps in frames.items():
            amplitude = eps[tuple(ordered.T)] @ inputs[ok]
            peaks[frame] = np.maximum(peaks[frame], np.abs(amplitude))
        count += 1
    return {
        'sector_size': len(triples), 'hopping_components': components,
        'reflections_of_triplet_trace_one': mirrors,
        'witness_index_for_each_lexicographic_port_pair': witnesses,
        'rotation_of_first_receiver': rotation, 'partitions': count,
        'max_per_preparation_singlet_weight': {
            frame: dict(zip(names, (np.round(amp**2, 11)).tolist(), strict=True))
            for frame, amp in peaks.items()},
        'common_ground_amplitudes_below_1e_minus_10': {
            name: bool(peaks['common'][i] < 1e-10)
            for i, name in enumerate(names) if name.startswith('ground')},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
