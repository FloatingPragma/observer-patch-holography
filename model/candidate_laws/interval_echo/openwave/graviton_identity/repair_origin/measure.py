"""Identify the finite quadrupolar representation and measure its native flight."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[5]))
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation, rotations
from model.executable_model.carrier import icosahedron_from_coordinates, isomorphism
from model.experiments.interval_echo_space.shell_source import prepare

SHOTS = 16384
TIMES = [0, 1, 2, 4, 8]


def angular_frame():
    car = cell()
    coords, seams = icosahedron_from_coordinates()
    coords = coords[isomorphism(car.seams, seams)]
    coords /= np.linalg.norm(coords, axis=1)[:, None]
    x, y, z = coords.T
    waves = np.column_stack([x*x - y*y, 2*x*y, 2*x*z, 2*y*z, 3*z*z - 1])
    waves /= np.linalg.norm(waves, axis=0)
    np.testing.assert_allclose(waves.T @ waves, np.eye(5), atol=1e-13)
    # The six antipodal-pair rows form a regular simplex. At most three
    # distinct gateway rows in its restriction are linearly independent.
    gram = waves @ waves.T
    identical = np.isclose(np.abs(coords @ coords.T), 1)
    np.testing.assert_allclose(gram, np.where(identical, 5/12, -1/12), atol=1e-13)
    traces = {}
    for g in rotations():
        rot = np.linalg.lstsq(coords, coords[list(g)], rcond=None)[0]
        rep = waves.T @ waves[list(g)]
        np.testing.assert_allclose(waves[list(g)], waves @ rep, atol=1e-13)
        # Symmetric square character minus the rotation-invariant trace.
        character = (np.trace(rot)**2 + np.trace(rot @ rot))/2 - 1
        np.testing.assert_allclose(np.trace(rep), character, atol=1e-13)
        tag = str(int(round(character)))
        traces[tag] = traces.get(tag, 0) + 1
    assert sum(int(k)**2 * n for k, n in traces.items()) == 60
    return car, waves, identical, traces


def source(car, target):
    return prepare(SimpleNamespace(car=car, active=car.seams), target)


def ramsey(car, waves, rng):
    uniform = np.ones(12) / np.sqrt(12)
    h = candidate.generator(car).toarray()
    energy = float(uniform @ h @ uniform)
    gaps = np.diag(waves.T @ h @ waves) - energy
    np.testing.assert_allclose(gaps, [3.]*5, atol=3e-14)
    np.testing.assert_allclose(h @ waves, waves * (energy + gaps), atol=3e-14)
    times = np.arange(17)/4
    output = {}
    for control, matrix in [('registered', h), ('hops_off', np.diag(h.diagonal()))]:
        records = []
        for j, q in enumerate(waves.T):
            probes = [source(car, (uniform+phase*q)/np.sqrt(2))[0] for phase in [1, -1, 1j, -1j]]
            all_counts, phasor = [], []
            for t in times:
                evolved = expm(-1j*t*matrix) @ np.stack(probes, axis=1)
                probabilities = np.abs(evolved.T)**2
                counts = np.array([rng.multinomial(SHOTS, p/p.sum()) for p in probabilities])
                all_counts.append(counts.tolist())
                contrasts = counts[::2] - counts[1::2]
                real, imag = np.sqrt(12)/2 * (contrasts @ q)/SHOTS
                phasor.append(real + 1j*imag)
            fit = float(np.polyfit(times, np.unwrap(np.angle(phasor)), 1)[0])
            records.append({'mode': j, 'frequency_from_complete_counts': fit, 'counts': all_counts})
        output[control] = records
    return {'exposures': times.tolist(), 'phase_order': ['0', 'pi', 'pi/2', '-pi/2'],
            'shots_per_setting': SHOTS, 'uniform_energy': energy,
            'quadrupole_energy': float(energy+gaps[0]), 'gap': 3, 'experiments': output}


def boundary_probe(car, waves, probes, receiver, rng):
    h = candidate.generator(car)
    # Subtract the common uniform-state phase before sparse propagation.
    h = h - 48*eye(car.n, format='csr')
    data = {}
    for name, projection in probes.items():
        prepared, gates = source(cell(), waves @ projection)
        state = np.zeros(car.n, complex)
        state[:12] = prepared
        flight = expm_multiply(-1j*h, state, start=0, stop=8, num=9)
        rows = []
        for t in TIMES:
            local = abs(flight[t, 12*receiver:12*receiver+12])**2
            arrival = float(local.sum())
            probability = np.append(local, max(0., 1 - arrival))
            np.testing.assert_allclose(probability.sum(), 1, atol=1e-13)
            counts = rng.multinomial(SHOTS, probability/probability.sum())
            rows.append({'time': t, 'arrival': arrival, 'counts_ports_then_quiet': counts.tolist()})
        data[name] = {'gates': gates, 'quadrupole_coefficients': projection.tolist(), 'readings': rows}
    assert max(row['arrival'] for row in data['dark']['readings']) < 1e-23
    assert max(row['arrival'] for row in data['bright']['readings']) > .001
    return data


def federation_probe(board, waves, axes, energy, rng):
    h = candidate.generator(board)
    census = {}
    for home in range(board.cells):
        block = np.zeros((board.n,5))
        block[12*home:12*home+12] = waves
        gateway_ports = np.where(board.degree[12*home:12*home+12] == 6)[0]
        assert len(gateway_ports) == 3
        distinct_axes = len({tuple(axes[p]) for p in gateway_ports})
        local = waves[gateway_ports]
        _, singular, right = np.linalg.svd(local, full_matrices=True)
        rank = int((singular > 1e-11).sum())
        assert rank == distinct_axes
        residue = h @ block - energy*block
        np.testing.assert_allclose(residue.T @ residue, .5*local.T @ local, atol=3e-13)
        np.testing.assert_allclose(residue @ right[rank:].T, 0, atol=3e-13)
        key = str(5-rank)
        census[key] = census.get(key, 0) + 1
        if home == 0:
            first = int(gateway_ports[0])
            receiver = next(q//12 for q in board.neighbours[first] if q//12 != home)
            bright = waves[first]/np.linalg.norm(waves[first])
            darkspace = right[rank:].T @ right[rank:]
            column = next(q for q in darkspace.T if np.linalg.norm(q) > .25)
            dark = column/np.linalg.norm(column)
            probe = boundary_probe(board, waves, {'dark': dark, 'bright': bright}, receiver, rng)
    return {'cells': board.cells, 'cells_by_dimension_of_trapped_quadrupoles': census,
            'first_source_cell': 0, 'first_contact': first, 'receiving_cell': receiver,
            'observer_experiments': probe}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    out = parser.parse_args().output
    car, waves, axes, traces = angular_frame()
    rng = np.random.default_rng(5364678)
    local = ramsey(car, waves, rng)
    boards = [federation_probe(federation(level), waves, axes, local['quadrupole_energy'], rng)
              for level in range(3)]
    with out.open('x') as f:
        json.dump({'proper_rotation_character_census': traces, 'local_ramsey': local,
                   'native_boards': boards, 'seed': 5364678}, f, indent=1, allow_nan=False)
        f.write('\n')


if __name__ == '__main__':
    main()
