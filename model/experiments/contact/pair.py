"""A source-on emitter and two source-on local pointers in the mobile two-load sector."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.replay import cell as load_cell
from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.cage import construction as cage
from .__main__ import HERE, check_sources, load_json, protocol

SOURCE_FILES = {'model/candidate_laws/contrast_memory_spherical_v1/replay.py', 'model/experiments/contact/pair.py',
                'model/experiments/contact/pair_protocol.json', 'model/experiments/contact/sources.json'}


def step(H, state, time):
    return expm_multiply(-1j * time * H, state, traceA=-1j * time * H.diagonal().sum())


def mobile(spec, g):
    """All 78 load configurations plus the emitter fuel, each with all 150 scalar modes."""
    car = cage.cell()
    D = cage.law(car)
    H, states, adjacent = load_cell(spec['kappa'])
    H -= 12 * spec['kappa'] * np.eye(len(H))
    pair = states.index(tuple(spec['emitter']))
    origin = car['index'][tuple(spec['scalar_seed_and_reader'])]
    occupancy = np.array([np.bincount(s, minlength=12) for s in states])
    a, b = np.asarray(car['edges']).T
    square = (occupancy[:, a] - occupancy[:, b]) ** 2
    potential = np.zeros((len(H) + 1, len(D)))
    potential[:-1, :len(a)] = square
    extended = sp.block_diag((sp.csr_matrix(H), np.array([[H[pair, pair]]])), format='csr')
    free = sp.kronsum(sp.csr_matrix(D), extended, format='csr') + sp.diags(g * potential.ravel())
    drive = sp.lil_matrix(extended.shape)
    drive[pair, -1] = drive[-1, pair] = spec['drive_strength']
    driven = free + sp.kron(drive.tocsr(), sp.eye(len(D)), format='csr')
    start = np.zeros(potential.size)
    start[-len(D) + origin] = 1
    pulse = math.pi / (2 * spec['drive_strength'])
    emitted = step(driven, start, pulse)
    states_at_read = np.column_stack([step(free, emitted, t) for t in spec['waits']])
    # H remains on throughout the terminal pulse. In the X basis of each
    # pointer the two controlled flips give four source Hamiltonians.
    p = np.zeros_like(potential)
    p[pair] = 1
    n = np.zeros_like(potential)
    n[:, origin] = 1
    branches = {}
    for x in (1, -1):
        for y in (1, -1):
            reading = free + sp.diags(spec['drive_strength'] * (x * p + y * n).ravel())
            branches[x, y] = step(reading, states_at_read, pulse)
    probabilities = []
    for i in (0, 1):
        for j in (0, 1):
            arm = sum(x**i * y**j * v for (x, y), v in branches.items()) / 4
            probabilities.append(np.sum(abs(arm) ** 2, axis=0))
    law = np.array(probabilities).T.reshape(-1, 2, 2)
    weights = abs(states_at_read.reshape(len(H) + 1, len(D), -1)) ** 2
    no_pair = weights[np.r_[~adjacent, True]].sum(axis=(0, 1))
    doubled = weights[np.r_[[a == b for a, b in states], False]].sum(axis=(0, 1))
    return law, no_pair, doubled


def binding_bound(spec, g):
    gap = 2 * spec['kappa'] - 1 - 3 * math.sqrt(6) - abs(g)
    require(gap > 0, 'cell block gap is closed')
    return min(1., (4 / gap + math.hypot(math.sqrt(2), g) / spec['drive_strength']) ** 2)


def distance_from_product(laws):
    product = laws.sum(axis=3, keepdims=True) * laws.sum(axis=2, keepdims=True)
    return .5 * abs(laws - product).sum(axis=(2, 3))


def calculate(spec):
    outcomes = [mobile(spec, g) for g in spec['couplings']]
    laws = np.array([row[0] for row in outcomes])
    return dict(pointer_laws=laws.tolist(), no_pair_or_fuel=[x[1].tolist() for x in outcomes],
                doubled=[x[2].tolist() for x in outcomes],
                all_time_bounds=[binding_bound(spec, g) for g in spec['couplings']],
                total_variation_from_product=distance_from_product(laws).tolist())


def assess(spec, result):
    menu = (len(spec['couplings']), len(spec['waits']))
    shapes = dict(pointer_laws=(*menu, 2, 2), no_pair_or_fuel=menu, doubled=menu,
                  all_time_bounds=menu[:1], total_variation_from_product=menu)
    require(type(result) is dict and result.keys() == shapes.keys(), 'mobile receipt keys differ')
    arrays = {}
    for name, shape in shapes.items():
        raw = np.asarray(result[name], dtype=object)
        require(raw.shape == shape and all(type(x) in (int, float) and 0 <= x <= 1 for x in raw.flat),
                f'invalid probability or shape: {name}')
        arrays[name] = raw.astype(float)
    laws = arrays['pointer_laws']
    require(np.max(abs(laws.sum(axis=(2, 3)) - 1)) < 1e-10, 'the mobile pointer record is incomplete')
    dependence = distance_from_product(laws)
    require(np.allclose(dependence, arrays['total_variation_from_product'], rtol=0, atol=2e-12),
            'dependence differs from the joint pointer probabilities')
    require(np.max(dependence[0]) < 1e-10, 'additive sources gave dependent records')
    require(np.all(np.max(dependence[1:], axis=1) > spec['minimum_joint_record_dependence']),
            'a signed contact lost its joint pointer signal')
    bounds = [binding_bound(spec, g) for g in spec['couplings']]
    require(np.allclose(arrays['all_time_bounds'], bounds, rtol=0, atol=2e-14), 'wrong all-time bounds')
    require(np.all(arrays['doubled'] == 0) and np.all(np.max(arrays['no_pair_or_fuel'], axis=1) <= bounds),
            'mobile sector bound failed')
    return arrays


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--record', type=Path)
    args = cli.parse_args()
    protocol()
    check_sources(HERE / 'pair_sources.json', SOURCE_FILES)
    spec = load_json(HERE / 'pair_protocol.json')
    result = calculate(spec)
    if args.record:
        with args.record.open('x', encoding='utf-8') as f:
            json.dump(result, f, indent=1, allow_nan=False)
            f.write('\n')
    else:
        saved = assess(spec, load_json(HERE / 'pair_records.json'))
        for k, v in result.items():
            require(np.allclose(saved[k], v, atol=2e-11, rtol=0), f'mobile receipt differs: {k}')
    assess(spec, result)
    print(json.dumps({'all_time_bounds': result['all_time_bounds'],
                      'max_record_dependence': np.max(result['total_variation_from_product'], axis=1).tolist()}))


if __name__ == '__main__':
    main()
