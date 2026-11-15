#!/usr/bin/env python3
"""Reproducible local integration smoke checks; never consumes frozen acceptance cases."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition: raise ValueError(message)


def run():
    sys.path.insert(0, str(ROOT / 'candidate'))
    from oph_openwave.geometry import flat_ring
    mass = .173
    construction = flat_ring(2, 3, mass)
    b = construction.incidence
    n, e = b.shape
    d = np.block([[np.zeros((n, n)), b], [b.T, np.zeros((e, e))]])
    gamma = np.diag(np.r_[np.ones(n), -np.ones(e)])
    mass = 0.173
    h = construction.hamiltonian.toarray() if hasattr(construction.hamiltonian, "toarray") else np.asarray(construction.hamiltonian)
    require(np.max(np.abs(h - (d + mass * gamma))) < 1e-12, 'local gate failed')
    square_error = float(np.max(np.abs(h @ h - (d @ d + mass**2 * np.eye(n + e)))))
    require(square_error < 1e-12, 'local gate failed')
    # Independent composition uses a Kronecker sum, checked on a small spectral subspace.
    indices = [int(construction.seams[0][0]), int(construction.seams[0][1]), n]
    small = h[np.ix_(indices, indices)]
    other = small + np.diag([0.1, -0.2, 0.3])
    joint = np.kron(small, np.eye(3)) + np.kron(np.eye(3), other)
    additive_error = float(np.max(np.abs(expm(-0.37j * joint) - np.kron(expm(-0.37j * small), expm(-0.37j * other)))))
    require(additive_error < 1e-12, 'local gate failed')
    rng = np.random.default_rng(20261114)
    initial = rng.integers(-5, 6, n)
    final, records = tuple(map(int, initial)), (0,)*e
    for _ in range(10000):
        final, records, _ = construction.committed_event(final, records, int(rng.integers(e)), int(rng.integers(2)))
        require(construction.gauss_register(final, records) == tuple(initial), 'local gate failed')
    final, records = np.array(final), np.array(records)
    gauss_error = float(np.max(np.abs(final + b @ records - initial)))
    require(gauss_error == 0, 'local gate failed')
    require(final.sum() == initial.sum(), 'local gate failed')
    packet = json.loads((ROOT / 'platform_packet/packet.json').read_text())
    rows = packet['rows']
    require(sorted(r['id'] for r in rows) == list(range(1,32)), 'incomplete row dispatch')
    return {'scope': 'local smoke checks; not platform validation or fresh acceptance',
            'construction_sha256': hashlib.sha256((ROOT / 'candidate/oph_openwave/construction.py').read_bytes()).hexdigest(),
            'geometry_sha256': hashlib.sha256((ROOT / 'candidate/oph_openwave/geometry.py').read_bytes()).hexdigest(),
            'seed': 20261114, 'ring_shape': [2, 3], 'mass': mass,
            'checks': {'mass_square_error': square_error, 'independent_composition_error': additive_error,
                       'native_record_gauss_error': gauss_error, 'official_untested_rows': len(rows)}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2) + '\n'
    if args.output:
        args.output.write_text(result)
    print(result, end='')
