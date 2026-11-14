"""Independent sector and measurement audit of the proposed local walk.

The matrix construction uses the graph and geometric matching definition.
It does not assume a compiled pulse implementation or a physical time unit.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'reverse-engineering-reality'))
sys.path.insert(0, str(HERE.parent / 'native_walk_bridge'))
from model.native_spatial_refinement import source as spatial
import common_walk as provider


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    paths = [Path(__file__), Path(provider.__file__), Path(spatial.__file__)]
    pins = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    g = spatial.graph(2)
    n = g.ports
    # Independently assign each internal edge to the declared matching and
    # each external edge to colour5; never use provider swap permutations.
    colours = {tuple(e): d for d, row in enumerate(provider.INTERNAL_MATCHINGS)
               for e in row}
    neighbour = np.full((n, 6), -1, dtype=int)
    for a, b in g.seams:
        d = colours[(a % 12, b % 12)] if a // 12 == b // 12 else 5
        assert neighbour[a, d] == neighbour[b, d] == -1
        neighbour[a, d], neighbour[b, d] = b, a
    assert (neighbour >= 0).all()
    swap = np.arange(n * 12)
    for v in range(n):
        for d in range(6):
            swap[12 * v + 2 * d + 1] = 12 * neighbour[v, d] + 2 * d + 1
    assert np.array_equal(swap[swap], np.arange(n * 12))
    A = np.kron(np.eye(n), np.ones((12, 1)) / np.sqrt(12))
    SA = A[swap]
    T = A.T @ SA
    laplacian = 6 * np.eye(n)
    for a, b in g.seams:
        laplacian[a, b] = laplacian[b, a] = -1
    compression_error = float(np.max(abs(T - (np.eye(n) - laplacian / 12))))
    assert compression_error < 1e-14
    assert np.max(abs(T - provider.native_transfer())) < 1e-14

    # The Gram matrix of [A,SA] exposes both the physical dispersive sector
    # and the omitted cycle sectors without diagonalizing the full walk.
    C = np.column_stack((A, SA))
    vals, vecs = np.linalg.eigh(C.T @ C)
    keep = vals > 1e-10
    Q = (C @ vecs[:, keep]) / np.sqrt(vals[keep])[None, :]
    assert Q.shape[1] == 2 * n - 1
    assert np.linalg.norm(Q.T @ Q - np.eye(2 * n - 1), 2) < 1e-12

    def walk(x):
        y = 2 * A @ (A.T @ x) - x
        return y[swap]

    wq = walk(Q)
    invariant_error = float(np.linalg.norm(wq - Q @ (Q.T @ wq), 2))
    assert invariant_error < 1e-12
    phase = np.linalg.eigvals(Q.T @ wq)
    lambdas = np.linalg.eigvalsh(T)
    target = np.r_[1 + 0j, np.exp(1j * np.arccos(lambdas[:-1])),
                   np.exp(-1j * np.arccos(lambdas[:-1]))]
    spectral_error = float(np.max(np.min(abs(phase[:, None] - target[None, :]), axis=1)))
    assert spectral_error < 1e-12

    # Generic vectors projected outside span(A,SA) obey W=-S. Splitting by
    # S parity verifies both flat sectors as full source states.
    rng = np.random.default_rng(1729)
    residual = rng.normal(size=(12 * n, 5)) + 1j * rng.normal(size=(12 * n, 5))
    residual -= Q @ (Q.T @ residual)
    plus = (residual - residual[swap]) / 2
    minus = (residual + residual[swap]) / 2
    flat_error = max(float(np.max(abs(walk(plus) - plus))),
                     float(np.max(abs(walk(minus) + minus))))
    assert flat_error < 1e-12
    multiplicities = {'physical_plus_one': 2 * n + 2,
                      'physical_minus_one': 8 * n,
                      'physical_dispersive_nonzero_phases': 2 * (n - 1),
                      'padded_plus_one': 2 * n + 2 + 32,
                      'padded_minus_one': 8 * n + 32 * 15 + n * 4}
    assert sum(multiplicities[k] for k in (
        'physical_plus_one', 'physical_minus_one', 'physical_dispersive_nonzero_phases')) == 12 * n
    assert multiplicities['padded_plus_one'] + multiplicities['padded_minus_one'] + 2 * (n - 1) == 2048

    # A complete position measurement can leave the dispersive sector.
    # Record every outcome, including zero weights, and use no conditioning.
    seed = np.zeros(n)
    seed[:12] = 1 / np.sqrt(12)
    state = A @ seed
    measurement = []
    for step in range(1, 5):
        state = walk(state)
        probabilities, flat_weights = [], []
        for v in range(n):
            branch = np.zeros(12 * n, dtype=complex)
            branch[12 * v:12 * v + 12] = state[12 * v:12 * v + 12]
            probabilities.append(float(np.vdot(branch, branch).real))
            flat_weights.append(float(np.vdot(branch, branch).real - np.vdot(Q.T @ branch, Q.T @ branch).real))
        assert abs(sum(probabilities) - 1) < 1e-12
        total = sum(flat_weights)
        assert total > 1e-3
        measurement.append({'step': step, 'outcome_probabilities': probabilities,
                            'unconditional_weight_outside_dispersive_sector': total,
                            'postselection_or_renormalization': False})

    # A source scalar sign is invisible in un-controlled position statistics
    # and observable in a probe that controls W.
    psi = A @ seed
    overlap = np.vdot(psi, walk(psi))
    sign_distance = float(abs(overlap))
    assert sign_distance > .8
    result = {
        'schema': 'oph.local12.kinematic_opposing.v1',
        'candidate_id': provider.MODEL_ID,
        'checks': {'geometric_involution': True, 'compression_independent': True,
                   'cyclic_sector_and_spectral_mapping': True, 'both_flat_sectors_retained': True,
                   'position_measurement_leaves_cyclic_sector': True,
                   'global_sign_has_controlled_readout': True, 'source_bytes_unchanged': True},
        'errors': {'compression': compression_error, 'cyclic_invariance': invariant_error,
                   'spectral_mapping': spectral_error, 'flat_eigenvectors': flat_error},
        'dimensions': multiplicities,
        'cyclic_sector_dimension': int(Q.shape[1]),
        'position_measurements': measurement,
        'controlled_sign_total_variation_distance': sign_distance,
        'scope': 'Exact kinematic sector theorem with numerical source controls. Native compiler, reset errors, autonomous field, physical clock and continuum interpretation are separate obligations.',
        'consequence': 'A continuum theorem restricted to span(A,SA) needs a preparation and instrument invariance proof. A complete sharp position instrument can populate the other physical sectors.',
        'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): sha(p) for p in paths}
    (HERE / 'KINEMATIC_REVIEW.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'errors', 'dimensions', 'controlled_sign_total_variation_distance')}))
    print(json.dumps({'position_measurement_flat_weights': [r['unconditional_weight_outside_dispersive_sector'] for r in measurement]}))


if __name__ == '__main__':
    run()
