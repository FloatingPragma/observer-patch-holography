"""Complete capped return records and head updates on the registered 80-cell board."""

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import measure

ROUTES = measure.reuse.ROUTES
SCORES = np.r_[np.arange(2, 13) / 4, 3.0]
DELAYS = (0, 1, 2, 4, 8)


def instruments(contrast):
    """Return cycle channels and the disjoint effects for ticks 2--12 and timeout."""
    station = measure.station
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    start = np.zeros((4 * board.n, 16))
    start[measure.bridge.beam.indices(ports)] = np.eye(16)
    channels = np.zeros((2, 12, 16, 16), complex)
    effects = np.zeros((2, 12, 12, 4, 4), complex)
    for b in range(2):
        for charge in (-1, 1):
            h = station.iso.ranging.generator(board, ports, charge, contrast, bool(b))
            op = -0.5j * h
            u = expm_multiply(
                op, start, start=0, stop=12, num=13, traceA=op.diagonal().sum()
            )
            u = u.reshape(13, board.n, 4, 4, 4)
            sampled = u[:, ports].transpose(0, 1, 3, 2, 4)
            hits = {r: station.iso.ranging.first_hits(sampled, *r) for r in ROUTES}
            survive = {}
            for i, j in ROUTES:
                series = []
                for t in range(13):
                    s = u[t, :, :, i].copy()
                    for a in range(1, t + 1):
                        s -= u[t - a, :, :, j] @ hits[i, j][a]
                    vectors = s.reshape(board.n, 16)
                    series.append(
                        (vectors.T @ vectors.conj())
                        .reshape(4, 4, 4, 4)
                        .transpose(0, 2, 1, 3)
                        .reshape(16, 16)
                    )
                survive[i, j] = series
            for k, (i, j) in enumerate(ROUTES):
                cycle = survive[i, j][12].copy()
                clocks = np.zeros((12, 4, 4), complex)
                for a in range(1, 13):
                    f = hits[i, j][a]
                    cycle += survive[j, i][12 - a] @ np.kron(f, f.conj())
                    for d in range(1, 13 - a):
                        ret = hits[j, i][d] @ f
                        cycle += np.kron(ret, ret.conj())
                        clocks[a + d - 2] += ret.conj().T @ ret
                clocks[-1] = np.eye(4) - clocks[:-1].sum(axis=0)
                effects[b, k] += clocks / 2
                channels[b, k] += cycle / 2
    np.testing.assert_allclose(
        np.eye(4).ravel() @ channels,
        np.broadcast_to(np.eye(4).ravel(), (2, 12, 16)),
        atol=2e-12,
    )
    np.testing.assert_allclose(
        effects.sum(axis=2), np.broadcast_to(np.eye(4), (2, 12, 4, 4)), atol=1e-13
    )
    return channels, effects


def evolve(channel, state, n):
    return np.stack(
        [
            (np.linalg.matrix_power(phi, n) @ rho.ravel()).reshape(4, 4)
            for phi, rho in zip(channel, state, strict=True)
        ]
    )


def census(state, effects):
    return np.einsum("brtij,bji->rt", effects, state).real
