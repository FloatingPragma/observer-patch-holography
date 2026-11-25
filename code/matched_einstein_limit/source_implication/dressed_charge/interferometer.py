"""Read the compact spatial coherences using isolated native blank-seam pulses."""

from collections import deque
from itertools import combinations

import numpy as np

from ..full_algebra import station


PORTS = (4, 5, 6, 7, 8, 9, 11)


def register(board):
    """Check both signs in all 99 contact rows and return the seven-port read circuit."""
    emission = station.iso.ranging.formation.emission
    h4 = emission.four_h(8)
    for sign in (-1, 1):
        ix = [emission.LAW.SRC.index(pair) for pair in ((1 + sign, 1), (1, 1 + sign))]
        columns = np.eye(len(h4), dtype=int)[:, ix]
        np.testing.assert_array_equal(
            h4 @ columns, columns @ np.array([[34, -2], [-2, 34]])
        )
    uniform = emission.LAW.SRC.index((1, 1))
    np.testing.assert_array_equal(h4[:, uniform], np.zeros(len(h4), int))
    return {
        i: [PORTS.index(q) for q in board.neighbours[p] if q in PORTS]
        for i, p in enumerate(PORTS)
    }


def gate(i, j, t):
    # The symmetric and antisymmetric blank modes have exact energies 8 and 9.
    u = np.eye(len(PORTS), dtype=complex)
    a, b = np.exp(-1j * 8 * t), np.exp(-1j * 9 * t)
    u[np.ix_([i, j], [i, j])] = np.array([[a + b, a - b], [a - b, a + b]]) / 2
    return u


def path(adj, start, end):
    todo = deque([[start]])
    seen = {start}
    while todo:
        trail = todo.popleft()
        if trail[-1] == end:
            return trail
        for nxt in adj[trail[-1]]:
            if nxt not in seen:
                seen.add(nxt)
                todo.append(trail + [nxt])
    raise ValueError("blank-seam read region disconnected")


def signed_y(adj, i, j):
    trail = path(adj, i, j)
    u = np.eye(len(PORTS), dtype=complex)
    for v, w in zip(trail[:-2], trail[1:-1]):
        u = gate(v, w, np.pi) @ u
    a, b = trail[-2:]
    u = gate(a, b, np.pi / 2) @ u
    signed = np.zeros(len(PORTS))
    signed[a], signed[b] = -1, 1
    return u.conj().T @ (signed[:, None] * u)


def terminal_settings(adj, spatial):
    """List of (weight, terminal detector effect) for the whole compiled spatial meter."""
    settings = [(1.0, np.diag(np.diag(spatial)))]
    for u, v in combinations(range(len(PORTS)), 2):
        weight = spatial[u, v]
        if weight == 0:
            continue
        w = next(x for x in adj[v] if x != u)
        y_uv, y_uw = signed_y(adj, u, v), signed_y(adj, u, w)
        for energy, sign in ((8, 1), (9, -1)):
            pre = gate(v, w, np.pi / (2 * energy))
            settings.append((-weight / 2, pre.conj().T @ y_uv @ pre))
            settings.append((-weight * sign / 2, pre.conj().T @ y_uw @ pre))
    return settings


def spatial_score(wave, rho, charge):
    kept = wave[list(PORTS)]
    return np.einsum(
        "ij,xai,ab,ybj->xy", rho.T, kept.conj(), charge, kept, optimize=True
    )
