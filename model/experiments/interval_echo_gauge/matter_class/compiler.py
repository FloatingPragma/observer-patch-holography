"""Universal load-address compilation for the twelve-dimensional B15 budget."""

from functools import cache

import networkx as nx
import numpy as np

from model.executable_model import federation
from ..scalar_grammar import measure as pulse


@cache
def to_gate(i, j):
    """Route indicated columns to local addresses 0,5 and record every native swap."""
    board = federation(0)
    graph = nx.Graph(board.seams).subgraph(range(12))
    labels = list(range(12))
    moves = []
    for src, target in [(i, 0), (j, 5)]:
        domain = graph if target == 0 else graph.subgraph(set(graph) - {0})
        path = nx.shortest_path(domain, labels.index(src), target)
        for a, b in zip(path, path[1:]):
            moves.append((a, b))
            labels[a], labels[b] = labels[b], labels[a]
    assert labels[0] == i and labels[5] == j
    return tuple(moves)


def compile_unitary(target):
    """Left elimination by diagonal phases and real rotations, then reverse the word."""
    u = np.asarray(target, dtype=complex).copy()
    if u.ndim != 2 or u.shape[0] != u.shape[1]:
        raise ValueError("The candidate action must be square")
    n = len(u)
    if not (1 <= n <= 12):
        raise ValueError("The candidate dimension must belong to the frozen budget")
    np.testing.assert_allclose(u.conj().T @ u, np.eye(n), rtol=0, atol=3e-13)
    reduction = []
    for column in range(n):
        for j in range(n - 1, column, -1):
            i = j - 1
            a, b = u[[i, j], column]
            if abs(b) < 1e-15:
                continue
            for row, entry in [(i, a), (j, b)]:
                angle = -np.angle(entry)
                reduction.append(("phase", row, row, angle))
                u[row] *= np.exp(1j * angle)
            theta = np.arctan2(abs(b), abs(a))
            c, s = np.cos(theta), np.sin(theta)
            u[[i, j]] = [[c, s], [-s, c]] @ u[[i, j]]
            reduction.append(("real", i, j, 2 * theta))
    for i in range(n):
        reduction.append(("phase", i, i, -np.angle(u[i, i])))
    return [(kind, i, j, -angle) for kind, i, j, angle in reversed(reduction)]


def execute(word):
    # π on an arbitrary bare 00/10/01 pair exchanges the two addresses
    # with unit amplitude. Verify against the same native 99-column service.
    swap = pulse.gate("contact", 0, 1, np.pi, pulse.candidate.KAPPA)
    np.testing.assert_allclose(swap, np.eye(3)[[1, 0, 2]], rtol=0, atol=2e-13)
    u = np.eye(12, dtype=complex)
    for kind, i, j, angle in word:
        parking = to_gate(i, j if kind == "real" else (i + 1) % 12)
        for a, b in parking:
            u[[a, b]] = u[[b, a]]
        program = pulse.ry(0, 1, angle) if kind == "real" else [pulse.phase(0, angle)]
        at_gate = pulse.execute(program)
        # Every off-edge address has the empty word 00 with exact energy zero.
        np.testing.assert_allclose(at_gate[2], [0, 0, 1], atol=1e-13)
        u[[0, 5]] = at_gate[:2, :2] @ u[[0, 5]]
        for a, b in reversed(parking):
            u[[a, b]] = u[[b, a]]
    return u
