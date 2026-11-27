"""Read the spatial Gram from two independently sampled arrival tables."""

from math import log, sqrt
import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo.openwave.row_29_maxwell.seam_symbol.measure import (
    spectrum,
)

SHOTS = 2**22
THRESHOLD = 1 / 200
RADIUS = sqrt(log(2 * 288 / 0.01) / (2 * SHOTS)) + 12 / 2**44


def arrival_counts():
    energies, modes = spectrum("cell", 1)
    flight = (modes * np.exp(-0.25j * energies)) @ modes.T
    prob = abs(flight) ** 2
    rng = np.random.Generator(np.random.PCG64(1652601))
    banks = []
    for _ in range(2):
        rows = []
        for row in prob:
            rounded = np.rint(row * 2**44).astype(np.int64)
            rounded[np.argmax(row)] += 2**44 - rounded.sum()
            rows.append(rng.multinomial(SHOTS, rounded / 2**44).tolist())
        banks.append(rows)
    return {"sampled": banks, "parked": (SHOTS * np.eye(12, dtype=int)).tolist()}


def graph_from_bank(counts):
    a = np.array(counts)
    if (
        a.shape != (12, 12)
        or a.dtype.kind not in "iu"
        or np.any(a < 0)
        or np.any(a.sum(axis=1) != SHOTS)
    ):
        raise ValueError("Incomplete arrival tables")
    p = a / SHOTS
    mask = ~np.eye(12, dtype=bool)
    if np.any(abs(p[mask] - THRESHOLD) <= RADIUS):
        raise ValueError("Unresolved contact")
    edges = (p > THRESHOLD) & mask
    if not np.array_equal(edges, edges.T) or np.any(edges.sum(axis=1) != 5):
        raise ValueError("Contact standard failed")
    graph = nx.from_numpy_array(edges)
    if not nx.is_isomorphic(graph, nx.icosahedral_graph()):
        raise ValueError("Receiver graph differs from the enrolled frame")
    return graph


def squared_separations(graph):
    # The unit-radius port Gram is R + S/sqrt(5). Its positive eigenvalue
    # is four and its columns belong to the adjacency eigenvalue sqrt(5).
    a = nx.to_numpy_array(graph, nodelist=range(12), dtype=int)
    path = np.asarray(nx.floyd_warshall_numpy(graph, nodelist=range(12)))
    r = (path == 0).astype(int) - (path == 3)
    s = (path == 1).astype(int) - (path == 2)
    if not (
        np.array_equal(a @ r, s)
        and np.array_equal(a @ s, 5 * r)
        and np.array_equal(5 * r @ r + s @ s, 20 * r)
        and np.array_equal(r @ s + s @ r, 4 * s)
        and r.trace() == 12
    ):
        raise ValueError("The recovered graph does not supply the rank-three slow Gram")
    edges = sorted(tuple(sorted(e)) for e in graph.edges)
    twice = np.zeros((len(edges), 12), dtype=int)
    for k, (i, j) in enumerate(edges):
        twice[k, [i, j]] = 1
    distances = []
    for g in (5 * r, s):
        h = twice @ g @ twice.T
        distances.append(h.diagonal()[:, None] + h.diagonal()[None, :] - 2 * h)
    # 20 times each squared midpoint distance is a+b*sqrt(5).
    return edges, [
        [(int(distances[0][i, j]), int(distances[1][i, j])) for j in range(30)]
        for i in range(30)
    ]
