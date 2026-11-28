"""Exact fastest-route distances and exterior release for several transmitting seams."""

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
import math

import numpy as np
from scipy.sparse.csgraph import dijkstra, shortest_path


@dataclass
class Survey:
    blank: object
    contacts: list
    center: int
    radii: object
    clearance: object
    starts: object
    index: object
    pairs: object
    neighbors: object
    degree: object


def survey(board, layout, external_only):
    marks = [seam for _, edges in layout for seam in edges]
    blank = board.adjacency().tolil()
    for a, b in marks:
        blank[a, b] = blank[b, a] = 0
    blank = blank.tocsr()
    source_ranges = shortest_path(
        blank, indices=np.unique(marks), unweighted=True, directed=False
    )
    assert np.isfinite(source_ranges).all()
    center = int(source_ranges.max(0).argmin())
    clearance = source_ranges.min(0).astype(int)
    radii = shortest_path(
        blank, indices=center, unweighted=True, directed=False
    ).astype(int)
    # Range and clearance use all three written seams of every source.
    eligible = np.ones(board.n, bool)
    if external_only:
        eligible[:] = False
        eligible[np.unique(board.external)] = True
    eligible &= clearance >= 2
    pairs = np.array(
        [
            p
            for r in np.unique(radii[eligible])
            for p in combinations(np.flatnonzero((radii == r) & eligible), 2)
        ],
        int,
    ).reshape(-1, 2)
    starts = np.unique(np.r_[center, pairs.ravel()])
    index = np.full(board.n, -1, int)
    index[starts] = np.arange(len(starts))
    return Survey(
        blank,
        [edges[-1] for _, edges in layout],
        center,
        radii,
        clearance,
        starts,
        index,
        pairs,
        board.table,
        board.degree,
    )


def nearest_sources(s, costs):
    if len(s.contacts) < 2:
        return []
    endpoints = np.array(s.contacts)

    def nearest(net):
        table = dijkstra(net, directed=False, indices=endpoints.ravel())[
            :, endpoints.ravel()
        ]
        distances = table.reshape(len(endpoints), 2, len(endpoints), 2).min((1, 3))
        np.fill_diagonal(distances, np.inf)
        return distances.min(1)

    network = s.blank.tolil()
    for (a, b), c in zip(s.contacts, costs, strict=True):
        network[a, b] = network[b, a] = float(c)
    # A route starting at either endpoint and stopping at its first
    # other source has no occasion to pay a transmitting-link toll.
    blank = nearest(s.blank)
    np.testing.assert_array_equal(nearest(network.tocsr()), blank)
    return blank.astype(int).tolist()


def encode_lengths(costs, diameter):
    # Every shortest route is simple. It can use a marked contact at most
    # once and no more blank hops than the length of the blank-only route.
    unique, counts = np.unique(np.asarray(costs, object), return_counts=True)
    exact = []
    for uses in product(*(range(int(n) + 1) for n in counts)):
        paid = sum((u * c for u, c in zip(uses, unique, strict=True)), Q(0))
        exact.extend((k + paid, k, uses) for k in range(diameter + 1))
    exact.sort()
    for bits in range(16, 48, 4):
        scale = 1 << bits
        price = [round(c * scale) for c in unique]
        codes = [
            k * scale + sum(n * p for n, p in zip(ns, price, strict=True))
            for _, k, ns in exact
        ]
        if max(codes) >= 2**52:
            break
        if all(
            (u < v and a < b) or (u == v and a == b)
            for u, v, (a, *_), (b, *_) in zip(
                codes[:-1], codes[1:], exact[:-1], exact[1:], strict=True
            )
        ):
            grid = {code: x[0] for code, x in zip(codes, exact, strict=True)}
            return scale, np.array([round(c * scale) for c in costs]), grid
    raise ArithmeticError("need an exact shortest path for closer candidate costs")


@dataclass
class Read:
    length: object
    curvature: object
    equal_radial: object


def measure(s, costs, witness=None):
    scale, prices, grid = encode_lengths(costs, 2 * int(s.radii.max()))
    network = (s.blank * scale).tolil()
    for (a, b), price in zip(s.contacts, prices, strict=True):
        network[a, b] = network[b, a] = price
    starts = s.starts if witness is None else np.unique([s.center, *witness])
    idx = (
        s.index if witness is None else np.searchsorted(starts, np.arange(len(s.index)))
    )
    raw = dijkstra(network.tocsr(), directed=False, indices=starts)
    # All edge additions and comparisons are exact integers in the binary64
    # domain. Recover their original rational length after the route choice.
    codes = np.array(list(grid), dtype=np.int64)
    rank = np.searchsorted(codes, raw)
    assert np.array_equal(codes[rank], raw)
    denominator = math.lcm(*(c.denominator for c in costs))
    lengths = np.array(list(grid.values()), object)
    ticks = np.array([int(x * denominator) for x in lengths], object)
    square = ticks**2
    if witness is not None:
        U, V = witness
        neighbors = [(V, k) for k in s.neighbors[U] if k >= 0] + [
            (U, k) for k in s.neighbors[V] if k >= 0
        ]

        def get(a, b):
            return lengths[rank[idx[a], b]]

        return get(U, V), [get(a, b) for a, b in neighbors], neighbors
    U, V = s.pairs.T
    chord = rank[idx[U], V]
    equal_radial = rank[idx[s.center], U] == rank[idx[s.center], V]
    cross = np.zeros(len(U), object)
    for start, opposite in ((U, V), (V, U)):
        for j in range(s.neighbors.shape[1]):
            dest = s.neighbors[start, j]
            use = dest >= 0
            cross[use] += square[rank[idx[opposite[use]], dest[use]]]
    initial = square[chord]
    numerator = cross - (s.degree[U] + s.degree[V]).astype(object) * (
        denominator**2 + initial
    )

    @lru_cache(maxsize=None)
    def ratio(a, b):
        return Q(a, b)

    return Read(
        lengths[chord],
        np.array(
            [ratio(a, b) for a, b in zip(numerator, initial, strict=True)], object
        ),
        equal_radial,
    )
