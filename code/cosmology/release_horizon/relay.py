"""Two sharp material registers, a local consultation service and its read cone."""

from functools import cache
from itertools import product
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from model.executable_model import CELL_SEAMS, federation
from model.candidate_laws.interval_echo import observations
from model.experiments.interval_echo_consensus.refinement_cooling.dimer_band import (
    measure as dimer,
)

CELL = nx.Graph(CELL_SEAMS)
STATES = list(product(range(2), repeat=2))


def transpose(a, b):
    path = nx.shortest_path(CELL, a, b)
    edges = list(zip(path, path[1:]))
    return edges + edges[:-1][::-1]


def align(slots, targets):
    labels = list(range(12))
    ops = []
    for slot, target in zip(slots, targets, strict=True):
        for a, b in transpose(labels.index(slot), target):
            labels[a], labels[b] = labels[b], labels[a]
            ops.append((a, b))
    assert [labels[t] for t in targets] == list(slots)
    return ops


def apply(word, ops):
    for a, b in ops:
        word[a], word[b] = word[b], word[a]


def sharp_state(x, d):
    # All five dimers are occupied. The other public and private rails are dark.
    w = [0] * 12
    for p in (x, 2 + d, 4, 6, 8):
        w[p] = 2
    return w


@cache
def transaction(g, h):
    """Transport the actual bright rail (or vacancy) through a single external seam."""
    ref = min(CELL[h])
    ret = align([10, 11], [h, ref])
    sender = [align([p], [g]) for p in [1, 10]]
    for x, d, y, f, live in product(range(2), repeat=5):
        a, b = sharp_state(x, d), sharp_state(y, f)
        old_a, old_b = a[:], b[:]
        apply(a, sender[d])
        apply(b, ret)
        if live:
            a[g], b[h] = b[h], a[g]
        assert b[ref] == 0
        received = int(b[h] == 2)
        if live:
            a[g], b[h] = b[h], a[g]
        apply(a, sender[d][::-1])
        apply(b, ret[::-1])
        assert a == old_a and b == old_b and received == live * x * (1 - d)
        if received and not y:
            apply(b, transpose(0, 1))
        assert b == sharp_state(y | received, f)
    return dict(sender=sender, receiver=ret, reference=ref)


def gates():
    law = observations.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + observations.KAPPA * contrast
    idx = [law.SRC.index(p) for p in [(2, 0), (0, 2)]]
    # This invariant two-rail block gives change probability sin(t/2)^2
    # for every exposure, including the branches of a general scheduler row.
    reduced = law.seam_term(observations.KAPPA)
    expected = np.zeros((law.DIM, 2))
    expected[: len(law.SRC)] = reduced[:, idx]
    np.testing.assert_allclose(h[:, idx], expected, atol=2e-13, rtol=0)
    np.testing.assert_array_equal(
        reduced[np.ix_(idx, idx)], [[32.5, -0.5], [-0.5, 32.5]]
    )
    flip = [law.SRC.index((b, a)) for a, b in law.SRC]
    for active, expected in [
        (h, flip),
        (observations.KAPPA * contrast, list(range(9))),
    ]:
        np.testing.assert_allclose(
            expm(-1j * np.pi * active)[:, :9],
            np.eye(law.DIM)[:, expected],
            atol=2e-13,
            rtol=0,
        )
    dimer.pointer()
    # The isolated decoder maps position labels to energies by inverse evolution.
    site, decode, _ = dimer.gates("native")
    levels, vectors = np.linalg.eigh(site[:2, :2])
    np.testing.assert_allclose(levels, [32, 33], atol=2e-13, rtol=0)
    np.testing.assert_allclose(
        abs(decode[:2, :2] @ vectors) ** 2, np.eye(2), atol=2e-13
    )
    for p in (0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 1):
        pulse = expm(-2j * np.arcsin(np.sqrt(p)) * h)[:, idx]
        target = np.zeros_like(abs(pulse))
        target[idx, :] = [[1 - p, p], [p, 1 - p]]
        np.testing.assert_allclose(abs(pulse) ** 2, target, atol=2e-13, rtol=0)
    # The selector starts on its current marker rail. Read the position,
    # restore it on the change branch, then advance it by a deterministic stroke.
    swap = expm(-1j * np.pi * h)[np.ix_(idx, idx)]
    for start, changed in product(range(2), repeat=2):
        after = start ^ changed
        restored = swap[:, after] if changed else np.eye(2)[:, after]
        np.testing.assert_allclose(
            abs(swap @ restored) ** 2, np.eye(2)[:, 1 - start], atol=2e-13
        )

    # Use a full controlled dimer swap, including both eigenvectors. Energy
    # states can be coherent superpositions of positions at checkpoints.
    for a, b in product(range(2), repeat=2):
        v = np.kron(vectors[:, a], vectors[:, b])
        route = transpose(4, 6) + transpose(5, 7)
        out = np.zeros(4, complex)
        for i, j in product(range(2), repeat=2):
            word = sharp_state(1, 1)
            word[4:8] = [0] * 4
            word[4 + i] = word[6 + j] = 2
            apply(word, route)
            assert (
                word[:4] == sharp_state(1, 1)[:4] and word[8:] == sharp_state(1, 1)[8:]
            )
            out[2 * (word[5] == 2) + (word[7] == 2)] += v[2 * i + j]
        np.testing.assert_allclose(
            out, np.kron(vectors[:, b], vectors[:, a]), atol=2e-13
        )
    # Any prepared occupancy can be addressed on the dimer/filter's path.
    for offset in (0, 2, 4, 6, 8):
        align([offset, offset + 1, 10, 11], [0, 1, 7, 6])
    costs = {}
    for x, d in STATES:
        q = (x, d)
        rset = {q, (1 - x, d), (x, 0)}
        assert rset == {(y, f) for y, f in STATES if f <= d and (x == y or d == f)}
        for y, f in rset:
            word = sharp_state(x, d)
            if y != x:
                apply(word, transpose(0, 1))
            if f != d:
                apply(word, transpose(2, 3))
            assert word == sharp_state(y, f)
        costs[str(q)] = [
            dict(output=r, shared=int(r[0] != x), private=int(r[1] != d))
            for r in sorted(rset)
        ]
    # Realize a three-output row on (0,1): publish, clear, idle.
    # A second local selector is drawn only on the first dark branch.
    weights = []
    remaining = 1.0
    for p in (1 / 5, 1 / 3):
        pulse = expm(-2j * np.arcsin(np.sqrt(p)) * h)[:, idx[0]]
        weights.append(remaining * abs(pulse[idx[1]]) ** 2)
        remaining *= abs(pulse[idx[0]]) ** 2
    np.testing.assert_allclose(
        weights + [remaining], [1 / 5, 4 / 15, 8 / 15], atol=2e-13
    )
    return costs


def topology(level):
    board = federation(level)
    g = nx.Graph()
    groups = []
    occupied = []
    for a, b in board.external:
        u, v = a // 12, b // 12
        g.add_edge(u, v)
        transaction(a % 12, b % 12)
        transaction(b % 12, a % 12)
        group = next(
            (i for i, s in enumerate(occupied) if u not in s and v not in s),
            len(groups),
        )
        if group == len(groups):
            groups.append([])
            occupied.append(set())
        groups[group].append((u, v))
        occupied[group].update((u, v))
    assert sorted(g) == list(range(board.cells)) and nx.is_connected(g)
    dist = shortest_path(
        nx.to_scipy_sparse_array(g, nodelist=range(board.cells), format="csr").astype(
            float
        ),
        directed=False,
        unweighted=True,
    ).astype(int)
    return board, groups, dist


def propagate(groups, n, root, depth, cut=False, stagger=True):
    # Shared writes are individually visible during consultation. The pending
    # bit gates forwarding until the subsequent separate acknowledgement.
    x = np.zeros(n, bool)
    d = np.ones(n, bool)
    x[root] = True
    d[root] = False
    arrival = np.full(n, -1, int)
    arrival[root] = 0
    clock = 0
    changes = 0
    for t in range(1, depth + 1):
        marker = 1 - clock
        for edges in groups:
            a, b = np.asarray(edges).T
            for u, v in ((a, b), (b, a)):
                seen = x[u] & ~d[u] if not cut else np.zeros(len(u), bool)
                fresh = v[seen & ~x[v]]
                x[fresh] = True
                arrival[fresh] = t
                if not stagger:
                    d[fresh] = False
        d[x] = False
        changes += int(clock != marker)
        clock = marker
    assert changes == depth
    return arrival
