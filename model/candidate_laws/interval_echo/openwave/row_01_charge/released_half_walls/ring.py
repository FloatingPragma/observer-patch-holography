"""Native blank seam sector on an induced ring with exterior zero spectators."""

from functools import cache
from itertools import product

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.special import jv

from model.executable_model import federation
from .... import candidate

HOLDS = np.array([0, 0.5, 1, 2, 3, 4, 6])
ANGLES = np.pi * np.array([0, 0.25, -0.25, 0.5, -0.5, 1])
K = float(candidate.KAPPA)


@cache
def geometry(level, length):
    board = federation(level)
    graph = nx.Graph(board.seams)
    gates = {
        (a // 12, b // 12): a for u, v in board.external for a, b in ((u, v), (v, u))
    }
    owners = nx.Graph(gates.keys())
    cycles = []

    def extend(path):
        for other in sorted(owners[path[-1]]):
            if other == path[0] and len(path) in (5, 6) and path[1] < path[-1]:
                if owners.subgraph(path).number_of_edges() == len(path):
                    cycles.append(tuple(path))
            if len(path) < 6 and other > path[0] and other not in path:
                extend(path + [other])

    for o in sorted(owners):
        extend([o])
    for cycle in sorted(cycles):
        menus = []
        for j, o in enumerate(cycle):
            a = gates[o, cycle[j - 1]]
            b = gates[o, cycle[(j + 1) % len(cycle)]]
            local = graph.subgraph(range(12 * o, 12 * o + 12))
            menus.append(
                sorted(
                    p
                    for p in nx.all_simple_paths(local, a, b, cutoff=3)
                    if local.subgraph(p).number_of_edges() == len(p) - 1
                )
            )
        # A later path of the same length cannot change feasibility in
        # another owner. Retain the first witness for each local length.
        unique = []
        for choices in menus:
            first = {}
            for path in choices:
                first.setdefault(len(path), path)
            unique.append(tuple(first.values()))
        menus = unique
        for segments in product(*menus):
            if sum(map(len, segments)) != length:
                continue
            ports = [i for p in segments for i in p]
            assert graph.subgraph(ports).number_of_edges() == length
            for j, site in enumerate(ports):
                active = {ports[j - 1], ports[(j + 1) % length]}
                assert active <= set(board.neighbours[site])
                assert (
                    len(set(board.neighbours[site]) - active) + 6 - board.degree[site]
                    == 4
                )
            return ports
    raise ValueError("Missing enrolled induced circuit")


@cache
def columns():
    law = candidate.native_law
    swap, contrast = law.local_term_parts()
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    part = blank.T @ swap @ blank
    c = blank.T @ contrast @ blank
    np.testing.assert_allclose(swap @ blank, blank @ part, atol=1e-13)
    np.testing.assert_allclose(contrast @ blank, blank @ c, atol=1e-13)
    for i, (a, b) in enumerate(law.SRC):
        expected = np.zeros(9)
        expected[i] += 0.5
        expected[law.SRC.index((b, a))] -= 0.5
        np.testing.assert_allclose(part[:, i], expected, atol=1e-13)
        target = np.zeros(9)
        target[i] = (a - b) ** 2
        np.testing.assert_allclose(c[:, i], target, atol=1e-13)
    return part, c


@cache
def basis(n):
    words = np.array([q for q in range(1 << n) if q.bit_count() == n // 2 + 1])
    sites = ((words[:, None] >> np.arange(n)) & 1).astype(np.int8)
    return words, sites


def launch(n):
    checker = sum(1 << j for j in range(1, n, 2))
    return checker | 1


def generator(n, amount, k):
    part, contrast = columns()
    words, sites = basis(n)
    idx = np.arange(len(words))
    diag = np.zeros(len(words))
    row = []
    col = []
    vals = []
    codes = np.array(
        [
            [candidate.native_law.SRC.index((amount * a, amount * b)) for b in (0, 1)]
            for a in (0, 1)
        ]
    )
    for i in range(n):
        j = (i + 1) % n
        a, b = sites[:, i], sites[:, j]
        local = codes[a, b]
        diag += part[local, local] + k * contrast[local, local]
        unequal = a != b
        row.extend(idx[unequal])
        col.extend(np.searchsorted(words, words[unequal] ^ (1 << i) ^ (1 << j)))
        vals.extend(part[codes[b[unequal], a[unequal]], local[unequal]])
    # Four physical fixed-zero contrast slots per ring site give a common
    # scalar in each sector; subtract that scalar and the diagonal centre.
    diag -= 0.5 * (max(diag) + min(diag))
    return sparse.diags(diag, format="csr") + sparse.csr_matrix(
        (vals, (row, col)), shape=(len(idx), len(idx))
    )


@cache
def released(n, amount, k):
    h = generator(n, amount, k)
    words, _ = basis(n)
    p = np.zeros(len(words))
    p[np.searchsorted(words, launch(n))] = 1
    radius = max(abs(h.diagonal())) + n / 2
    h = h / radius
    z = HOLDS * radius
    order = int(np.ceil(z.max() + 18 * np.cbrt(z.max() + 1) + 40))
    powers = np.arange(order + 1)
    coeff = 2 * (-1j) ** (powers % 4) * jv(powers[None, :], z[:, None])
    coeff[:, 0] /= 2
    assert abs(coeff[:, -3:]).max() < 1e-20
    q = h @ p
    out = coeff[:, 0, None] * p + coeff[:, 1, None] * q
    for j in range(2, order + 1):
        p, q = q, 2 * (h @ q) - p
        out += coeff[:, j, None] * q
    result = abs(out) ** 2
    np.testing.assert_allclose(result.sum(axis=1), 1, atol=9e-12, rtol=0)
    return result


def regional(n, sign, amount):
    _, s = basis(n)
    if sign < 0:
        s = 1 - s
    bonds = s + np.roll(s, -1, axis=1) - 1
    left = np.array([-5, -4, -3]) % n
    right = np.array([2, 3, 4])
    two = (bonds != 0).sum(axis=1) == 2
    gates = (
        two
        & ((bonds[:, left] == sign).sum(axis=1) == 1)
        & ((bonds[:, right] == sign).sum(axis=1) == 1)
    )
    tapered = np.stack(
        [bonds[:, arm].sum(axis=1) * amount for arm in (left, right)], axis=1
    )  # doubled units
    flat = np.stack(
        [
            (s[:, np.r_[arm, arm[-1] + 1]].sum(axis=1) - 2) * amount
            for arm in (left, right)
        ],
        axis=1,
    )
    return gates, two, tapered, flat
