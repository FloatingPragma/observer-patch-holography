"""Native translation blocks for odd alternating load clusters on the physical route."""

from functools import cache
from itertools import combinations
from math import pi

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

from ..conjugate_flight.measure import LAW, edge, layout, seam_data


def translate(word, j):
    j %= len(word)
    return word[-j:] + word[:-j] if j else word


def core_word(n, d):
    return tuple((0 if a % 2 == 0 else 2) if a < d else 1 for a in range(n))


def prepare(n_defects):
    """Write every unused local contact with equal loads, returning the sharp core each time."""
    board, _, ring, cuts, _, donor = layout()
    n = len(ring)
    core = core_word(n, n_defects)
    board_word = np.ones(board.n, dtype=int)
    board_word[ring] = core
    board_word[donor] = 2
    assert sum(board_word - 1) == 0
    before = board_word.copy()
    writes = []

    def write(kind, endpoints, allowed):
        for shift in range(n):
            # A right rotation is n-1 adjacent swaps from the far end toward zero.
            candidate = board_word.copy()
            candidate[ring] = translate(core, shift)
            if allowed(candidate):
                break
        else:
            raise AssertionError("No sharp write with the correct loads")
        swaps = list(zip(ring[:-1], ring[1:]))[::-1] * shift
        for a, b in swaps:
            assert edge(a, b) in board.seams
            board_word[a], board_word[b] = board_word[b], board_word[a]
        assert np.array_equal(board_word, candidate)
        for port in endpoints:
            assert board_word[port] == (0 if kind == "private" else 1)
        writes.append([kind, endpoints, shift])
        for a, b in reversed(swaps):
            board_word[a], board_word[b] = board_word[b], board_word[a]
        assert np.array_equal(board_word, before)

    for port in ring:
        if board.degree[port] == 5:
            write("private", [port], lambda w: w[port] == 0)
    for a, b in cuts:
        write("physical", [a, b], lambda w: w[a] == w[b] == 1)
    live = set(board.seams) - set(cuts)
    incident = {e for e in live if set(e).intersection(ring)}
    assert incident == {edge(ring[a - 1], ring[a]) for a in range(n)}
    return dict(route=ring, compensator=donor, initial_route=core, writes=writes)


@cache
def orbits(n, d):
    holes, doubles = (d + 1) // 2, (d - 1) // 2
    representatives = []
    address = {}

    def enter(word):
        if word in address:
            return
        i = len(representatives)
        representatives.append(word)
        for shift in range(n):
            moved = translate(word, shift)
            assert moved not in address
            address[moved] = (i, shift)

    enter(core_word(n, d))
    for vacant in combinations(range(n), holes):
        for excess in combinations(sorted(set(range(n)) - set(vacant)), doubles):
            word = [1] * n
            for a in vacant:
                word[a] = 0
            for b in excess:
                word[b] = 2
            enter(tuple(word))
    return representatives, address


@cache
def columns(n, d, kappa):
    """Compress the blank block checked on all nine native inputs; every inactive write is dark."""
    native, _ = seam_data(kappa)
    h = np.rint(2 * native).astype(int)
    np.testing.assert_allclose(h / 2, native, atol=1e-13)
    words, address = orbits(n, d)
    rows, cols, shifts, values, diagonal = [], [], [], [], []
    outside_bound = 0
    for j, word in enumerate(words):
        cost = tail = 0
        for a in range(n):
            b = (a + 1) % n
            p = LAW.SRC.index((word[a], word[b]))
            cost += h[p, p]
            if word[a] == word[b]:
                continue
            other = list(word)
            other[a], other[b] = other[b], other[a]
            i, s = address[tuple(other)]
            element = h[LAW.SRC.index((other[a], other[b])), p]
            assert element == -1
            rows.append(i)
            cols.append(j)
            shifts.append(s)
            values.append(element / 2)
            if i and j:
                tail += abs(element)
        diagonal.append(cost / 2)
        if j:
            outside_bound = max(outside_bound, (cost + tail) / 2)
    if d > 1:
        assert all(i or j for i, j in zip(rows, cols))
        assert cols.count(0) == d + 1
        # The translated columns of QHP have disjoint support: no noncore word meets two cores.
        neighbors = [j for i, j in zip(rows, cols) if i == 0]
        assert len(neighbors) == len(set(neighbors)) == d + 1
    return (
        np.array(rows),
        np.array(cols),
        np.array(shifts),
        np.array(values),
        np.array(diagonal),
        outside_bound,
    )


def sector(n, d, kappa, mode):
    rows, cols, shifts, values, diagonal, upper = columns(n, d, kappa)
    phase = np.exp(-2j * pi * mode * shifts / n)
    size = len(diagonal)
    offset = diagonal[0]
    h = coo_matrix((values * phase, (rows, cols)), shape=(size, size)).tocsr()
    h = h + diags(diagonal - offset)
    return h, offset, upper


@cache
def eigenpair(n, d, kappa, m):
    h, center, other_upper = sector(n, d, kappa, m)
    if h.shape == (1, 1):
        return float(h[0, 0].real), np.array([1.0 + 0j]), 0.0, center
    start = np.zeros(h.shape[0])
    start[0] = 1
    vals, vecs = eigsh(h, k=1, which="LA", v0=start, tol=1e-14, ncv=14)
    v = vecs[:, 0]
    # Recenter before the Rayleigh quotient so the narrow recoil survives subtraction.
    energy = float(np.vdot(v, h @ v).real / np.vdot(v, v).real)
    residual = np.linalg.norm(h @ v - energy * v)
    assert center + energy - residual > other_upper and residual < 2e-12
    return energy, v, float(residual), center


def line(n, d, kappa, m):
    energy, v, residual, center = eigenpair(n, d, kappa, m)
    return dict(
        offset=energy,
        core_weight=float(abs(v[0]) ** 2),
        residual=float(residual),
        center=center,
    )
