"""Two real paths with a logical address in a common launch and receiving cell."""

from itertools import pairwise
from functools import cache

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight.measure import (
    LAW,
    closed_walk,
    edge,
    seam_data,
)
from model.executable_model import federation

# Port addresses in the twenty-owner carrier. The paths have disjoint interiors.
RAILS = (
    (0, 10, 20, 15, 14, 24, 25, 32, 46, 36, 37, 50, 51),
    (5, 1, 82, 78, 80, 192, 202, 194, 85, 84, 94, 56, 54),
)


@cache
def preparation():
    board = federation(0)
    graph = nx.Graph(board.seams)
    nodes = set(RAILS[0] + RAILS[1])
    assert len(nodes) == sum(map(len, RAILS))
    active = {edge(a, b) for rail in RAILS for a, b in pairwise(rail)}
    assert all(graph.has_edge(a, b) for a, b in active)
    cuts = sorted(e for e in board.seams if nodes.intersection(e) and e not in active)
    buffers = [r[i] for i in (0, -1) for r in RAILS]
    assert all(board.degree[b] == 5 for b in buffers)
    assert all(edge(*buffers[i : i + 2]) in cuts for i in (0, 2))
    private = sorted(n for n in nodes.difference(buffers) if board.degree[n] == 5)
    dist = nx.multi_source_dijkstra_path_length(graph, nodes)
    donor = max(range(board.n), key=lambda n: (dist[n], -n))
    # At the private write the mobile load and degree reference both read zero.
    word = np.ones(board.n, int)
    start = buffers[0]
    word[start], word[donor] = 0, 2
    for n in private:
        route = closed_walk(nx.shortest_path(graph, start, n))
        for a, b in route:
            word[a], word[b] = word[b], word[a]
        assert word[n] == 0
        for a, b in reversed(route):
            word[a], word[b] = word[b], word[a]
    assert word[start] == 0 and word[donor] == 2

    # For either charge, park it at rail 0's buffer before preparing its spin input.
    # Except for the launch shutter, the cuts are written on two filled loads.
    for cut in cuts:
        occupied = RAILS[0][1] if start in cut else start
        assert occupied not in cut and donor not in cut
    # Both rails may be populated on decoding. Reopen the receiving triangle
    # and cell by reversing the writes while shuttling any touched buffer.
    receiver = [(r[-1], r[-2]) for r in RAILS]
    home = receiver[0][0] // 12
    local_cuts = [c for c in cuts if all(x // 12 == home for x in c)]
    for cut in local_cuts:
        parked = [y if x in cut else x for x, y in receiver]
        assert not set(parked).intersection(cut)
    return dict(
        active=sorted(active),
        cuts=cuts,
        silent_private=private,
        buffers=buffers,
        compensator=donor,
        receiver_cuts=local_cuts,
    )


@cache
def pair(k, sign):
    blank, _ = seam_data(k)
    ids = [LAW.SRC.index(x) for x in ((1 + sign, 1), (1, 1 + sign))]
    return blank[np.ix_(ids, ids)]


def native_matrix(k, sign, cut=False):
    """Sum active physical seams and each unsilenced degree term above a filled rail."""
    prep = preparation()
    board = federation(0)
    ports = list(RAILS[0] + RAILS[1])
    indices = {p: i for i, p in enumerate(ports)}
    n = len(ports)
    h = np.zeros((n, n))
    for a, b in prep["active"]:
        if cut and edge(a, b) == edge(RAILS[1][1], RAILS[1][2]):
            continue
        ix = [indices[a], indices[b]]
        h[np.ix_(ix, ix)] += pair(k, sign)
    for i, p in enumerate(ports):
        if board.degree[p] == 5 and p not in prep["silent_private"]:
            h[i, i] += k * ((1 + sign) ** 2 - 1)
    return h


def spectral_word(target, k, sign):
    """Compile a real orbital on one rail from its buffer with native contacts."""
    n = len(target)
    state = np.zeros(n, complex)
    state[0] = 1
    phase_rate = -k * ((1 + sign) ** 2 - 1)
    contacts = 0

    def phase(site, angle):
        # Full adjacent pi swaps bring one address to the retained blank private
        # buffer, then restore all ports. The differential clock rate is physical.
        nonlocal contacts
        dwell = (angle * np.sign(phase_rate) % (2 * np.pi)) / abs(phase_rate)
        for j in range(site, 0, -1):
            state[[j, j - 1]] = state[[j - 1, j]]
        state[0] *= np.exp(1j * phase_rate * dwell)
        for j in range(site):
            state[[j, j + 1]] = state[[j + 1, j]]
        contacts += 2 * site + 1

    for a in range(n - 1):
        t = 2 * np.arctan2(np.linalg.norm(target[a + 1 :]), abs(target[a]))
        phase(a + 1, np.pi / 2)
        state[a : a + 2] = expm(-1j * t * pair(k, sign)) @ state[a : a + 2]
        phase(a, (k + 0.5) * t)
        phase(a + 1, (k + 0.5) * t - np.pi / 2)
        contacts += 1
    for i, amplitude in enumerate(target):
        phase(i, np.angle(amplitude))
    return state, contacts


def exact_lines(k, sign):
    """Simple Jacobi levels have nonzero weight at both measured gateways."""
    from fractions import Fraction

    n = len(RAILS[0])
    h2 = np.rint(2 * native_matrix(k, sign)[:n, :n]).astype(int)
    np.testing.assert_array_equal(np.diag(h2, 1), -np.ones(n - 1, int))
    # With p_j = det(d0 I - H_j) for the twice-energy open path, evaluate
    # the continuant exactly. A zero at the gateway would force E=d0/2.
    d0 = int(h2[0, 0])
    prev, cur = 1, 0
    for d in np.diag(h2)[1:]:
        prev, cur = cur, (d0 - int(d)) * cur - prev
    assert cur != 0
    assert Fraction(d0, 2) == Fraction(int(h2[-1, -1]), 2)
    return cur
