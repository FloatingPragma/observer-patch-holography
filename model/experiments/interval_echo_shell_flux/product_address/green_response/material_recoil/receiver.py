"""Native contact phases and motion of one independently stocked test load."""

from collections import deque
from functools import cache
from itertools import product

import numpy as np
from scipy.linalg import expm

from ... import measure as rails

LEFT, RIGHT = 49, 53
SENSE, REFERENCE = 55, 52
CLOCK = np.pi / 32.5
FLIGHT = np.pi / np.hypot(1, 32.5)
MASKS = tuple(product((0, 1), repeat=3))


@cache
def pad_routes():
    board, axes, _ = rails.layout()
    assert board.degree[LEFT] == board.degree[RIGHT]
    assert RIGHT in board.neighbours[LEFT]
    assert SENSE in board.neighbours[LEFT] and SENSE not in board.neighbours[RIGHT]
    assert (
        REFERENCE in board.neighbours[RIGHT] and REFERENCE not in board.neighbours[LEFT]
    )
    paths = {SENSE: (SENSE,)}
    queue = deque([SENSE])
    while queue:
        u = queue.popleft()
        for v in sorted(map(int, board.neighbours[u])):
            if v not in paths and v not in (LEFT, RIGHT, REFERENCE):
                paths[v] = (*paths[u], v)
                queue.append(v)
    routes = {}
    for axis in axes:
        for port in axis:
            path = paths[port]
            steps = list(zip(path, path[1:]))
            swaps = steps + steps[-2::-1] if len(steps) > 1 else steps
            labels = list(range(board.n))
            for a, b in swaps:
                labels[a], labels[b] = labels[b], labels[a]
            target = list(range(board.n))
            target[SENSE], target[port] = target[port], target[SENSE]
            assert labels == target
            routes[port] = path
    return routes


@cache
def terms(q, k):
    law = rails.candidate.native_law
    seam = law.seam_term(k)
    exchange, contrast = law.local_term_parts()
    np.testing.assert_allclose(
        (exchange + k * contrast)[:, :9],
        np.vstack((seam, np.zeros((90, 9)))),
        atol=1e-14,
    )
    select = [law.SRC.index((q, 0)), law.SRC.index((0, q))]
    diag = seam.diagonal().reshape(3, 3)
    flight = seam[np.ix_(select, select)]
    assert flight[0, 1] == flight[1, 0] == -0.5
    return diag, flight


def pad_energies(q, k, left_pad, right_pad):
    diag, _ = terms(q, k)
    return np.array(
        [diag[q, left_pad] + diag[0, right_pad], diag[0, left_pad] + diag[q, right_pad]]
    )


@cache
def probabilities(q, k, sigma, matches, axes, time, parked):
    """Right-port read after each of eight physical phase schedules."""
    _, kinetic = terms(q, k)
    pads = (2, 0) if sigma < 0 else (0, 2)
    h = kinetic + np.diag(pad_energies(q, k, *pads))
    if parked:
        h = np.diag(h.diagonal())
    propagator = expm(-1j * time * (h - h.trace() * np.eye(2) / 2))
    probs = []
    for mask in MASKS:
        psi = np.ones(2, dtype=complex) / np.sqrt(2)
        for i in axes:
            if mask[i]:
                e = pad_energies(q, k, 2 * matches[i], 2)
                if sigma < 0:
                    e = e[::-1]
                psi *= np.exp(-1j * CLOCK * (e - e[0]))
        probs.append(float(abs((propagator @ psi)[1]) ** 2))
    return tuple(probs)


def response_table(q, k, axes, time, parked):
    return np.array(
        [
            [
                np.mean(probabilities(q, k, sign, word, axes, time, parked))
                for word in product((0, 1), repeat=3)
            ]
            for sign in (-1, 1)
        ]
    )


def exact_coefficient(q, k, time):
    d = (4 * k + 0.5) * q
    return d * (1 - np.cos(np.hypot(1, d) * time)) / (2 * (1 + d * d))
