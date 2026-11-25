"""Local operator certificates and actual transport on blank tree contacts."""

from collections import deque
from itertools import pairwise

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law as law
from model.experiments.interval_echo_consensus.asynchronous import verify_native_service


def certify():
    verify_native_service()
    w = law._write_pattern()
    rev, contrast = law._label_ops()
    identity = np.eye(90, dtype=int)
    source = np.array(law.SRC)
    for g, native in zip(
        (np.eye(10, dtype=int) - rev, 2 * contrast), law.local_term_parts(), strict=True
    ):
        g = np.kron(g, np.eye(9, dtype=int))
        complement = 2 * identity - w @ w.T
        np.testing.assert_array_equal(g @ complement, complement @ g)
        four_h = np.zeros((99, 99), dtype=int)
        four_h[:9, :9] = w.T @ g @ w
        four_h[9:, 9:] = g @ complement
        np.testing.assert_allclose(four_h, 4 * native, atol=8e-14)
        for statistic in (
            np.tile((source == 0).sum(axis=1), 11),
            np.tile((source == 2).sum(axis=1), 11),
            np.repeat(np.arange(11) != 0, 9).astype(int),
        ):
            np.testing.assert_array_equal(
                four_h * statistic, statistic[:, None] * four_h
            )
    swap = np.eye(9)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    blank = law.seam_term(KAPPA)
    np.testing.assert_allclose(expm(-1j * np.pi * blank), swap, atol=4e-14)
    # The subtraction here preserves the contrast Hamiltonian and all nine load words.
    disabled = np.diag(KAPPA * (source[:, 0] - source[:, 1]) ** 2)
    np.testing.assert_allclose(expm(-1j * np.pi * disabled), np.eye(9), atol=4e-14)
    return {
        "local_dimension": 99,
        "free_generator_parts": 2,
        "commuting_census_projectors": 3,
        "routing_angle_in_pi": 1,
    }


class Routing:
    def __init__(self, board):
        self.n = board.n
        self.parent = {0: None}
        frontier = deque([0])
        while frontier:
            a = frontier.popleft()
            for b in sorted(board.neighbours[a]):
                if b not in self.parent:
                    self.parent[b] = a
                    frontier.append(b)
        tree = {tuple(sorted((v, p))) for v, p in self.parent.items() if p is not None}
        self.chords = sorted(set(board.seams) - tree)
        if len(tree) != self.n - 1 or len(self.chords) < self.n // 2:
            raise ValueError("Insufficient routing or write capacity")

    def ancestry(self, v):
        path = [v]
        while self.parent[v] is not None:
            v = self.parent[v]
            path.append(v)
        return path

    def transpose(self, a, b):
        if a == b:
            return []
        left, right = self.ancestry(a), self.ancestry(b)
        while left and right and left[-1] == right[-1]:
            lca = left.pop()
            right.pop()
        path = left + [lca] + right[::-1]
        moves = list(pairwise(path))
        return moves + moves[:-1][::-1]

    def legs(self, pair, chord):
        i, j = pair
        a, b = chord
        if j == a:
            j = i
        return self.transpose(i, a) + self.transpose(j, b)


def attempt(loads, routing, pair, chord, disabled=False):
    route = routing.legs(pair, chord)
    if not disabled:
        for i, j in route:
            loads[i], loads[j] = loads[j], loads[i]
    a, b = chord
    written = abs(loads[a] - loads[b]) == 2
    if written:
        loads[a] = loads[b] = 1
    if not disabled:
        for i, j in reversed(route):
            loads[i], loads[j] = loads[j], loads[i]
    return bool(written), 2 * len(route)
