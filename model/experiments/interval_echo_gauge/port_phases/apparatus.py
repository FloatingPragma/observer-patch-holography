"""Native ternary writer, rail routes and complete terminal measurements."""

from collections import deque
from functools import cache
from itertools import permutations
from math import pi

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

LAW = candidate.native_law
K = candidate.KAPPA
SERVICE, GUARD, OUTSIDE = 13, 19, 24
SIZE = 14  # twelve data rails, external reference, guard
ANGLES = tuple(2 * pi * ((7 * p + 3) % 29) / 29 for p in range(12))


@cache
def native_contact(k):
    x, c = LAW.local_term_parts()
    h = x + k * c
    expected = LAW.seam_term(k)
    np.testing.assert_allclose(
        h[:, :9], np.vstack([expected, np.zeros((90, 9))]), atol=4e-14, rtol=0
    )
    return expected


@cache
def contact(k, duration):
    return expm(-1j * duration * native_contact(k))


def reference(k, duration):
    return np.diag(np.exp(-1j * k * duration * np.repeat([0.0, 1.0, 4.0], 3)))


@cache
def writer(angle, menu):
    k = 0 if menu == "exchange_only" else K
    if menu == "square":
        return reference(k, angle / K)
    plus = reference(k, (angle / 2) / K)
    parity = reference(k, pi / K)
    inverse_angle = (-angle / 2) % (2 * pi)
    half = contact(k, inverse_angle / (2 * (K + 0.5)))
    return parity @ half @ parity @ half @ plus


def path(board, start, end, avoid):
    queue = deque([start])
    parent = {start: None}
    while queue:
        node = queue.popleft()
        if node == end:
            steps = []
            while parent[node] is not None:
                steps.append((parent[node], node))
                node = parent[node]
            return steps[::-1]
        for nxt in board.neighbours[node]:
            if nxt not in avoid and nxt not in parent:
                parent[nxt] = node
                queue.append(nxt)
    raise ValueError("No route between the declared ports")


def swap_word(labels, edges):
    for u, v in edges:
        labels[u], labels[v] = labels[v], labels[u]


def route_pair(board, starts, ends):
    labels = list(range(board.n))
    edges = []
    fixed = {GUARD}
    for start, end in zip(starts, ends, strict=True):
        steps = path(board, labels.index(start), end, fixed)
        swap_word(labels, steps)
        edges.extend(steps)
        fixed.add(end)
    if tuple(labels[x] for x in ends) != starts or labels[GUARD] != GUARD:
        raise ValueError("Bad two-rail transport")
    return edges


def stations(board):
    def pair(owner):
        return next(
            (a, b)
            for a, b in board.seams
            if a // 12 == b // 12 == owner and board.degree[a] == board.degree[b] == 5
        )

    return pair(0), pair(13)


def layout():
    board = federation(0)
    if GUARD not in board.neighbours[SERVICE] or board.degree[SERVICE] != 5:
        raise ValueError("Inaccessible nonlinear reference")
    writers = []
    for p in range(12):
        edges = path(board, p, SERVICE, {GUARD, OUTSIDE})
        labels = list(range(board.n))
        swap_word(labels, edges)
        if labels[SERVICE] != p:
            raise ValueError("Writer reached wrong rail")
        swap_word(labels, edges[::-1])
        if labels != list(range(board.n)):
            raise ValueError("Writer left displaced ports")
        writers.append(edges)
    benches = stations(board)
    deliveries = []
    for source, receiver in (benches, benches[::-1]):
        arms = []
        for p in range(12):
            before = route_pair(board, source, (p, OUTSIDE))
            after = route_pair(board, (p, OUTSIDE), receiver)
            labels = list(range(board.n))
            swap_word(labels, after)
            varied = {*range(12), OUTSIDE, GUARD}
            read = []
            protected = {*receiver, GUARD}
            for target in (*receiver, GUARD):
                neighbor = next(
                    q for q in board.neighbours[target] if q not in protected
                )
                # Fetch a certified unit even if the changed-law writer disturbed all
                # twelve data ports. Undo each witness route before the next reading.
                if labels[neighbor] in varied:
                    stock = next(
                        q
                        for q in range(board.n)
                        if labels[q] not in varied and q not in protected
                    )
                    fetch = path(board, stock, neighbor, protected)
                else:
                    fetch = []
                installed = labels.copy()
                swap_word(installed, fetch)
                if installed[neighbor] in varied or any(
                    installed[q] != labels[q] for q in protected
                ):
                    raise ValueError("Missing unit witness")
                read.append(dict(target=target, witness=neighbor, fetch=fetch))
            arms.append(dict(before=before, after=after, read=read))
        deliveries.append(dict(source=source, receiver=receiver, arms=arms))
    return dict(
        level=0,
        writer_site=SERVICE,
        guard=GUARD,
        reference=OUTSIDE,
        writer_routes=writers,
        stations=deliveries,
    )


class Chamber:
    """Complete multiset of a ternary pair among unit spectators, including the guard."""

    def __init__(self, pair):
        exc = [x for x in pair if x != 1]
        words = []
        for pos in permutations(range(SIZE), len(exc)):
            word = [1] * SIZE
            for j, value in zip(pos, exc, strict=True):
                word[j] = value
            words.append(tuple(word))
        self.index = {w: j for j, w in enumerate(words)}
        self.words = np.array(words, dtype=np.int8)
        self.pair = pair

    @cache
    def support(self, i, j):
        flip = self.words.copy()
        flip[:, [i, j]] = flip[:, [j, i]]
        return (
            3 * self.words[:, i] + self.words[:, j],
            3 * self.words[:, j] + self.words[:, i],
            np.array([self.index[tuple(w)] for w in flip]),
        )

    def apply(self, vector, gate, i, j):
        xy, yx, dest = self.support(i, j)
        return gate[xy, xy] * vector + (xy != yx) * gate[xy, yx] * vector[dest]

    def receive(self, p, circuit, menu, quadrature):
        # Physical pi-swaps in layout transport exactly these labeled rails.
        reference_slot, guard_slot = 12, 13
        start = [1] * SIZE
        a, b = self.pair
        start[p], start[reference_slot] = a, b
        vector = np.zeros(len(self.words), complex)
        vector[self.index[tuple(start)]] = 1
        k = 0 if menu == "exchange_only" else K
        vector = self.apply(vector, contact(k, pi / 2), p, reference_slot)
        for site, gate in enumerate(circuit):
            vector = self.apply(vector, gate, site, guard_slot)
        dwell = quadrature * pi / (2 * K * (b * b - a * a))
        vector = np.exp(-1j * k * dwell * self.words[:, p] ** 2) * vector
        vector = self.apply(vector, contact(k, 3 * pi / 2), p, reference_slot)
        # Each remote contrast head compared to a unit witness reads n-1.
        bins = (
            9 * self.words[:, p]
            + 3 * self.words[:, reference_slot]
            + self.words[:, guard_slot]
        )
        return np.bincount(bins, weights=abs(vector) ** 2, minlength=27)


@cache
def chamber(a, b):
    return Chamber((a, b))


def free_records(theta, t, exchange):
    board = federation(0)
    # One vacancy above the unit-filled reference; remove its common energy.
    h = native_contact(K)
    diag_pair = h[LAW.SRC.index((0, 1)), LAW.SRC.index((0, 1))]
    move = h[LAW.SRC.index((0, 1)), LAW.SRC.index((1, 0))]
    diagonal = diag_pair * board.degree - K * (6 - board.degree)
    generator = diags(diagonal) + (move * board.adjacency() if exchange else 0)
    basis = np.zeros(board.n, complex)
    basis[:2] = contact(K, pi / 2)[[1, 3], 1]
    # The background's phase is common; the vacancy at 0 lacks one unit.
    gate = writer(theta, "native")
    basis[0] *= gate[1, 1]
    basis[1] *= gate[4, 4]
    out = expm_multiply(-1j * t * generator, basis)
    # One contrast reading on the physical 0--1 edge separates all three bins.
    return np.array(
        [abs(out[0]) ** 2, abs(out[1]) ** 2, np.vdot(out[2:], out[2:]).real]
    )
