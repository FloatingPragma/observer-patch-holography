"""Two-coin rotations from a coherent mobile rail and the actual right head contact."""

from collections import deque
from itertools import pairwise
from math import pi

import numpy as np
from scipy.linalg import expm

from ....full_geometry.station_balls.reciprocal_screen.live_reference import (
    apparatus as rail,
)

ID2, X, Y, Z = rail.IDENTITY, rail.X, rail.Y, rail.Z
H = (X + Z) / np.sqrt(2)
S_INV = np.diag([1, -1j])
ORDER = [0, 1, 3, 2]


def shuttle():
    board = rail.screen.camera.front.optics.pairing.formation.federation(1)
    a, b, c, d = rail.screen.camera.front.optics.pairing.formation.region(board)
    forbidden = {frozenset((a, b)), frozenset((c, d))}

    def path(source, target, block):
        parents = {source: None}
        queue = deque([source])
        while target not in parents:
            x = queue.popleft()
            for y in sorted(board.neighbours[x]):
                if (
                    y != block
                    and y not in parents
                    and frozenset((x, y)) not in forbidden
                ):
                    parents[y] = x
                    queue.append(y)
        steps = [target]
        while steps[-1] != source:
            steps.append(parents[steps[-1]])
        return steps[::-1]

    paths = [path(a, c, b), path(b, d, c)]
    # These are permutations of physical ports on one excess rail; the other
    # logical rail position is protected even when its amplitude is unknown.
    positions = [a, b]
    for steps in paths:
        for i, j in pairwise(steps):
            positions = [j if x == i else i if x == j else x for x in positions]
    assert positions == [c, d]
    return paths


def native_cz():
    g = rail.contact(False)
    half = expm(-1j * pi / 16 * g)
    pulse = np.kron(Z, ID2)
    both = np.kron(S_INV, S_INV) @ pulse @ half @ pulse @ half
    return both / both[0, 0]


def rot(axis, angle):
    return np.cos(angle / 2) * ID2 - 1j * np.sin(angle / 2) * axis


def controlled(u, control, bit, cz):
    target = 1 - control

    def one(site, v):
        return np.kron(v, ID2) if site == 0 else np.kron(ID2, v)

    cx = one(target, H) @ cz @ one(target, H)
    alpha = np.angle(np.linalg.det(u)) / 2
    v = np.exp(-1j * alpha) * u
    s, p = np.angle(v[1, 0]), np.angle(v[0, 0])
    beta, delta = s - p, -s - p
    gamma = 2 * np.arctan2(abs(v[1, 0]), abs(v[0, 0]))
    a = rot(Z, beta) @ rot(Y, gamma / 2)
    b = rot(Y, -gamma / 2) @ rot(Z, -(delta + beta) / 2)
    c = rot(Z, (delta - beta) / 2)
    gate = (
        one(control, np.diag([1, np.exp(1j * alpha)]))
        @ one(target, a)
        @ cx
        @ one(target, b)
        @ cx
        @ one(target, c)
    )
    if bit == 0:
        gate = one(control, X) @ gate @ one(control, X)
    return gate


def two_level(u, i, j, cz):
    key = i ^ j
    assert key in (1, 2)
    control = 0 if key == 1 else 1
    bit = (i >> (1 - control)) & 1
    if i & key:
        u = X @ u @ X
    return controlled(u, control, bit, cz)


def compile_gate(u, cz):
    """Gray-ordered QR uses only two-level rotations with one common bit."""
    m = u.copy()
    factors = []
    for pos, col in enumerate(ORDER[:-1]):
        for k in range(3, pos, -1):
            a, b = ORDER[k - 1 : k + 1]
            x, y = m[[a, b], col]
            n = np.hypot(abs(x), abs(y))
            if n < 1e-14:
                continue
            g = np.array([[x.conjugate(), y.conjugate()], [-y, x]]) / n
            m[[a, b]] = g @ m[[a, b]]
            factors.append((g.conj().T, a, b))
    total = two_level(np.diag(m.diagonal()[:2]), 0, 1, cz) @ two_level(
        np.diag(m.diagonal()[2:]), 2, 3, cz
    )
    for g, a, b in reversed(factors):
        total = two_level(g, a, b, cz) @ total
    return total
