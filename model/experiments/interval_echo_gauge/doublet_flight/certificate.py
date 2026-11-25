"""Certify analytic rank two without floating point, using an exact finite-field witness."""

from math import factorial

import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    pair_index,
)
from model.candidate_laws.interval_echo import currents
from model.executable_model import cell

from .measure import fixtures

PRIME = 101
ROOT = 69
ORDER = 64


def modes(axis):
    car = cell()
    sector = PairSector(car)
    anti = currents._geometry(car)[0]
    turn = next(p for p in currents._one_edge_turns(car) if p[axis] == axis)
    start = min(car.neighbours[axis])
    orbit = [start, turn[start], turn[turn[start]]]
    stars = np.zeros((sector.P, 3), dtype=np.int64)
    for i, port in enumerate(orbit):
        for aligned, sign in ((port, 1), (anti[port], -1)):
            for neighbor in car.neighbours[aligned]:
                stars[sector.index(aligned, neighbor), i] = sign
    raw = np.diff(stars)
    k = (2 * sector.hamiltonian(currents.KAPPA)).astype(np.int64)
    f = (1, 34, -128, -6904, 41936)
    # F vanishes on this cyclic source, as checked over the integers.
    remainder = raw.copy()
    for c in f[1:]:
        remainder = k @ remainder + c * raw
    assert np.all(remainder == 0)
    # Q(K;r) for the monic quotient Q(X;r) = (F(X)-F(r))/(X-r).
    vec = raw.copy()
    coefficient = 1
    for c in f[1:-1]:
        coefficient = (coefficient * ROOT + c) % PRIME
        vec = (k @ vec + coefficient * raw) % PRIME
    return vec


def witness(board, destination, p, q, kappa):
    sector = PairSector(board)
    small = PairSector(cell())
    source_index = pair_index(board.n, small.aa, small.bb)
    receiver_index = pair_index(
        board.n, 12 * destination + small.aa, 12 * destination + small.bb
    )
    v = np.zeros((sector.P, 2), dtype=np.int64)
    src = modes(p)
    v[source_index] = src
    w = modes(q)
    gs = src.T @ src % PRIME
    gd = w.T @ w % PRIME
    adj = np.array([[gd[1, 1], -gd[0, 1]], [-gd[1, 0], gd[0, 0]]])
    assert (gs[0, 0] * gs[1, 1] - gs[0, 1] ** 2) % PRIME
    assert (gd[0, 0] * gd[1, 1] - gd[0, 1] ** 2) % PRIME
    k = (2 * sector.hamiltonian(kappa)).astype(np.int64)
    moments = []
    for n in range(ORDER + 1):
        den = pow(factorial(n) % PRIME, -1, PRIME)
        moments.append((w.T @ v[receiver_index]) * den % PRIME)
        v = k @ v % PRIME
    determinant = []
    unequal_axes = []
    for n in range(ORDER + 1):
        total = sum(
            a[0, 0] * b[1, 1] - a[0, 1] * b[1, 0]
            for a, b in zip(moments[: n + 1], reversed(moments[: n + 1]), strict=True)
        )
        if total % PRIME:
            determinant.append((n, int(total % PRIME)))
        gram = sum(
            (-1) ** j * moments[j].T @ adj @ moments[n - j] for j in range(n + 1)
        )
        unequal = (gram[0, 0] * gs[1, 1] - gram[1, 1] * gs[0, 0]) % PRIME
        if unequal:
            unequal_axes.append((n, int(unequal)))
    return determinant, unequal_axes


def run():
    x = sp.Symbol("t")
    f = sp.Poly(x**4 + 34 * x**3 - 128 * x**2 - 6904 * x + 41936, x, domain=sp.ZZ)
    assert f.is_irreducible and f.eval(ROOT) % PRIME == 0
    assert f.diff().eval(ROOT) % PRIME != 0
    setups, p, q = fixtures()
    for board, destination in setups:
        for kappa in (currents.KAPPA, 0):
            data = witness(board, destination, p, q, kappa)
            print(board.name, kappa, [part[:2] for part in data])


if __name__ == "__main__":
    run()
