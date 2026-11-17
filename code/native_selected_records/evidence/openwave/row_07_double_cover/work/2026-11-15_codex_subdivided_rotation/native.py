"""Finite triangle native endpoint and full current instruments; stdlib only."""
from fractions import Fraction as F
from math import log10

CALLS = {"native": 0, "swap": 0}
EDGES = ((0, 1), (1, 2), (0, 2))
I = tuple(tuple(F(i == j) for j in range(3)) for i in range(3))


def endpoint(a, b, coin, law):
    CALLS[law] += 1
    if law == "swap":
        return b, a
    assert law == "native"
    lo = (a + b) // 2
    hi = a + b - lo
    return (lo, hi) if coin == 0 else (hi, lo)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return tuple(dot(row, v) for row in a)


def mm(a, b):
    return tuple(tuple(dot(row, col) for col in zip(*b)) for row in a)


def add(a, b, s=1):
    return tuple(tuple(x + s * y for x, y in zip(r, q)) for r, q in zip(a, b))


def scale(a, s):
    return tuple(tuple(s * x for x in row) for row in a)


def fine(edges, law, vacuum=False):
    out = []
    for a, b in edges:
        for coin in (0, 1):
            flags = [[[0] * 3 for _ in range(3)] for _ in range(3)]
            for x in range(3):
                state = [int(j == x) for j in range(3)]
                if vacuum:
                    endpoint(0, 0, coin, law)
                    current, y = 0, x
                else:
                    old = state[a]
                    state[a], state[b] = endpoint(state[a], state[b], coin, law)
                    current, y = old - state[a], state.index(1)
                flags[current + 1][y][x] = 1
            out.extend(tuple(tuple(row) for row in a) for a in flags)
    return tuple(out)


def eraser(branches):
    prefix, rows = scale(I, 0), []
    count = len(branches) // 3
    for k, a in enumerate(branches):
        if k:
            rows.append((add(prefix, a, -k), count * k * (k + 1)))
        prefix = add(prefix, a)
    return ((prefix, count * len(branches)), *rows)


def routing(slot, subdivisions):
    vertex, k = divmod(slot % (3 * subdivisions), subdivisions)
    other = (vertex + 1) % 3
    return (tuple(e for e in EDGES if vertex in e) * (subdivisions - k)
            + tuple(e for e in EDGES if other in e) * k)


def source(slot, subdivisions, law, full=False):
    branches = fine(routing(slot, subdivisions), law)
    numerator = scale(I, 0)
    for a in branches:
        numerator = add(numerator, a)
    mean = scale(numerator, F(1, 4 * subdivisions))
    return mean, eraser(branches) if full else None


def logprob(p):
    return log10(p.numerator) - log10(p.denominator)
