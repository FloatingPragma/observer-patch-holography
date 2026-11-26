"""A signed vacancy coupled by two native triangles to distinct receiving observers."""

from itertools import product
import numpy as np

from model.executable_model import federation
from ...row_09_neutral.passive_decay.material import edge, commission
from ...row_09_neutral.passive_decay.measure import physical

X = np.array([[0.0, 1.0], [1.0, 0.0]])
Z = np.diag([1.0, -1.0])
IDENTITY = np.eye(2)
SOURCE = np.diag([0.0, 2.0]) - X / 2


def triangle(board, port):
    nearby = sorted(q for q in board.neighbours[port] if q // 12 == port // 12)
    for x, y in product(nearby, repeat=2):
        if y in board.neighbours[x] and board.degree[x] == board.degree[y]:
            yield x, y


def mount():
    for level in range(3):
        board = federation(level)
        for p, q in sorted(edge(*s) for s in board.seams):
            if p // 12 == q // 12 or board.degree[p] != board.degree[q]:
                continue
            for b, c in product(triangle(board, p), triangle(board, q)):
                dimers = ((p, q), b, c)
                sites = set(sum(dimers, ()))
                if len(sites) != 6:
                    continue
                heads = {
                    edge(p, b[0]): 33,
                    edge(p, b[1]): 32,
                    edge(q, c[0]): 33,
                    edge(q, c[1]): 32,
                }
                for r, cost in [(q, 33), (p, 32), (q, 33), (p, 32)]:
                    choices = sorted(
                        s
                        for s in board.neighbours[r]
                        if s not in sites and edge(r, s) not in heads
                    )
                    if not choices:
                        break
                    heads[edge(r, choices[0])] = cost
                else:
                    try:
                        stock = commission(board, dimers, heads)
                    except ValueError:
                        continue
                    return board, dimers, heads, stock
    raise ValueError("No compatible fork on the three federations")


def generator(board, dimers, heads, cuts=(), changed=False):
    heads = dict(heads)
    for arm in cuts:
        heads[edge(dimers[0][arm - 1], dimers[arm][0])] = 32
    if changed:
        heads = {key: 32 for key in heads}
    h = np.diag(physical(board, dimers, heads)).astype(float)
    words = np.arange(8)
    for bit in (1, 2, 4):
        h[words, words ^ bit] -= 0.5
    a, b, c = ((words >> k) & 1 for k in (2, 1, 0))
    expected = (
        np.zeros(8)
        if changed
        else 2 * a - (1 not in cuts) * a * b - (2 not in cuts) * (1 - a) * c
    )
    np.testing.assert_array_equal(np.diag(h) - h[0, 0], expected)
    return h - h[0, 0] * np.eye(8)


def formation(board, dimers, heads):
    """At both weak cuts, subtract a source filled to 11 on the same material."""
    heads = dict(heads)
    for arm in (1, 2):
        heads[edge(dimers[0][arm - 1], dimers[arm][0])] = 32
    # Turn both receiver sharp primary bits off for the fixed energy reference.
    rails = np.ones(board.n, dtype=int)
    for first, _ in dimers[1:]:
        rails[first] = 0
    filled = float((8 * rails**2) @ (6 - board.degree))
    filled += sum(cost * (1 - rails[p] * rails[q]) for (p, q), cost in heads.items())
    # Source unequal blank pair contributes 8 and the diagonal 1/2 of its
    # exchange projector. It contributed zero on 11. Receiver blank terms cancel.
    delta = physical(board, dimers, heads)[::4] - filled + 8.5
    return np.diag(delta) - 0.5 * X
