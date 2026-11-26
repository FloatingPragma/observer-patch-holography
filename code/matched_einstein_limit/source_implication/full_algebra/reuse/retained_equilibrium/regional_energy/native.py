"""Enclose cell-energy effects on the original two-head input after four native ticks."""

import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

import numpy as np

from .....directional_stress import bridge
from .....directional_stress.translation import measure as region

HERE = Path(__file__).parent
U, C, Q, EFFECT = 2**55, 2**46, 2**43, 2**42
K = 180


def scaled_product(x, c):
    # floor(x*c/2^46) without a wide multiply. |x|<2^56 and |c|<3*2^46.
    low = (1 << 23) - 1
    hi, lo = x >> 23, x & low
    cj, cl = c >> 23, c & low
    return hi * cj + ((hi * cl + lo * cj + ((lo * cl) >> 23)) >> 23)


def flight(h, starts, weights):
    row = int(abs(h).sum(axis=1).max())
    assert row < 64
    word = np.zeros((h.shape[0], len(starts)), dtype=np.int64)
    word[starts] = U * np.eye(len(starts), dtype=np.int64)
    prior = np.zeros_like(word)
    output = np.zeros((2, *word.shape), dtype=np.int64)
    eps = F(0)
    root = isqrt(h.shape[0]) + 1
    for k in range(K + 1):
        assert int(np.max(abs(word))) < 2**56
        c = int(weights[k]) * C // bridge.native.prior.FRAC
        assert abs(c) < 3 * C
        sign = 1 if k % 4 in (0, 3) else -1
        output[k % 2] += sign * scaled_product(word, c)
        # Contractive scaled Hamiltonian; Chebyshev U_j has operator norm <=j+1.
        # Each downward rounded state recurrence injects norm <=root/U.
        eps += F(abs(c) * k * (k + 1) * root, 2 * C * U) + F(1, C) + F(root, U)
        if k < K:
            word, prior = (h @ word >> (6 if k == 0 else 5)) - prior, word
    # Integer safety, independently of cancellation in the accumulator.
    assert (
        sum(abs(int(x) * C // bridge.native.prior.FRAC) for x in weights[: K + 1])
        < 15 * C
    )
    eps += F(sum(map(abs, weights[K + 1 :])), bridge.native.prior.FRAC)
    eps += F(len(weights), bridge.native.prior.FRAC) + F(1, 10**59) + F(2 * root, Q)
    assert eps < F("0.0000000003")
    assert np.max(abs(output)) < 2**56
    return (output // (U // Q)).reshape(2, -1, 4, len(starts))


def product(x, h):
    # Ten launches are mixed independently. Never interfere distinct launches.
    a, b = x
    ha, hb = h @ a, h @ b
    out = np.zeros((2, 4, 4), dtype=object)
    for start in range(0, a.shape[1], 4):
        s = slice(start, start + 4)
        ar, ai = a[:, s].astype(object), b[:, s].astype(object)
        out[0] += ar.T @ ha[:, s] + ai.T @ hb[:, s]
        out[1] += ar.T @ hb[:, s] - ai.T @ ha[:, s]
    return out


def cells_in_charge(board, route, flight_state, q, odd):
    # Inputs to product are scaled integer columns. Local matrices represent
    # four times mobile energy: twice the twice-internal, and half each bridge.
    blocks = bridge.operators(board, route, q, odd)
    bridge4 = np.kron(np.array([[17, -1], [-1, 17]]), np.eye(4, dtype=int))
    out = []
    for c in range(board.cells):
        x = flight_state[:, 12 * c : 12 * c + 12].reshape(2, 48, -1)
        assert int(abs(2 * blocks[c]).sum(axis=1).max()) + 3 * 18 < 4 * 200
        out.append(product(x, 2 * blocks[c]))
    for u, v in board.external:
        b = product(flight_state[:, [u, v]].reshape(2, 8, -1), bridge4)
        out[u // 12] += b
        out[v // 12] += b
    rows, error = [], []
    for i, mat in enumerate(out):
        values = np.empty((2, 4, 4), dtype=np.int64)
        for a, x, y in np.ndindex(2, 4, 4):
            values[a, x, y] = round(F(int(mat[a, x, y]) * EFFECT, 40 * Q**2))
        assert np.array_equal(values[0].T, values[0]) and np.array_equal(
            values[1].T, -values[1]
        )
        rows.append(values.tolist())
        supported = set(range(12 * i, 12 * i + 12))
        for a, b in board.external:
            if a // 12 == i or b // 12 == i:
                supported.update((a, b))
        # At most 15 supporting ports. A four-column block has Frobenius
        # norm <= sqrt(2*15*4*4) * max(grid amplitude)/Q (real + imaginary).
        assert len(supported) <= 15
        d = F("0.0000000006")  # <= 2 times any single column error
        amp = min(
            1 + d, F(22 * int(np.max(abs(flight_state[:, sorted(supported)]))), Q)
        )
        err = 200 * (2 * amp * d + d * d) + F(4, EFFECT)
        error.append(-((-err * 10**12) // 1))
    return rows, error


def build():
    atlas = region.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    starts = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    weights = bridge.native.prior.weights()[7]
    rows, errors = {}, {}
    for q in (-1, 1):
        by_sector, bounds = [], []
        for odd in (False, True):
            h = bridge.native.generator(board, route, q, odd)
            w = flight(h, starts, weights)
            calculated, error = cells_in_charge(board, route, w, q, odd)
            by_sector.append(calculated)
            bounds.append(error)
        rows[str(q)] = by_sector
        errors[str(q)] = [max(row) for row in zip(*bounds, strict=True)]
    return {
        "scale": EFFECT,
        "cell_error_trillionths": errors,
        "marks": [route[0] // 12, route[2] // 12],
        "chart": region.chosen(atlas)["cells"],
        "rows": rows,
    }


if __name__ == "__main__":
    data = build()
    (HERE / "effects.json").write_text(json.dumps(data, separators=(",", ":")) + "\n")
