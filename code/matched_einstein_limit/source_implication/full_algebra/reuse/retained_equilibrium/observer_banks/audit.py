"""Read the archived integer banks without a source simulation, eigensolver or RNG."""

from itertools import combinations, permutations, product
import json
from math import ceil, floor, isqrt
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
iv = mp.iv
IDENTITY = (1, 0, 0, 1)
X = (0, 1, 1, 0)
Y = (0, -1j, 1j, 0)
Z = (1, 0, 0, -1)
PAULIS = (IDENTITY, X, Y, Z)
GRID = 1 << 40
N = 1 << 66
Q = 1 << 36
EDGES = list(combinations(range(4), 2))
ROUTES = list(permutations(range(4), 2))


def marginal_check(archive):
    assert (
        int(archive["trials_per_setting"]) == N and archive["tickets_per_setting"] == Q
    )
    largest = 0
    for endpoint, tickets in zip(archive["counts"], archive["tickets"], strict=True):
        for name, width, settings in [("pauli", 4, 15), ("clock", 12, 12)]:
            assert len(endpoint[name]) == settings
            for raw, prob in zip(endpoint[name], tickets[name], strict=True):
                row = list(map(int, raw))
                assert (
                    len(row) == len(prob) == width and min(row) >= 0 and min(prob) >= 0
                )
                assert sum(row) == N and sum(prob) == Q
                for c, k in zip(row, prob, strict=True):
                    assert k > 0 or c == 0
                    largest = max(largest, abs(c * Q - k * N))
    # Over all 408 bins Hoeffding and a union bound give failure probability <=816 exp(-18).
    assert largest <= 3 * isqrt(N) * Q
    return largest


def head_block(records, b):
    counts = [list(map(int, row)) for row in records]
    weights = [sum(counts[0][2 * b : 2 * b + 2])] + [
        c[2 * b] - c[2 * b + 1] for c in counts
    ]
    panel = [[None] * 4 for _ in range(4)]
    for i, j in product(range(4), repeat=2):
        rr = ii = 0
        for w, (p, q) in zip(weights, product(PAULIS, repeat=2), strict=True):
            z = p[2 * (i // 2) + j // 2] * q[2 * (i % 2) + j % 2]
            rr += int(z.real) * w
            ii += int(z.imag) * w
        panel[i][j] = (rr, ii)
    return panel


def characteristic(matrix, point):
    """Determinant of 4N times the Hermitian Pauli estimate, shifted by a dyadic point."""
    shift = point * (4 * N // GRID)
    det = 0
    imag = 0
    for perm in permutations(range(4)):
        re, im = 1, 0
        for i, j in enumerate(perm):
            a, b = matrix[i][j]
            a -= shift * (i == j)
            re, im = re * a - im * b, im * a + re * b
        sign = (-1) ** sum(i > j for k, i in enumerate(perm) for j in perm[k + 1 :])
        det += sign * re
        imag += sign * im
    assert imag == 0
    return det


def entropy_box(raw, root_boxes, head_error):
    ans = iv.mpf(0)
    for b, points in enumerate(root_boxes):
        matrix = head_block(raw, b)
        assert len(points) == 4 and all(lo > 0 and lo < hi for lo, hi in points)
        assert all(hi < lo for (_, hi), (lo, _) in zip(points, points[1:]))
        for lo, hi in points:
            assert characteristic(matrix, lo) * characteristic(matrix, hi) < 0
            x = iv.mpf([lo, hi]) / GRID
            radius = 2 * head_error
            full = x + iv.mpf([-1, 1]) * radius
            ans += -x * iv.ln(x) + iv.mpf([-1, 1]) * radius * abs(1 + iv.ln(full))
    return ans


def geometry_box(raw, width):
    weights = [*range(2, 13), 12]
    ranges = [
        sum(int(c) * k for c, k in zip(row, weights, strict=True)) * iv.mpf(1) / (4 * N)
        for row in raw
    ]
    powers = [
        [
            int((i, j) == e or (j, i) == e) - (int(i in e) - int(j in e)) / 4
            for i, j in ROUTES
        ]
        for e in EDGES
    ]

    def geometry(rs):
        log_r = [iv.ln(x) for x in rs]
        edges = [
            iv.exp(sum(g * v for g, v in zip(row, log_r, strict=True)))
            for row in powers
        ]
        grad = [iv.mpf(0) for _ in edges]
        area = iv.mpf(0)
        for face in combinations(range(4), 3):
            ids = [EDGES.index(e) for e in combinations(face, 2)]
            a, b, c = [edges[i] for i in ids]
            root = iv.sqrt(2 * (a * b + a * c + b * c) - a * a - b * b - c * c)
            area += root / 4
            for i, side, others in zip(
                ids, (a, b, c), (b + c, a + c, a + b), strict=True
            ):
                grad[i] += side * (others - side) / (4 * root)
        a, b, c, d, e, f = edges
        p, q, r = (a + b - d) / 2, (a + c - e) / 2, (b + c - f) / 2
        det = a * b * c + 2 * p * q * r - a * r * r - b * q * q - c * p * p
        dp, dq, dr = 2 * (q * r - c * p), 2 * (p * r - b * q), 2 * (p * q - a * r)
        jac = [
            b * c - r * r + (dp + dq) / 2,
            a * c - q * q + (dp + dr) / 2,
            a * b - p * p + (dq + dr) / 2,
            -dp / 2,
            -dq / 2,
            -dr / 2,
        ]
        return (area, iv.sqrt(det) / 6), (
            grad,
            [
                der * edge / (12 * iv.sqrt(det))
                for der, edge in zip(jac, edges, strict=True)
            ],
        )

    centers, _ = geometry(ranges)
    whole = [x + iv.mpf([-1, 1]) * width for x in ranges]
    _, grad = geometry(whole)
    out = []
    for val, deriv in zip(centers, grad, strict=True):
        err = 0
        for i, r in enumerate(whole):
            err += (
                abs(sum(row[i] * d for row, d in zip(powers, deriv, strict=True)) / r)
                * width
            )
        out.append(val + iv.mpf([-1, 1]) * err)
    return out


def rounded(value):
    # Decimal interval endpoints in integer millionths, outward rounded.
    scaled = value * 10**6
    return [floor(float(scaled.a)) - 1, ceil(float(scaled.b)) + 1]


def _audit(archive, frames):
    marginal_check(archive)
    # Outward-rounded coverage for the continuum target, including tickets, finite settling,
    # and the explicitly registered numerical allowance (rather than RNG identity).
    tail = iv.sqrt(iv.ln(17600) / (2 * N))
    settle = 2 * (iv.mpf(7) / 8 + iv.mpf("0.00001")) ** 256
    head = 2 * (tail + iv.mpf(4) / Q + settle) + iv.mpf("1e-10")
    clock = iv.mpf("2.5") * (tail + iv.mpf(12) / Q + settle) + iv.mpf("1e-10")
    data = []
    for arm, boxes in zip(archive["counts"], frames, strict=True):
        data.append(
            [entropy_box(arm["pauli"], boxes, head), *geometry_box(arm["clock"], clock)]
        )
    ds, da, dv = [(p - m) * 50000 for m, p in zip(*data, strict=True)]
    eta = iv.mpf(["23.10423617469", "23.10423621780"])
    return {
        k: rounded(v)
        for k, v in zip(
            ("entropy", "area", "volume", "ordinary_residual", "double_gain_control"),
            (ds, da, dv, ds + eta * da, ds + 2 * eta * da),
            strict=True,
        )
    }


def audit(archive, frames):
    previous = iv.dps
    try:
        iv.dps = 45
        return _audit(archive, frames)
    finally:
        iv.dps = previous


def make_frames(data):
    # Approximate eigenvalues propose brackets only. Exact determinants certify the brackets.
    import numpy as np

    out = []
    for arm in data["counts"]:
        block = []
        for b in range(2):
            m = head_block(arm["pauli"], b)
            v = np.array([[complex(a, d) / (4 * N) for a, d in row] for row in m])
            block.append(
                [
                    [floor(x * GRID) - 1, ceil(x * GRID) + 1]
                    for x in np.linalg.eigvalsh(v)
                ]
            )
        out.append(block)
    return out


if __name__ == "__main__":
    data = json.loads((HERE / "records.json").read_text())
    archive = json.loads((HERE / "certificate.json").read_text())
    assert (
        audit(data, archive["dyadic_root_boxes"])
        == archive["secant_intervals_millionths"]
    )
    print(
        "Verified the integer terminal counts, eigenvalue signs and simultaneous response intervals."
    )
