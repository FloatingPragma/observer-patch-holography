"""Dyadic enclosure of the two-clock score for every static loading profile."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from . import port_clocks as clocks

BASE = 2**24
GRID = BASE**2
ORDER = 16
SCORE = (8357, -91643)


def square(real, imag):
    # Exact signed integer products in two limbs. Rounding is applied only
    # after summing both real products that feed each output component.
    rh, rl = np.divmod(real, BASE)
    ih, il = np.divmod(imag, BASE)
    size = len(real)
    limb = max(BASE, int(abs(rh).max()), int(abs(ih).max()))
    assert 4 * size * limb**2 + 2 * GRID < 2**63

    def finish(high, cross, low):
        q, rem = np.divmod(cross, BASE)
        return high + q + (rem * BASE + low + GRID // 2) // GRID

    return (
        finish(
            rh @ rh - ih @ ih, rh @ rl + rl @ rh - ih @ il - il @ ih, rl @ rl - il @ il
        ),
        finish(
            rh @ ih + ih @ rh, rh @ il + rl @ ih + ih @ rl + il @ rh, rl @ il + il @ rl
        ),
    )


def block_exponential(h4, center):
    n = len(h4)
    # At time twelve, tH=3(4H). Remove the same center for every parity.
    h = sparse.csr_matrix(3 * (h4 - center * np.eye(n, dtype=np.int64)))
    norm = int(abs(h).sum(axis=1).max())
    depth = max(0, (2 * norm - 1).bit_length()) if norm else 0
    steps = 2**depth
    x = F(norm, steps)
    assert x <= F(1, 2)
    real = GRID * np.eye(n, dtype=np.int64)
    imag = np.zeros_like(real)
    tr, ti = real.copy(), imag.copy()
    for k in range(1, ORDER + 1):
        den = k * steps
        assert norm * max(int(abs(tr).max()), int(abs(ti).max())) + den < 2**63
        tr, ti = (h @ ti + den // 2) // den, (-(h @ tr) + den // 2) // den
        assert (
            max(int(abs(real).max()), int(abs(imag).max()))
            + max(int(abs(tr).max()), int(abs(ti).max()))
            < 2**63
        )
        real += tr
        imag += ti
    # Each rounded complex matrix has operator error <= Frobenius <= n/GRID.
    # Taylor term error <= (n/GRID) sum_j (1/2)^j <= 2n/GRID.
    rounding = F(n, GRID)
    error = 2 * ORDER * rounding + 2 * x ** (ORDER + 1) / factorial(ORDER + 1)
    for _ in range(depth):
        real, imag = square(real, imag)
        # For U unitary, ||(U+D)^2-U^2|| <= 2e+e². Round the
        # scalar bound upward to the grid to keep its denominator small.
        up = (2 * error + error**2 + rounding) * GRID
        error = F(-(-up.numerator // up.denominator), GRID)
    return real, imag, error


def sign_banks(n):
    model = clocks.joint.source.timing.model
    mats = [model(n, p)[0] for p in range(8)]
    center = (
        min(a.diagonal().min() for a in mats) + max(a.diagonal().max() for a in mats)
    ) // 2
    banks = []
    worst = F(0)
    for a in mats:
        count, labels = connected_components(sparse.csr_matrix(a), directed=False)
        real, imag = np.zeros_like(a), np.zeros_like(a)
        for j in range(count):
            ids = np.flatnonzero(labels == j)
            ix = np.ix_(ids, ids)
            real[ix], imag[ix], err = block_exponential(a[ix], int(center))
            worst = max(worst, err)
        banks.append((real, imag))
    # The direct sum norm is the maximum over invariant blocks. The four
    # sign-sector transforms have coefficients ±1/4 and hence error <= worst.
    return banks, worst


def score_populations(n):
    banks, error = sign_banks(n)
    board, _, _, _ = clocks.joint.source.matched.fp.make_board()
    starts = clocks.joint.source.timing.model(n, 0)[1]
    # Source density numerators in coin,endpoint order, with denominator 8.
    paired = np.diag([0, 4, 4, 0])
    f = paired.copy()
    f[1, 2] = f[2, 1] = -2
    g = np.diag([3, 1, 1, 3])
    g[0, 3] = g[3, 0] = -3
    g[1, 2] = g[2, 1] = 1
    scores = {
        k: F(0) for k in (["empty"] if not n else [*range(board.n), "paired", "f", "g"])
    }
    for source_sign in range(2):
        for near, weight in enumerate(SCORE):
            # The far sign is always + in this score.
            r, s = [
                sum(
                    (-1) ** (i * near) * banks[source_sign + 2 * i + 4 * j][k]
                    for i, j in np.ndindex(2, 2)
                ).astype(object)
                for k in (0, 1)
            ]
            diag = np.sum(r * r + s * s, axis=0)
            factor = F(weight, 100000 * 2 * (4 * GRID) ** 2)
            if not n:
                scores["empty"] += factor * F(sum(diag), 8)
                continue
            for p, row in enumerate(diag.reshape(board.n, 8)):
                scores[p] += factor * F(sum(row), 8)
            r, s = r[:, starts], s[:, starts]
            gram = (r.T @ r + s.T @ s).reshape(4, 4, 4, 4)
            eff = np.trace(gram, axis1=1, axis2=3)
            for name, rho in (("paired", paired), ("f", f), ("g", g)):
                scores[name] += factor * F(int(np.sum(rho * eff)), 32)
    return scores, 2 * error + error**2


def run():
    vac, e0 = score_populations(0)
    one, e1 = score_populations(1)
    scale = 10**10

    def interval(value, err):
        lo, hi = scale * (value - err), scale * (value + err)
        return [lo.numerator // lo.denominator, -(-hi.numerator // hi.denominator)]

    endpoint = {k: interval((one[k] - one["paired"]) / 2, e1) for k in ("f", "g")}
    ports = {
        str(k): interval((v - vac["empty"]) / 2, (e0 + e1) / 2)
        for k, v in one.items()
        if isinstance(k, int)
    }
    return dict(
        protocol_sha256=sha256(
            (clocks.HERE / "PORT_ENCLOSURE.md").read_bytes()
        ).hexdigest(),
        exposure=12,
        coefficients=list(SCORE),
        denominator=100000,
        grid=GRID,
        interval_denominator=scale,
        endpoints=endpoint,
        ports=ports,
        source_lower=min(pair[0] for pair in endpoint.values()),
        maximum_port_upper=max(pair[1] for pair in ports.values()),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = run()
    path = clocks.HERE / "port_enclosure.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        assert out == json.loads(path.read_bytes())
    print(out["source_lower"], out["maximum_port_upper"], out["interval_denominator"])
