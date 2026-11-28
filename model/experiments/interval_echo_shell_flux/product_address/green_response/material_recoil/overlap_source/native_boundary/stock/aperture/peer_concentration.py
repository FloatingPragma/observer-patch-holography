"""A finite type-count certificate for material released by a single agreeing peer."""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from sympy import Matrix

from . import peer_fibers as fibers


def logarithm(q):
    # Decimal.ln is correctly rounded. For these positive integers each
    # individual logarithm is below 1000, hence its absolute rounding
    # error is below 1e-55 at 60 significant digits.
    assert 0 < q.numerator < 10**100 and 0 < q.denominator < 10**100
    with localcontext() as ctx:
        ctx.prec = 60
        return Q(Decimal(q.numerator).ln()) - Q(Decimal(q.denominator).ln())


def witness(B, reference, hold, membership, product, size):
    r = np.array(reference, float)
    h = np.array(hold, float)
    u = np.array(product, float)
    t = np.zeros(len(B))
    for _ in range(30):
        v = r + t @ B
        table = h + membership @ v
        grad = np.log(v / r) + membership.T @ np.log(table / u)
        H = (B / v) @ B.T + ((B @ membership.T) / table) @ (membership @ B.T)
        step = np.linalg.solve(H, B @ grad)
        if max(abs(step @ B)) < 1e-16:
            break
        while min(v - step @ B) <= 0:
            step *= 0.5
        t -= step
    counts = (
        np.array([int(x * size) for x in reference])
        + np.rint(size * t).astype(np.int64) @ B
    )
    assert min(counts) > 0
    return np.array([Q(int(n), size) for n in counts], object)


def run():
    words, cells, _, cat, types, n, protected, _, B = fibers.apparatus()
    size = 56 * 10**12
    original = fibers.source.records()[1]["original"]
    ref = np.array([Q(int(x), len(words)) for x in n], object)
    hold = np.array([Q(int(x), len(words)) for x in protected], object)
    codes = np.array([t[0] for t in types])
    L = np.equal.outer(np.arange(4), codes).astype(int)
    initial = (hold + L @ ref).reshape(2, 2)
    product = np.outer(initial.sum(1), initial.sum(0)).ravel()
    v = witness(B, ref, hold, L, product, size)
    table = hold + L @ v
    assert np.array_equal(table.reshape(2, 2).sum(0), initial.sum(0))
    assert np.array_equal(table.reshape(2, 2).sum(1), initial.sum(1))
    log_v = np.array([logarithm(a / b) for a, b in zip(v, ref, strict=True)], object)
    log_t = np.array(
        [logarithm(a / b) for a, b in zip(table, product, strict=True)], object
    )
    Bq = Matrix(B.tolist())
    proj = Bq.T * (Bq * Bq.T).inv() * Bq
    projection = np.array([[Q(x) for x in row] for row in proj.tolist()], object)
    residual = projection @ (log_v + L.T @ log_t)
    rho = max(abs(residual)) + 4 * Q(1, 10**55) * max(abs(projection).sum(1))
    assert np.array_equal(projection @ (v - ref), v - ref)
    assert sum(v) == sum(ref) and all((x * size).denominator == 1 for x in v)
    # Probability of any type and one agreeing shuffle has exponential part
    # exp(N*F), where F = -D(v||ref)-D(table||product). Fixed protected
    # weights contribute zero to the first divergence. Exposure includes
    # drawing all free members at native frequencies before the ordinary gate.
    score = -v @ log_v - table @ log_t
    score_err = 4 * Q(1, 10**55)
    eps = Q(1, 50000)
    # Count at most (N+1)^16 types. A given type probability has lower
    # polynomial factor (N+1)^-18 and upper factor (N+1). Pinsker gives
    # F(p)-F(v) <= 2*rho - 2*TV(p,v)^2, including the four fixed cells.
    log_upper = (size + 1).bit_length()  # natural log of N+1 is smaller
    log_tail = 35 * log_upper + size * (2 * rho - 2 * eps**2)
    assert log_tail < -100
    tv = eps + Q(1, 2**100)
    log_accept = 17 * log_upper + size * (score + score_err + 2 * rho)
    assert log_accept < 0
    z = original.copy()
    for i, num in enumerate(n):
        z[cat == i] = v[i] / int(num)
    assert sum(z) == 1
    blank = (z - original) @ (words[:, -1] == 0)
    singles = (words[:, :4] == 1).sum(1)
    dn = (z - original) @ singles
    b = original @ np.isin(words[:, -1], fibers.source.PARTIAL)
    scale = 4 / b

    def bounds(x, width):
        return list(map(str, (x - width * tv, x + width * tv)))

    assert blank > tv
    delta0 = (z - original) @ (cells == 0)
    assert delta0 > tv
    return dict(
        sha256={
            p: sha256((fibers.source.tidal.HERE / p).read_bytes()).hexdigest()
            for p in ("PEER_RETURN.md", "peer_population.json")
        },
        register_size=size,
        free_counts=[int(p * size) for p in v],
        protected_cell_counts=[int(p * size) for p in hold],
        projected_gradient_upper=str(rho),
        maximum_conditional_mean_tv=str(tv),
        tail_log_upper=math.ceil(log_tail),
        # e > 2, so a negative natural-log bound is also a valid base-2 upper exponent.
        joint_acceptance_upper_power_of_two=math.ceil(log_accept),
        blank_change_bounds=bounds(blank, 1),
        variance_change_bounds=bounds(delta0, 1),
        single_count_change_bounds=bounds(dn, 4),
        ordinary_per_port_bounds=bounds(scale * blank, scale),
        witness_ordinary_coefficients=list(map(str, (scale * blank, dn, -dn))),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = fibers.source.tidal.HERE / "peer_concentration.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(
        [float(Q(x)) for x in out["blank_change_bounds"]],
        out["tail_log_upper"],
        out["joint_acceptance_upper_power_of_two"],
    )
