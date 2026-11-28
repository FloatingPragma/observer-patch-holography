"""Predict sigma-eight for independently active cubic pools at the required mean density."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np

from . import pooled_ordinary as pooled


def overlap(x):
    # Intersection of two equal balls / single ball volume; x is d/R.
    x = np.minimum(x, 2.0)
    return (2 - x) ** 2 * (4 + x) / 16


def window_factor(side, radius, order):
    # Differences of two independent uniform cube positions have triangular
    # coordinate densities. Integrate over the positive octant by symmetry.
    node, weight = np.polynomial.legendre.leggauss(order)
    t = (1 + node) / 2
    weight = weight * (1 - t)
    d = np.sqrt(t[:, None, None] ** 2 + t[None, :, None] ** 2 + t[None, None, :] ** 2)
    integral = np.einsum(
        "i,j,k,ijk->", weight, weight, weight, overlap(d * side / radius)
    )
    return float(side**3 * integral / (4 * np.pi * radius**3 / 3))


def count_variance(n, w, a, b, c):
    d = a * pooled.PORTS * n + c
    active = -math.expm1(n * math.log1p(-w))
    count = n * w
    mean = d * active + b * count
    var = (
        d**2 * active * (1 - active)
        + 2 * d * b * count * (1 - active)
        + b * b * count * (1 - w)
    )
    return active, mean, var


def sensitivity(p, z, coeff, mu, acceptance=1.0):
    # The BBN value with empirical nuclear rate carries statistical and
    # rate errors. Bound the effect of every point in the indicated box.
    da = z * math.hypot(0.00016, 0.00033) / 0.02235
    pa, pb = p * (1 - da), p * (1 + da)

    def rule(x):
        return 6 * x * x * (1 - x) ** 2 * 448 / 2110 * acceptance

    wa, wb = rule(pa), rule(pb)
    assert 0 < pa <= pb < 0.5
    cold_fraction_min = (0.1200 - z * 0.0012) / (0.1430 + z * 0.0011)
    ratio_max = (0.1200 + z * 0.0012) / (0.02237 - z * 0.00015)
    a, b, c = coeff
    assert a > 0 and b < 0 < c
    # A matching mean imposes an upper activation probability. S increases
    # with both pool size and w, giving a limit for every feasible size.
    smax = (4 * pb * ratio_max / mu - b * wb) / (pooled.PORTS * a)
    nmax = math.ceil(math.log1p(-smax) / math.log1p(-wa)) + 1
    cell = pooled.WIDTH / 1000
    Lmax = cell * nmax ** (1 / 3)
    rmax = 8 / (0.6736 - z * 0.0054)
    rmin = 8 / (0.6736 + z * 0.0054)
    # CV² >= (1-S)/S and S<=n*w. Jensen on the equal-ball intersection:
    # E[d] <= sqrt(E[d²]) = side/sqrt(2); the window is convex decreasing.
    raw = max(0, 1 / wb - nmax) * cell**3 / (4 * np.pi / 3)
    jensen = raw * overlap(Lmax / (np.sqrt(2) * rmin)) / rmax**3
    # For any grid offset all intersecting cells are inside the ball
    # enlarged by a cube diagonal. Cauchy bounds their squared weights.
    packing = raw / (rmax + np.sqrt(3) * Lmax) ** 3
    return dict(
        pool_tiles_upper=nmax,
        sigma8_lower=float(cold_fraction_min * np.sqrt(max(jensen, packing))),
        planck_sigma8_upper=0.8111 + z * 0.0060,
    )


def forecasts(cos, coeff, mu, acceptance=1.0):
    a, b, c = coeff
    out = {}
    for label, rec in cos["tests"].items():
        p = rec["filling"]
        n = rec["inferred_tiles"]
        w = 6 * p**2 * (1 - p) ** 2 * 448 / 2110 * acceptance
        active, mean, var = count_variance(n, w, a, b, c)
        ratio = mean * float(mu) / (4 * p * n)
        volume = rec["inferred_side_kpc"] / 1000
        gain = ratio * 0.02237 / 0.1430 * np.sqrt(var) / mean
        coarse, fine = [
            gain * np.sqrt(window_factor(volume, 8 / 0.6736, m)) for m in (20, 40)
        ]
        out[label] = dict(
            activation_probability=active,
            pooled_relative_rms=np.sqrt(var) / mean,
            predicted_sigma8=float(fine),
            quadrature_difference=float(abs(fine - coarse)),
            sensitivity={
                str(z): sensitivity(p, z, (a, b, c), float(mu), acceptance)
                for z in (0, 1, 2)
            },
        )
    return out


def run():
    ref = pooled.HERE / "pooled_ordinary.json"
    native = pooled.scales.budget.LAB / "collective_relay_lab.json"
    coeff = [
        float(Q(x)) for x in json.loads(native.read_bytes())["ordinary_coefficients"]
    ]
    _, mu = pooled.scales.source_data()
    cos = json.loads(ref.read_bytes())["homogeneous"]
    out = forecasts(cos, coeff, mu)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (ref, native, pooled.HERE / "POOLED_FLUCTUATIONS.md")
        },
        planck_h=0.6736,
        planck_omega_m=0.1430,
        planck_sigma8=0.8111,
        planck_sigma8_error=0.0060,
        tests=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = pooled.HERE / "pooled_fluctuations.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        pooled.sky.neighbor.inventory.rp.galaxy.check(
            json.loads(dst.read_bytes()), result
        )
    print(min(v["predicted_sigma8"] for v in result["tests"].values()))
