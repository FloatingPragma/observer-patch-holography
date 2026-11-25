"""Exact four-tick integer wave amplitudes and interval bounds on five formed-source cameras."""

import json

import numpy as np
import sympy as sp
from mpmath import iv

from ....directional_stress.traction.certificate import four_tick_columns
from ....full_algebra import native_precision as native
from . import bounds, measure

DEN = 2**80
ERR = "2e-18"
PROBABILITY_DEN = 8 * DEN**2 * native.prior.FRAC**2


def quantize(c):
    return np.array(
        [
            [[int(round(part(c[i, j]) * DEN)) for j in range(4)] for i in range(4)]
            for part in (sp.re, sp.im)
        ],
        dtype=object,
    )


def optical(block, coin):
    a, b = block
    c, d = coin
    e = a @ c - b @ d
    f = a @ d + b @ c
    sr = e @ a.T + f @ b.T
    si = f @ a.T - e @ b.T
    z = sum((c * sr - d * si).flat)
    return int(z)


def product(rho):
    left = sp.Matrix(
        2, 2, lambda i, j: sum(rho[2 * i + b, 2 * j + b] for b in range(2))
    )
    right = sp.Matrix(
        2, 2, lambda i, j: sum(rho[2 * b + i, 2 * b + j] for b in range(2))
    )
    return sp.kronecker_product(left, right)


def entropy_bounds(rho, tau):
    delta = rho - tau
    answers = []
    for c, err, lower in (
        (rho, sp.Rational(1, 10**18), sp.Rational(6, 10**9)),
        (tau, sp.Rational(2, 10**18), sp.Rational(1, 500)),
    ):
        vals, norm = native.prior.entropy_coordinates(c, err, lower, [c, delta])
        # The normalized state has trace norm one, the difference at most two.
        # The true source-minus-product direction differs by <3e-18 in trace norm.
        tol = 2 * norm + 21 * sp.Rational(3, 10**18)
        answers.append(
            [
                bounds.pair(
                    iv.mpf(int(v.p)) / int(v.q)
                    + iv.mpf([-1, 1]) * int(tol.p) / int(tol.q)
                )
                for v in vals
            ]
        )
    return dict(
        formed_entropy=answers[0][0],
        formed_modular_slope=answers[0][1],
        product_entropy=answers[1][0],
        product_modular_slope_exact=0,
    )


def probabilities(rho, tau):
    stocks = [quantize(s) for s in (tau, (rho + tau) / 2, rho)]
    board = measure.camera.front.optics.pairing.formation.federation(1)
    route = measure.camera.front.optics.pairing.formation.region(board)
    mounts = measure.mounts(board, route)
    src = sorted({i for s, _ in mounts.values() for i in s})
    dst = sorted({i for _, s in mounts.values() for i in s})
    ix = measure.camera.indices(src)
    iy = measure.camera.indices(dst)
    coeff = native.prior.weights()[7]
    records = {}
    for q in (-1, 1):
        panels = {}
        for odd in (False, True):
            h = native.generator(board, route, q, odd)
            w = (
                four_tick_columns(h, iy, ix, coeff)
                .reshape(2, len(dst), 4, len(src), 4)
                .transpose(1, 3, 0, 2, 4)
            )
            for edge, (s, t) in mounts.items():
                if not odd:
                    panels[str(edge)] = np.zeros((5, 5, 5), dtype=object)
                table = panels[str(edge)]
                for j, v in enumerate(t):
                    for i, u in enumerate(s):
                        block = w[dst.index(v), src.index(u)]
                        if odd:
                            # All tables use twice the single-sector denominator.
                            table[3, j, i] = table[0, j, i] // 2 + optical(
                                block, stocks[0]
                            )
                            table[4, j, i] = table[2, j, i] // 2 + optical(
                                block, stocks[2]
                            )
                        else:
                            for k, z in enumerate(stocks):
                                table[k, j, i] = 2 * optical(block, z)
        records[str(q)] = panels
    return records


def encode(p):
    # Round to the nearest attoprobability for compact exact reproducibility.
    return (
        np.array(
            [
                (int(z) * 10**18 + PROBABILITY_DEN // 2) // PROBABILITY_DEN
                for z in p.flat
            ],
            dtype="int64",
        )
        .reshape(p.shape)
        .tolist()
    )


def decode(serial):
    iv.dps = 40
    bands = []
    for a in np.asarray(serial).flat:
        bands.append(iv.mpf(int(a)) / 10**18 + iv.mpf([-1, 1]) * iv.mpf("4e-18"))
    return np.asarray(bands, dtype=object).reshape(np.asarray(serial).shape)


def report(exact=None):
    iv.dps = 40
    rho = native.formed_state()
    tau = product(rho)
    if exact is None:
        exact = probabilities(rho, tau)
    summary = {}
    for charge, cameras in exact.items():
        summary[charge] = {}
        for g, p in cameras.items():
            values = encode(p)
            table = decode(values)
            # Verify the published rounded center encloses the full analytic
            # wave/formation error, including decimal/binary conversions.
            for z, outer in zip(p.flat, table.flat, strict=True):
                mid = iv.mpf(int(z)) / PROBABILITY_DEN
                assert outer.a < mid.a - iv.mpf(ERR) and mid.b + iv.mpf(ERR) < outer.b
            summary[charge][g] = dict(
                probabilities_times_1e18=values,
                response=bounds.analyze(table),
                formed_source=bounds.analyze(table, formed=True),
            )
    pooled = {}
    for g in summary["-1"]:
        table = (
            decode(summary["-1"][g]["probabilities_times_1e18"])
            + decode(summary["1"][g]["probabilities_times_1e18"])
        ) / 2
        pooled[g] = dict(
            response=bounds.analyze(table),
            formed_source=bounds.analyze(table, formed=True),
        )
    return dict(
        native_probability_radius=4e-18,
        material=entropy_bounds(rho, tau),
        charges=summary,
        pooled=pooled,
    )


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
