"""Outward enclosures for the two outer stationary sources and their finite readout differences."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import measure
from .. import invisible
from ... import certificate as c
from ....directional_stress.regional_modular import energy_odd

SCALE = 2**38
R = sp.Rational


def grid(x):
    return sp.Matrix(np.rint(x.real * SCALE).astype(np.int64)) + sp.I * sp.Matrix(
        np.rint(x.imag * SCALE).astype(np.int64)
    )


def norm(m):
    return c.iv.sqrt(c.ri(sum(sp.re(x) ** 2 + sp.im(x) ** 2 for x in m)))


def center(channel, weights, proposed, basis):
    # Inputs are charge-averaged binary64 entries rounded to a dyadic grid,
    # then they are combined with exact rational lottery probabilities.
    phi = (
        sum(
            (R(float(w)) * m for w, m in zip(weights, channel, strict=True)),
            sp.zeros(16),
        )
        / SCALE
    )
    jac = (basis.H * phi * basis).applyfunc(sp.re)
    shift = (basis.H * phi * sp.eye(4).reshape(16, 1) / 4).applyfunc(sp.re)
    upper = R(".919")
    _, pivots = (upper**2 * sp.eye(15) - jac.T * jac).LDLdecomposition(hermitian=False)
    if min(pivots.diagonal()) <= 0:
        raise ArithmeticError("contraction certificate failed")
    # One unit trace input has HS <=1. The Hermitian trace-norm perturbation
    # uses ||X||1<=2||X||HS on the traceless sector. A full 16x16 complex
    # grid contributes <32/SCALE in operator norm, including charge averaging.
    grid_error = R(32, SCALE)
    error = R(invisible.stacked_error())
    if upper + 2 * error + grid_error >= R(".92"):
        raise ArithmeticError("true channel contraction failed")
    coordinates = np.einsum("kij,ji->k", measure.BASIS, proposed).real
    coords = sp.Matrix([R(round(float(z) * SCALE), SCALE) for z in coordinates])
    rho = sp.eye(4) / 4 + (basis * coords).reshape(4, 4)
    deviation = norm(rho - sp.eye(4) / 4)
    if deviation.b >= c.ri(".13"):
        raise ArithmeticError("source center positivity failed")
    residual = norm(jac * coords + shift - coords)
    radius = (residual + c.ri(error + grid_error)) / c.ri(".08")
    return rho, radius, deviation


def entropy_correction(rho, radius, dev):
    # Entropy of each normalized 4x4 block minus log(4), by the series
    # of (1+x)log(1+x) for X=4*rho-I. Integrate the centered gradient
    # along its certified trace-preserving fixed-point correction.
    x = 4 * rho - sp.eye(4)
    # All factors live on one dyadic grid. Integer products avoid symbolic expression growth.
    d = 2 * SCALE
    xint = (d * x).applyfunc(sp.expand)
    power = xint
    ans = R(0)
    order = 16
    for n in range(2, order + 1):
        power = (power * xint).applyfunc(sp.expand)
        ans -= (-1) ** n * power.trace() / (4 * d**n * n * (n - 1))
    z = 4 * dev
    tail = z ** (order + 1) / ((order + 1) * order * (1 - z))
    lip = (dev + radius) / (c.ri(".25") - dev - radius)
    return c.ri(ans) + c.iv.mpf([-1, 1]) * (tail + lip * radius)


def read(ranges):
    logs = np.full((4, 4), c.ri(0), dtype=object)
    for value, (a, b) in zip(ranges, measure.reuse.ROUTES, strict=True):
        logs[a, b] = c.iv.ln(value)
    gain = [sum(logs[i, j] - logs[j, i] for j in range(4)) / 4 for i in range(4)]
    lengths = [
        c.iv.exp((logs[a, b] + logs[b, a] - gain[a] - gain[b]) / 2)
        for a, b in combinations(range(4), 2)
    ]
    area, volume = c.enc.geometry_jet(np.array(ranges, dtype=object))
    return lengths, area.v, volume.v


def verify(bank, effects, receipt):
    c.iv.dps = 32
    basis = sp.Matrix.hstack(
        *(c.prior.sym_matrix(m).reshape(16, 1) for m in measure.BASIS)
    )
    channels = [[grid(m) for m in sector] for sector in bank]
    # The integer radar estimate has operator error below 1e-14 before
    # conversion. Four complex rows at 1/SCALE dominate the extra rounding.
    effects = [[grid(a) / SCALE for a in pair] for pair in effects]
    odds = [R(8657, 10000), R(1343, 10000)]
    # Ordinary energy current on the central native seam, reconstructed
    # independently by the integer four-tick kernel on the twelve-port collar.
    slopes, _ = energy_odd.energy_rows()
    fair = (slopes[0, :] + slopes[1, :]) / 2
    current_ops = [
        (basis * sp.Matrix(fair[0, 15 * b : 15 * (b + 1)]).T).reshape(4, 4)
        for b in range(2)
    ]
    out = {}
    for mode, row in receipt["comparisons"].items():
        data = []
        for item in row["endpoints"][-2:]:
            w = item["routes"]
            if sum(R(float(v)) for v in w) != 1 or min(w) <= 0:
                raise ArithmeticError(
                    "lottery must be normalized and strictly positive"
                )
            guess, _ = measure.prepare(bank, np.array(w))
            states = [
                center(ch, w, guess[b] / float(odds[b]), basis)
                for b, ch in enumerate(channels)
            ]
            s = sum(
                c.ri(o) * entropy_correction(*state)
                for o, state in zip(odds, states, strict=True)
            )
            ranges = []
            for pair in effects:
                mean = sum(
                    o * c.inner(h, rho)
                    for o, h, (rho, _, _) in zip(odds, pair, states, strict=True)
                )
                width = sum(
                    c.ri(o) * norm(h - h.trace() * sp.eye(4) / 4) * r
                    for o, h, (_, r, _) in zip(odds, pair, states, strict=True)
                )
                ranges.append(
                    c.ri(mean) + c.iv.mpf([-1, 1]) * (width + c.ri(R(8, SCALE)))
                )
            cur = sum(
                c.ri(o * c.inner(h, rho))
                for o, h, (rho, _, _) in zip(odds, current_ops, states, strict=True)
            )
            cur_error = sum(
                c.ri(o) * norm(h) * radius
                for o, h, (_, radius, _) in zip(odds, current_ops, states, strict=True)
            )
            data.append(
                (
                    s,
                    *read(ranges),
                    cur + c.iv.mpf([-1, 1]) * (cur_error + c.ri("2e-15")),
                )
            )
        low, high = data
        ds = high[0] - low[0]
        if ds.a <= 0:
            raise ArithmeticError("entropy gap unresolved")
        dj = high[4] - low[4]
        if dj.b >= 0:
            raise ArithmeticError("seam energy current change unresolved")
        out[mode] = dict(
            entropy_increase_nats=str(ds),
            central_energy_current_change=str(dj),
            maximum_absolute_edge_difference=str(
                max(abs(a - b).b for a, b in zip(low[1], high[1], strict=True))
            ),
            absolute_area_difference_upper=str(abs(high[2] - low[2]).b),
            absolute_volume_difference_upper=str(abs(high[3] - low[3]).b),
        )
    return out


def run():
    bank, effects = measure.load_native()
    return verify(
        bank, effects, json.loads((measure.HERE / "results.json").read_text())
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(output)
    else:
        print(output, end="")
