"""Short-pulse dark-port and heterodyne statistics on the complete connected generator."""

import argparse
from fractions import Fraction
import json
from math import factorial

import numpy as np
from mpmath import iv
from scipy.sparse.linalg import expm_multiply

from . import measure

NATIVE = measure.action.native
PORTS = (0, 1, 38, 39)
TIMES = (1 / 128, 1 / 64, 1 / 32)
PHASES = (1, 1j, -1, -1j)
ROOT = measure.HERE
SHOTS = 2**74
CHUNK = 2**59
# Two charges, three time bins, five source states, three heterodyne reads
# with four phases, and two dark ports, with two-sided error.
TAIL = 2 * 2 * 3 * 5 * (3 * 4 + 2)


def draw(rng, p):
    if not 0 <= p <= 1:
        raise ValueError(f"invalid complete probability {p}")
    return sum(map(int, rng.binomial(CHUNK, p, size=SHOTS // CHUNK)))


def kernels(board, route, q, t):
    idx = np.array([4 * p + j for p in PORTS for j in range(4)])
    sectors = []
    for odd in (False, True):
        h = NATIVE.generator(board, route, q, odd)
        a = -0.5j * t * h
        x = np.zeros((h.shape[0], 16))
        x[idx] = np.eye(16)
        u = (
            expm_multiply(a, x, traceA=a.diagonal().sum())[idx]
            .reshape(4, 4, 4, 4)
            .transpose(0, 2, 1, 3)
        )
        sectors.append(u)
    return np.array(sectors)


def mean(k, x):
    return np.einsum("bri,bij,brj->", k, x, k.conj()).real


def estimate(kernels, t, x, rng, tail):
    def read(k):
        n = draw(rng, mean(k, x))
        w = (iv.sqrt(2 * n * tail) + 2 * tail) / SHOTS
        return [n, SHOTS - n], iv.mpf(n) / SHOTS + iv.mpf([-1, 1]) * w

    bank = {}
    values = {}
    for name, a, b, order in (
        ("left", 0, 1, 1),
        ("right", 2, 3, 1),
        ("third", 0, 3, 3),
    ):
        epsilon = t**order
        observed = []
        raw = []
        for z in PHASES:
            k = (kernels[:, b, a] + epsilon * z * kernels[:, b, b]) / np.sqrt(
                1 + epsilon**2
            )
            counts, p = read(k)
            observed.append(p)
            raw.append(counts)
        signal = (
            iv.mpc(observed[0] - observed[2], observed[1] - observed[3])
            * (1 + epsilon**2)
            / (4 * epsilon)
        )
        values[name] = iv.mpc(0, 1) ** order * factorial(order) * signal / (t**order)
        bank[name] = raw
    dark = []
    for a, b in ((0, 3), (1, 2)):
        k = (
            kernels[:, a, a] + kernels[:, a, b] - kernels[:, b, a] - kernels[:, b, b]
        ) / 2
        counts, p = read(k)
        dark.append(p / t**2)
        bank["dark_" + str(PORTS[a])] = counts
    values["dark"] = (dark[0] + dark[1]) / 2
    return bank, values


def derive(means, angles):
    out = {}
    for key in means["equilibrium"]:
        ref = means["equilibrium"][key]
        rotation = (means["plus"][key] - means["minus"][key]) / 2
        out[key] = dict(
            reference=ref,
            **{
                n: means[n][key] - ref + alpha * rotation for n, alpha in angles.items()
            },
        )
    left, right = out["left"], out["right"]
    lr = measure.product(left, right)
    out["memory"] = {k: out["third"][k] + lr[k] / 2 for k in lr}

    def real(v):
        return {k: str(a.real) for k, a in v.items()}

    return {name: real(out[name]) for name in ("third", "memory", "dark")}


def run():
    iv.dps = 30
    rng = np.random.default_rng(381663)
    tail = iv.ln(TAIL / iv.mpf("0.01"))
    rows = json.loads((ROOT.parent.parent / "signed_transfers.json").read_text())[
        "tests"
    ]
    angles = {
        ("ordinary" if r["source"] == "ordinary" else "pairing"): iv.mpf(
            r["angle_per_unit"]
        )
        for r in rows
    }
    states = measure.action.preparations()
    board = NATIVE.iso.pairing.formation.federation(1)
    route = NATIVE.iso.pairing.formation.region(board)
    charges = {}
    for q in (-1, 1):
        counts = {}
        shots = {}
        for t in TIMES:
            kernels_t = kernels(board, route, q, t)
            c = {}
            v = {}
            for label, x in states.items():
                c[label], v[label] = estimate(kernels_t, t, x, rng, tail)
            counts[str(t)] = c
            shots[t] = v
        # Positive-time quadratic extrapolation of the divided signals.
        # Intervals here concern the actual finite three-dwell statistic.
        # The exact zero-dwell coefficients are enclosed in results.json.
        weights = [
            iv.mpf(f.numerator) / f.denominator
            for f in (Fraction(8, 3), Fraction(-2), Fraction(1, 3))
        ]
        extrapolated = {
            n: {
                k: sum(w * shots[t][n][k] for w, t in zip(weights, TIMES))
                for k in shots[TIMES[0]][n]
            }
            for n in states
        }
        charges[str(q)] = dict(
            counts=counts,
            three_dwell_nominal_intervals_99pct=derive(extrapolated, angles),
            single_dwell_nominal_intervals_99pct={
                str(t): derive(shots[t], angles) for t in TIMES
            },
        )
    return dict(
        seed=381663,
        shots_per_setting=SHOTS,
        phase_order=["0", "pi/2", "pi", "3pi/2"],
        times=list(TIMES),
        charges=charges,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (ROOT / "counts.json").write_text(out)
    else:
        print(out, end="")
