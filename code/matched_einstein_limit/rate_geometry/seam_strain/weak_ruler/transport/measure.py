"""Intermediate site clocks on a freely propagating, endpoint-counted mobile.

Times are given in native half ticks. The late occupation and the transported
early occupation need not commute. The real ordered product is measurable
using either a sharp intermediate site mark or a tunably weak pointer.
"""

import argparse
from fractions import Fraction
import json
from math import log, sqrt
from pathlib import Path

import numpy as np
from mpmath import iv

from .. import flight, hinge
from ... import live_bounds, live_cadence, live_radar
from .....source_implication.isometric_collar import certify as fixed

HERE = Path(__file__).resolve().parent
PORTS = (0, 1, 38, 39)
DIRECTIONS = ((0, 1, 3), (3, 2, 0))
# Two settings selected by surveying the reference half-tick appointments.
APPOINTMENTS = ((8, 11), (4, 9))
ARMS = ("reference", "ordinary", "decorrelated")
SAMPLES = 2**40
SEED = 48210593


def multiply(a, b):
    return np.stack((a[0] @ b[0] - a[1] @ b[1], a[0] @ b[1] + a[1] @ b[0]))


def adj(a):
    return np.stack((a[0].T, -a[1].T))


def effects(raw, s, t, a, m, b):
    early = raw[:, s - 1, m, a] * fixed.FRAC
    through = multiply(raw[:, t - s - 1, b, m], raw[:, s - 1, m, a])
    free = raw[:, t - 1, b, a] * fixed.FRAC
    mixed = multiply(adj(free), through)

    def gram(x):
        return 2 * multiply(adj(x), x)

    return {
        "arrival": gram(free),
        "site": gram(early),
        "marked_arrival": gram(through),
        "clock": mixed + adj(mixed),
    }


def exact_centers():
    from sympy import ceiling

    states = hinge.sources()
    out = {}
    for tag in (*ARMS, "right"):
        center, radius = states[tag]
        # Tr(E rho) pairs real parts and imaginary parts of Hermitian matrices.
        parts = [list(map(Fraction, p)) for p in center.as_real_imag()]
        erho = Fraction(int(ceiling(radius * 2**90)), 2**90)
        out[tag] = np.array(parts, dtype=object).reshape(2, 8, 8), erho
    return out


def bounds(operators, states, e):
    # Every raw block is within e of a contraction. The through block error
    # is <=2e+e² <=3e; its Gram error is <=6e+9e² <=8e. The mixed effect
    # has smaller error. Every approximate effect has norm <2. The physical
    # source directions have trace norm <=2, including the rail difference.
    assert e < Fraction(1, 10)
    out = {}
    for k, operator in operators.items():
        out[k] = {}
        for tag, (x, erho) in states.items():
            mean = np.sum(operator * x) / (2 * fixed.FRAC**4)
            radius = 2 * erho + 16 * e
            out[k][tag] = live_bounds.rational_interval(mean) + iv.mpf(
                [-1, 1]
            ) * live_bounds.rational_interval(radius)
    return out


def population(values):
    p, site, through, w = (
        values[v] for v in ("arrival", "site", "marked_arrival", "clock")
    )
    avoid = p + through - 2 * w
    ep = live_bounds.endpoints
    return {
        "unmarked": list(map(ep, (p, 1 - p))),
        "sharp": list(map(ep, (through, site - through, avoid, 1 - site - avoid))),
        "joint_clock": ep(w),
        "conditional_clock": ep(w / p),
        "minimum_relative_arrival_change_under_dephasing": ep(-2 * w / p)
        if w.b < 0
        else None,
    }


def transfer(v, modular, volume):
    p, w = v["arrival"], v["clock"]
    z = {
        tag: w[tag] / p["reference"] - w["reference"] * p[tag] / p["reference"] ** 2
        for tag in ("ordinary", "decorrelated", "right")
    }

    # Right rotation holds the actual filled five-station radar volume.
    def matched(x, tag):
        return x[tag] - volume[tag] * x["right"] / volume["right"]

    o, c = (matched(z, tag) for tag in ("ordinary", "decorrelated"))
    ko, kc = (matched(modular, tag) for tag in ("ordinary", "decorrelated"))
    predicted = o * kc / ko
    return {
        k: live_bounds.endpoints(x)
        for k, x in dict(
            ordinary=o,
            decorrelated=c,
            predicted_from_retained_modular_charge=predicted,
            measured_minus_predicted=c - predicted,
        ).items()
    }


def bands(rng, bank):
    """Independent positive event draws at the interval midpoints."""
    out = {}
    for tag in ("unmarked", "sharp"):
        p = np.array(bank[tag]).mean(axis=1)
        p /= p.sum()
        if max(abs(p[:, None] - bank[tag]).ravel()) > 1e-14:
            raise ArithmeticError(
                "simulated propensities exceed the rounding allowance"
            )
        out[tag] = rng.multinomial(SAMPLES, p).tolist()
    return out


def count_clock(counts):
    base, mark = (
        np.asarray(counts[k], dtype=float) / SAMPLES for k in ("unmarked", "sharp")
    )
    # Bernstein for a Bernoulli mean gives |x-p| <= sqrt(2ap)+2a/3,
    # a=log(2/delta)/n. Solve the upper inequality first, then bound both
    # tails. Union-bound 3 event marginals per read in 24 settings.
    a = log(2 * 3 * 24 / 0.01) / SAMPLES
    b = 2 * a / 3

    def half_width(z):
        upper = z + b + a + sqrt(a * a + 2 * a * (z + b))
        return b + sqrt(2 * a * upper) + 1e-14

    eps = sum(map(half_width, (base[0], mark[0], mark[2]))) / 2
    x = (base[0] + mark[0] - mark[2]) / 2
    return [x - eps, x + eps]


def native(charge):
    board = live_radar.source.formation.federation(1)
    return live_radar._fixed_kernels(
        board, live_radar.source.formation.region(board), charge, True, PORTS
    )


def run():
    state = exact_centers()
    modular = {
        t: iv.mpf(v)
        for t, v in json.loads((hinge.HERE.parent / "live_modular.json").read_text())[
            "modular_slopes"
        ].items()
    }
    volume = {
        t: iv.mpf(v)
        for t, v in json.loads((hinge.HERE / "hinge.json").read_text())[
            "finite_appointment_certificate"
        ]["raw_volume_slopes"].items()
    }
    rng = np.random.default_rng(SEED)
    rows = []
    for charge in (-1, 1):
        raw, err = native(charge)
        h, addressed = flight.hardware(charge)
        for a, m, b in DIRECTIONS:
            # For a projector P, ||[H,P]||=||PH(1-P)||. The native integer
            # leakage bound controls its square. Cauchy-Schwarz then gives
            # |dW/ds| <= sqrt(p_arrival)*||[H,P]|| for every positive source.
            j2 = live_cadence.leakage_bound(h, addressed[m])
            for s, t in APPOINTMENTS:
                v = bounds(effects(raw, s, t, a, m, b), state, err)
                arms = {}
                for tag in ARMS:
                    vals = {
                        key: u["reference"] + (0 if tag == "reference" else u[tag])
                        for key, u in v.items()
                    }
                    bank = population(vals)
                    if vals["clock"].b < 0:
                        width = -vals["clock"] / (
                            2
                            * iv.sqrt(
                                live_bounds.rational_interval(j2) * vals["arrival"]
                            )
                        )
                        bank["negative_shutter_half_width_lower"] = (
                            live_bounds.endpoints(width)[0]
                        )
                    counts = bands(rng, bank)
                    bank["counts"], bank["joint_from_counts_99_family"] = (
                        counts,
                        count_clock(counts),
                    )
                    arms[tag] = bank
                rows.append(
                    dict(
                        charge=charge,
                        ports=[PORTS[i] for i in (a, m, b)],
                        mark_time=s / 2,
                        arrival_time=t / 2,
                        arms=arms,
                        fixed_volume_transfer=transfer(v, modular, volume),
                    )
                )
    return dict(
        ports=list(PORTS),
        shots_per_arm_and_instrument=SAMPLES,
        seed=SEED,
        mark_order=[
            "site and arrive",
            "site and miss",
            "elsewhere and arrive",
            "elsewhere and miss",
        ],
        setting_selection="Reference half-appointment search, then apply the same two settings to both charges, directions and held preparations.",
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(report)
    else:
        print(report, end="")
