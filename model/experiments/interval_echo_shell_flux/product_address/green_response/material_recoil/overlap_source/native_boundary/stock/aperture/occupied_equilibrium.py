"""Certified shell work and internal energy for a defect in an occupied chart."""

import argparse
from hashlib import sha256
import json
import math

from mpmath import iv, mp
import numpy as np
from scipy import sparse

from . import occupied_release as flight

HERE = flight.HERE
BITS = 40
BASE = 1 << BITS
SCALE = BASE**3
ORDER = 112
SHELLS = tuple(range(7))
DECIMALS = 56


def next_power(m, digits, shift):
    """Floor each exact entry of M times a fixed-point matrix divided by 2^shift."""
    mask = BASE - 1
    scaled = [m @ x for x in digits]
    for i in (0, 1):
        scaled[i + 1] += scaled[i] >> BITS
        scaled[i] &= mask
    # The top digit can represent unity in the normalized fixed-point matrix.
    return [
        (
            (scaled[j] >> shift)
            | ((scaled[j + 1] & ((1 << shift) - 1)) << (BITS - shift))
        )
        for j in (0, 1)
    ] + [scaled[2] >> shift]


def heat_shells(twice, radii, order=ORDER, *, include_trace=False):
    # m=2*(C I-H) is nonnegative. L is a power of two whose row norm
    # dominates m. This makes floor-propagation contract in the max norm.
    n = len(twice)
    offset = int(max(twice.diagonal()))
    m = sparse.csr_matrix(offset * np.eye(n, dtype=np.int64) - twice)
    assert min(m.data) >= 0
    shift = int(np.asarray(m.sum(axis=1)).max()).bit_length()
    length = 2**shift
    assert shift < BITS and (length + 1) * BASE < 2**62
    z = np.zeros((n, n), dtype=np.int64)
    digits = [z, z, np.eye(n, dtype=np.int64) * BASE]
    beta = iv.log(iv.mpf(8) / 7)
    t = beta * length / 2
    masks = [np.ones(n, bool), *(radii == r for r in SHELLS)]
    coeff = iv.mpf(1)
    total = [iv.mpf(0) for _ in masks]
    derivative = iv.mpf(0)
    sizes = [int(mask.sum()) for mask in masks]
    for k in range(order + 1):
        for index, sel in enumerate(masks):
            numerator = sum(
                int(v.diagonal()[sel].sum()) * BASE**j for j, v in enumerate(digits)
            )
            # Repeated flooring loses at most k/SCALE per entry.
            v = iv.mpf([numerator, numerator + k * sizes[index]]) / SCALE
            total[index] += coeff * v
            if index == 0:
                derivative += k * coeff * v
        if k < order:
            digits = next_power(m, digits, shift)
            coeff *= t / (k + 1)
    ratio = t / (order + 2)
    assert ratio < 1 and order > t
    first = coeff * t / (order + 1)
    tail = first / (1 - ratio)
    for j, size in enumerate(sizes):
        total[j] += iv.mpf([0, size]) * tail
    derivative += (
        iv.mpf([0, n]) * first * ((order + 1) / (1 - ratio) + ratio / (1 - ratio) ** 2)
    )
    energy = iv.mpf(offset) / 2 - derivative / (total[0] * beta)
    probabilities = [v / total[0] for v in total[1:]]
    # Independent spectral evaluation of all retained moments.
    energies, basis = np.linalg.eigh(twice / 2)
    gibbs = np.exp(-float(beta.mid) * (energies - min(energies)))
    gibbs /= gibbs.sum()
    density = basis**2 @ gibbs
    np.testing.assert_allclose(float(energy.mid), energies @ gibbs, atol=4e-12)
    np.testing.assert_allclose(
        [float(q.mid) for q in probabilities],
        [sum(density[s]) for s in masks[1:]],
        atol=4e-14,
    )
    if include_trace:
        return probabilities, energy, total[0] * iv.exp(-beta * offset / 2)
    return probabilities, energy


def mixture(vals, weights):
    return [
        sum(w * v for w, v in zip(weights, cell, strict=True))
        for cell in zip(*vals, strict=True)
    ]


def scaled(value):
    factor = 10**DECIMALS
    y = value * factor
    with mp.workdps(90):
        return [
            str(int(mp.floor(mp.make_mpf(y._mpi_[0])))),
            str(int(mp.ceil(mp.make_mpf(y._mpi_[1])))),
        ]


def summarize(twice, marks, rad, p, delta, rng):
    vals = []
    for a, b in flight.SECTORS:
        prob, en = heat_shells(twice + np.diag(2 * (a * marks[0] + b * marks[1])), rad)
        vals.append([*prob, en])
    sharp = [
        mixture(vals, [iv.mpf(int(round(w * 16))) / 16 for w in row])
        for row in flight.CONDITIONAL
    ]
    # Native counts (denominator 448) and the signed fixed-marginal direction.
    native = mixture(sharp, [iv.mpf(int(round(x * 448))) / 448 for x in p])
    mixed = mixture(sharp, flight.MIX.tolist())
    beta = iv.log(iv.mpf(8) / 7)
    gain = [m / n for m, n in zip(mixed[:-1], native[:-1], strict=True)]
    measured = iv.mpf(delta.tolist())

    def work(i, j):
        # Cancel the common unknown head displacement before interval division.
        return (
            -iv.log1p(measured * (gain[j] - gain[i]) / (1 + measured * gain[i])) / beta
        )

    changes = {
        f"{i}_{j}": work(x, y)
        for x, i in enumerate(SHELLS)
        for y, j in enumerate(SHELLS)
        if x < y
    }
    shape = abs(changes["4_6"] / changes["2_4"])
    assert math.isfinite(float(shape.b))
    # Direct simulated position measurements read complete shell frequencies.
    data = []
    for arm in (native, [n + m / 14 for n, m in zip(native, mixed, strict=True)]):
        probs = [float(x.mid) for x in arm[:-1]]
        data.append(flight.counts(rng, [*probs, 1 - sum(probs)]))
    out = flight.lab.interval_output
    return dict(
        work={k: out(v) for k, v in changes.items()},
        shells=[*SHELLS, "other"],
        decimal_places=DECIMALS,
        native_probability_bounds=[scaled(v) for v in native[:-1]],
        mixed_probability_bounds=[scaled(v) for v in mixed[:-1]],
        counts=data,
        absolute_outer_to_inner=out(shape),
        internal_energy_change=out(measured * mixed[-1]),
        internal_energy_per_joint_increment=out(mixed[-1]),
        per_sector_trace_bins="individually normalized at beta=log(8/7)",
    )


def run():
    board, old, edges, rad, _ = flight.chart()
    p, _ = flight.sector_prior(old)
    delta = flight.lab.probe.log_factor(old["source"])
    iv.dps = 65
    rng = np.random.default_rng(58942633)
    panels = {}
    for load, tag in ((0, "vacancy"), (2, "excess")):
        twice, marks = flight.apparatus(board, edges, load, -1, 1)
        panels[tag] = summarize(twice, marks, rad, p, delta, rng)
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("OCCUPIED_RELEASE.md", "pair_contact.json")
        },
        order=ORDER,
        precision_bits=3 * BITS,
        shots_per_arm=flight.N,
        results=panels,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "occupied_equilibrium.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    for label, row in out["results"].items():
        print(label, row["absolute_outer_to_inner"], row["internal_energy_change"])
