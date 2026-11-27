"""Acquire two-clock, seven-receiver spectra after the fixed coin schedule."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import binom

from .. import measure as flight
from . import material, readback

HERE = Path(__file__).parent
SEED = 625013
STEP = 1 / 8
SAMPLES = 4096
SHOTS = 1 << 19
CLOCK_SHOTS = 1 << 24
GRID = 1 << 26


def integer_draws(rng, pmf, shots):
    """Sequential conditional binomial quantiles from uniform 52-bit midpoints."""
    weights = np.rint(pmf * GRID).astype(np.int64)
    largest = np.argmax(weights, axis=-1)[..., None]
    residual = GRID - weights.sum(axis=-1, keepdims=True)
    np.put_along_axis(
        weights,
        largest,
        np.take_along_axis(weights, largest, axis=-1) + residual,
        axis=-1,
    )
    if np.min(weights) < 0:
        raise ArithmeticError("signed sampling weight")
    out = np.empty(weights.shape, dtype=np.int32)
    remaining = np.full(weights.shape[:-1], shots, dtype=np.int32)
    mass = np.full(weights.shape[:-1], GRID, dtype=np.int64)
    for i in range(weights.shape[-1] - 1):
        w = weights[..., i]
        part = np.where(w == mass, remaining, 0)
        use = (w > 0) & (w < mass) & (remaining > 0)
        uniforms = (rng.integers(0, 1 << 52, size=np.count_nonzero(use)) + 0.5) / (
            1 << 52
        )
        part[use] = binom.ppf(uniforms, remaining[use], w[use] / mass[use]).astype(
            np.int32
        )
        out[..., i] = part
        remaining -= part
        mass -= w
    out[..., -1] = remaining
    return out


def acquire():
    flight.inherited.strand.native_checks()
    route, source, reference, head, first = flight.inherited.geometry()
    sites = tuple(range(first, len(route)))
    waits = np.arange(SAMPLES) * STEP
    times = np.r_[0, 1, 2 + waits]
    rng = np.random.default_rng(SEED)
    phase = (1 + np.stack([np.cos(times), np.sin(times)], axis=-1)) / 2
    clocks = integer_draws(rng, np.stack([phase, 1 - phase], axis=-1), CLOCK_SHOTS)[
        ..., 0
    ]
    original = json.loads((HERE.parent / "records.json").read_text())
    if original["route"] != [list(pair) for pair in route]:
        raise ValueError("direct-propagation record uses a different strand")
    output = []
    for kappa in (8, 0):
        joint = material.joint(
            kappa, len(route), head, sites, (source, reference), waits
        )
        bank = joint.reshape(*joint.shape[:3], 8)
        # The previous full-orientation calculation measured the first receiver
        # on this 14-dimer board at all three of its terminal appointments.
        for row in original["banks"]:
            if row["kappa"] == kappa and row["menu"] == "spaced":
                trial = round((row["appointment"] - 2) / STEP)
                origin = (source, reference).index(row["source"])
                np.testing.assert_allclose(
                    bank[trial, 0, origin], row["probabilities"], rtol=0, atol=3e-10
                )
        output.append([integer_draws(rng, bank, SHOTS) for _ in range(2)])
    return dict(
        counts=np.array(output),
        clocks=clocks,
        sites=np.array(sites),
        route=np.array(route),
        sources=np.array([source, reference]),
        shots=np.array(SHOTS),
        hand_shots=np.array(CLOCK_SHOTS),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = acquire()
    decoded = readback.decode(**raw)
    if args.write:
        np.savez_compressed(HERE / "counts.npz", **raw)
        (HERE / "results.json").write_text(json.dumps(decoded, indent=2) + "\n")
    else:
        with np.load(HERE / "counts.npz") as original:
            for key in raw:
                np.testing.assert_array_equal(raw[key], original[key])
        readback.verify(decoded, json.loads((HERE / "results.json").read_text()))
    print(json.dumps(decoded["summary"], indent=2))


if __name__ == "__main__":
    main()
