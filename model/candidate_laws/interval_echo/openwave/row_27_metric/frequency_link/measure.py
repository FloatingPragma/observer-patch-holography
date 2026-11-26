"""Replicate the two-station spectra and held frequency links on the same stationary source."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from . import native
from .readout import D, GRID, N, T, curve, errors, ratio, read

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "results.json.gz"
CONFIGS = [(8, "coin_even"), (8, "coin_odd"), (8, "sham"), (0, "coin_even")]
OUT = [3, 5]


def samples(p, rng):
    p = np.atleast_2d(p)
    if np.min(p) < -1e-13 or np.max(abs(p.sum(axis=1) - 1)) > 1e-12:
        raise ArithmeticError("Incomplete probability bank")
    p = np.maximum(p, 0)
    weights = np.rint(p * GRID).astype(np.int64)
    indices = np.argmax(p, axis=1)
    weights[np.arange(len(p)), indices] += GRID - weights.sum(axis=1)
    if np.min(weights) < 0:
        raise ArithmeticError("Invalid finite sampler")
    hist = np.array([rng.multinomial(N, w / GRID) for w in weights])
    return dict(native_probability=p.tolist(), counts=hist.tolist())


def two_port(p, rng):
    return samples(np.stack([1 - p, p], axis=-1), rng)


def read_clock(h, rng):
    raw = two_port(native.hand(h, T), rng)
    raw["fit"] = read(np.asarray(raw["counts"])[:, 1])
    return raw


def forecast(source, flights, t):
    model = source["fit"]
    coeff = model["coefficients"]
    probabilities = np.asarray(flights["counts"]) / N
    predictions = curve(model["frequency"], coeff, t[:, None] - D)
    return predictions @ probabilities / len(D)


def branch(k, arm, q, blocks, inserts, rng):
    h, clocks = native.histories(k, q, arm, blocks, inserts)
    refs = [read_clock(clock, rng) for clock in clocks]
    links = []
    for sender, (end, clock) in enumerate(zip(OUT, clocks, strict=True)):
        receiver = 1 - sender
        dest = OUT[receiver]
        for parked in (False, True):
            flights = samples(native.impulse(h, end, D, parked), rng)
            raw = samples(native.output(clock, h, end, T, D, parked), rng)
            com = np.asarray(raw["counts"])[:, 1 + dest]
            f = read(com)
            r = ratio(f, refs[receiver]["fit"])
            predicted = forecast(refs[sender], flights, T[512:])[:, dest]
            observed = com[512:] / N
            # Independent sharp clock probabilities at the two emission epochs and
            # flight calibration predict the entire held address record directly.
            known = np.asarray(refs[sender]["counts"])[:, 1] / N
            estimates = np.column_stack(
                [known[512 - len_shift : 1024 - len_shift] for len_shift in (32, 36)]
            )
            logged = np.asarray(flights["counts"])[:, dest] / N
            direct = estimates @ logged / 2
            direct_bound = (
                errors(observed)
                + (
                    (estimates + errors(estimates)) @ errors(logged)
                    + errors(estimates) @ logged
                )
                / 2
            )
            max_fraction = float(np.max(abs(direct - observed) / direct_bound))
            if max_fraction > 1:
                raise ArithmeticError(
                    "Independent two-delay forecast failed its simultaneous interval"
                )
            links.append(
                dict(
                    sender=sender,
                    receiver=receiver,
                    parked=parked,
                    impulse=flights,
                    observation=raw,
                    fit=f,
                    ratio=r,
                    independent_forecast_held=dict(
                        rms=float(np.sqrt(np.mean((predicted - observed) ** 2))),
                        max_abs=float(np.max(abs(predicted - observed))),
                        max_observed_count_error=float(np.max(errors(observed))),
                        independent_delayed_clock_max_error_fraction=max_fraction,
                        bounds_from_source_prefix="local sinusoid and impulse fits only",
                    ),
                )
            )
    return dict(contrast=k, source=arm, tag=q, local=refs, links=links)


def run():
    rng = np.random.default_rng(270526)
    blocks, inserts, chart = native.contact_certificate()
    # Measure one appointment step against the separate unit blank dimer.
    step = native.hand(blocks[8, 1], np.array([1 / 16]))
    standard = two_port(step, rng)
    y = standard["counts"][0][1] / N
    err = float(errors(y))
    standard.update(
        reading=float(2 * np.arcsin(np.sqrt(y))),
        interval=[
            float(2 * np.arcsin(np.sqrt(z))) for z in (max(0, y - err), min(1, y + err))
        ],
    )
    return dict(
        seed=270526,
        shots=N,
        sample_probability_denominator=GRID,
        stochastic_law="PCG64 multinomial draws on a fixed probability grid from the calculated finite Born law",
        positions=["no launch", *chart["perimeter"]],
        geometry=chart,
        appointment_step=standard,
        source_census=inserts,
        source_preparation=native.formed_source(),
        storage_routes=native.preparation_routes(chart),
        terminal_tree=native.terminal_tree(chart["perimeter"]),
        banks=[
            branch(k, arm, q, blocks, inserts, rng)
            for k, arm in CONFIGS
            for q in (1, 2)
        ],
    )


def compare(observed, stored):
    if isinstance(stored, dict):
        if stored.keys() != observed.keys():
            raise ValueError("Different receipt fields")
        for k in stored:
            compare(observed[k], stored[k])
    elif isinstance(stored, list):
        if len(stored) != len(observed):
            raise ValueError("Different receipt lengths")
        if stored and (
            isinstance(stored[0], (float, int))
            or (
                isinstance(stored[0], list)
                and all(len(v) == len(stored[0]) for v in stored)
            )
        ):
            a, b = np.asarray(observed), np.asarray(stored)
            if np.issubdtype(b.dtype, np.number):
                good = (
                    np.array_equal(a, b)
                    if b.dtype.kind in "iu"
                    else np.allclose(a, b, atol=5e-9, rtol=2e-8)
                )
                if not good:
                    raise ValueError("Different numerical bank")
                return
        for x, y in zip(observed, stored, strict=True):
            compare(x, y)
    elif isinstance(stored, float):
        if not np.isclose(observed, stored, atol=5e-9, rtol=2e-8):
            raise ValueError("Different reading")
    elif observed != stored:
        raise ValueError("Different setting or count")


if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--output", type=Path)
    target = args.parse_args().output
    receipt = run()
    if target:
        with target.open("xb") as f:
            f.write(
                gzip.compress(
                    json.dumps(
                        receipt, separators=(",", ":"), allow_nan=False
                    ).encode(),
                    mtime=0,
                )
            )
    else:
        compare(receipt, json.loads(gzip.decompress(ARCHIVE.read_bytes())))
        print("Reciprocal material frequency-link receipt verified")
