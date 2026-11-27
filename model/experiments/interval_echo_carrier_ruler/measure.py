"""Acquire all three exit bins; infer comparison lengths independently of source labels."""

import argparse
import gzip
import hashlib
import json
from inspect import getfile, unwrap
from math import acos, ceil, floor, inf, isfinite, log, pi, sqrt
from pathlib import Path

import numpy as np

from . import apparatus as app

HERE = Path(__file__).resolve().parent
SHOTS = 2**24
GRID = 2**28
SEED = 731126
BUDGET = 4096
ARMS = {
    "nominal": (1.5, 1, "native"),
    "dilated": (2, 1, "native"),
    "gain": (1.5, 0.8, "native"),
    "folded": (8.5, 1, "native"),
    "blind": (1.5, 1, "exchange_off"),
}


def draw(rng, p):
    np.testing.assert_allclose(sum(p), 1, atol=4e-13, rtol=0)
    scaled = np.asarray(p[:2]) * GRID
    if np.any(abs(scaled - np.floor(scaled) - 0.5) < 1e-6):
        raise ArithmeticError("Unstable probability bin")
    n = np.rint(scaled).astype(np.int64)
    # A full rod state has exactly zero support outside its two resolved exits.
    if p[2] < 1e-23:
        n[1] = GRID - n[0]
    total = np.r_[n, GRID - sum(n)]
    if any(total < 0):
        raise ValueError("Invalid exit distribution")
    return rng.multinomial(SHOTS, total / GRID).tolist()


def interval(x, n, conditional):
    if not n:
        return 0, 1
    delta = sqrt(log(2 * BUDGET / 0.01) / (2 * n)) + (64 if conditional else 2) / GRID
    return max(0, x / n - delta), min(1, x / n + delta)


def phase(counts):
    p = interval(counts[1], sum(counts[:2]), True)
    return tuple(acos(1 - 2 * x) for x in p)


def length_interval(sample, rod):
    lo, hi = phase(sample)
    lowref, hiref = phase(rod)
    return (lo / hiref if hiref else 0, hi / lowref if lowref else inf)


def bounds_pair(lo, hi):
    # Round away from the confidence set on a reproducible decimal grid.
    return [floor(lo * 1e7) / 1e7, ceil(hi * 1e7) / 1e7 if isfinite(hi) else None]


def analysis(calibration, held):
    scales = {}
    for arm, rows in calibration.items():
        bounds = []
        by_wave = []
        for record in rows:
            block = [
                length_interval(counts, record["reference"][station])
                for station, out in enumerate(record["edges"])
                for counts in out
            ]
            bounds += block
            by_wave.append(
                bounds_pair(min(a for a, _ in block), max(b for _, b in block))
            )
        responsive = all(
            phase(rows[0]["reference"][i])[0] > phase(rows[1]["reference"][i])[1]
            for i in (0, 1)
        )
        lo, hi = max(a for a, _ in bounds), min(b for _, b in bounds)
        if not responsive:
            status = "unresolved_reference"
        elif lo > hi:
            status = "no_common_branch"
        else:
            status = "calibrated"
        scales[arm] = dict(
            status=status,
            joint_intersection=bounds_pair(lo, hi) if status == "calibrated" else None,
            calibration_hulls=by_wave,
        )
    energies = {}
    for arm, axes in held.items():
        energies[arm] = {}
        for name, banks in axes.items():
            intervals = [
                tuple(
                    sum(values)
                    for values in zip(
                        *(interval(c[1], sum(c), False) for c in v), strict=True
                    )
                )
                for v in banks
            ]
            energies[arm][name] = bounds_pair(
                max(a for a, _ in intervals),
                min(b for _, b in intervals),
            )
    return scales, energies


def acquire():
    instrument = app.Ruler()
    rng = np.random.Generator(np.random.PCG64(SEED))
    points = instrument.points
    calibrations = {}
    for arm, (a, gain, receiver) in ARMS.items():
        rows = []
        for n, g in enumerate((pi / 5, pi / 7)):
            wave = gain * g
            rod = instrument.prepare(f"rod/{arm}/{n}", app.ROD, [0, wave])
            reference = [
                draw(rng, instrument.read(rod, app.ROD, i, receiver)) for i in (0, 1)
            ]
            banks = [[], []]
            for j, (u, v) in enumerate(app.EDGES):
                axis = points[v] - points[u]
                np.testing.assert_allclose(axis @ axis, 1, atol=2e-14, rtol=0)
                psi = instrument.prepare(
                    f"cell/{arm}/{n}/{j}", range(12), wave * a * (points @ axis)
                )
                for i in (0, 1):
                    banks[i].append(
                        draw(rng, instrument.read(psi, (u, v), i, receiver))
                    )
            rows.append(dict(reference=reference, edges=banks))
        calibrations[arm] = rows
    held = {}
    for arm in ("nominal", "dilated"):
        a = ARMS[arm][0]
        by_axis = {}
        for name, axis in instrument.axes.items():
            psi = instrument.prepare(
                f"held/{arm}/{name}", range(12), (8 / 3) * a * (points @ axis)
            )
            by_axis[name] = [
                [draw(rng, instrument.read(psi, e, i, "native")) for e in app.EDGES]
                for i in (0, 1)
            ]
        held[arm] = by_axis
    scales, energies = analysis(calibrations, held)
    sources = [
        HERE / "PROTOCOL.md",
        Path(__file__),
        Path(app.__file__),
        Path(app.symbol.__file__),
        Path(app.symbol.native_law.__file__),
        Path(app.symbol.carrier.__file__),
        *(
            Path(getfile(unwrap(f)))
            for f in (
                app.symbol.gates,
                app.exchange_word,
                app.symbol.path,
                app.symbol.geometry,
            )
        ),
    ]
    return dict(
        provenance="Simulated multinomial draws from complete native routed sign and off-seam probabilities.",
        reference="One independently enrolled two-port rod defines L. SI conversion is external.",
        scale_units="L",
        nominal="3/2",
        seed=SEED,
        shots=SHOTS,
        grid=GRID,
        budget=BUDGET,
        layout=instrument.layout(),
        preparations=instrument.calendars,
        calibration=calibrations,
        held=held,
        intervals=scales,
        energies=energies,
        source_sha256={
            str(f.resolve().relative_to(HERE.parents[2])): hashlib.sha256(
                f.read_bytes()
            ).hexdigest()
            for f in sources
        },
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    packet = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    if args.output:
        with args.output.open("xb") as out:
            out.write(gzip.compress(packet, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != packet:
        raise ValueError("Ruler replay differs")
