"""Spectroscopy and public position readout of a coherently written occupied core."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from . import core

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results.json.gz"
SHOTS = 2**22
SEED = 310526
GRID = 2**35
BINS = 4096
STEP = 0.25
INDEX = np.arange(BINS)
PERIOD = 2 * np.pi / STEP
ENERGY = (INDEX * PERIOD / BINS + 18) % PERIOD - 18


def sample(prob, rng):
    p = np.asarray(prob)
    if not np.all(np.isfinite(p)) or np.min(p) < 0 or abs(np.sum(p) - 1) > 2e-10:
        raise ArithmeticError("Incomplete outcome law")
    weights = np.rint(p * GRID).astype(np.int64)
    weights[np.argmax(p)] += GRID - weights.sum()
    if min(weights) < 0:
        raise ArithmeticError("Invalid sampling grid")
    return dict(
        counts=rng.multinomial(SHOTS, weights / GRID).tolist(),
        native_probability=p.tolist(),
        l1_grid_error=float(np.abs(p - weights / GRID).sum()),
    )


def timing(e, source):
    bank = np.zeros((2, BINS))
    for start in range(0, len(e), 128):
        end = start + 128
        cycles = e[start:end, None] * BINS / PERIOD
        delta = (cycles - INDEX + BINS / 2) % BINS - BINS / 2
        p = (np.sinc(delta) / np.sinc(delta / BINS)) ** 2
        bank += (abs(source[:, start:end]) ** 2) @ p
    return bank


def calibration(rng, z):
    # Supplied bracket [0,64] removes winding for this coarse quadrature.
    wait = np.pi / 128
    raw = []
    for part in (z.real, z.imag):
        bright = (1 + part) / 2
        raw.append(sample([bright, 1 - bright], rng)["counts"])
    x, y = [2 * c[0] / SHOTS - 1 for c in raw]
    angle = np.arctan2(y, x)
    radius = np.arcsin(2 * np.sqrt(np.log(2 * 64 / 0.01) / SHOTS)) / wait
    return dict(
        real_and_imag_counts=raw, reading=float(angle / wait), radius=float(radius)
    )


def panel(board, seam, k, kind, odd, rng):
    e, src, positions = core.standing_modes(board, seam, k, kind, odd)
    if (
        e.min() < -18
        or e.max() > 7
        or np.max(abs(np.sum(abs(src) ** 2, axis=1) - 1)) > 1e-10
    ):
        raise ArithmeticError("Energy alias or incomplete source")
    spectra = timing(e, src)
    neg = np.flatnonzero(e < -1e-9)
    lines = sorted(
        (
            [float(e[i]), float(abs(src[0, i]) ** 2), float(abs(src[1, i]) ** 2)]
            for i in neg
        ),
        key=lambda x: x[0],
    )
    banks = []
    for sign, state, p in zip(
        ("ordinary", "one_endpoint_phase"), src, spectra, strict=True
    ):
        clock = sample(p, rng)
        peak = int(np.argmax(clock["counts"]))
        flights = []
        for t in (0, 1, 8, 32):
            spatial = positions(state * np.exp(-1j * t * e))
            observed = sample(spatial, rng)
            observed.update(
                time=t,
                endpoint_probability=float(spatial[list(seam)].sum()),
                endpoint_counts=sum(observed["counts"][i] for i in seam),
            )
            flights.append(observed)
        banks.append(
            dict(
                source=sign,
                clock=clock,
                peak_bin=peak,
                peak_energy=float(ENERGY[peak]),
                flight=flights,
            )
        )
    # Calibrate binding only from the ordinary source's maximum clock bin.
    # Use the opposite-phase bank for the held comparison, with all bins retained.
    anchor = banks[0]["peak_energy"]
    windows = []
    if anchor < 0:
        for principal in (2, 3):
            target = anchor / principal**2
            select = np.abs(ENERGY - target) <= 0.1
            counts = np.array(banks[1]["clock"]["counts"])
            windows.append(
                dict(
                    n=principal,
                    target=target,
                    halfwidth=0.1,
                    probability=float(spectra[1, select].sum()),
                    counts=int(counts[select].sum()),
                )
            )
    return dict(
        core=kind,
        contrast=k,
        directional_parity="odd" if odd else "even",
        below_minus_1e_9=lines,
        range=[float(e.min()), float(e.max())],
        banks=banks,
        held_clock_windows=windows,
    )


def run():
    rng = np.random.default_rng(SEED)
    cert = core.integer_blocks()
    references = core.separated_standards(np.pi / 128)
    free_control = references.pop("uniform_blank_0")
    standards = {k: calibration(rng, z) for k, z in references.items()}
    result = dict(
        starts=SHOTS,
        seed=SEED,
        bins=BINS,
        step=STEP,
        grid=GRID,
        certificate=cert,
        separated_energy_standards=standards,
        boards={},
    )
    for level in range(3):
        board, seam, chart = core.chart(level)
        data = dict(chart=chart, terminal=core.terminal_reader(board, seam), panels=[])
        result["boards"][str(board.cells)] = data
        if level == 2:
            continue
        for kind, k in (
            ("occupied", 8),
            ("occupied", 0),
            ("vacuum_write", 8),
            ("blank", 8),
        ):
            for odd in (False, True):
                data["panels"].append(panel(board, seam, k, kind, odd, rng))
    standards["uniform_blank_0"] = calibration(rng, free_control)
    return result


def verify(actual, archive):
    if isinstance(archive, dict):
        if actual.keys() != archive.keys():
            raise ValueError("Receipt keys differ")
        for key in archive:
            verify(actual[key], archive[key])
    elif isinstance(archive, list):
        for a, b in zip(actual, archive, strict=True):
            verify(a, b)
    elif isinstance(archive, float):
        if not np.isclose(actual, archive, rtol=1e-8, atol=5e-10):
            raise ValueError("Computed observable differs")
    elif actual != archive:
        raise ValueError("Count or setting differs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    target = parser.parse_args().output
    results = run()
    if target:
        payload = json.dumps(results, allow_nan=False, separators=(",", ":")).encode()
        with target.open("xb") as fh:
            fh.write(gzip.compress(payload, mtime=0))
    else:
        verify(results, json.loads(gzip.decompress(OUTPUT.read_bytes())))
        print("Occupied-core clock and position receipts verified")
