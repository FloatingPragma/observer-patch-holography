"""Acquire and replay all programmed endpoint words and keep every start."""

import argparse
from collections import Counter
import gzip
import json
from math import ceil, log, sqrt
from pathlib import Path

import numpy as np

from ...signed_effects.sequential import grid_counts
from . import apparatus as a
from .classify import partitions


HERE = Path(__file__).resolve().parent
SHOTS = 2**16
SEED = 261130
RECEIPT = HERE / "records.json.gz"


def menus():
    for level in range(3):
        cells = a.geometry(level)["cells"]
        for menu, writer in a.writers().items():
            values = a.level_probabilities(writer)
            for j, pair in enumerate(a.INPUTS):
                yield f"{cells}/{menu}/{pair[0]}{pair[1]}", values[:, :, j]
    for kappa in (8, 0):
        for arm in (False, True):
            for t, values in a.free_wait(kappa, arm).items():
                for j, pair in enumerate(a.INPUTS):
                    yield (
                        f"wait/{kappa}/{int(arm)}/{t + 1}/{pair[0]}{pair[1]}",
                        values[:, :, j],
                    )
    for i, values in enumerate(a.clock()):
        yield f"clock/{i}", values
    for turns in (1, 3):
        for sign in (False, True):
            yield f"phase/{turns}/{int(sign)}", a.phase_return(turns * np.pi / 2, sign)


def acquire():
    a.verify_marginal_channels()
    rng = np.random.default_rng(SEED)
    banks = {}
    for name, levels in menus():
        probabilities = a.read_distribution(levels)
        if name in banks:
            raise ValueError("duplicated detector bank")
        draws = grid_counts(rng, SHOTS, probabilities)
        banks[name] = {
            "counts": draws.tolist(),
            # Record probabilities rounded for display/replay; decisions include rounding margin.
            "probabilities": [round(float(p), 11) for p in probabilities.ravel()],
        }
    return {
        "boards": [a.geometry(i) for i in range(3)],
        "words": list(a.WORDS),
        "calibrated_load_words": list(a.detector_words()),
        "shots": SHOTS,
        "seed": SEED,
        "banks": banks,
        "source_partitions": partitions(),
    }


def mismatch(counts):
    grid = np.asarray(counts).reshape(4, 4)
    return grid.sum() - np.trace(grid)


def expected_names():
    pairs = [f"{i}{j}" for i, j in a.INPUTS]
    immediate = {
        f"{size}/{menu}/{pair}"
        for size in (20, 80, 320)
        for menu in ("repair", "unchanged_loads", "parked")
        for pair in pairs
    }
    waits = {
        f"wait/{k}/{arm}/{t}/{pair}"
        for k in (0, 8)
        for arm in (0, 1)
        for t in (1, 2)
        for pair in pairs
    }
    return (
        immediate
        | waits
        | {f"clock/{i}" for i in range(3)}
        | {f"phase/{t}/{s}" for t in (1, 3) for s in (0, 1)}
    )


def decode(data):
    banks = data["banks"]
    if data["shots"] != SHOTS or set(banks) != expected_names():
        raise ValueError("Missing or unregistered acquisition")
    if data["words"] != ["00", "01", "10", "11"] or data["calibrated_load_words"] != [
        "01",
        "00",
        "10",
    ]:
        raise ValueError("Unknown detector-word calibration")
    n = data["shots"]
    radius = ceil(1e9 * (sqrt(log(2 * 17 * len(banks) / 0.01) / (2 * n)) + 1e-10)) / 1e9
    reads = {}
    for name, b in banks.items():
        c = np.asarray(b["counts"])
        if (
            c.shape != (16,)
            or c.dtype.kind not in "iu"
            or np.any(c < 0)
            or c.sum() != n
        ):
            raise ValueError(f"Missing or malformed outcomes in {name}")
        reads[name] = c.reshape(4, 4)
    # Judge the incoming label forecasts from the terminal words alone.
    for size in (20, 80, 320):
        where = (1, 0, 2)
        for j, pair in enumerate(a.INPUTS):
            low, high = a.possible(pair)
            forecast = np.zeros((4, 4))
            forecast[where[low], where[high]] += 0.5
            forecast[where[high], where[low]] += 0.5
            obs = reads[f"{size}/repair/{pair[0]}{pair[1]}"]
            if abs(obs / n - forecast).max() > radius or mismatch(obs) != n * (
                low != high
            ):
                raise ValueError(
                    f"Native repair forecast failed for board {size}, input {j}"
                )
    return reads, radius


def summarize(data):
    banks = data["banks"]
    n = data["shots"]
    reads, radius = decode(data)
    for name, b in banks.items():
        c = reads[name].ravel()
        p = np.asarray(b["probabilities"])
        if max(abs(c / n - p).max(), abs(mismatch(c) / n - mismatch(p))) > radius:
            raise ValueError(f"Forecast outside enrolled interval in {name}")
    hist = Counter(len(parts) for _, parts in data["source_partitions"])
    return {
        "banks": len(banks),
        "starts": n * len(banks),
        "joint_99_half_width": radius,
        "supports_by_number_of_sharp_labels": dict(sorted(hist.items())),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = acquire()
    summary = summarize(data)
    if args.write:
        RECEIPT.write_bytes(
            gzip.compress(
                (json.dumps(data, separators=(",", ":")) + "\n").encode(), mtime=0
            )
        )
    elif json.loads(gzip.decompress(RECEIPT.read_bytes())) != data:
        raise SystemExit("Saved enrollment, probabilities or counts differ")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
