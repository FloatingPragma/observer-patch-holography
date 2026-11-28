"""Repeat the independent survey on the same six directed routes and two clocks."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json

import numpy as np

from . import measure

HERE = measure.HERE
RECEIPT = HERE / "radar_precision.json"
SHOTS = 2**80
SEED = 190826


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def acquire(primary):
    rng = np.random.default_rng(SEED)
    records = {}
    for probe in measure.PROBES:
        records[probe] = {}
        for stock in measure.STOCKS:
            p = primary["banks"]["8"]["population"][probe][stock]
            clocks = [measure.camera.precision_draw(rng, SHOTS, v) for v in p["clock"]]
            ranges = []
            for row in p["range"]:
                probs = np.array(row[:-1])[:, None]
                # Sampling smaller blocks limits temporary memory. Every
                # independent virtual start contributes to exactly one bin.
                bins = [0] * len(row)
                block = 2**76
                for _ in range(SHOTS // block):
                    draw = measure.camera.precision_draw(rng, block, probs)[0]
                    bins = [n + m for n, m in zip(bins, draw, strict=True)]
                ranges.append(bins)
            records[probe][stock] = {"clock": clocks, "range": ranges}
    data = {
        "seed": SEED,
        "shots_per_ruler_setting": SHOTS,
        "counts": records,
        "protocol_sha256": digest(HERE / "RADAR_PRECISION.md"),
        "primary_sha256": digest(measure.RECEIPT),
    }
    data["decoded"] = decode(primary, data)
    return data


def decode(primary, data):
    if data["seed"] != SEED or data["shots_per_ruler_setting"] != SHOTS:
        raise ValueError("Unregistered radar acquisition")
    if data["protocol_sha256"] != digest(HERE / "RADAR_PRECISION.md") or data[
        "primary_sha256"
    ] != digest(measure.RECEIPT):
        raise ValueError("Radar precision inputs differ from archive")
    counts = deepcopy(primary["banks"]["8"]["counts"])
    for probe, stocks in counts.items():
        for stock, panel in stocks.items():
            new = data["counts"][probe][stock]
            if set(new) != {"clock", "range"}:
                raise ValueError("Precision repeat changes an unregistered instrument")
            panel.update(new)
    return measure.decode(counts, primary["half_steps"], True, ruler_shots=SHOTS)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = p.parse_args()
    primary = json.loads(measure.RECEIPT.read_text())
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        if decode(primary, data) != data["decoded"]:
            raise ValueError("Decoded measures differ from archive")
    else:
        data = acquire(primary)
        if args.write:
            with RECEIPT.open("x") as handle:
                json.dump(data, handle, indent=2, allow_nan=False)
                handle.write("\n")
        elif data != json.loads(RECEIPT.read_text()):
            raise ValueError("Radar acquisition differs from archive")
    print(json.dumps(data["decoded"]["tests"]["radar"], indent=2))


if __name__ == "__main__":
    main()
