"""Infer the five directional energies solely from labeled clicks and misses."""

import json
from math import log, sqrt
from pathlib import Path

DIRECTIONS = ("zero", "port", "edge", "face", "oblique")
ENROLLMENT = 2**34
TOLERANCE = sqrt(2 * log(2 * 302 / 0.01) / ENROLLMENT)
HALF_BIN = 0.5 / 2**32


def read_bank(batches):
    if len(batches) != 8:
        raise ValueError("Missing batch")
    counts = [0] * 4
    for batch in batches:
        if len(batch) != 4 or any(len(arm) != 2 for arm in batch):
            raise ValueError("Missing arm or outcome")
        if any(type(n) is not int or n < 0 for arm in batch for n in arm):
            raise ValueError("Invalid count")
        if sum(sum(arm) for arm in batch) != ENROLLMENT // 8:
            raise ValueError("Incorrect batch enrollment")
        for arm, (click, _) in enumerate(batch):
            counts[arm] += click

    dark, bright, signal, half = counts
    span = bright - dark
    allowance = ENROLLMENT * TOLERANCE
    half_passes = abs(half - dark - span / 2) <= allowance
    if span <= 4 * allowance:
        return "unresolved", (0.0, 1.0), half_passes
    lo = max(0.0, (signal - dark - allowance) / span - HALF_BIN)
    hi = min(1.0, (signal - dark + allowance) / span + HALF_BIN)
    return ("resolved" if lo <= hi else "inconsistent"), (lo, hi), half_passes


def report(records):
    if (
        records["arms"] != ["dark", "bright", "signal", "half"]
        or records["outcomes"] != ["click", "miss"]
        or records["events_per_bank"] != ENROLLMENT
    ):
        raise ValueError("Unrecognized record format")
    enrolled = {(d, e, False) for d in DIRECTIONS for e in range(30)}
    enrolled.add(("port", 0, True))
    found = {}
    for bank in records["banks"]:
        key = (bank["axis"], bank["seam"], bank["blind"])
        if key in found or key not in enrolled:
            raise ValueError("Duplicate or unenrolled bank")
        found[key] = read_bank(bank["batches"])
    if found.keys() != enrolled:
        raise ValueError("Incomplete receiver enrollment")

    energies = {}
    for direction in DIRECTIONS:
        results = [found[direction, i, False] for i in range(30)]
        if any(status == "inconsistent" or not passed for status, _, passed in results):
            energies[direction] = None
        else:
            energies[direction] = [
                sum(bounds[k] for _, bounds, _ in results) for k in (0, 1)
            ]
    lower = energies["face"]
    upper = energies["port"]
    contrast = (
        None
        if lower is None or upper is None
        else [lower[0] - upper[1], lower[1] - upper[0]]
    )
    return dict(
        energies=energies,
        face_minus_port=contrast,
        blind=found["port", 0, True],
        half_checks_passing=sum(passed for _, _, passed in found.values()),
    )


if __name__ == "__main__":
    print(
        json.dumps(
            report(json.loads(Path(__file__).with_name("records.json").read_text())),
            indent=2,
        )
    )
