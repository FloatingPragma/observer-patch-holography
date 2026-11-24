"""Infer the rest beat from receiver and hand counts, without a time or energy label."""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def estimate(records):
    # The four phase settings use separate copies, ordered +1, +i, -1, -i.
    # The last axis is (hand0 quiet, hand0 hit, hand1 quiet, hand1 hit).
    totals = records.sum(axis=-1)
    hits = records[..., [1, 3]].sum(axis=-1) / totals
    hand = records[..., 2:].sum(axis=(-1, -3)) / totals.sum(axis=-2)
    times = 2 * np.arcsin(np.sqrt(hand))
    amplitude = (hits[..., 0, :] - hits[..., 2, :]) + 1j * (
        hits[..., 1, :] - hits[..., 3, :]
    )
    phase = np.unwrap(np.angle(amplitude), axis=-1)
    centered = times - times.mean(axis=-1, keepdims=True)
    return (centered * phase).sum(axis=-1) / (centered**2).sum(axis=-1)


def decode(panel, control, draws=1024):
    groups = {}
    for row in panel["interferometry"]:
        if row["internal_exchange"] == control:
            groups.setdefault(tuple(row["phase"]), []).append(row["counts"])
    records = np.array([groups[key] for key in [(1, 0), (0, 1), (-1, 0), (0, -1)]])
    # Each resampled trial stays in the same joint hand/receiver outcome space.
    rng = np.random.default_rng([303441, panel["cells"], control])
    samples = np.empty((draws, *records.shape), dtype=int)
    for index in np.ndindex(*records.shape[:2]):
        row = records[index]
        samples[(slice(None), *index)] = rng.multinomial(
            row.sum(), row / row.sum(), draws
        )
    interval = np.quantile(estimate(samples), [0.025, 0.975])
    return float(estimate(records)), interval.tolist()


if __name__ == "__main__":
    for panel in json.loads((HERE / "results.json").read_text())["panels"]:
        for control in (1, 2):
            print(panel["cells"], control, decode(panel, control))
