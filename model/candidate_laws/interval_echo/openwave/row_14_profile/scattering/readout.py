"""Calibrated radial signs and central port populations from detector counts."""

import json
from pathlib import Path

import numpy as np

STARTS = 2**40
GRID = 2**40
COMPARISONS = 1000


def score_interval(counts, scores):
    """Hoeffding interval for the measured plus-minus score on admitted starts."""
    n = sum(counts)
    if n == 0:
        return [-1.0, 1.0]
    est = sum(c * v for c, v in zip(counts, scores, strict=True)) / n
    stat = np.sqrt(2 * np.log(2 * COMPARISONS / 0.01) / n)
    # Multiplicative Chernoff inversion gives p_grid >= n/(2*STARTS)
    # except with chance exp[-n(log 2 - 1/2)] at this fixed setting.
    if n * (np.log(2) - 0.5) < np.log(COMPARISONS / 0.01):
        return [-1.0, 1.0]
    acceptance_lower = n / (2 * STARTS)
    err = stat + 2 * (7 / GRID + 2e-11) / acceptance_lower
    return [est - err, est + err]


def decode(bank, standard, size, topology=None):
    bins = np.array(bank["counts"][:-1]).reshape(2, 3)
    counts = bins.sum(axis=0) if topology is None else bins[topology]
    signed = score_interval(counts.tolist(), (1, -1, 0))
    calibration_interval = score_interval(standard["counts"][:-1], (1, -1, 0))
    if calibration_interval[0] <= 0:
        return None
    ends = [size * a / b for a in signed for b in calibration_interval]
    return [min(ends), max(ends)]


def center_composition(fringe, elastic, standard, topology):
    """Separate the two charge signs using first and second moments at one site."""
    signed = decode(fringe, standard, 1, topology)
    unit_power = score_interval(standard["counts"][:-1], (1, 1, -1))
    if unit_power[0] <= 0 or signed is None:
        return None
    bins = np.array(elastic["counts"][:-1]).reshape(2, 2)
    counts = bins.sum(axis=0) if topology is None else bins[topology]
    click = [(1 + side) / 2 for side in score_interval(counts.tolist(), (1, -1))]
    second = [a / b for a in click for b in unit_power]
    lo, hi = min(second), max(second)
    return {
        "excess": [(lo + signed[0]) / 2, (hi + signed[1]) / 2],
        "vacancy": [(lo - signed[1]) / 2, (hi - signed[0]) / 2],
    }


if __name__ == "__main__":
    result = json.loads((Path(__file__).parent / "results.json").read_text())
    calibration = result["standards"]["native"]["2"]
    for source, panel in result["panels"].items():
        for center in ("intrinsic", "gateway_1"):
            banks = panel["settings"][f"{center}/center"]
            for name, group in (("all admitted", None), ("connected", 1)):
                estimates = center_composition(
                    banks["native_fringe"], banks["native_elastic"], calibration, group
                )
                print(source, center, name, json.dumps(estimates))
