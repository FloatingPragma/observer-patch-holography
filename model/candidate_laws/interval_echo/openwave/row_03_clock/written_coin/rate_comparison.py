"""Compare two clock species under a change of contrast using only recorded counts."""

import json
import math
from pathlib import Path

import numpy as np
from model.experiments.count_intervals import clopper_pearson

HERE = Path(__file__).resolve().parent
FREQUENCY_LIMIT = 46.0  # Upper search limit of the original clock protocol.
CONFIDENCE = 0.99


def frequency_intervals(times, flip):
    """Invert sin²(wt/2) over all phase branches in the enrolled frequency range."""
    t_low, t_high = times
    a, b = (2 * math.asin(math.sqrt(p)) for p in flip)
    period = 2 * math.pi
    possible = []
    for cycle in range(math.floor(t_high * FREQUENCY_LIMIT / period) + 1):
        for left, right in ((a, b), (period - b, period - a)):
            low = (cycle * period + left) / t_high
            high = (cycle * period + right) / t_low if t_low > 0 else math.inf
            if low <= FREQUENCY_LIMIT:
                possible.append((low, min(FREQUENCY_LIMIT, high)))
    return possible


def merge_overlaps(intervals):
    merged = []
    for low, high in sorted(intervals):
        if merged and low <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], high))
        else:
            merged.append((low, high))
    return merged


def frequency_from_records(rows, starts, timer_starts, tail):
    admitted = [(0.0, FREQUENCY_LIMIT)]
    for timer, same, flip, blank in rows:
        if same + flip + blank != starts:
            raise ValueError("Clock outcomes must account for every start")
        timer_band = clopper_pearson(timer, timer_starts, tail)
        times = [16 * math.asin(math.sqrt(p)) for p in timer_band]
        gate = frequency_intervals(times, clopper_pearson(flip, starts, tail))
        admitted = merge_overlaps(
            (max(low, a), min(high, b))
            for low, high in admitted
            for a, b in gate
            if max(low, a) <= min(high, b)
        )
    if not admitted:
        raise ValueError("The record excludes a single free sinusoidal clock")
    # Take the hull so reported tests retain every unresolved phase branch.
    return np.array([admitted[0][0], admitted[-1][1]])


def compare(records):
    panels = records["panels"]
    samples = {
        key: panels[key]["counts"][selection]
        for key, selection in (
            ("8_even", slice(0, 97)),
            ("8_odd", slice(0, 97)),
            ("4_even", slice(0, 97)),
            ("4_odd", slice(97, None)),
        )
    }
    # Each appointment contributes two two-sided binomial intervals.
    tail = (1 - CONFIDENCE) / (4 * sum(map(len, samples.values())))
    bands = {
        key: frequency_from_records(
            rows, records["starts"], records["timer_starts"], tail
        )
        for key, rows in samples.items()
    }
    native_even, native_odd = bands["8_even"], bands["8_odd"]
    target_even, target_odd = bands["4_even"], bands["4_odd"]
    common = native_odd * target_even / native_even[::-1]
    additive = native_odd + target_even - native_even[::-1]

    def overlaps(left, right):
        return bool(max(left[0], right[0]) <= min(left[1], right[1]))

    return {
        "joint_coverage": CONFIDENCE,
        "appointments": {k: len(v) for k, v in samples.items()},
        "frequency_bands": {k: v.tolist() for k, v in bands.items()},
        "predicted_4_odd_common_rate": common.tolist(),
        "predicted_4_odd_equal_shift": additive.tolist(),
        "common_rate_compatible": overlaps(common, target_odd),
        "equal_shift_compatible": overlaps(additive, target_odd),
    }


if __name__ == "__main__":
    result = compare(json.loads((HERE / "records.json").read_text()))
    print(json.dumps(result, indent=2))
