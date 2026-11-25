"""Predict one-way terminal records using only successful two-leg echo counts."""
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta


def limits(counts, shots, tail):
    counts = np.asarray(counts)
    lower = np.zeros(counts.shape)
    upper = np.ones(counts.shape)
    pos = counts > 0
    below = counts < shots
    lower[pos] = beta.ppf(tail, counts[pos], shots - counts[pos] + 1)
    upper[below] = beta.isf(tail, counts[below] + 1, shots - counts[below])
    return lower, upper


def check(data):
    n = data["shots"]
    settings = {(cells, gain, t, bit) for cells in (20, 320) for gain in (1, 0.5)
                for t in (0, 0.25, 1, 2, 4) for bit in (0, 1)}
    if n != 2**22 or len(data["panels"]) != len(settings):
        raise ValueError("incomplete calibration")
    # 12 resolved echoes, 13 challenge bins, two aggregate echoes, one local arrival.
    tail = 0.01 / (2 * len(settings) * 28)
    scores = {}
    misses = [0, 0]
    checked = 0
    for row in data["panels"]:
        key = (row["cells"], row["gain"], row["exposure"], row["bit"])
        settings.remove(key)
        banks = [np.array(row[k]) for k in ("resolved", "coarse", "challenge")]
        for bank, size in zip(banks, (25, 3, 13), strict=True):
            if bank.shape != (size,) or bank.dtype.kind not in "iu" or bank.min() < 0 or bank.sum() != n:
                raise ValueError("incomplete enrolled bank")
        sharp, coarse, direct = banks
        echo = sharp[:12]
        lo, hi = np.sqrt(limits(echo, n, tail))
        lo, hi = np.r_[lo, max(0., 1 - hi.sum())], np.r_[hi, min(1., 1 - lo.sum())]
        lower, upper = limits(direct, n, tail)
        if (lo > upper).any() or (hi < lower).any():
            raise ValueError(f"echo forecast misses a one-way detector bin on {key}")
        checked += direct.size
        # A coarse return preserves phase between the possible midpoints. Ignoring
        # recorded midpoint labels performs a different experiment.
        _, echo_upper = limits([coarse[0], echo.sum()], n, tail)
        arrival_lower, _ = limits([direct[:12].sum()], n, tail)
        failures = np.sqrt(echo_upper) < arrival_lower[0]
        misses = [a + int(b) for a, b in zip(misses, failures, strict=True)]
        estimate = np.sqrt(echo / n)
        quiet = 1 - estimate.sum()
        if quiet < 0:
            raise ValueError("echo point estimate exceeds a complete probability law")
        scores.setdefault(key[:3], {})[key[3]] = (np.r_[estimate, quiet], direct)
        if key == (20, 1, 4, 0):
            print("20 cells, bit zero, native exposure 4")
            print("Direct local arrival, resolved echo forecast, coarse square root, unlabelled square root:",
                  direct[:12].sum() / n, estimate.sum(), *np.sqrt(np.array([coarse[0], echo.sum()]) / n))
    for key, terms in scores.items():
        zero, one = terms[0][0], terms[1][0]
        table = one > zero
        hits = terms[0][1][~table].sum() + terms[1][1][table].sum()
        if key[0] == 20 and key[2] == 4:
            print(f"20 cells, gain {key[1]}, exposure 4, correct decisions: {hits}/{2*n}")
    print(f"{checked} one-way bin intervals overlap their echo forecasts.")
    print("Coarse and unlabelled square-root forecasts excluded in", *misses, "of 40 panels respectively.")


if __name__ == "__main__":
    check(json.loads(Path(__file__).with_name("records.json").read_text()))
