"""Resolve all-launch direction contrasts after both loads start in a single owner."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta

from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .precision import HERE, SHOTS, run


def direction_interval(counts, n, tail):
    # For two launch signs, the oriented contrast is (++ plus --) minus (+- plus -+).
    arr = np.asarray(counts)[:, :9].reshape(2, 3, 3).sum(axis=1)
    successes = [arr[0, 0], arr[1, 2], arr[0, 2], arr[1, 0]]
    lows = [beta.ppf(tail, k, n - k + 1) if k else 0.0 for k in successes]
    highs = [beta.isf(tail, k + 1, n - k) if k < n else 1.0 for k in successes]
    return [sum(lows[:2]) - sum(highs[2:]), sum(highs[:2]) - sum(lows[2:])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run((20, 23, 24, 27, 28), seed=291126)
    tail = 0.01 / (2 * 4 * len(data["times"]) * len(data["panels"]))
    data["individual_tail_probability"] = tail
    for panel in data["panels"]:
        for read in panel["reads"]:
            read["exact_binomial_contrast_interval"] = direction_interval(
                read["counts"], SHOTS, tail
            )
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "flight.json").read_text()),
            "local_flight",
        )
        print(
            "Verified independent localized trials and simultaneous binomial intervals."
        )


if __name__ == "__main__":
    main()
