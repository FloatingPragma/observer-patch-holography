"""Confidence intervals from the archived integer records, independent of the law."""

import json
from pathlib import Path

import numpy as np

from model.experiments.count_intervals import clopper_pearson

HERE = Path(__file__).resolve().parent
TAIL = 0.01 / (2 * 2500)


def interval(values, counts, low, high):
    """E[X] = low + sum of the probabilities P[X >= j] for integer X."""
    if min(values) < low or max(values) > high:
        raise ValueError("Count outside instrument range")
    tails = [int(counts[values >= j].sum()) for j in range(low + 1, high + 1)]
    bounds = clopper_pearson(tails, int(counts.sum()), TAIL)
    return [float(low + side.sum()) for side in bounds]


def read(profile):
    counts = np.array(profile["counts"])
    labels = np.array(profile["keys"])
    if int(counts.sum()) != profile["accepted"]:
        raise ValueError("Profile counts do not exhaust accepted shots")
    return dict(
        shells=[
            interval(labels[:, i], counts, *([-1, 1] if i == 0 else [-2, 2]))
            for i in range(4)
        ],
        distant=interval(labels[:, 2] + labels[:, 3], counts, -2, 2),
        connected=interval(labels[:, 4], counts, 0, 1),
    )


def run(records):
    answer = {}
    for name, sectors in records["panels"].items():
        answer[name] = {}
        for tag, panel in sectors.items():
            profiles = panel if tag == "unfiltered_profile" else panel["profile"]
            answer[name][tag] = {key: read(row) for key, row in profiles.items()}
    return answer


if __name__ == "__main__":
    for branch, panels in run(json.loads((HERE / "records.json").read_text())).items():
        for mode, labels in panels.items():
            for centre, result in labels.items():
                print(branch, mode, centre, result)
