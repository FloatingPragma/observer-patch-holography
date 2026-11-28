"""Compare two spatial windings in the same nine pairs of binary analyzer banks."""

import argparse
from hashlib import sha256
from itertools import product
import json
from math import log, sqrt
from pathlib import Path
from random import Random

import numpy as np

from ....fan_recoil.self_start.measure import DRAW_PAD, draw
from ....replay import compare
from ..measure import PAULI


HERE = Path(__file__).resolve().parent
TRIALS = 2**36
VECTORS = [[-6 + 4j, -2 + 7j, 5j, 8], [8, 3 + 5j, 4 + 5j, -5 - 2j]]


def coefficients(v):
    v = np.array(v)
    return np.array(
        [[np.vdot(v, np.kron(a, b.T) @ v).real for b in PAULI] for a in PAULI]
    ) / (4 * np.vdot(v, v).real)


def run():
    file = HERE / "records.json"
    receipt = json.loads(file.read_text())
    settings = [(0, 0.5), (1, 0.5), (1, 0)]
    coeff = [coefficients(v) for v in VECTORS]
    x, y = np.array(list(product((1, -1), repeat=2))).T
    rng = Random(10451053)
    cases = []
    for f, g in settings:
        flight = receipt["flights"][f]
        case = next(
            c
            for c in flight["counter_readings"]
            if c["ring_rate"] == g and c["wait_in_pi_over_8"] == 8
        )
        t = [
            np.array(case[label]["return_word_source_pauli_moments_per_original"])
            for label in ("positive_return", "negative_return")
        ]
        p = [mom[0, 0] for mom in t]
        pre = [
            flight[k]
            for k in (
                "rejected_before_terminal_fan",
                "terminal_fan_away_from_hub",
                "wrong_hub_word",
            )
        ]
        other = flight["passed_on_original_launch"] - sum(p)
        score = [c[0, 0] for c in coeff]
        ranges = np.zeros(2)
        banks = []
        for i, j in product(range(1, 4), repeat=2):
            outcomes = [
                (mom[0, 0] + x * mom[i, 0] + y * mom[0, j] + x * y * mom[i, j]) / 4
                for mom in t
            ]
            counts = draw([*outcomes[0], *outcomes[1], *pre, other], TRIALS, rng)
            assert sum(counts) == TRIALS
            for m, c in enumerate(coeff):
                bins = counts[4 * m : 4 * m + 4]
                weights = c[i, j] * x * y
                if i == j:
                    weights = weights + c[i, 0] * x + c[0, j] * y
                score[m] += bins @ weights / sum(bins)
                ranges[m] += np.ptp(weights) ** 2 / sum(bins)
            banks.append(
                {"receiver_source_axes": "IXYZ"[i] + "IXYZ"[j], "all_counts": counts}
            )
        witnesses = []
        for m, c in enumerate(coeff):
            # Bounds the unnormalized trace-distance error as well as the
            # denominator error. The grid changes four individual label bins.
            pad = (
                4 * flight["input_pauli_error_upper"]
                + 2 * receipt["counter_transfer_operator_error_upper"]
            )
            sys = (2 * pad + 16 * DRAW_PAD * (1 + abs(c).sum())) / (
                p[m] - pad - 4 * DRAW_PAD
            )
            radius = sqrt(log(2 * 6 / 0.01) * ranges[m] / 2) + sys
            interval = [score[m] - radius, score[m] + radius]
            expected = (t[m] * c).sum() / p[m]
            assert interval[0] < expected < interval[1]
            assert (interval[1] < 0) if (f, g) == (1, 0.5) else (interval[0] > 0)
            witnesses.append(
                {
                    "predicted_conditional_score": expected,
                    "score_from_counts": score[m],
                    "simultaneous_99pct_interval": interval,
                    "numerical_allowance": sys,
                }
            )
        cases.append(
            {
                "source_flight": flight["source_flight"],
                "ring_rate": g,
                "banks": banks,
                "positive_negative_witnesses": witnesses,
            }
        )
    return {
        "count_protocol_sha256": sha256(
            (HERE / "COUNT_TEST.md").read_bytes()
        ).hexdigest(),
        "material_receipt_sha256": sha256(file.read_bytes()).hexdigest(),
        "simulated_original_launches_per_setting": TRIALS,
        "count_order": [
            "+1:++",
            "+1:+-",
            "+1:-+",
            "+1:--",
            "-1:++",
            "-1:+-",
            "-1:-+",
            "-1:--",
            "outside_fan",
            "outside_hub",
            "other_word",
            "other_material_result",
        ],
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(data, indent=1) + "\n")
    else:
        compare(data, json.loads((HERE / "counts.json").read_text()), "counted_return")
        print(
            "Verified both single-circuit count witnesses and all four control intervals."
        )


if __name__ == "__main__":
    main()
