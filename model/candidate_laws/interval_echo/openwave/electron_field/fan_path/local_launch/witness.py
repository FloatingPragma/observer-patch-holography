"""Translate a local entanglement witness into nine pairs of binary detector records."""

import argparse
from hashlib import sha256
from itertools import product
import json
from math import log, sqrt
from pathlib import Path

import numpy as np

from .measure import HERE, PAULI

TRIALS = 2**23
LABEL = [1045, 0, 334, 365]
SETTINGS = [
    (4, 27.0625, False),
    (4, 30.8125, False),
    (2, 30.8125, False),
    (4, 30.8125, True),
]


def run():
    source = HERE / "word_reads.json"
    data = json.loads(source.read_text())
    flight = next(
        s["flights"]
        for s in data["cases"][-1]["stock_arms"]
        if s["preparation"] == "clock_bypassed"
    )
    v = np.array([2, 1, 1, -1 - 1j])
    coef = np.array(
        [[np.vdot(v, np.kron(a, b.T) @ v).real for b in PAULI] for a in PAULI]
    )
    np.testing.assert_array_equal(
        coef, [[8, 2, 2, 2], [2, -2, 4, 6], [-2, -4, -6, 2], [2, 6, -2, 4]]
    )
    denom = float(4 * np.vdot(v, v).real)
    signs = np.array(list(product((1, -1), repeat=2)))
    x, y = signs.T
    rng = np.random.default_rng(981045)
    cohorts = []
    for rate, wait, erased in SETTINGS:
        record = next(
            f
            for f in flight
            if (f["material_bridge_rate"], f["appointment"]) == (rate, wait)
        )
        pair = next(
            pair
            for pair in record["pair_readings"]
            if pair["port_face_and_words"] == LABEL
        )
        t = np.asarray(pair["receiver_reference_pauli_moments"])
        if erased:
            t[:, 1:3] = 0
        p = t[0, 0]
        exact_mean = float(np.sum(coef * t) / (denom * p))
        bound = record["observable_absolute_error_upper_per_pair"]
        sys_error = (np.abs(coef).sum() / denom + abs(exact_mean)) * bound / (p - bound)
        banks = []
        w = coef[0, 0] / denom
        sum_range = 0.0
        for i, j in product((1, 2, 3), repeat=2):
            probs = (p + x * t[i, 0] + y * t[0, j] + x * y * t[i, j]) / 4
            np.testing.assert_allclose(probs.sum(), p, rtol=0, atol=2e-16)
            assert np.all(probs >= 0)
            counts = rng.multinomial(TRIALS, [*probs, 1 - p])
            assert sum(counts) == TRIALS
            accepted = sum(counts[:4])
            score = coef[i, j] * x * y
            if i == j:
                score = score + coef[i, 0] * x + coef[0, j] * y
            score = score / denom
            w += np.dot(counts[:4], score) / accepted
            sum_range += np.ptp(score) ** 2 / accepted
            banks.append({"axes": "IXYZ"[i] + "IXYZ"[j], "counts": counts.tolist()})
        radius = sqrt(log(2 * len(SETTINGS) / 0.01) * sum_range / 2) + sys_error
        cohorts.append(
            {
                "rate": rate,
                "wait": wait,
                "reference_dephased": erased,
                "acceptance_per_original": float(p),
                "predicted_conditional_score": exact_mean,
                "measured_conditional_score": float(w),
                "simultaneous_99pct_interval": [float(w - radius), float(w + radius)],
                "numerical_score_allowance": float(sys_error),
                "settings": banks,
            }
        )
    return {
        "pre_count_protocol_sha256": sha256(
            (HERE / "WITNESS.md").read_bytes()
        ).hexdigest(),
        "source_receipt_sha256": sha256(source.read_bytes()).hexdigest(),
        "counter_order": ["++", "+-", "-+", "--", "refused"],
        "preparations_per_setting": TRIALS,
        "cohorts": cohorts,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        if result != json.loads((HERE / "witness.json").read_text()):
            raise ValueError("The local analyzer counts or reported bounds changed.")
        for cohort in result["cohorts"]:
            print(
                cohort["rate"],
                cohort["wait"],
                cohort["reference_dephased"],
                cohort["simultaneous_99pct_interval"],
            )


if __name__ == "__main__":
    main()
