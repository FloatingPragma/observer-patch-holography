"""Finite launch-count predictions of separately surveyed four-station volumes."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as cert
from . import force_forecast, metric_forecast

HERE = Path(__file__).resolve().parent


def volume_band(intervals):
    # log det is concave. Its minimum on the squared-length box is
    # attained at a corner; its supporting plane at the center bounds the maximum.
    squared = [
        (sp.Rational(str(lo)) ** 2, sp.Rational(str(hi)) ** 2) for lo, hi in intervals
    ]
    corners = [metric_forecast.geometry.gram(x) for x in product(*squared)]
    for g in corners:
        assert cert.global_survey.positive_definite(g)
    lower = cert.iv.sqrt(cert.ri(min(g.det() for g in corners))) / 6
    center = metric_forecast.geometry.gram([(lo + hi) / 2 for lo, hi in squared])
    inverse = center.inv()
    excursion = max((inverse * (g - center)).trace() for g in corners)
    upper = cert.iv.sqrt(cert.ri(center.det()) * cert.iv.exp(cert.ri(excursion))) / 6
    return [
        float(np.nextafter(float(lower.a), -np.inf)),
        float(np.nextafter(float(upper.b), np.inf)),
    ]


def run():
    cert.iv.dps = 32
    frozen = json.loads((HERE / "metric_forecast.json").read_text())
    effects, _ = metric_forecast.direct_effects(frozen["joint_events"])
    trials = 2**72
    batch = 2**58
    rng = np.random.default_rng(191327)
    archive = {}
    for name, state in force_forecast.preparations().items():
        probabilities = np.einsum("qkbij,bji->qk", effects, state).real
        held = frozen["held_population_comparisons"][name]
        arms = []
        for q, (p, panel) in enumerate(
            zip(probabilities, frozen["absolute_readers"], strict=True)
        ):
            rates = np.r_[p, 1 - p.sum()]
            assert min(rates) >= 0
            rows = rng.multinomial(batch, rates, size=trials // batch)
            counts = [int(v) for v in rows.astype(object).sum(axis=0)]
            assert sum(counts) == trials
            frequency = np.array([v / trials for v in counts[:-1]])
            intervals = []
            for k, reader in enumerate(panel):
                unit = int(reader["denominator"])
                w = np.array([int(v) / unit for v in reader["score_numerators"]])
                intercept = int(reader["intercept_numerator"]) / unit
                estimate = float(intercept + w @ frequency)
                width = np.ptp(np.r_[w, 0]) * np.sqrt(
                    np.log(2 * 36 / 0.01) / (2 * trials)
                )
                # Correct for the independent numerical simulator's probability
                # error and three-way comparison with the separately surveyed target.
                width += (
                    held["float_probability_and_roundoff_budgets"][q][k]
                    + reader["uniform_error_upper"]
                    + 1e-6
                )
                intervals.append(
                    [
                        float(np.nextafter(estimate - width, -np.inf)),
                        float(np.nextafter(estimate + width, np.inf)),
                    ]
                )
            band = volume_band(intervals)
            assert band[0] < held["separate_volume"] < band[1]
            arms.append(
                dict(
                    charge=frozen["species"][q],
                    counts=counts,
                    edge_bands=intervals,
                    volume_band=band,
                )
            )
        archive[name] = arms
    for q in range(2):
        assert (
            archive["dephased_Xplus_Yplus"][q]["volume_band"][1]
            < archive["kept_Zplus_Zminus"][q]["volume_band"][0]
        )
    return dict(
        synthetic_trials_per_preparation_per_probe=str(trials),
        batch_size=batch,
        simultaneous_coverage=0.99,
        labels=frozen["joint_events"] + [["other"]],
        records=archive,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(content)
    else:
        print(content, end="")
