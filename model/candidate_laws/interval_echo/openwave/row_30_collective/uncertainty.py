"""Perturb the held-out shape forecasts with a common bootstrap of the bench meters."""

import argparse
import json
from pathlib import Path

import numpy as np

from .coupling_meter import HAND_TRIALS, TRIALS, measure
from .experiment import HERE, verify


def run():
    meter = json.loads((HERE / "coupling_meter.json").read_text())
    probs = np.array(meter["all_counts"]) / TRIALS
    cal = np.array(meter["hand_calibration"]) / HAND_TRIALS
    rng = np.random.default_rng(371093)
    samples = np.array(
        [
            measure(
                rng.multinomial(TRIALS, probs),
                rng.binomial(HAND_TRIALS, cal) / HAND_TRIALS,
            )[0][2:4]
            for _ in range(400)
        ]
    )
    mass, speed2 = samples.T
    boards = []
    for b in json.loads((HERE / "holdout.json").read_text())["boards"]:
        k = np.array([r["k"] for r in b["observer"]["modes"]])
        observed = np.array([r["omega"] for r in b["observer"]["modes"]])
        q2 = 2 - 2 * np.cos(k)
        lattice = np.sqrt(mass[:, None] ** 2 + speed2[:, None] * q2)
        par = mass[:, None] + speed2[:, None] * q2 / (2 * mass[:, None])

        def rms(model):
            err = model[:, 1:] - model[:, [0]] - (observed[1:] - observed[0])
            return np.sqrt(np.mean(err**2, axis=1))

        ratio = rms(lattice) / rms(par)
        absolute = np.max(abs(lattice - observed), axis=1)
        boards.append(
            dict(
                level=b["level"],
                dimers=b["dimers"],
                rms_ratio_percentile95=np.quantile(ratio, [0.025, 0.975]).tolist(),
                fraction_below_half=float(np.mean(ratio < 0.5)),
                maximum_absolute_error_percentile95=np.quantile(
                    absolute, [0.025, 0.975]
                ).tolist(),
            )
        )
    return dict(
        bootstrap_replicates=len(samples),
        seed=371093,
        interpretation="Conditional on the held observer records; common bench clock and both measured couplings resampled jointly.",
        boards=boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run()
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
    else:
        verify(out, json.loads((HERE / "uncertainty.json").read_text()))
    print(json.dumps(out, indent=2))
