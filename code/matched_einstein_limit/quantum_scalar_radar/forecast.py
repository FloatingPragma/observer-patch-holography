"""Forecast later radar echoes from independent quantum-field and clock records."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .. import observer_field as classical
from . import measure


def rotated_records(preparation, rng):
    source = measure.propagate(preparation)
    unit_frequency = math.sqrt(1 + classical.MASS_SQUARED)
    mapping = source(measure.TIMES[0])[6:18].reshape(3, 2, 2)
    rows = mapping.sum(axis=1) / math.sqrt(2)
    widths = np.sqrt(
        rows[:, 0] ** 2 / unit_frequency + rows[:, 1] ** 2 * unit_frequency
    )
    return [
        [measure.meter(width, mix, rng) for _ in range(2)]
        for width, mix in zip(widths, preparation, strict=True)
    ]


def initialize(tape, rotated):
    start, arms, rates, *_ = classical.local_state(tape["radar"])
    moments = []
    for axis, angles in zip(tape["modes"], rotated, strict=True):
        qq, pp = [
            sum(measure.second_moment(row[key]) for row in axis) for key in ("q", "p")
        ]
        qp = sum(map(measure.second_moment, angles)) - (qq + pp) / 2
        moments.append([qq, qp, pp])
    return start, arms, rates, np.array(moments)


def forecast_tapes(initial, mass2, coupling, vacuum, permute=False):
    start, lengths, rates, moments = initial
    spacing = measure.SPACINGS[-1]
    if permute:
        gradients = moments[:, 0] / lengths**2
        hi, lo = np.argmax(gradients), np.argmin(gradients)
        moments = moments.copy()
        moments[[hi, lo]] = moments[[lo, hi]]

    def slope(_, state):
        arms = state[:3]
        expansion = state[3:6]
        qq, qp, pp = state[6:15].reshape(3, 3).T
        volume = arms.prod()
        frequency_squared = mass2 + arms ** (-2)
        potential = mass2 * qq.sum() / (2 * measure.VOL0)
        directional = qq / (measure.VOL0 * arms**2)
        acceleration = (
            -expansion.sum() * expansion + vacuum + coupling * (potential + directional)
        )
        covariance = np.column_stack(
            (
                2 * qp / volume,
                pp / volume - volume * frequency_squared * qq,
                -2 * volume * frequency_squared * qp,
            )
        )
        return np.r_[arms * expansion, acceleration, covariance.ravel(), 1 / arms]

    solution = solve_ivp(
        slope,
        (start, measure.TIMES[-1] + 5 * spacing),
        np.r_[lengths, rates, moments.ravel(), np.zeros(3)],
        method="DOP853",
        atol=2e-13,
        rtol=2e-13,
        dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    def radar_state(time):
        state = solution.sol(time)
        return np.r_[state[:3], np.zeros(5), state[-3:]]

    return [classical.records(radar_state, t, spacing) for t in measure.TIMES[1:]]


def errors(predictions, observed):
    result = []
    for predicted, truth in zip(predictions, observed, strict=True):
        scale = truth["tick_unit"] / (2 * truth["comoving_arm_spacing"])

        def flight_times(tape):
            return np.diff(
                [r["radar_send_return_ticks"] for r in tape["readings"]], axis=-1
            )

        result.append(
            float(np.max(abs(flight_times(truth) - flight_times(predicted))) * scale)
        )
    return result


def run(reference):
    survey = reference["surveys"][-1]
    result = {}
    for index, name in enumerate(("y", "z", "even_x")):
        angles = rotated_records(
            measure.SOURCES[name], np.random.default_rng([measure.SEED, 2, index])
        )
        tapes = survey["records"][name]
        initial = initialize(tapes[0], angles)
        arms = {}
        for label, permute in (("measured", False), ("permuted", True)):
            predicted = forecast_tapes(
                initial,
                reference["recovered_mass2"],
                **survey["calibration"],
                permute=permute,
            )
            arms[label] = dict(
                predicted=predicted,
                max_errors_by_time=errors(predicted, [t["radar"] for t in tapes[1:]]),
            )
        result[name] = dict(rotated_quadratures=angles, comparisons=arms)
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "FORECAST_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        source_record_sha256=hashlib.sha256(
            (measure.HERE / "results.json").read_bytes()
        ).hexdigest(),
        surveys=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(json.loads((measure.HERE / "results.json").read_text()))
    with args.output.open("x") as out:
        json.dump(report, out, allow_nan=False, separators=(",", ":"))
        out.write("\n")
    for label, survey in report["surveys"].items():
        print(
            label,
            {k: v["max_errors_by_time"] for k, v in survey["comparisons"].items()},
        )
