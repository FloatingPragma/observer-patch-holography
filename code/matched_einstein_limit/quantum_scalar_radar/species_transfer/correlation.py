"""Read gravity sourced by a cross-field interaction at fixed single-field marginal states."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import measure as transfer
from .. import balls, measure as scalar
from ... import observer_field

HERE = transfer.HERE
SEED = 298903
TIMES = (0.60, 0.85, 1.10)
DIAGONAL_MASS = observer_field.MASS_SQUARED
REFERENCE = math.sqrt(1 + DIAGONAL_MASS)
WAVES = np.ones(2)


def field_source(sign, interaction):
    number = transfer.prepare((0, None) if sign == 1 else (None, 0))
    masses = DIAGONAL_MASS + interaction * np.array([1, -1])
    params = dict(masses=masses, waves=WAVES, reference_frequency=np.full(2, REFERENCE))
    before = transfer.geometry_source(number, 1, end=-0.6, **params)
    after = transfer.geometry_source(number, 1, end=0.95, **params)

    def source(times):
        delta = np.asarray(times) - TIMES[0]
        if delta.ndim == 0:
            return (after if delta >= 0 else before)(float(delta))
        data = np.empty((len(after(0)), len(delta)))
        positive = delta >= 0
        if any(positive):
            data[:, positive] = after(delta[positive])
        if any(~positive):
            data[:, ~positive] = before(delta[~positive])
        return data

    return number, source


def acquire(source, number, rng):
    def radar_source(time):
        state = source(time)
        return np.r_[np.exp(state[:3]), np.zeros(5), state[-3:]]

    tapes = []
    for time in TIMES:
        maps = source(time)[6:-3].reshape(2, 3, 2, 2)
        variances = np.array(
            [scalar.field_widths(matrix, REFERENCE) for matrix in maps]
        )
        # A reduced physical quadrature reads the occupied component of its canonical
        # superposition. Tracing the orthogonal oscillator gives a vacuum/one mixture.
        average = variances.mean(axis=0)
        component = (variances * number[:, :, None]).sum(axis=0) / (2 * average)
        singles = []
        for _ in range(2):
            readings = []
            for scale, odd in zip(np.sqrt(average), component, strict=True):
                axis = []
                for _ in range(2):
                    axis.append(
                        {
                            key: scalar.meter(s, (1 - n, n, 0.0), rng)
                            for key, s, n in zip(("q", "p"), scale, odd, strict=True)
                        }
                    )
                readings.append(axis)
            singles.append(readings)
        rotations = [
            [
                [
                    scalar.meter(
                        math.sqrt(variances[i, axis, 0]),
                        scalar.ODD if number[i, axis] else scalar.GROUND,
                        rng,
                    )
                    for _ in range(2)
                ]
                for axis in range(3)
            ]
            for i in range(2)
        ]
        radar = [observer_field.records(radar_source, time, h) for h in scalar.SPACINGS]
        tapes.append(dict(field=singles, rotations=rotations, radar=radar))
    return tapes


def read(tape, level, mass, interaction):
    radar = tape["radar"][level]
    _, lengths, *_ = observer_field.local_state(radar)
    g00, spatial, *_ = observer_field.reconstruct(radar)
    moments = np.array(
        [
            [
                [sum(scalar.second_moment(r[m]) for r in axis) for m in ("q", "p")]
                for axis in species
            ]
            for species in tape["field"]
        ]
    )
    qq, pp = moments.sum(axis=0).T
    kinetic = pp.sum() / (scalar.VOL0 * lengths.prod() ** 2)
    potential = mass * qq.sum() / scalar.VOL0
    gradient = qq / (scalar.VOL0 * lengths**2)
    marginal = np.r_[
        (kinetic + potential + gradient.sum()) / 2,
        (kinetic - potential - gradient.sum()) / 2 + gradient,
    ]
    combination = [
        sum(scalar.second_moment(row) for axis in basis for row in axis)
        for basis in tape["rotations"]
    ]
    covariance = (combination[0] - combination[1]) / (2 * scalar.VOL0)
    cross = interaction * covariance * np.array([1, -1, -1, -1])
    return np.r_[g00, spatial], marginal, cross, moments, covariance


def flat_spectra(interaction, rng):
    counts, masses = [], []
    for sign in (1, -1):
        omega = math.sqrt(1 + DIAGONAL_MASS + sign * interaction)
        count = int(rng.binomial(scalar.SHOTS, (1 + math.cos(omega / 2)) / 2))
        counts.append(count)
        masses.append((2 * math.acos(2 * count / scalar.SHOTS - 1)) ** 2 - 1)
    return counts, (masses[0] - masses[1]) / 2


def analyze(records, area, cal, mass, observed_interaction):
    refinements = []
    for level, calibration in enumerate(cal):
        kappa = calibration["coupling"]
        levels = []
        for plus, minus in zip(records["plus"], records["minus"], strict=True):
            a = read(plus, level, mass, observed_interaction)
            b = read(minus, level, mass, observed_interaction)
            delta_g = a[0] - b[0]
            local = kappa * (a[1] - b[1])
            complete = local + kappa * (a[2] - b[2])
            levels.append(
                dict(
                    geometry=delta_g.tolist(),
                    marginal_forecast=local.tolist(),
                    complete_forecast=complete.tolist(),
                    max_complete_error=float(max(abs(complete - delta_g))),
                    max_marginal_moment_difference=float(np.max(abs(a[3] - b[3]))),
                    measured_delta_cross_field_covariance=float(a[4] - b[4]),
                    delta_axial_null=(delta_g[0] + delta_g[1:]).tolist(),
                )
            )
        refinements.append(dict(spacing=scalar.SPACINGS[level], times=levels))
    fixed_volume = [
        dict(radius=a["reference"], delta_Guu=a["estimated_Guu"] - b["estimated_Guu"])
        for a, b in zip(area["plus"], area["minus"], strict=True)
    ]
    return dict(radar=refinements, area=fixed_volume)


def run():
    calibration_path = scalar.HERE / "results.json"
    primary = json.loads(calibration_path.read_text())
    calibrations = transfer.calibrate(primary)
    mass = primary["recovered_mass2"]
    arms = {}
    for i, (name, interaction) in enumerate((("mixed", 0.3), ("separate", 0.0))):
        bright, calibrated = flat_spectra(
            interaction, np.random.default_rng([SEED, 0, i])
        )
        records, area = {}, {}
        for j, (label, sign) in enumerate((("plus", 1), ("minus", -1))):
            number, source = field_source(sign, interaction)
            records[label] = acquire(
                source, number, np.random.default_rng([SEED, 1, i, j])
            )
            area[label] = balls.survey_source(source)
        arms[name] = dict(
            spectrometer_bright=bright,
            calibrated_interaction=calibrated,
            records=records,
            apertures=area,
            result=analyze(records, area, calibrations, mass, calibrated),
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "CORRELATION_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        calibration_receipt_sha256=hashlib.sha256(
            calibration_path.read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots=scalar.SHOTS,
        diagonal_mass_from_spectrometer=mass,
        calibrations=calibrations,
        arms=arms,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet = run()
    with args.output.open("x") as out:
        json.dump(packet, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
    for name, row in packet["arms"].items():
        print(
            name, row["calibrated_interaction"], json.dumps(row["result"]["radar"][-1])
        )
        print(row["result"]["area"])
