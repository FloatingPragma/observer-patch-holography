"""Read two independently metered scalar stresses on one evolving semiclassical metric."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ... import observer_field
from .. import balls
from .. import measure as scalar

HERE = Path(__file__).resolve().parent
SEED = 298902
MASSES_SQUARED = np.array([observer_field.MASS_SQUARED, 2.1])
WAVE_NUMBERS = np.array([1.0, 2.0])
REFERENCES = {
    "ground": (None, None),
    "b_y": (None, 1),
    "a_x": (0, None),
    "a_x_b_y": (0, 1),
    "a_x_b_z": (0, 2),
}
ADDITIONS = (("ground", "b_y"), ("a_x", "a_x_b_y"), ("a_x", "a_x_b_z"))
MECHANISMS = {"common": 1.0, "selective": 0.5}


def prepare(axes):
    number = np.zeros((2, 3), dtype=int)
    for species, axis in enumerate(axes):
        if axis is not None:
            number[species, axis] = 1
    return number


def geometry_source(
    number,
    b_weight,
    *,
    masses=MASSES_SQUARED,
    waves=WAVE_NUMBERS,
    reference_frequency=None,
    end=1.5,
):
    """Hamilton equations for log lengths, ADM momenta and scalar canonical maps."""
    kappa, vacuum = observer_field.COUPLING, observer_field.VACUUM
    multiplicities = number + 0.5
    weights = np.array([1.0, b_weight])
    natural = np.sqrt(masses + waves**2)
    if reference_frequency is None:
        omega0 = natural
        energy_per_pair = 2 * omega0
    else:
        omega0 = reference_frequency
        energy_per_pair = omega0 + natural**2 / omega0
    rho = energy_per_pair * multiplicities.sum(axis=1) / scalar.VOL0
    shear = np.array([0.040, -0.025, -0.015])
    expansion = math.sqrt((vacuum + kappa * weights @ rho + 0.5 * shear @ shear) / 3)
    rates0 = expansion + shear
    momenta = scalar.VOL0 * (rates0 - rates0.sum()) / kappa
    initial_maps = np.tile(np.eye(2), (2, 3, 1, 1))
    initial = np.r_[np.zeros(3), momenta, initial_maps.ravel(), np.zeros(3)]

    def derivative(_, state):
        lengths = np.exp(state[:3])
        volume = lengths.prod()
        momentum = state[3:6]
        maps = state[6:-3].reshape(2, 3, 2, 2)
        gradients = waves[:, None] ** 2 / lengths[None, :] ** 2
        omega2 = masses[:, None] + gradients
        rates = kappa * (momentum - 0.5 * momentum.sum()) / (scalar.VOL0 * volume)
        widths = maps[:, :, :, 0] ** 2 / omega0[:, None, None]
        widths += maps[:, :, :, 1] ** 2 * omega0[:, None, None]
        qq, pp = np.moveaxis(multiplicities[:, :, None] * widths, -1, 0)
        # Sum -d<H_species>/d(log length) for sine and cosine partners together.
        force = (pp / volume - volume * omega2 * qq).sum(axis=1)[:, None]
        force = force + 2 * volume * gradients * qq
        pdot = 0.5 * momentum @ rates - vacuum * scalar.VOL0 * volume / kappa
        pdot = pdot + weights @ force
        map_rate = np.empty_like(maps)
        map_rate[:, :, 0] = maps[:, :, 1] / volume
        map_rate[:, :, 1] = -volume * omega2[:, :, None] * maps[:, :, 0]
        return np.r_[rates, pdot, map_rate.ravel(), 1 / lengths]

    source = solve_ivp(
        derivative,
        (0, end),
        initial,
        method="DOP853",
        dense_output=True,
        atol=2e-13,
        rtol=2e-13,
    )
    if not source.success:
        raise RuntimeError(source.message)
    return source.sol


def acquire(source, number, rng):
    omega0 = np.sqrt(MASSES_SQUARED + WAVE_NUMBERS**2)
    samples = []

    def radar_geometry(time):
        state = source(time)
        return np.r_[np.exp(state[:3]), np.zeros(5), state[-3:]]

    for time in scalar.TIMES:
        maps = source(time)[6:-3].reshape(2, 3, 2, 2)
        np.testing.assert_allclose(np.linalg.det(maps), 1, rtol=0, atol=3e-11)
        fields = []
        for species, matrices in enumerate(maps):
            widths = np.sqrt(scalar.field_widths(matrices, omega0[species]))
            readings = []
            for axis, (position, momentum) in enumerate(widths):
                mix = scalar.ODD if number[species, axis] else scalar.GROUND
                readings.append(
                    [
                        dict(
                            q=scalar.meter(position, mix, rng),
                            p=scalar.meter(momentum, mix, rng),
                        )
                        for _ in range(2)
                    ]
                )
            fields.append(readings)
        radar = [
            observer_field.records(radar_geometry, time, h) for h in scalar.SPACINGS
        ]
        samples.append(dict(field=fields, radar=radar))
    return samples


def calibrate(primary):
    rows = []
    for S in primary["surveys"]:
        ground, excited = [
            scalar.read(S["records"][arm][0], primary["recovered_mass2"])
            for arm in ("ground", "x")
        ]
        coupling, vacuum = scalar.inferred(ground, excited)
        rows.append(dict(coupling=float(coupling), vacuum=float(vacuum)))
    return rows


def read(tape, level, masses):
    radar = tape["radar"][level]
    _, lengths, *_ = observer_field.local_state(radar)
    temporal, spatial, *_ = observer_field.reconstruct(radar)
    matter = []
    for meter, mass2, k in zip(tape["field"], masses, WAVE_NUMBERS, strict=True):
        qq, pp = [
            np.array(
                [sum(scalar.second_moment(m[key]) for m in axis) for axis in meter]
            )
            for key in ("q", "p")
        ]
        kinetic = pp.sum() / (scalar.VOL0 * lengths.prod() ** 2)
        potential = mass2 * qq.sum() / scalar.VOL0
        directional = k**2 * qq / (scalar.VOL0 * lengths**2)
        rho = 0.5 * (kinetic + potential + directional.sum())
        pressures = 0.5 * (kinetic - potential - directional.sum()) + directional
        matter.append(np.r_[rho, pressures])
    return np.r_[temporal, spatial], np.array(matter)


def infer_fraction(delta_g, stress_delta, coupling):
    delta_a, delta_b = stress_delta
    return float((delta_g / coupling - delta_a) / delta_b)


def analyze(rows, apertures, calibrations, masses):
    levels = []
    for i, cal in enumerate(calibrations):
        kappa, vacuum = cal["coupling"], cal["vacuum"]
        measurements = {
            name: [read(t, i, masses) for t in tapes] for name, tapes in rows.items()
        }
        residuals = []
        for measurements_by_time in measurements.values():
            for geometry, stress in measurements_by_time:
                forecast = kappa * stress.sum(axis=0) + vacuum * np.array(
                    [1, -1, -1, -1]
                )
                residuals.append(abs(geometry - forecast))
        fractions = {}
        for base, added in ADDITIONS:
            name = base + "_to_" + added
            fractions[name] = [
                infer_fraction(g1[0] - g0[0], s1[:, 0] - s0[:, 0], kappa)
                for (g0, s0), (g1, s1) in zip(
                    measurements[base], measurements[added], strict=True
                )
            ]
        levels.append(
            dict(
                spacing=scalar.SPACINGS[i],
                max_diagonal_errors=np.max(residuals, axis=0).tolist(),
                inferred_b_coupling_ratios=fractions,
            )
        )
    # Both optical and area comparisons use the same measured local rest frame.
    added = measurements["a_x_b_y"][1]
    base = measurements["a_x"][1]
    delta = added[1][:, 0] - base[1][:, 0]
    area = []
    for sphere, ground in zip(apertures["a_x_b_y"], apertures["a_x"], strict=True):
        curvature = sphere["estimated_Guu"] - ground["estimated_Guu"]
        area.append(
            dict(
                radius=sphere["reference"],
                delta_curvature=float(curvature),
                inferred_b_coupling_ratio=infer_fraction(
                    curvature, delta, calibrations[-1]["coupling"]
                ),
            )
        )
    return dict(radar=levels, equal_volume=area)


def run():
    calibration_path = scalar.HERE / "results.json"
    primary = json.loads(calibration_path.read_text())
    cal = calibrate(primary)
    omega_b = math.sqrt(MASSES_SQUARED[1] + WAVE_NUMBERS[1] ** 2)
    bright = int(
        np.random.default_rng([SEED, 0]).binomial(
            scalar.SHOTS, (1 + math.cos(omega_b / 2)) / 2
        )
    )
    b_mass = (2 * math.acos(2 * bright / scalar.SHOTS - 1)) ** 2 - WAVE_NUMBERS[1] ** 2
    masses = [primary["recovered_mass2"], float(b_mass)]
    arms = {}
    for j, (mechanism, b_weight) in enumerate(MECHANISMS.items()):
        records = {}
        apertures = {}
        for i, (name, axes) in enumerate(REFERENCES.items()):
            number = prepare(axes)
            source = geometry_source(number, b_weight)
            records[name] = acquire(
                source, number, np.random.default_rng([SEED, 1, j, i])
            )
            if name in ("a_x", "a_x_b_y"):
                apertures[name] = balls.survey_source(source)
        arms[mechanism] = dict(
            records=records,
            apertures=apertures,
            result=analyze(records, apertures, cal, masses),
        )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        calibration_receipt_sha256=hashlib.sha256(
            calibration_path.read_bytes()
        ).hexdigest(),
        shots=scalar.SHOTS,
        seed=SEED,
        b_spectrometer_bright=bright,
        measured_masses_squared=masses,
        calibrations=cal,
        arms=arms,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"), allow_nan=False)
        stream.write("\n")
    for arm, row in result["arms"].items():
        print(arm, json.dumps(row["result"]))
