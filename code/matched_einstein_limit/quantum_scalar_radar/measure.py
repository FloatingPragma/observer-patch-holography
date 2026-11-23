"""Read quantum stress and radar curvature on a controlled semiclassical box."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erf

from .. import observer_field as classical


HERE = Path(__file__).resolve().parent
VOL0 = (2 * math.pi) ** 3
SPACINGS = (0.08, 0.04, 0.02)
TIMES = (0.35, 0.60, 0.85)
SHOTS = 2**27
SEED = 762011
BINS = np.arange(769) / 64
GROUND, ODD, EVEN = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.5, 0.0, 0.5)
SOURCES = dict(
    ground=(GROUND,) * 3,
    x=(ODD, GROUND, GROUND),
    y=(GROUND, ODD, GROUND),
    z=(GROUND, GROUND, ODD),
    even_x=(EVEN, GROUND, GROUND),
)


def field_widths(mapping, omega0):
    return mapping[:, :, 0] ** 2 / omega0 + mapping[:, :, 1] ** 2 * omega0


def propagate(preparation, decoupled=False):
    coupling, vacuum, mass2 = (
        classical.COUPLING,
        classical.VACUUM,
        classical.MASS_SQUARED,
    )
    occupation_factor = np.asarray(preparation) @ np.array([0.5, 1.5, 2.5])
    omega0 = math.sqrt(1 + mass2)
    # One oscillator of each spatial phase per axis. Its summed initial energy
    # is homogeneous in space even for inputs with non-Gaussian higher moments.
    rho0 = 2 * omega0 * occupation_factor.sum() / VOL0
    shear = np.array([0.040, -0.025, -0.015])
    hubble = math.sqrt(
        (vacuum + (0 if decoupled else coupling * rho0) + 0.5 * float(shear @ shear))
        / 3
    )
    rates0 = hubble + shear
    gravit_momentum = VOL0 / coupling * (rates0 - rates0.sum())
    initial = np.r_[
        np.zeros(3), gravit_momentum, np.tile(np.eye(2).ravel(), 3), np.zeros(3)
    ]

    def derivative(_, state):
        lengths = np.exp(state[:3])
        volume = lengths.prod()
        momentum = state[3:6]
        symplectic = state[6:18].reshape(3, 2, 2)
        omega2 = mass2 + 1 / lengths**2
        rates = coupling / (VOL0 * volume) * (momentum - 0.5 * momentum.sum())
        gravitational_kinetic = 0.5 * float(momentum @ rates)
        squares = occupation_factor[:, None] * field_widths(symplectic, omega0)
        q2, p2 = squares.T
        # -d<H_field>/d(log a_i), holding the canonical field operators fixed.
        force = (
            np.full(3, np.sum(p2 / volume - volume * omega2 * q2))
            + 2 * volume * q2 / lengths**2
        )
        pdot = (
            gravitational_kinetic
            - vacuum * VOL0 * volume / coupling
            + (0 * force if decoupled else force)
        )
        udot = np.empty_like(symplectic)
        udot[:, 0, :] = symplectic[:, 1, :] / volume
        udot[:, 1, :] = -volume * omega2[:, None] * symplectic[:, 0, :]
        return np.r_[rates, pdot, udot.ravel(), 1 / lengths]

    solution = solve_ivp(
        derivative,
        (0, 1.5),
        initial,
        dense_output=True,
        method="DOP853",
        atol=2e-13,
        rtol=2e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.sol


def folded_cdf(scale, probabilities):
    """Absolute-value quadrature distribution for the three declared Fock levels."""
    x = BINS / scale
    vacuum = erf(x)
    term = 2 * np.exp(-(x**2)) / math.sqrt(math.pi)
    cdf = (
        probabilities[0] * vacuum
        + probabilities[1] * (vacuum - term * x)
        + probabilities[2] * (vacuum - term * (x**3 + 0.5 * x))
    )
    return np.maximum.accumulate(np.clip(cdf, 0, 1))


def meter(scale, probabilities, rng):
    cdf = folded_cdf(scale, probabilities)
    probs = np.r_[np.diff(cdf), 1 - cdf[-1]]
    bins = rng.multinomial(SHOTS, probs / probs.sum())
    if bins[-1] or probs[-1] > 1e-12:
        raise ValueError(
            "Instrument range failed to contain the quadrature distribution"
        )
    return bins[: np.flatnonzero(bins)[-1] + 1].tolist()


def acquire(solution, preparation, spacing, rng):
    omega0 = math.sqrt(1 + classical.MASS_SQUARED)

    def radar_source(time):
        s = solution(time)
        # Adapter for the comoving three-arm null-reflection instrument.
        return np.r_[np.exp(s[:3]), np.zeros(3), 0.0, 0.0, s[-3:]]

    tapes = []
    for time in TIMES:
        state = solution(time)
        mapping = state[6:18].reshape(3, 2, 2)
        np.testing.assert_allclose(np.linalg.det(mapping), 1, rtol=0, atol=2e-11)
        scales = np.sqrt(field_widths(mapping, omega0))
        records = []
        for mix, (q, p) in zip(preparation, scales, strict=True):
            # Read q, p and number parity on separate copies of each partner mode.
            records.append(
                [
                    dict(
                        q=meter(q, mix, rng),
                        p=meter(p, mix, rng),
                        odd=int(rng.binomial(SHOTS, mix[1])),
                    )
                    for _ in range(2)
                ]
            )
        tapes.append(
            dict(radar=classical.records(radar_source, time, spacing), modes=records)
        )
    return tapes


def second_moment(counts):
    if sum(counts) != SHOTS:
        raise ValueError("Incomplete quadrature record")
    centers = (BINS[: len(counts)] + BINS[1 : len(counts) + 1]) / 2
    return float(np.asarray(counts) @ centers**2 / SHOTS)


def read(tape, mass2):
    _, lengths, rates, acc, _, _, _ = classical.local_state(tape["radar"])
    volume = lengths.prod()
    q2 = np.array([sum(second_moment(x["q"]) for x in axis) for axis in tape["modes"]])
    p2 = np.array([sum(second_moment(x["p"]) for x in axis) for axis in tape["modes"]])
    kinetic = p2.sum() / (VOL0 * volume**2)
    potential = mass2 * q2.sum() / VOL0
    gradient = q2 / (VOL0 * lengths**2)
    energy = 0.5 * (kinetic + potential + gradient.sum())
    pressure = 0.5 * (kinetic - potential - gradient.sum()) + gradient
    temporal = sum(rates[i] * rates[j] for i, j in ((0, 1), (1, 2), (2, 0)))
    spatial = np.array(
        [-acc[i] - acc[j] - rates[i] * rates[j] for i, j in ((1, 2), (2, 0), (0, 1))]
    )
    return dict(energy=energy, pressure=pressure, temporal=temporal, spatial=spatial)


def inferred(calibration, excited):
    coupling = (excited["temporal"] - calibration["temporal"]) / (
        excited["energy"] - calibration["energy"]
    )
    return coupling, calibration["temporal"] - coupling * calibration["energy"]


def compare(rows, coupling, vacuum):
    errors, anisotropy = [], []
    for obs in rows:
        pressure = obs["pressure"]
        spatial = obs["spatial"]
        errors.extend(abs(spatial - (coupling * pressure - vacuum)).tolist())
        errors.append(abs(obs["temporal"] - vacuum - coupling * obs["energy"]))
        high, low = np.argmax(pressure), np.argmin(pressure)
        if high != low:
            predicted = coupling * (pressure[high] - pressure[low])
            anisotropy.append(abs((spatial[high] - spatial[low]) / predicted - 1))
    return float(max(errors)), float(max(anisotropy, default=0.0))


def run():
    mass_rng = np.random.default_rng([SEED, 0])
    expected_bright = (1 + math.cos(math.sqrt(1 + classical.MASS_SQUARED))) / 2
    bright = int(mass_rng.binomial(SHOTS, expected_bright))
    mass2 = math.acos(2 * bright / SHOTS - 1) ** 2 - 1
    sources = {label: propagate(prep) for label, prep in SOURCES.items()}
    sources["decoupled_y"] = propagate(SOURCES["y"], decoupled=True)
    results = []
    for level, spacing in enumerate(SPACINGS):
        tapes = {
            name: acquire(
                source,
                SOURCES["y" if name == "decoupled_y" else name],
                spacing,
                np.random.default_rng([SEED, 1, level, i]),
            )
            for i, (name, source) in enumerate(sources.items())
        }
        observations = {
            name: [read(t, mass2) for t in arm] for name, arm in tapes.items()
        }
        coupling, vacuum = inferred(observations["ground"][0], observations["x"][0])
        held = (
            observations["y"]
            + observations["z"]
            + observations["even_x"]
            + observations["x"][1:]
        )
        error, angular = compare(held, coupling, vacuum)
        no_backreaction = compare(observations["decoupled_y"], coupling, vacuum)[0]
        pooled = max(
            float(
                np.max(
                    abs(obs["spatial"] - (coupling * obs["pressure"].mean() - vacuum))
                )
            )
            for obs in held
        )
        # All spectrum-sensitive comparisons use member-labelled matter tapes.
        parity_matched_radar = all(
            a["radar"] == b["radar"]
            for a, b in zip(tapes["x"], tapes["even_x"], strict=True)
        )
        null_residual = [
            float(
                np.max(
                    abs(
                        (o["temporal"] + o["spatial"])
                        - coupling * (o["energy"] + o["pressure"])
                    )
                )
            )
            for o in held
        ]
        results.append(
            dict(
                spacing=spacing,
                calibration=dict(coupling=coupling, vacuum=vacuum),
                max_held_error=error,
                max_held_axial_null_error=max(null_residual),
                max_relative_pressure_anisotropy_error=angular,
                decoupled_error=no_backreaction,
                isotropic_stress_error=pooled,
                identical_radar_opposite_parity=parity_matched_radar,
                records=tapes,
            )
        )
    return dict(
        source="semiclassical ADM Bianchi I box with six quantized scalar modes",
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        seed=SEED,
        shots_per_meter=SHOTS,
        maximum_quadrature=12,
        bin_width=1 / 64,
        isolated_frequency_bright=bright,
        recovered_mass2=mass2,
        surveys=results,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
    for row in result["surveys"]:
        print({k: v for k, v in row.items() if k != "records"})
