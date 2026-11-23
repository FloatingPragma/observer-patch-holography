"""Independent Hamiltonian, area, sampling and cross-calibration checks."""

import hashlib
import json
import math

import numpy as np

from . import measure as trial


def test_two_field_hamiltonian_constraint_and_one_field_reduction():
    # Use the energy of each transported number state computed in canonical variables.
    numbers = trial.prepare((0, 1))
    for weight in (0.0, 0.5, 1.0):
        solution = trial.geometry_source(numbers, weight)
        for time in np.linspace(0, 1.3, 11):
            state = solution(time)
            lengths, p = np.exp(state[:3]), state[3:6]
            volume = lengths.prod()
            physical_volume = trial.scalar.VOL0 * volume
            gravit_energy = (
                trial.observer_field.COUPLING
                / (2 * physical_volume)
                * (p @ p - p.sum() ** 2 / 2)
                + trial.observer_field.VACUUM
                * physical_volume
                / trial.observer_field.COUPLING
            )
            matter_energies = []
            for i, maps in enumerate(state[6:-3].reshape(2, 3, 2, 2)):
                initial_omega = math.sqrt(
                    trial.MASSES_SQUARED[i] + trial.WAVE_NUMBERS[i] ** 2
                )
                squares = (
                    trial.scalar.field_widths(maps, initial_omega)
                    * (numbers[i] + 0.5)[:, None]
                )
                omega2 = (
                    trial.MASSES_SQUARED[i] + trial.WAVE_NUMBERS[i] ** 2 / lengths**2
                )
                matter_energies.append(
                    np.sum(squares[:, 1] / volume + volume * omega2 * squares[:, 0])
                )
            assert (
                abs(gravit_energy + matter_energies[0] + weight * matter_energies[1])
                < 2e-10
            )
        if weight == 0:
            reference = trial.scalar.propagate(trial.scalar.SOURCES["x"])
            times = np.linspace(0, 1.3, 11)
            for current, expected in (
                (solution(times)[:18], reference(times)[:18]),
                (solution(times)[-3:], reference(times)[-3:]),
            ):
                np.testing.assert_allclose(current, expected, rtol=1e-11, atol=1e-11)


def test_separate_radar_and_rim_readers_on_the_held_species():
    packet = json.loads((trial.HERE / "results.json").read_text())
    calibration = trial.scalar.HERE / "results.json"
    assert (
        packet["calibration_receipt_sha256"]
        == hashlib.sha256(calibration.read_bytes()).hexdigest()
    )
    assert (
        packet["protocol_sha256"]
        == hashlib.sha256((trial.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert abs(packet["measured_masses_squared"][1] - trial.MASSES_SQUARED[1]) < 0.002
    _, radial_weights = np.polynomial.legendre.leggauss(6)
    interior = trial.balls.grid(4)[-1]
    rim = trial.balls.grid(5)[-1]
    for name, arm in packet["arms"].items():
        target = trial.MECHANISMS[name]
        for preparation in arm["records"].values():
            for appointment in preparation:
                assert all(
                    sum(pair[meter]) == packet["shots"]
                    for species in appointment["field"]
                    for axis in species
                    for pair in axis
                    for meter in ("q", "p")
                )
        measured = trial.analyze(
            arm["records"],
            arm["apertures"],
            packet["calibrations"],
            packet["measured_masses_squared"],
        )
        finest = measured["radar"][-1]
        ratios = np.array(list(finest["inferred_b_coupling_ratios"].values()))
        np.testing.assert_allclose(ratios, target, rtol=0, atol=0.02)
        for i in range(len(trial.scalar.SPACINGS)):
            np.testing.assert_allclose(
                measured["radar"][i]["max_diagonal_errors"],
                arm["result"]["radar"][i]["max_diagonal_errors"],
                rtol=1e-6,
                atol=3e-11,
            )
        if target == 1:
            assert max(finest["max_diagonal_errors"]) < 0.001
        else:
            assert min(abs(1 - ratios).ravel()) > 0.25
            assert max(finest["max_diagonal_errors"]) > 0.005
        area_fractions = np.array(
            [r["inferred_b_coupling_ratio"] for r in measured["equal_volume"]]
        )
        np.testing.assert_allclose(area_fractions[-1], target, rtol=0, atol=0.03)
        assert all(np.diff(abs(area_fractions - target)) < 0)
        for spheres in arm["apertures"].values():
            for sphere in spheres:
                radius = sphere["reference"]
                v = (
                    0.5
                    * sphere["matched_radius"]
                    * trial.balls.TICK
                    * radial_weights
                    @ (np.array(sphere["interior_ticks"]) @ interior)
                )
                assert abs(v - 4 * math.pi * radius**3 / 3) < 2e-13
                area = trial.balls.TICK * (np.array(sphere["boundary_ticks"]) @ rim)
                np.testing.assert_allclose(
                    15 * (1 - area / (4 * math.pi * radius**2)) / radius**2,
                    sphere["estimated_Guu"],
                    atol=3e-10,
                )
