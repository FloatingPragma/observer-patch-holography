"""Independent source-constraint and held-observation checks."""

import json
import math

import numpy as np

from . import forecast, measure
from .. import observer_field as classical


def test_source_satisfies_constraint_and_retains_distinct_quantum_preparations():
    sources = {name: measure.propagate(prep) for name, prep in measure.SOURCES.items()}
    # Equal total number and rotated wave vectors fix identical initial geometry.
    for name in ("y", "z", "even_x"):
        np.testing.assert_allclose(
            sources[name](0)[:6], sources["x"](0)[:6], atol=1e-14
        )
    for name, source in sources.items():
        number = np.asarray(measure.SOURCES[name]) @ [0.5, 1.5, 2.5]
        frequency = math.sqrt(1 + classical.MASS_SQUARED)
        for time in np.linspace(0, 1.4, 15):
            state = source(time)
            arms, momenta = np.exp(state[:3]), state[3:6]
            volume = arms.prod()
            integral = measure.VOL0 * volume
            maps = state[6:18].reshape(3, 2, 2)
            np.testing.assert_allclose(np.linalg.det(maps), 1, atol=2e-11, rtol=0)
            qq, pp = (number[:, None] * measure.field_widths(maps, frequency)).T
            field = sum(
                pp / volume + volume * (arms ** (-2) + classical.MASS_SQUARED) * qq
            )
            kinetic = (
                classical.COUPLING
                * (momenta @ momenta - 0.5 * momenta.sum() ** 2)
                / (2 * integral)
            )
            vacuum = classical.VACUUM * integral / classical.COUPLING
            assert abs(kinetic + vacuum + field) < 1e-10


def test_observer_records_and_controls_reproduce():
    receipt = json.loads((measure.HERE / "results.json").read_text())
    assert measure.run() == receipt
    resolutions = [s["max_held_error"] for s in receipt["surveys"]]
    assert resolutions[0] > resolutions[1] > resolutions[2]
    fine = receipt["surveys"][-1]
    assert fine["max_held_error"] < 0.001
    assert fine["max_relative_pressure_anisotropy_error"] < 0.05
    assert fine["decoupled_error"] > 0.005
    assert fine["isotropic_stress_error"] > 0.005
    assert fine["calibration"]["coupling"] > 0
    for odd, even in zip(fine["records"]["x"], fine["records"]["even_x"], strict=True):
        assert odd["radar"] == even["radar"]
        assert all(row["odd"] == measure.SHOTS for row in odd["modes"][0])
        assert all(row["odd"] == 0 for row in even["modes"][0])
    saved = json.loads((measure.HERE / "forecast.json").read_text())
    assert forecast.run(receipt) == saved
    for survey in saved["surveys"].values():
        fitted = max(survey["comparisons"]["measured"]["max_errors_by_time"])
        permuted = max(survey["comparisons"]["permuted"]["max_errors_by_time"])
        assert fitted < 0.0001 and permuted > 5 * fitted


def test_volume_and_fixed_volume_area_readers():
    from . import balls, volume

    primary = json.loads((measure.HERE / "results.json").read_text())
    packet = json.loads((measure.HERE / "volume.json").read_text())
    assert volume.run(primary) == packet
    errors = [cohort["max_held_error"] for cohort in packet["refinements"]]
    assert errors[0] > errors[1] > errors[2] and errors[2] < 0.0001
    for cohort in packet["refinements"]:
        assert cohort["flat"]["volume_rate"] < 0
        assert abs(cohort["flat"]["ricci"]) < 1e-7
        assert abs(cohort["flat"]["raw"]) > 14

    packet = json.loads((measure.HERE / "balls.json").read_text())
    # Recompute the reported inference directly from independent interior and rim records.
    _, _, _, inside_weights = balls.grid(4)
    _, _, _, rim_weights = balls.grid(5)
    _, radial_weights = np.polynomial.legendre.leggauss(6)
    for label, cohort in packet["records"].items():
        residuals = []
        for sphere in cohort["spheres"]:
            ell = sphere["reference"]
            normalization = sphere["matched_radius"] * balls.TICK / 2
            observed = (
                normalization
                * radial_weights
                @ (np.array(sphere["interior_ticks"]) @ inside_weights)
            )
            assert abs(observed - 4 * math.pi * ell**3 / 3) < 1e-13
            area = np.array(sphere["boundary_ticks"]) @ rim_weights * balls.TICK
            guu = 15 * (1 - area / (4 * math.pi * ell**2)) / ell**2
            assert abs(guu - sphere["estimated_Guu"]) < 1e-10
            residuals.append(abs(guu - cohort["stress_prediction"]))
        assert residuals[0] > residuals[1] > residuals[2] and residuals[2] < 0.001
    reproduced = balls.acquire(measure.SOURCES["y"])
    np.testing.assert_allclose(
        [s["estimated_Guu"] for s in reproduced],
        [s["estimated_Guu"] for s in packet["records"]["y"]["spheres"]],
        atol=1e-8,
        rtol=0,
    )
    assert all(abs(row["Guu"]) < 1e-8 for row in packet["flat_control"])


def test_reserved_matter_meter_cross_check():
    from . import partial_meter

    primary, geometry, saved = [
        json.loads((measure.HERE / name).read_text())
        for name in ("results.json", "balls.json", "partial_meter.json")
    ]
    assert partial_meter.run(primary, geometry) == saved
    for row in saved.values():
        assert abs(row["comparisons"][-1]["discrepancy"]) < 0.00005
    assert saved["z"]["reserved_change"] > 0 > saved["y"]["reserved_change"]
