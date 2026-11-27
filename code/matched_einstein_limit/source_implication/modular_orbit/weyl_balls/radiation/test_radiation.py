import json
from math import pi
from pathlib import Path

import numpy as np

from .measure import audit, clock_record, conformal_ball, geodesic_map, proper_ball, volume
from .core_modes import optical_clearance


def test_flat_exponential_map_and_independent_volume_stops():
    mu = np.linspace(-0.95, 0.95, 35)
    for speed in [0, 4 / 3, np.sqrt(8)]:
        w, jac = geodesic_map(0.09, mu, 0, speed)
        assert np.allclose(w, 0.09)
        assert np.allclose(jac, 0.09**2)
        radius, area = proper_ball(0.09, 0, speed)
        assert np.isclose(area, 4 * pi * radius**2)
        assert np.isclose(volume(radius, 0, speed), 4 * pi * 0.09**3 / 3)
        assert np.allclose(conformal_ball(0.09, 0, speed), [radius, area])


def test_comoving_geodesic_rest_rim_and_clock_normalization():
    h, r = 0.42, 0.24
    w, jac = geodesic_map(r, np.array([0, 0.5, 0.9]), h, 0)
    assert np.all(w == w[0])
    assert np.allclose(jac, (w * np.cos(h * w))**2)
    assert np.allclose(w / 2 + np.sin(2 * h * w) / (4 * h), r)
    assert np.isclose(clock_record(h)["inferred_expansion_rate"], h)


def test_one_eigenword_sources_all_frames_and_a_changed_aperture_is_visible():
    data = audit()
    assert data["flipped_band_pairs"] == 56 and data["radiation_calibration"]["flipped_band_pairs"] == 8
    assert np.isclose(data["ordinary_coupling"], 6 / data["kinetic_stiffness"])
    assert np.isclose(3 * data["held_clock"]["inferred_expansion_rate"]**2,
                      data["ordinary_coupling"] * data["comoving_energy_density"])
    errors, discrepancies = [], []
    for radius in [0.12, 0.06, 0.03, 0.015, 0.0075]:
        rows = [row for row in data["proper_balls"] if row["radius"] == radius]
        assert len(rows) == 10
        assert all(row["proper_area"] < 4 * pi * radius**2 for row in rows)
        assert all(row["volume_at_fixed_rod_radius"] < 4 * pi * radius**3 / 3 for row in rows)
        errors.append(max(abs(row["area_coefficient"] - row["stress_forecast"]) for row in rows))
        discrepancies.append(max(abs(row["diamond_minus_true_area_over_R4"]) for row in rows))
    assert all(later < 0.26 * earlier for earlier, later in zip(errors, errors[1:]))
    assert all(later < 0.26 * earlier for earlier, later in zip(discrepancies, discrepancies[1:]))
    assert 0 < errors[-1] < 0.003
    first, second = data["finite_rest_entropy"]
    assert first["same_sites"] and not second["same_sites"]
    assert 0 < second["registered_core_entropy_change"] < 0.15
    assert second["all_reconstructed_sites_entropy_change"] > 2
    controls = data["fixed_occupation_controls"]["area_coefficients"]
    assert abs(controls["flat_clock_and_rods"]) < 1e-8
    assert np.isclose(2 * controls["doubled_stiffness"], controls["declared"], rtol=0.002)


def test_boosted_optical_core_clearance_and_archived_finite_modes():
    for speed in (0., 0.5, np.sqrt(8)):
        assert np.isclose(optical_clearance(0.1, 0., speed), 0.1)
    data = json.loads(Path(__file__).with_name("core_modes.json").read_text())
    assert len(data["frames"]) == 10 and data["mode_count"] == 8
    for row in data["frames"]:
        assert data["packet_support_radius"] < row["sampled_optical_interior_clearance"]
        assert 0 < row["registered_entropy_change"] < row["weighted_stress_forecast"]
        assert abs(row["vacuum_transport_change_at_finite_cutoff"]) < 0.001
    assert data["frames"][1]["registered_entropy_change"] > data["frames"][0]["registered_entropy_change"]
