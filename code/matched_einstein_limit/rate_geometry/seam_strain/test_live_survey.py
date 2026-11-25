import json

import numpy as np
import pytest
from scipy.spatial.distance import cdist

from . import live_bounds, live_matching, live_survey


def test_intrinsic_cone_geometry_and_smooth_flat_control():
    # A Euclidean tetrahedron built around its actual centroid is flat.
    points = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1], [0, 0, 0.0]])
    flat = cdist(points, points)
    np.testing.assert_allclose(
        live_survey.star(flat), [8 * np.sqrt(3), 8 / 3, 4 * np.pi], atol=1e-12
    )
    # A regular 4-simplex gives four regular 3-cells and a singular apex.
    d = np.ones((5, 5)) - np.eye(5)
    area, volume, angle = live_survey.star(d)
    np.testing.assert_allclose(
        [area, volume, angle],
        [np.sqrt(3), 4 / (6 * np.sqrt(2)), 4 * (3 * np.arccos(1 / 3) - np.pi)],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        live_survey.slopes(d, d), [2 * area, 3 * volume, 0], atol=1e-12
    )
    assert live_survey.admission(d)["safe_cone_radius"] == pytest.approx(np.sqrt(2 / 3))


def test_replay_all_prep_station_and_detector_panels():
    saved = json.loads((live_survey.HERE / "live_survey.json").read_text())
    fresh = live_survey.run()
    assert fresh["stations"] == saved["stations"] == [0, 1, 38, 39, 5, 7, 40, 42]
    assert len(fresh["panels"]) == len(saved["panels"]) == 18
    for archived, panel in zip(saved["panels"], fresh["panels"], strict=True):
        assert list(archived) == list(panel)
        assert (archived["read"], archived["coupled"], archived["deadline"]) == (
            panel["read"],
            panel["coupled"],
            panel["deadline"],
        )
        np.testing.assert_allclose(archived["lengths"], panel["lengths"], atol=1e-11)
        for key in archived["length_slopes"]:
            np.testing.assert_allclose(
                archived["length_slopes"][key], panel["length_slopes"][key], atol=2e-12
            )
        for old, star in zip(archived["centers"], panel["centers"], strict=True):
            assert star["solid_angle"] == pytest.approx(old["solid_angle"], abs=1e-10)
            assert 0 < star["solid_angle"] < 4 * np.pi
            for side in ("left", "right"):
                for arm in ("ordinary", "decorrelated", "hidden"):
                    row = star["matched_boundary"][side][arm]
                    assert row["area_slope"] == pytest.approx(
                        old["matched_boundary"][side][arm]["area_slope"],
                        rel=2e-6,
                        abs=1e-9,
                    )
        assert (
            panel["round_comparison"]["all_radii"]["maximum_distance_change_lower"]
            > 0.2
        )
    primary = next(
        p
        for p in fresh["panels"]
        if p["coupled"] and p["read"] == "pooled" and p["deadline"] == 12
    )
    saved_match = json.loads((live_survey.HERE / "live_matching.json").read_text())
    finite, fixed = live_matching.finite_pairs(primary)
    for row, old in zip(finite, saved_match["finite_pairs"], strict=True):
        assert row["area_slope"] == pytest.approx(
            old["area_slope"], rel=1e-6, abs=2e-11
        )
        assert abs(row["volume_change"]) < 1e-18
    for row, old in zip(fixed, saved_match["matched_area_volume"], strict=True):
        assert max(abs(row[k]) for k in ("area_change", "volume_change")) < 1e-18
        np.testing.assert_allclose(row["angles"], old["angles"], rtol=2e-6, atol=1e-10)


def test_enclosed_same_source_geometry_and_modular_changes():
    old = json.loads((live_survey.HERE / "live_bounds.json").read_text())
    new = live_bounds.certificate()
    for record in (new, old):
        assert record["all_radius_distance_error_lower"] > 0.302
        assert record["two_rotation_jacobian"][1] < -9.7e-13
        assert record["charge_determinant_at_fixed_geometry"][1] < -0.004
        assert (
            record["fixed_geometry"]["decorrelated"]["fixed_outer_area_and_volume"][
                "modular"
            ]["retained"][0]
            > 0.0245
        )
        assert (
            record["fixed_geometry"]["decorrelated"]["fixed_cone_angle"]["modular"][
                "retained"
            ][0]
            > 0.025
        )
        assert record["methods"]["right"]["retained"]["ordinary_eta"][0] > 69
        assert record["methods"]["right"]["retained"]["predicted_area"]["hidden"][0] > 0
        assert record["methods"]["right"]["signed_area"]["hidden"][1] < 0
        assert record["methods"]["left"]["retained"]["ordinary_eta"][1] < 0
    for key in (
        "outer_area",
        "interior_volume",
        "solid_angle",
        "two_rotation_jacobian",
    ):
        assert max(old[key][0], new[key][0]) <= min(old[key][1], new[key][1])
