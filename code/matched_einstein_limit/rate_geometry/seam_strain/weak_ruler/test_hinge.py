"""Euclidean controls, exact physical inequalities and finite positive-source transfer."""

from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
import pytest
from scipy.spatial.distance import cdist

from . import hinge


RECEIPT = hinge.HERE / "hinge.json"
SHA = "3f788c73b5b0eef3069ef5175e6c26e281ab79900bcab15021428368d2094dec"


@pytest.fixture(scope="module")
def saved():
    raw = RECEIPT.read_bytes()
    assert sha256(raw).hexdigest() == SHA
    return json.loads(raw)


def test_flat_subdivision_regular_tetrahedra_and_schlaefli_identity():
    corners = np.array([[0, 0, 0], [2.0, 0, 0], [0.3, 2.4, 0], [0.1, 0.6, 1.8]])
    flat = np.vstack((corners, np.array([0.1, 0.2, 0.3, 0.4]) @ corners))
    x = cdist(flat, flat)
    v, deficits, _ = hinge.bend(x)
    np.testing.assert_allclose(deficits, 0, atol=8e-15, rtol=0)
    assert abs(v[0]) < 1e-13 and abs(v[2] - 4 * np.pi) < 8e-15
    assert abs(v[1] - 0.8 * 1.8) < 2e-15
    regular = np.ones((5, 5)) - np.eye(5)
    vals, deficits, boundary = hinge.bend(regular)
    theta = np.arccos(1 / 3)
    np.testing.assert_allclose(deficits, 2 * np.pi - 3 * theta, atol=8e-15, rtol=0)
    np.testing.assert_allclose(boundary, np.pi - 2 * theta, atol=8e-15, rtol=0)
    assert abs(vals[0] - 8 * (2 * np.pi - 3 * theta)) < 1e-13
    assert abs(vals[1] - np.sqrt(2) / 3) < 1e-15
    assert abs(vals[2] + sum(deficits) - 4 * np.pi) < 8e-15
    assert abs(vals[3] - 12 * (np.pi - 2 * theta)) < 1e-13
    var = np.triu(np.arange(25).reshape(5, 5) / 17, 1)
    var = var + var.T
    dh = hinge.derive(regular, var)
    direct = 2 * (
        var[4, :4] @ deficits
        + np.array([var[e] for e in hinge.EDGES if 4 not in e]) @ boundary
    )
    assert abs(dh[0] + dh[3] - direct) < 3e-14
    scaled = hinge.derive(regular, regular)
    np.testing.assert_allclose(
        scaled, vals * np.array([1, 3, 0, 1, 2]), atol=3e-14, rtol=0
    )


def test_certified_retained_curvature_transfers_and_positive_preparations(saved):
    rebuilt = hinge.run()
    cert = rebuilt["finite_appointment_certificate"]
    old = saved["finite_appointment_certificate"]
    assert (
        cert["pooled_range_effect_operator_radius"]
        == old["pooled_range_effect_operator_radius"]
    )
    p = rebuilt["panels"][0]
    # Independent recurrence detector vs native exact conditional flights.
    assert (
        cert["reference_bulk_scalar"][0]
        < p["reference"]["bulk_scalar"]
        < cert["reference_bulk_scalar"][1]
    )
    assert (
        cert["filled_volume"][0]
        < p["reference"]["filled_volume"]
        < cert["filled_volume"][1]
    )
    assert 28 < cert["primary_cone_mean_times_r_squared"][0] < 29
    for tag in ("ordinary", "decorrelated", "hidden", "right", "left"):
        lo, hi = cert["raw_volume_slopes"][tag]
        assert lo < p["source_slopes"][tag]["filled_volume"] < hi
    centers, den, _ = hinge.exact_ranges()
    operators = (
        np.array(centers[:, 0], float) + 1j * np.array(centers[:, 1], float)
    ) / float(den)
    ref, changes, _ = hinge.live_survey.retained(True)
    radar = json.loads((hinge.HERE / "flight.json").read_text())
    clock = hinge.flight.pooled(radar, 16)
    for name, state in {"reference": ref, **changes}.items():
        actual = hinge.live_survey.read(operators, state)
        expected = [clock[name][e] for e in hinge.EDGES]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=3e-12)
    right = cert["matched"]["right"]
    for tag in ("decorrelated", "hidden"):
        row = right["sources"][tag]
        lo, hi = row["residual"]
        assert lo < p["matched"]["right"]["held"][tag]["measured_minus_predicted"] < hi
        oldlo, oldhi = old["matched"]["right"]["sources"][tag]["residual"]
        assert max(lo, oldlo) < min(hi, oldhi)
    for exp in (saved, rebuilt):
        geometry = exp["finite_appointment_certificate"]["matched"]["right"]
        collar, hidden = (
            geometry["sources"][key] for key in ("decorrelated", "hidden")
        )
        assert collar["residual"][1] < -2.60e-5 and hidden["residual"][0] > 2.20e-4
        assert collar["total_ordinary_residual"][1] < -2.2e-6
        assert hidden["total_ordinary_residual"][0] > 6.7e-5
        assert geometry["determinant_bulk_total_modular"][1] < -2.4e-11
        bending = geometry["implicit_geometry_constraints"]
        assert bending["volume_and_two_curvatures_control_minor"][0] > 3.8e-16
        left = bending["witnesses"]["left"]
        assert left["decorrelation_coefficient"][1] < -0.031
        assert left["modular_derivative_at_fixed_geometry"][1] < -0.00035
        area = geometry["equal_volume_outer_area_bulk"]
        assert area["control_minor"][0] > 1.3e-16
        assert area["partial_repairing_per_left_radian"][1] < -0.026
        assert area["retained_entropy_per_left_radian"][1] < -0.00027
        assert area["twelve_port_entropy_per_left_radian"][1] < -0.0038
        assert area["regional_ports"] == [0, 1, 5, 7, 8, 9, 38, 39, 40, 42, 46, 47]
        for record, third, slope in (
            (
                "finite_equal_volume_outer_area_and_bulk",
                "outer_area",
                area["retained_entropy_per_left_radian"],
            ),
            (
                "finite_equal_volume_bulk_and_boundary",
                "bulk_plus_boundary",
                left["modular_derivative_at_fixed_geometry"],
            ),
        ):
            trials = exp[record]
            assert [trial["left_rotation"] for trial in trials] == [-1e-4, -1e-5, -1e-6]
            for trial in trials:
                assert 0 < trial["partial_repairing"] < 1
                assert all(0 < probability < 1 for probability in trial["rail_lottery"])
                assert abs(sum(trial["rail_lottery"]) - 1) < 2e-14
                assert trial["minimum_retained_eigenvalue"] > 0.0014
                assert set(trial["geometry_error"]) == {"volume", "bulk_scalar", third}
                assert max(map(abs, trial["geometry_error"].values())) < 2e-17
                assert trial["entropy_change"] > 0 and trial["modular_change"] > 0
            assert (
                slope[0]
                < trials[-1]["entropy_change"] / trials[-1]["left_rotation"]
                < slope[1]
            )
        nominal = exp["panels"][0]["matched"]["right"]["held"]["decorrelated"][
            "measured_minus_predicted"
        ]
        rows = exp["finite_positive_right_matches"]
        last = rows[-1]
        assert all(
            abs(a["volume_error"]) < 1e-17 for r in rows for a in r["arms"].values()
        )
        assert abs(last["held_residual"] - nominal) < 1e-6
        assert (
            abs(
                last["arms"]["ordinary"]["modular_secant"]
                - p["matched"]["right"]["sources"]["ordinary"]["modular"]
            )
            < 1e-6
        )


def test_continuous_meshes_and_the_all_source_held_cone(saved):
    f = hinge.uniform()
    assert f == saved["uniform_held_star"] and f["stations"] == [0, 39, 5, 7, 40]
    for mode, row in f["bounds"].items():
        assert Fraction(row["cone_mean_scalar_times_r_squared_lower"]) > 13
        assert Fraction(row["safe_cone_radius_squared_lower"]) > 1
        assert Fraction(row["dihedral_cosine_lower"]) > 0
    fine = saved["panels"][-1]
    coarse = saved["panels"][-2]
    assert fine["continuous"] and coarse["continuous"]
    assert max(fine["maximum_real_derivative_discrepancy"].values()) < 1e-11
    assert (
        abs(
            fine["cone_mean_scalar_times_r_squared"]
            - coarse["cone_mean_scalar_times_r_squared"]
        )
        < 0.001
    )
    assert fine["spherical_link_identity_error"] < 1e-13
    for arm in ("ordinary", "decorrelated", "hidden"):
        a = fine["cone_volume_matched"][arm]
        w = fine["reference"]["spherical_link"]
        assert (
            abs(a["dI_per_r"] + 4 * (2 * np.pi + w) / w * a["dA_per_r_squared"]) < 1e-15
        )
        if arm != "ordinary":
            x = fine["matched"]["right"]["held"][arm]["measured_minus_predicted"]
            y = coarse["matched"]["right"]["held"][arm]["measured_minus_predicted"]
            assert abs(x - y) < 1e-6 and abs(x) > 1e-5
