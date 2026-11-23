"""Regressions for both velocity bounds, enclosure selection and the lensing witness."""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.outer_orbit_transfer import linear_profile, measure


def assert_numeric_receipt(actual, archived):
    """Keep source identities and integers exact across BLAS implementations."""
    if isinstance(archived, dict):
        assert actual.keys() == archived.keys()
        for key in archived:
            assert_numeric_receipt(actual[key], archived[key])
    elif isinstance(archived, list):
        assert len(actual) == len(archived)
        for result, expected in zip(actual, archived):
            assert_numeric_receipt(result, expected)
    elif isinstance(archived, float):
        assert math.isclose(actual, archived, rel_tol=3e-14, abs_tol=1e-25), (
            actual,
            archived,
        )
    else:
        assert type(actual) is type(archived) and actual == archived


def test_sample_is_photometric_and_checked_against_the_catalog():
    inputs = json.loads((measure.HERE / "inputs.json").read_text())["galaxies"]
    result = measure.run()
    for row in inputs:
        assert abs(row["loggstar_check_dex"]) < 0.005
        assert abs(row["loggobs_check_dex"]) < 0.005
        assert measure.enclosed_fraction(row) <= 1
    panels = result["panels"]
    assert len(panels["primary"]["names"]) == 9
    assert len(panels["secondary"]["names"]) == 14
    assert panels["secondary_deep"]["names"] == ["NGC2685", "NGC3838", "NGC6798"]
    assert_numeric_receipt(
        result, json.loads((measure.HERE / "results.json").read_text())
    )


def test_shape_witness_and_deletion_bounds():
    result = measure.run()["panels"]["primary"]
    shape = result["orbital_shape"]
    assert shape["min_speed_error_multiplier"] == pytest.approx(4.4434276444)
    assert shape["one_speed_error"]["worst_excess"] > 1.45
    knots = shape["witness"]["knots"].values()
    for i in knots:
        f = max(k["lower"] * min(1, i["x"] / k["x"]) for k in knots)
        assert i["lower"] - 1e-12 <= f <= i["upper"] + 1e-12
    for name, stack in result["cases"]["central"]["surveys"].items():
        full = stack["all_bins"]
        assert all(
            full["weighted_upper"] <= upper < full["observed"]
            for upper in full["leave_one_out"].values()
        )
        wide = result["cases"]["two_error_box"]["surveys"][name]
        assert wide["all_bins"]["weighted_upper"] > full["observed"]
        assert stack["exterior_scaled_bins"]["bins"] == 6
        assert wide["exterior_scaled_bins"]["bins"] == 5


def test_constructive_witness_with_independent_forward_velocities_and_shears():
    receipt = linear_profile.run()
    assert_numeric_receipt(
        receipt, json.loads((measure.HERE / "linear.json").read_text())
    )
    assert receipt["lens_only_forecast_at_central_distance_kms"][
        "NGC3626"
    ] == pytest.approx([235.79, 244.02], abs=0.01)
    assert (
        receipt["panels"]["secondary_deep"]["distance_and_speed"]["error_multiplier"]
        < 0.575
    )
    assert receipt["panels"]["primary"]["distance_and_speed"]["error_multiplier"] > 3
    rows = {
        r["name"]: r
        for r in json.loads((measure.HERE / "inputs.json").read_text())["galaxies"]
    }
    survey = json.loads((measure.HERE.parent / "finite_halo/results.json").read_text())[
        "surveys"
    ]
    inputs = json.loads(measure.lensing_profile.RECEIPT.read_text())["profiles"]
    for panel in receipt["panels"].values():
        for kind, witness in panel.items():
            scale, k = witness["a0_m_s2"], witness["error_multiplier"]
            for name, record in witness["witness"].items():
                row = rows[name]
                s = record["distance_mpc"] / row["distance_mpc"]
                stellar_sq = (
                    0.8 * (row["vdisk_unit_kms"] ** 2 + row["vbul_unit_kms"] ** 2) * s
                )
                gm = (
                    measure.GM_SUN
                    * (0.8 * row["stellar_unit_msun"] + row["helium_hi_msun"])
                    * s**2
                )
                speed = math.sqrt(stellar_sq + math.sqrt(gm * scale) / 1e6)
                assert speed == pytest.approx(record["predicted_speed_kms"])
                assert abs(speed - row["speed_kms"]) <= (
                    (k if kind == "distance_and_speed" else 0) * row["speed_error_kms"]
                    + 1e-9
                )
                assert (
                    abs(record["distance_mpc"] - row["distance_mpc"])
                    <= k * row["distance_error_mpc"] + 1e-10
                )
            for name, stack in survey.items():
                p = inputs[name]
                mask = np.array(p["lower_bin_edge_m_s2"]) >= 1e-13
                g = np.array(p["g_bar_m_s2"])[mask]
                mean = np.array(stack["positive_weights"]) @ (
                    (4 / math.pi * g + np.sqrt(g * scale))
                    / measure.lensing_profile.FOUR_G
                )
                assert (
                    abs(mean - stack["weighted_signal"])
                    <= 2 * stack["statistical_se"] + 1e-10
                )
