"""Independent metric, shear, and public cross-data checks."""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.finite_halo import measure, transfer


def test_gravity_and_clock_instruments():
    halo = measure.Halo(0.002, 0.006, 30.0)
    for r in (0.2, 1.0, 8.0, 31.0, 50.0):
        step = r * 1e-5
        beta = (
            r
            * (halo.log_lapse_squared(r + step) - halo.log_lapse_squared(r - step))
            / (4 * step)
        )
        np.testing.assert_allclose(beta, halo.orbital_beta2(r), rtol=3e-9)
        if r < halo.edge:
            # Obtain stress by differentiating the metric, independently of the radial closure.
            def data(x):
                m = halo.mass(x)
                derivative = (
                    halo.log_lapse_squared(x + step) - halo.log_lapse_squared(x - step)
                ) / (4 * step)
                return derivative, -0.5 * math.log1p(-2 * m / x)

            phi, lam = data(r)
            pp = (data(r + step)[0] - data(r - step)[0]) / (2 * step)
            lp = (data(r + step)[1] - data(r - step)[1]) / (2 * step)
            radial = math.exp(-2 * lam) * (2 * phi / r + 1 / r**2) - 1 / r**2
            tangential = math.exp(-2 * lam) * (pp + phi**2 - phi * lp + (phi - lp) / r)
            assert abs(radial) * r * r < 2e-8
            assert abs(tangential - halo.slope / r**2 * beta) * r * r < 2e-7
    # Matching conditions on the timelike edge.
    r = halo.edge
    assert (
        abs(halo.log_lapse_squared(r - 1e-8) - halo.log_lapse_squared(r + 1e-8)) < 1e-10
    )
    assert abs(halo.orbital_beta2(r - 1e-8) - halo.orbital_beta2(r + 1e-8)) < 1e-10
    vacuum = measure.Halo(1e-7, 0.0, 1.0)
    assert math.isclose(vacuum.deflection(2.0), 4e-7 / 2, rel_tol=5e-6)
    halo = measure.Halo(0.0, 1e-4, 1e6)
    x = halo.orbital_beta2(1.0)
    limit = math.pi * (math.sqrt(1 + 2 * x) / (1 - x) - 1)
    assert abs(halo.deflection(1.0) - limit) < 3e-10


def test_truncated_shear_by_projecting_shells():
    # Set A=R=1. Cap area and exterior shell contributions are integrated separately.
    for turn_over_edge in (0.01, 0.49, 0.9, 0.999, 1.0, 1.1, 5.0):
        edge = 1 / turn_over_edge
        enclosed = 4 * math.pi * min(1.0, edge)
        column = 0.0
        if edge > 1:

            def shell_cap(r):
                c = math.sqrt(1 - 1 / r**2)
                return 4 * math.pi / (r * r * (1 + c))

            enclosed += quad(shell_cap, 1.0, edge, epsabs=1e-12, epsrel=2e-12)[0]
            column = quad(lambda z: 2 / (1 + z * z), 0, math.sqrt(edge**2 - 1))[0]
        ratio = (enclosed / math.pi - column) / math.pi
        np.testing.assert_allclose(
            ratio, measure.shear_factor(turn_over_edge), atol=3e-11
        )
        assert ratio <= measure.MAX_SHEAR_FACTOR + 1e-13
    # A positive thin shell outside the aperture lowers its differential surface density.
    c = np.linspace(0.001, 0.999, 81)
    np.testing.assert_allclose(
        (1 - c) - (1 - c * c) / (2 * c), -((1 - c) ** 2) / (2 * c), atol=1e-12
    )


def test_data_and_profile_envelopes():
    raw = json.loads(measure.lensing_profile.RECEIPT.read_text())
    saved = json.loads((measure.HERE / "results.json").read_text())
    for label, summary in saved["surveys"].items():
        checked = measure.shear_stack(
            raw["profiles"][label],
            raw["fitted_sparc_a0_m_s2"],
            raw["sparc_cluster_bootstrap_95pct_m_s2"][1],
        )
        assert checked == summary
        # Dense independent edge/radius enumeration verifies the finite-knot supremum.
        p = raw["profiles"][label]
        low = np.asarray(p["lower_bin_edge_m_s2"])
        select = low >= 1e-13
        high = 2 * np.asarray(p["g_bar_m_s2"])[select] - low[select]
        low = low[select]
        limits = checked["common_edge_arbitrary_bin_weights"]["nominal"]
        weights = np.array(checked["positive_weights"])
        for edge in np.geomspace(low.min() / 2, high.max() * 2, 120):
            b = np.geomspace(low, high, 49)
            projected = (
                4 / math.pi * b
                + np.sqrt(raw["fitted_sparc_a0_m_s2"] * b)
                * np.vectorize(measure.shear_factor)(np.sqrt(edge / b))
            ) / measure.lensing_profile.FOUR_G
            assert (
                np.max(projected, axis=0) @ weights <= limits["weighted_upper"] + 1e-12
            )
    observed = transfer.run()
    saved = json.loads((measure.HERE / "transfer.json").read_text())
    for label, case in observed["cases"].items():
        for survey_name, row in case["surveys"].items():
            prior = saved["cases"][label]["surveys"][survey_name]
            for key in ("weighted_upper", "gap_over_statistical_se"):
                assert math.isclose(row[key], prior[key], rel_tol=2e-12)
                assert math.isclose(
                    row["sublinear_shared_profile"][key],
                    prior["sublinear_shared_profile"][key],
                    rel_tol=2e-12,
                )
            assert row["largest_allowed_gbar_m_s2"] < case["min_gt_at_low_scale"]


def test_reverse_calibration_units_and_roundtrip():
    from cosmology.rar_deep_regime.finite_halo import calibration

    value = calibration.run()
    saved = json.loads((measure.HERE / "calibration.json").read_text())
    assert json.loads(json.dumps(value)) == saved
    rings = json.loads(
        (measure.HERE.parent / "udg_transfer" / "rings.json").read_text()
    )
    outer = rings["panels"]["figure8"]["rings"][-1]
    radius = outer["radius_kpc"] * transfer.KPC
    gbar = (1000 * outer["baryonic_kms"]) ** 2 / radius
    expected = gbar + np.sqrt(gbar * value["dwarf"]["measured_disk_force_inversion"])
    assert np.isclose(
        np.sqrt(radius * expected) / 1000, outer["circular_kms"], atol=1e-10
    )
    reader = json.loads((measure.HERE / "results.json").read_text())
    for name, row in reader["surveys"].items():
        w, b = np.array(row["positive_weights"]), np.array(row["centres"])
        observed = value["shear_inversions"][name]
        reconstructed = (
            w
            @ (4 / math.pi * b + np.sqrt(b * observed["centre_a0"]))
            / measure.lensing_profile.FOUR_G
        )
        assert math.isclose(reconstructed, row["weighted_signal"], rel_tol=1e-12)
        assert (
            observed["generous_shear_floor"]
            > value["dwarf"]["spherical_generous_upper"]
        )
