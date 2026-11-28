"""Checks independent of the profiled distance calculation and the vacuum decoder."""

from math import sqrt

import numpy as np
from scipy.linalg import solve_triangular

from . import capacity_distances as bao
from .dust_or_vacuum import run as source_readings


def test_dust_distance_integral_and_exact_gaussian_contrasts():
    data = bao.DistanceData()
    a = 0.0004
    dm = 2 * (1 - 1 / np.sqrt(1 + data.z)) / sqrt(a)
    dh = (1 + data.z) ** -1.5 / sqrt(a)
    reference = np.where(
        data.kind == "DH_over_rs",
        dh,
        np.where(data.kind == "DM_over_rs", dm, np.cbrt(data.z * dh * dm**2)),
    )
    np.testing.assert_allclose(bao.distance_vector(data, a, 0), reference, rtol=2e-14)
    # Different redshifts are uncorrelated in this release. This sum is
    # independent of the production matrix projection and uses per-tracer 2x2 blocks.
    chisq = 0.0
    for z in data.z[data.kind == "DH_over_rs"]:
        (i,) = np.where((data.z == z) & (data.kind == "DM_over_rs"))[0]
        (j,) = np.where((data.z == z) & (data.kind == "DH_over_rs"))[0]
        f = dm[i] / dh[j]
        var = data.cov[i, i] - 2 * f * data.cov[i, j] + f**2 * data.cov[j, j]
        chisq += (data.obs[i] - f * data.obs[j]) ** 2 / var
    assert abs(bao.dust_contrasts(data)["chi2"] - chisq) < 1e-9


def test_capacity_prediction_and_held_out_calibration():
    data = bao.DistanceData()
    b = bao.capacity_B(bao.R_DRAG, bao.L_PLANCK)
    result = bao.prediction(data, b)
    assert abs(bao.X - 280.9654644692834) < 1e-12
    assert abs(bao.X * (bao.PRICE_BUDGET - bao.X) - 6 * bao.pi) < 1e-10
    i = result["anchor_row"]
    assert result["rows"][i]["kind"] == "DH_over_rs" and data.z[i] == 2.33
    a = result["all_data_A"]
    prediction = bao.distance_vector(data, a, b)
    w = solve_triangular(data.full, data.obs - prediction, lower=True)
    assert abs(w @ w - result["all_data_profile_chi2"]) < 1e-8
    assert 17.26 < w @ w < 17.28
    data.obs[0] += 1
    changed = bao.prediction(data, b)
    np.testing.assert_equal(
        [r["forecast"] for r in result["rows"]],
        [r["forecast"] for r in changed["rows"]],
    )


def test_measured_dust_and_vacuum_separate_under_refinement():
    cases = source_readings()
    expected = dict(
        ordinary=0, held=0, omit_counted_collar=4 / 11, half_collar_force=-2 / 11
    )
    for case, exact in expected.items():
        readings = cases[case]
        for name, limit in (
            ("dust_coefficient", exact),
            ("vacuum_intercept", 0),
            ("max_other_frame_error", 0),
        ):
            errors = [abs(r[name] - limit) for r in readings]
            assert all(new < old / 14 for old, new in zip(errors, errors[1:]))
            assert errors[-1] < 2e-7
    held = cases["omit_counted_collar"][-1]["residual_range"]
    assert held[0] < 0.37 and held[1] > 3.27


def test_joint_matter_coverage_and_finite_gaussian_ellipsoid():
    data = bao.DistanceData()
    b = bao.capacity_B(bao.R_DRAG, bao.L_PLANCK)
    result = bao.prediction(data, b)
    lo, hi = result["shared_A_in_all_13_simultaneous_bands"]
    assert lo < result["A_from_anchor"] < hi
    sigma = np.sqrt(np.diag(data.cov))

    def max_residual(a):
        return np.max(np.abs(data.obs - bao.distance_vector(data, a, b)) / sigma)

    k = result["simultaneous_box_k"]
    np.testing.assert_allclose([max_residual(lo), max_residual(hi)], [k, k], atol=1e-8)
    assert max_residual(lo - 1e-8) > k and max_residual(hi + 1e-8) > k
    assert max_residual((lo + hi) / 2) < k
    assert bao.prediction(data, 0)["shared_A_in_all_13_simultaneous_bands"] is None
    assert result["all_data_profile_chi2"] < result["gaussian_13D_95_threshold"]


def test_time_dilation_inequalities_and_measured_volume_radial_contrast():
    data = bao.DistanceData()
    result = bao.nonacceleration_test(data)
    assert result["family_size"] == 72
    strongest = result["strongest"]
    assert (strongest["kind"], strongest["z1"], strongest["z2"]) == (
        "volume_to_later_slope",
        0.295,
        0.934,
    )
    (i,) = np.flatnonzero(data.kind == "DV_over_rs")
    (j,) = np.flatnonzero((data.kind == "DH_over_rs") & (data.z == 0.934))
    size = np.cbrt((1 + data.z[i]) / (data.z[i] * np.log1p(data.z[i]) ** 2))
    speed = 1 + data.z[j]
    diff = size * data.obs[i] - speed * data.obs[j]
    error = sqrt(
        size**2 * data.cov[i, i]
        + speed**2 * data.cov[j, j]
        - 2 * size * speed * data.cov[i, j]
    )
    np.testing.assert_allclose(
        [strongest["contrast"], strongest["sigma"]], [diff, error]
    )
    assert 0.00101 < result["gaussian_family_p_upper"] < 0.00102
    # Populate the same observed distance coordinates from a coasting spacetime.
    # Then check pressureless flat matter, which strictly decelerates.
    coasting_dm = np.log1p(data.z)
    coasting_dh = 1 / (1 + data.z)
    masks = (data.kind == "DM_over_rs", data.kind == "DH_over_rs")
    data.obs = np.select(
        masks,
        (coasting_dm, coasting_dh),
        np.cbrt(data.z * coasting_dh * coasting_dm**2),
    )
    assert all(
        abs(r["contrast"]) < 1e-13 for r in bao.nonacceleration_test(data)["rows"]
    )
    data.obs = bao.distance_vector(data, 0.0004, 0)
    assert all(r["contrast"] > 0 for r in bao.nonacceleration_test(data)["rows"])


def test_acceleration_bound_with_unfixed_spatial_curvature():
    data = bao.DistanceData()
    result = bao.curved_nonacceleration_test(data)
    positive = [r for r in result["rows"] if r["gap"] > 0]
    assert len(positive) == 1 and positive[0]["z_radial"] == 0.934
    witness = positive[0]
    assert 0.0074 < witness["gap"] < 0.0075
    # Invert the constraint the other way. The maximum transverse distance
    # inferred from the observed near D_V lies below every allowed geometry.
    z = result["z_volume"]
    dm_upper = sqrt(
        result["measured_DV_upper"] ** 3 * (1 + z) / (z * witness["u_lower"])
    )
    assert dm_upper < witness["near_transverse_lower"]
    # Coasting direct images for three curvatures all satisfy the bound.
    # The closed example places the far marker past the equator, before its antipode.
    chi = 30 * np.log1p(data.z)
    dh = 30 / (1 + data.z)
    for dm in (chi, 15 * np.sin(chi / 15), 30 * np.sinh(chi / 30)):
        data.obs = np.select(
            (data.kind == "DM_over_rs", data.kind == "DH_over_rs"),
            (dm, dh),
            np.cbrt(data.z * dm**2 * dh),
        )
        assert all(r["gap"] < 0 for r in bao.curved_nonacceleration_test(data)["rows"])
