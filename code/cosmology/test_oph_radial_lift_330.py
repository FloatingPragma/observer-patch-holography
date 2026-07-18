from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import spherical_jn

from oph_radial_lift_330 import (
    approximate_dilation_shape_bound,
    build_radial_receipt,
    derivative_mellin_norm,
    dilation_intertwiner_receipt,
    finite_window_stability_bound,
    forward_residual,
    mellin_spherical_bessel_square,
    minimum_prior_continuation,
    primordial_amplitude_from_screen,
    primordial_amplitude_log_sensitivity,
    radial_null_space_report,
    radial_projection_matrix,
    screen_amplitude_from_primordial,
    screen_gamma_ratio_cl,
    source_powerlaw,
    thin_shell_powerlaw_cl,
    window_powerlaw_cl_quadrature,
)


def test_scale_invariant_mellin_identity() -> None:
    for ell in range(2, 12):
        expected = 1.0 / (2.0 * ell * (ell + 1.0))
        assert mellin_spherical_bessel_square(ell, 0.0) == pytest.approx(expected, rel=2e-14)


def test_mellin_formula_against_direct_quadrature() -> None:
    ell = 2
    theta = 0.8
    integrand = lambda x: x ** (-theta - 1.0) * spherical_jn(ell, x) ** 2
    numeric = quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-11, limit=500)[0]
    numeric += quad(integrand, 1.0, 2000.0, epsabs=1e-12, epsrel=1e-10, limit=1200)[0]
    exact = mellin_spherical_bessel_square(ell, theta)
    # The omitted positive tail is O(X^{-theta-2}) and is below this tolerance.
    assert numeric == pytest.approx(exact, rel=5e-8)


def test_source_native_conversion_factor() -> None:
    theta = 1.630968209403959 / 48.0
    result = primordial_amplitude_from_screen(
        1.0, theta, Z_q=1.0, R_star=1.0, k_pivot=1.0
    )
    assert result.conversion_factor == pytest.approx(
        0.16080676040273598, rel=2e-14
    )


def test_derivative_norm_is_positive() -> None:
    for ell in (2, 3, 8):
        for theta in (0.1, 0.5, 1.5):
            assert derivative_mellin_norm(ell, theta) > 0.0


def test_exact_amplitude_round_trip_and_general_pivot() -> None:
    theta = 0.03397850436258248
    A_zeta = 2.2e-9
    Z_q = 0.87
    R_star = 14_000.0
    k_pivot = 0.05
    A_q = screen_amplitude_from_primordial(
        A_zeta, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    lifted = primordial_amplitude_from_screen(
        A_q, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    assert lifted.A_zeta == pytest.approx(A_zeta, rel=5e-15)
    assert lifted.kR_pivot == pytest.approx(k_pivot * R_star)


def test_thin_shell_formula_matches_gamma_ratio() -> None:
    theta = 0.31
    A_zeta = 3.1e-8
    Z_q = 1.2
    R_star = 2.7
    k_pivot = 0.8
    ell = np.arange(2, 20)
    A_q = screen_amplitude_from_primordial(
        A_zeta, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    expected = screen_gamma_ratio_cl(ell, A_q, theta)
    actual = thin_shell_powerlaw_cl(
        ell, A_zeta, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=0.0)


def test_amplitude_sensitivity_is_finite() -> None:
    value = primordial_amplitude_log_sensitivity(0.03397850436258248, 1.0)
    assert math.isfinite(value)
    assert 0.0 < value < 1.0


def test_finite_window_bound_contains_numerical_difference() -> None:
    ell = 2
    theta = 0.8
    A_zeta = 1.7
    Z_q = 1.0
    k_pivot = 1.0
    R_star = 1.0
    radii = np.array([0.97, 1.0, 1.03])
    radial_weights = np.array([0.25, 0.5, 0.25])

    # Logarithmic quadrature.  The chosen theta makes both tails small on this grid.
    logk = np.linspace(-16.0, 16.0, 180_001)
    k = np.exp(logk)
    dlnk = np.gradient(logk)
    window_cl = window_powerlaw_cl_quadrature(
        ell,
        A_zeta,
        theta,
        Z_q=Z_q,
        k_pivot=k_pivot,
        k=k,
        dlnk_weights=dlnk,
        radii=radii,
        radial_weights=radial_weights,
    )
    shell_cl = thin_shell_powerlaw_cl(
        ell,
        A_zeta,
        theta,
        Z_q=Z_q,
        R_star=R_star,
        k_pivot=k_pivot,
    )
    cert = finite_window_stability_bound(
        ell,
        theta,
        A_zeta=A_zeta,
        Z_q=Z_q,
        k_pivot=k_pivot,
        R_star=R_star,
        radii=radii,
        radial_weights=radial_weights,
    )
    # Allow a tiny quadrature error in the comparison.
    assert abs(window_cl - shell_cl) <= cert.absolute_cl_bound + 2e-6 * shell_cl


def test_finite_radial_operator_has_reported_nullity() -> None:
    ell = np.arange(2, 9)
    k = np.geomspace(1e-2, 20.0, 24)
    logk = np.log(k)
    A = radial_projection_matrix(
        ell,
        k,
        np.gradient(logk),
        Z_q=1.0,
        radii=[1.0],
        radial_weights=[1.0],
    )
    report = radial_null_space_report(A)
    assert report.rank <= len(ell)
    assert report.nullity == len(k) - report.rank
    assert len(report.null_basis) == report.nullity


def test_minimum_prior_continuation_and_resolution() -> None:
    # Full-row-rank toy map with a two-dimensional null space.
    A = np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]])
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    C = A @ truth
    p0 = np.zeros(4)
    Q = np.diag([1.0, 2.0, 3.0, 4.0])
    result = minimum_prior_continuation(
        A, C, prior_center=p0, prior_precision=Q
    )
    p = np.asarray(result.p)
    R = np.asarray(result.resolution)
    N = np.asarray(result.null_projector)
    np.testing.assert_allclose(A @ p, C, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(R @ R, R, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(A @ N, np.zeros_like(A), rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(R + N, np.eye(4), rtol=0.0, atol=1e-12)


def test_source_powerlaw_forward_residual_is_not_a_fit() -> None:
    ell = np.arange(2, 10)
    k = np.geomspace(1e-4, 1e4, 4000)
    logk = np.log(k)
    weights = np.gradient(logk)
    theta = 0.7
    A_zeta = 0.8
    A = radial_projection_matrix(
        ell,
        k,
        weights,
        Z_q=1.0,
        radii=[1.0],
        radial_weights=[1.0],
    )
    p = source_powerlaw(k, A_zeta, theta, 1.0)
    C_exact = np.asarray(
        thin_shell_powerlaw_cl(
            ell, A_zeta, theta, Z_q=1.0, R_star=1.0, k_pivot=1.0
        )
    )
    residual = forward_residual(A, p, C_exact)
    assert residual["relative_l2_residual"] < 2e-3


def test_dilation_intertwiner_receipt_accepts_powerlaw_and_rejects_wiggle() -> None:
    k = np.geomspace(1e-4, 1e4, 5000)
    theta = 0.23
    p = source_powerlaw(k, 2.0, theta, 1.0)
    good = dilation_intertwiner_receipt(
        k, p, theta, scale_ratios=[1.1, 1.7, 2.0], tolerance=2e-9
    )
    assert good.passed
    wiggle = p * np.exp(0.03 * np.sin(3.0 * np.log(k)))
    bad = dilation_intertwiner_receipt(
        k, wiggle, theta, scale_ratios=[1.7, 2.0], tolerance=1e-3
    )
    assert not bad.passed


def test_approximate_dilation_bound() -> None:
    k = np.geomspace(0.1, 10.0, 101)
    eps = np.full_like(k, 0.02)
    bound = approximate_dilation_shape_bound(k, eps, k_pivot=1.0)
    expected = 0.02 * np.abs(np.log(k))
    np.testing.assert_allclose(bound, expected, rtol=2e-3, atol=2e-5)


def test_receipts_fail_closed() -> None:
    clean_dag = {"nodes": [{"id": "collar_source"}, {"id": "mode_intertwiner"}]}
    clean = build_radial_receipt(
        receipt="SCR330_RADIAL_DILATION_INTERTWINER_RECEIPT",
        passed=True,
        claim_tier="E4",
        source_dag=clean_dag,
    )
    assert clean["passed"]

    dirty = build_radial_receipt(
        receipt="SCR330_RADIAL_PROMOTION_RECEIPT",
        passed=True,
        claim_tier="E4",
        source_dag={"nodes": [{"id": "planck", "measurement": True}]},
    )
    assert not dirty["passed"]
    assert "measurement_fit_or_likelihood_ancestor" in dirty["blockers"]

    early_transfer = build_radial_receipt(
        receipt="SCR330_TRANSFER_FIREWALL_RECEIPT",
        passed=True,
        claim_tier="E4",
        source_dag=clean_dag,
        physical_tt_te_ee_claim=True,
    )
    assert not early_transfer["passed"]
    assert "tt_te_ee_claim_before_E5" in early_transfer["blockers"]
