"""Astronomical consequences of a level-202 horizon and a flat dust + Lambda background."""

import json
from math import sqrt

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2

from discrete_spherical_refinement import report as screen_report

# Loubser (2025), equations 4–6 and the following observed correlation coefficients.
Z = np.array([0.46, 0.67, 0.83])
H = np.array([88.48, 119.45, 108.28])
STAT = np.array([0.57, 6.39, 10.07])
SYST = np.array([12.32, 16.64, 15.08])
CORR = np.array([[1.0, 0.932, 0.830], [0.932, 1.0, 0.776], [0.830, 0.776, 1.0]])
MPC_KM = 3.0856775814913673e19
C_CM_S = 29979245800
TEN_YEARS_S = 10 * 365.25 * 86400


def forecast(vacuum, anchor):
    matter = (anchor**2 - vacuum**2) / (1 + Z[0]) ** 3
    if matter <= 0:
        raise ValueError("The clock anchor must leave positive pressureless matter")
    return {"anchor_H": anchor, **background(vacuum, matter)}


def background(vacuum, matter):
    present = sqrt(vacuum**2 + matter)

    def expansion(z):
        return sqrt(vacuum**2 + matter * (1 + z) ** 3)

    return dict(
        matter_coefficient=matter,
        H0=present,
        q0=(matter / 2 - vacuum**2) / present**2,
        redshift_at_zero_acceleration=(2 * vacuum**2 / matter) ** (1 / 3) - 1,
        H_at_other_redshifts=[expansion(z) for z in Z[1:]],
        drift_cm_s_per_ten_julian_years={
            str(z): C_CM_S * TEN_YEARS_S / MPC_KM * (present - expansion(z) / (1 + z))
            for z in (1, 4)
        },
    )


def profile(vacuum, covariance):
    cholesky = np.linalg.cholesky(covariance)

    def residual(params):
        predicted = np.sqrt(vacuum**2 + params[0] * (1 + Z) ** 3)
        return solve_triangular(cholesky, H - predicted, lower=True)

    fit = least_squares(residual, [1500.0], bounds=(0, np.inf), ftol=1e-13, xtol=1e-13)
    if not fit.success:
        raise ValueError(fit.message)
    statistic = float(fit.fun @ fit.fun)
    return dict(
        matter_coefficient=fit.x[0],
        chi_squared=statistic,
        nominal_chi_squared_2_tail=chi2.sf(statistic, 2),
    )


def measure():
    vacuum = float(screen_report()["level202"]["H_dS_s_minus_1"]) * MPC_KM
    sigma = np.hypot(STAT, SYST)
    covariance = CORR * np.outer(sigma, sigma)
    central = forecast(vacuum, H[0])
    predicted = np.array(central["H_at_other_redshifts"])
    anchor_derivative = ((1 + Z[1:]) / (1 + Z[0])) ** 3 * H[0] / predicted
    jacobian = np.column_stack([-anchor_derivative, np.eye(2)])
    residual_cov = jacobian @ covariance @ jacobian.T
    return dict(
        inputs=dict(
            chronometer_source="https://arxiv.org/html/2511.02730v1#S3.SS3",
            redshift_drift_source="https://arxiv.org/abs/astro-ph/9802122",
            screen_source="discrete_spherical_refinement.py",
            redshifts=Z.tolist(),
            H_observed=H.tolist(),
            stat_sigma=STAT.tolist(),
            syst_sigma=SYST.tolist(),
            correlations=CORR.tolist(),
            units_H="km/s/Mpc",
        ),
        assumptions="Observed horizon benchmark selects L=202. The proposed physical cell area and de Sitter entropy set the vacuum contribution. Spatially flat Friedmann evolution, conserved pressureless matter, constant vacuum density, negligible radiation, and comoving clocks and spectral tracers.",
        vacuum_H=vacuum,
        central_forecast_using_lowest_z=central,
        forecasts_at_minus_plus_one_anchor_sigma=[
            forecast(vacuum, H[0] + sign * sigma[0]) for sign in (-1, 1)
        ],
        measured_minus_anchored_prediction=(H[1:] - predicted).tolist(),
        linearized_residual_covariance=residual_cov.tolist(),
        residual_sigma=np.sqrt(residual_cov.diagonal()).tolist(),
        all_bin_profiles={
            name: profile(rate, covariance)
            for name, rate in (
                ("level201", 2 * vacuum),
                ("level202", vacuum),
                ("level203", vacuum / 2),
                ("zero_vacuum", 0.0),
            )
        },
        inference="Chronometer anchor and tests use different redshift bins with correlated uncertainties. The residual covariance propagates the anchor to first order. The separate profiles fit one nonnegative matter coefficient to all three H values with the published total covariance held fixed. Tail areas use an approximate Gaussian chi-squared reference with two degrees of freedom.",
    )


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
