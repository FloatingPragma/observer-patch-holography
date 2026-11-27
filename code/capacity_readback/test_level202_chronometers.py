"""Check the clock extrapolation and correlated profile by independent formulas."""

from math import exp
import numpy as np
import pytest
from scipy.optimize import brentq

from level202_chronometers import CORR, H, STAT, SYST, Z, measure


def test_correlated_profile_and_clock_kinematics():
    r = measure()
    h = r["vacuum_H"]
    f = r["central_forecast_using_lowest_z"]
    x = (1 + Z) ** 3
    measured_matter = (H[0] ** 2 - h * h) / x[0]
    assert np.allclose(
        np.sqrt(h * h + measured_matter * x[1:]), f["H_at_other_redshifts"]
    )
    t = f["redshift_at_zero_acceleration"]
    assert measured_matter * (1 + t) ** 3 == pytest.approx(2 * h * h)
    assert f["drift_cm_s_per_ten_julian_years"]["4"] == pytest.approx(
        -6.296454, abs=1e-6
    )
    cov = np.diag(np.hypot(STAT, SYST)) @ CORR @ np.diag(np.hypot(STAT, SYST))

    def gradient(a):
        pred = np.sqrt(h * h + a * x)
        return (x / pred) @ np.linalg.solve(cov, pred - H)

    fit = brentq(gradient, 0.0, 1e4)
    residual = np.sqrt(h * h + fit * x) - H
    stat = residual @ np.linalg.solve(cov, residual)
    assert r["all_bin_profiles"]["level202"][
        "nominal_chi_squared_2_tail"
    ] == pytest.approx(exp(-stat / 2))
    assert stat == pytest.approx(8.202964, abs=1e-6)


def test_window_rate_against_integrated_clock_and_independent_spline_basis():
    from scipy.integrate import quad
    from level202_chronometer_windows import measure as window_report, SUB_Z, BREAKS

    result = window_report()
    m = result["background"]["matter_coefficient"]
    h = measure()["vacuum_H"]
    ages = [
        -quad(lambda x: 1 / ((1 + x) * np.sqrt(h * h + m * (1 + x) ** 3)), 0.30, z)[0]
        for z in SUB_Z
    ]
    # Columns give an intercept plus the length in each distinct redshift segment.
    b0, b1 = BREAKS
    design = np.column_stack(
        [
            np.ones_like(SUB_Z),
            np.minimum(SUB_Z, b0) - 0.30,
            np.clip(SUB_Z - b0, 0, b1 - b0),
            np.maximum(0, SUB_Z - b1),
        ]
    )
    s = np.linalg.lstsq(design, ages, rcond=None)[0][1:]
    assert -1 / ((1 + Z) * s) == pytest.approx(
        [88.48, 106.520897, 125.964595], abs=1e-6
    )
    assert result["marginal_standardized_residuals"] == pytest.approx(
        [1.98753, -1.42242], abs=1e-5
    )
    assert result["approximate_tail_two_dof"] == pytest.approx(0.03957186, abs=1e-7)


def test_observer_zero_periods_and_chronometer_coverage():
    from level202_chronometer_windows import observer_zero_pulse_readout

    result = observer_zero_pulse_readout()
    assert result["receiver"] == [0, 0]
    assert result["chronometer_H_km_s_Mpc"] is None
    panels = {p["query_tick"]: p for p in result["panels"]}
    assert set(panels) == {0, 4, 16, 96}
    assert panels[0]["pulse_shift_counts"] == {"0": 92, "1/4": 5}
    assert panels[96]["pulse_shift_counts"] == {"0": 96}
    assert panels[96]["complete_pulse_and_echo"] == 95
    for panel in panels.values():
        received = sum(panel["pulse_shift_counts"].values())
        assert received + panel["missing_pulse_pairs"] == result["enrolled_beacons_per_epoch"]
        assert panel["pulse_pairs_in_galaxy_fit_range"] == 0
    assert result["enrolled_beacons_per_epoch"] == 274
