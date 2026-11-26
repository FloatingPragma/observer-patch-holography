"""Check the projected angular inequalities against direct line integration."""

import json
from math import exp, hypot, log, sin, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog

from . import measure
from . import strong_radial


def integrated_kernel(R, s, beta, weight):
    Z = sqrt(s * s - R * R)
    return quad(
        lambda z: (1 - beta(hypot(R, z)) * R**2 / (R**2 + z * z))
        * weight(hypot(R, z), s)
        / s,
        0,
        Z,
    )[0]


def test_general_sign_changing_jeans_bounds():
    b = 0.17

    def beta(r):
        return b * sin(log(r))

    def weight(r, s):
        return exp(2 * b * (np.cos(log(r)) - np.cos(log(s))))

    for s in (0.35, 2.07, 9, 110):
        for f in (0.001, 0.25, 0.9, 0.99):
            R = f * s
            val = integrated_kernel(R, s, beta, weight) / sqrt(1 - f * f)
            assert (1 - b) / (1 + 2 * b) <= val <= (1 + b) / (1 - 2 * b)


def test_outer_angular_response_for_radial_stars():
    a = 0.85
    # Smooth radial profiles with an analytic integrating factor.
    for power in (0.0, 0.5, 1.0):

        def beta(r):
            return power * r * r / (r * r + a * a)

        def w(r, s):
            return ((a * a + s * s) / (a * a + r * r)) ** power

        for s in (2.07, 4, 12):
            axis = quad(lambda z: w(z, s), 0, s)[0] / s
            for R in (0.02, 0.7, 1.98):
                value = integrated_kernel(R, s, beta, w)
                lo = sqrt(1 - (R / measure.S) ** 2) / (1 + (R / a) ** 2) ** 2
                assert lo - 1e-13 <= value / axis <= 1 + 1e-13
                if power == 1:
                    np.testing.assert_allclose(
                        value, strong_radial.projected_kernel(R, s, a), rtol=1e-10
                    )

    # Deliberately nonmonotone beta within the same radius envelope.
    def beta(r):
        return (0.5 + 0.4 * sin(r)) * r * r / (r * r + a * a)

    def w(r, s):
        return exp(quad(lambda t: 2 * beta(t) / t, r, s)[0])

    for s in (2.07, 4):
        axis = quad(lambda z: w(z, s), 0, s)[0] / s
        R = 0.91
        assert (
            sqrt(1 - (R / measure.S) ** 2) / (1 + (R / a) ** 2) ** 2
            <= integrated_kernel(R, s, beta, w) / axis
            <= 1
        )


def test_moment_extremes_with_independent_lp():
    survey = object.__new__(measure.Survey)
    survey.outer = np.linspace(0.64, 1.37, 9)
    x = survey.outer[1:]
    a = np.array([20, 19, 14, 7, 10, 16, 21, 22]) * 1e3
    k = np.array([1, 4, 7, 2, 5, 3, 8, 6.0])
    # With z=p/<k> the original fractional program becomes a linear program.
    c = measure.muse.OUTER_UPPER**2 - a
    ans = linprog(
        -c, A_ub=[0.85 - x], b_ub=[0.0], A_eq=[k], b_eq=[1.0], bounds=(0, None)
    )
    assert ans.success
    np.testing.assert_allclose(-ans.fun, survey.exterior_bound(a, k, 0.85), rtol=1e-12)


def test_counted_apertures_and_refinement():
    out = measure.Survey().run()
    assert out == json.loads((measure.HERE / "results.json").read_text())
    assert out["uniform_anisotropy_envelope"]["central_threshold_abs_beta"] > 0.018
    assert out["radial_anisotropy_envelope"]["central_threshold_radius_arcsec"] < 3.11
    original = out["strong_radial_jeans_control"]
    refined = strong_radial.evaluate(2701, 1900)
    for key in (
        "keck_central_km_s",
        "muse_outer_at_luminosity_center_km_s",
        "held_muse_central_footprint_lower_km_s",
    ):
        assert abs(original[key] - refined[key]) < 0.05
    assert original["held_muse_central_footprint_lower_km_s"] > 290
