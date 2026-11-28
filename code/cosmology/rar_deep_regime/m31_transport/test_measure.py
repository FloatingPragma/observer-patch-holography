"""Independent checks of the profile, positive-source bounds and data ancestry."""

import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import roots_legendre

from matched_einstein_limit.collar_identifiability.refinement import transport
from cosmology.rar_deep_regime.m31_transport import aperture, measure as m


def test_finite_mark_chart_converges_to_declared_green_function():
    radii = np.array([6, 10, 14])
    mu = 0.15
    expected = np.exp(-mu * radii) / radii
    max_errors = []
    for spacing in (1, 0.5, 0.25):
        survival = 1 / (1 + (mu * spacing) ** 2 / 8)
        axis = transport.axial_visits(int(60 / spacing), radii / spacing, survival)
        field = np.pi * survival * np.array(axis) / (2 * spacing)
        max_errors.append(max(abs(field / expected - 1)))
    assert max_errors[1] < 0.26 * max_errors[0]
    assert max_errors[2] < 0.26 * max_errors[1]
    assert max_errors[2] < 0.0009


def test_yukawa_kernel_and_enclosed_total_from_independent_angular_integrals():
    ds = 0.002
    edge = np.arange(10001) * ds
    r = (edge[1:] + edge[:-1]) / 2
    src = SimpleNamespace(edge=edge, r=r)
    ix = [140, 897, 1600]
    weights = np.zeros_like(r)
    weights[ix] = [0.21, 0.32, 0.47]
    a = r[ix]
    p = weights[ix]
    mu = 0.19
    for radius in (0.05, 1.0, 4.2, 22.8):
        direct = 0
        for y, w in zip(a, p):

            def kernel(c):
                d = np.sqrt(radius**2 + y**2 - 2 * radius * y * c)
                return np.exp(-mu * d) / d

            direct += w * quad(kernel, -1, 1, epsabs=1e-12)[0] / 2
        found = m.ru(src, weights, mu, np.array([radius]))[0] / radius
        assert found == pytest.approx(direct, rel=1e-11)
        assert (
            0
            < found
            <= m.ru(src, weights, 0, np.array([radius]))[0] / radius
            <= 1 / radius
        )
    # Integral of the product of two Yukawa kernels over R^3 equals
    # (2 pi / mu) exp(-mu |y-z|). Average the absolute separation directly.
    expected = 0
    for y, wy in zip(a, p):
        for z, wz in zip(a, p):
            value = quad(
                lambda c: np.exp(-mu * np.sqrt(y * y + z * z - 2 * y * z * c)), -1, 1
            )[0]
            expected += wy * wz * value / (4 * mu)
    assert m.total_length(src, weights, mu) == pytest.approx(expected, abs=2e-7)


def test_bound_includes_high_loss_and_extended_shell_controls():
    shells = np.array([0.5, 3.0, 8.0, 17.5, 23.0])
    weights = np.array([0.40, 0.25, 0.15, 0.10, 0.10])
    src = SimpleNamespace(r=shells)
    r1, r2 = 15.0, 23.0
    bound = aperture.aperture_bounds(src, weights, r1, r2)

    # Integrate (2 mu r u)^2 directly with finite exponential differences, no
    # hyperbolic-function overflow, and cut quadrature at each emitting shell.
    def f(r, mu):
        contributions = np.exp(-mu * np.abs(r - shells)) * -np.expm1(
            -2 * mu * np.minimum(r, shells)
        )
        return (np.dot(weights / shells, contributions)) ** 2

    for mu in (0.0001, 0.01, 0.25, 0.99, 1.0, 2.33, 41.0):

        def integral(lo, hi):
            points = shells[(shells > lo) & (shells < hi)]
            return quad(f, lo, hi, points=points, args=(mu,), epsabs=1e-17)[0]

        ratio = 1 + integral(r1, r2) / integral(0, r1)
        assert ratio <= bound["upper_mass_ratio_all_lifetimes"]
        if mu >= 1:
            assert ratio <= bound["upper_mass_ratio_high_loss"]
    # Concentrating injections at an outer shell can produce a steep rise in
    # the strong-loss limit. The source-distribution condition matters.
    outer = weights.copy()
    outer[[0, 4]] = [0.01, 0.49]
    other = aperture.aperture_bounds(src, outer, r1, r2)
    assert (
        other["upper_mass_ratio_all_lifetimes"]
        > bound["upper_mass_ratio_all_lifetimes"]
    )


def test_radius_only_forecast_and_spatial_refinement():
    inputs = json.loads(m.INPUT.read_text())
    original = json.loads(m.RESULT.read_text())
    changed = copy.deepcopy(inputs)
    for rows in changed[
        "outer_radius_kpc_speed_error_kms_and_inclination_deg"
    ].values():
        for row in rows:
            row[1:] = [1000, 1, 30]
    altered = m.run(changed)
    fine = m.run(inputs, dr=0.025)
    for name, survey in original["surveys"].items():
        for orig, swap, refined in zip(
            survey["rows"],
            altered["surveys"][name]["rows"],
            fine["surveys"][name]["rows"],
        ):
            for key in (
                "ordinary_speed",
                "cosmic_share",
                "lossless",
                "all_radial_injectors_speed_ceiling",
            ):
                assert orig[key] == pytest.approx(swap[key], abs=1e-9)
                assert orig[key] == pytest.approx(refined[key], abs=0.012)
    fine_src = m.source(inputs["HI_radius_kpc_and_surface_msun_pc2"], 0.025)
    mass = 7.1e10 * fine_src.disk + 6.4e10 * fine_src.bulge + fine_src.gas
    fine_bound = aperture.aperture_bounds(fine_src, mass / mass.sum(), 28.17, 38.09)
    reference = json.loads((m.HERE / "aperture.json").read_text())["surveys"]["Chemin"][
        "source_cases"
    ]["disk_7.1_bulge_6.4_disk"]
    assert fine_bound["upper_mass_ratio_all_lifetimes"] == pytest.approx(
        reference["upper_mass_ratio_all_lifetimes"], abs=0.0002
    )


def test_gas_potential_independent_angular_quadrature_and_height():
    inputs = json.loads(m.INPUT.read_text())
    src = m.source(inputs["HI_radius_kpc_and_surface_msun_pc2"], 0.05)
    radii = np.array([28.17, 38.09])
    angles, w = roots_legendre(180)
    cos = np.cos(np.pi * (angles + 1))
    for rad in radii:

        def potential(x):
            dist2 = (
                (x - src.r[:, None]) ** 2 + 2 * x * src.r[:, None] * (1 - cos) + 0.5**2
            )
            return -m.G * src.gas @ (1 / np.sqrt(dist2) @ (w / 2))

        fd = rad * (potential(rad + 0.001) - potential(rad - 0.001)) / 0.002
        assert fd == pytest.approx(m.gas_force(src, np.array([rad]))[0], rel=3e-7)
    speed = m.predict(src, radii, 7.1e10, 6.4e10, 8.750402247718735e-11)["cosmic_share"]
    for h in (0.25, 1.0):
        modified = np.sqrt(
            speed**2 + m.gas_force(src, radii, h) - m.gas_force(src, radii)
        )
        assert max(abs(modified - speed)) < 0.2
