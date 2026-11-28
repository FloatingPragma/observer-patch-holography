import json
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import beta as beta_function

from cosmology.rar_deep_regime.m31_transport import tracer as t


def test_jeans_projector_recovers_flat_potential_for_radial_and_tangential_orbits():
    gamma = 3.2
    v = 211.0
    for beta in (-0.5, 0, 0.5):
        measured = t.jeans(lambda r: np.full_like(r, v**2), 33, gamma, beta)
        analytic = v**2 * (1 - beta * (gamma - 1) / gamma) / (gamma - 2 * beta)
        assert measured == pytest.approx(analytic, rel=3e-7)


def isotropic_direct(v2, R, gamma=3.2, breaks=()):
    # Interchange the order of radial-pressure and line-of-sight integration.
    # This one-dimensional formula is independent of the nested quadrature.
    points = [R / r for r in breaks if r > R]
    num = quad(
        lambda x: x ** (gamma - 2) * np.sqrt(1 - x * x) * v2(R / x),
        0,
        1,
        points=points,
        epsabs=0.002,
        limit=200,
    )[0]
    den = beta_function((gamma - 1) / 2, 0.5) / 2
    return num / den


def test_shell_floor_and_extended_force_against_isotropic_single_integral():
    R = t.R_M31[1]
    for name, v2 in t.fields().items():
        assert t.jeans(v2, R) == pytest.approx(isotropic_direct(v2, R), rel=7e-7), name
    rows = json.loads(t.m.RESULT.read_text())["surveys"]["Chemin"]["rows"]
    rings, incs = t.mass_floor(rows, 1)

    def v2(r):
        return incs[rings <= r].sum() / r

    assert t.floor_extra(rings, incs, R, 3.2, 0) == pytest.approx(
        isotropic_direct(v2, R, breaks=rings), rel=8e-9
    )


def test_quadrature_converges_and_measured_curve_extends_analytically():
    R = t.R_M31[1]
    fields = t.fields()
    for beta in (-0.5, 0, 0.5):
        a = t.jeans(fields["finite"], R, 2.9, beta)
        b = t.jeans(fields["finite"], R, 2.9, beta, n=144)
        assert np.sqrt(a) == pytest.approx(np.sqrt(b), abs=0.002)
    setting = json.loads(t.m.RESULT.read_text())
    r = 1e7
    assert fields["finite"](r) * r / t.m.G == pytest.approx(
        (1 + setting["external_halo_dark_to_baryon_mass"])
        * setting["total_baryonic_injector_mass_msun"]
        + 1.4e8,
        rel=1e-10,
    )
    assert fields["lossless"](r) == pytest.approx(
        t.m.G * setting["lossless_outer_mass_slope_msun_kpc"], rel=3e-6
    )
