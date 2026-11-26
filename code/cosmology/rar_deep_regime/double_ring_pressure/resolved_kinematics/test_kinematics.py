"""Independent Jeans projection, seeing-kernel and archived data checks."""

import json
from math import hypot, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_hermitenorm

from . import measure as m
from .enclosure import added_mass_floor, upper_limit


def test_independent_jeans_projection():
    scale = m.lens.geometry()[0]
    logs = m.lens.LOG_MASSES["Chabrier"]
    accel = m.lens.acceleration()[0]

    def radial_pressure(r):
        return quad(
            lambda s: m.density(s)
            * m.circular_speed_squared(s, logs, scale, accel)
            / s,
            r,
            np.inf,
            epsabs=1e-8,
            epsrel=1e-7,
        )[0]

    for x in (0.17, 1.04):
        brightness = quad(lambda z: m.density(hypot(x, z)), 0, np.inf)[0]
        projected = quad(lambda z: radial_pressure(hypot(x, z)), 0, np.inf)[0]
        via_shell = m.los_moment(np.array([x]), logs, scale, accel, None, 128)[0]
        np.testing.assert_allclose(
            [brightness, projected], [m.intensity(x), via_shell], rtol=1e-7
        )


def test_exterior_shell_after_cartesian_seeing():
    nodes, weights = roots_hermitenorm(50)
    weights = np.outer(weights, weights) / (2 * np.pi)
    for x in (0.0, 0.64, 1.37):
        r = np.hypot(x + m.PSF * nodes[:, None], m.PSF * nodes[None, :])
        for s in (m.RING, 2 * m.RING, 100 * m.RING):
            kernel = np.sum(weights * np.sqrt(np.maximum(0, 1 - (r / s) ** 2)))
            assert 1 - (x**2 + 2 * m.PSF**2) / m.RING**2 <= kernel <= 1.000000000001


def test_archive_and_monotone_enclosure():
    bounds = m.run()
    assert bounds == json.loads((m.HERE / "results.json").read_text())
    value = upper_limit()
    assert (
        bounds["models"]["upper_Chabrier"][
            "central_upper_arbitrary_positive_exterior_km_s"
        ]
        < value
    )
    assert (
        value
        < bounds["outward_upper_km_s"]
        < bounds["central_lower_graphical_bar_km_s"]
    )


def test_conservative_equivalent_mass_floor():
    # Independent adaptive evaluation of the PSF-capped integral and Bessel brightness.
    kernel_integral = quad(
        lambda s: m.density(s) * min(1, s**2 / (2 * m.PSF**2)) / s,
        0,
        m.RING,
        points=[sqrt(2) * m.PSF],
    )[0]
    sky = np.linspace(0, 6, 3001)
    brightness = m.blurred(m.CENTRAL_END, sky, m.intensity(sky))
    limit = (
        (270**2 - 262**2)
        * brightness
        * m.lens.geometry()[0]
        / (m.lens.G * kernel_integral)
    )
    assert 8.4e8 < added_mass_floor(270, 262) < limit < 8.5e8
