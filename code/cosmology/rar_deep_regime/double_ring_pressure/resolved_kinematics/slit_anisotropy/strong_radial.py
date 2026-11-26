"""Replay one orbital second-moment model against three independent apertures."""

import json
from math import sqrt, log

import numpy as np
from scipy.integrate import simpson
from scipy.special import i0e, roots_legendre

from .. import measure as m

ANISOTROPY_RADIUS = 0.4
SHELL_MASS = 8e10


def projected_kernel(sky, s, a):
    z = np.sqrt(np.maximum(0, s**2 - sky**2))
    den = a**2 + sky**2
    # Integral of (z²+a²)/(z²+R²+a²)² along the positive half sightline.
    return (
        (s * s + a * a)
        / s
        * (
            (sky**2 + 2 * a * a) / (2 * den**1.5) * np.arctan(z / np.sqrt(den))
            - sky**2 * z / (2 * den * (s * s + a * a))
        )
    )


def fields(sky, n, extent):
    edge = np.r_[np.geomspace(1e-5, m.RING, n), np.geomspace(m.RING, extent, n)[1:]]
    s = np.sqrt(edge[:-1] * edge[1:])
    scale = m.lens.geometry()[0]
    vc2 = m.circular_speed_squared(
        s, (10.94, 11.32), scale, m.lens.acceleration()[1][1]
    )
    vc2 += m.lens.G * SHELL_MASS * (s >= m.RING) / (scale * s)
    weights = m.density(s) * vc2 * np.diff(edge)
    out = np.empty_like(sky)
    for i in range(0, len(sky), 64):
        rr = sky[i : i + 64, None]
        out[i : i + len(rr)] = projected_kernel(rr, s, ANISOTROPY_RADIUS) @ weights
    return out


def smooth(sky, field, offsets, seeing):
    kernel = (
        sky[None, :]
        * np.exp(-0.5 * ((sky[None, :] - offsets[:, None]) / seeing) ** 2)
        * i0e(sky[None, :] * offsets[:, None] / seeing**2)
    )
    return simpson(kernel * field, x=sky, axis=1) / seeing**2


def evaluate(n=2001, extent=1000):
    sky = np.linspace(0, 5, 1501)
    moment = fields(sky, n, extent)
    light = m.intensity(sky)
    u, w = roots_legendre(12)
    offsets = np.hypot(0.5 * u[:, None], 0.105 * u).ravel()
    weights = np.outer(w, w).ravel()
    seeing = 0.8 / sqrt(8 * log(2))
    keck = sqrt(
        (weights @ smooth(sky, moment, offsets, seeing))
        / (weights @ smooth(sky, light, offsets, seeing))
    )
    source = json.loads((m.HERE / "muse.json").read_text())
    center_pixel = source["central_and_nearest_to_one_arcsec"][1]["marker_pixel"][0]
    tick = source["tick_x_pixels"]
    outer = (center_pixel - tick[0]) * 2.5 / (tick[-1] - tick[0])
    central = np.linspace(0, m.CENTRAL_END, 65)
    x = np.r_[central, outer]
    j = smooth(sky, moment, x, m.PSF)
    i = smooth(sky, light, x, m.PSF)
    # Both profiles decrease. Compare the right numerator with the left denominator.
    lower = np.sqrt(np.min(j[1 : len(central)] / i[: len(central) - 1]))
    return dict(
        anisotropy_radius_arcsec=ANISOTROPY_RADIUS,
        exterior_shell_Msun=SHELL_MASS,
        keck_central_km_s=round(keck, 3),
        muse_outer_at_luminosity_center_km_s=round(sqrt(j[-1] / i[-1]), 3),
        held_muse_central_at_origin_km_s=round(sqrt(j[0] / i[0]), 3),
        held_muse_central_footprint_lower_km_s=round(float(lower), 3),
    )
