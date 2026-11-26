"""Monotone shell and annulus sums for the all-upper interior force case."""

from math import pi, sqrt

import numpy as np
from scipy.special import chndtr

from . import measure as m


def annulus(x, radial, field, upper):
    # The squared radial offset after isotropic Gaussian seeing is noncentral chi-square.
    cdf = chndtr((radial / m.PSF) ** 2, 2, (x / m.PSF) ** 2)
    weights = np.diff(cdf)
    if upper:
        return np.dot(weights, field[:-1]) + (1 - cdf[-1]) * field[-1]
    return np.dot(weights, field[1:])


def upper_limit():
    scale, _, _ = m.lens.geometry()
    acceleration = m.lens.acceleration()[1][1]
    masses = 10 ** np.array([10.94, 11.32])
    shell = np.r_[0, np.geomspace(0.001, m.RING, 3001)]
    low, high = shell[:-1], shell[1:]
    mass_hi = m.lens.fractions(high[:, None])[0] @ masses
    # The central spherical density bounds the mass enclosed by the first shell.
    central_mass_coefficient = (2 / (3 * pi)) * sum(
        masses * (m.lens.EDGE + m.lens.CORE) / (m.lens.EDGE * m.lens.CORE) ** 2
    )
    mass_over_radius = mass_hi / np.where(low > 0, low, 1.0)
    mass_over_radius[0] = central_mass_coefficient * high[0] ** 2
    mass_hi[0] = central_mass_coefficient * high[0] ** 3
    speed_hi = m.lens.G * mass_over_radius / scale + np.sqrt(
        m.lens.G * acceleration * m.lens.KPC_KM / 1000 * mass_hi
    )
    coeff = m.density(low) * speed_hi * np.diff(shell)
    radial = np.linspace(0, 6, 6001)
    inner_j = np.zeros_like(radial)
    for start in range(0, np.searchsorted(radial, m.RING), 128):
        r = radial[start : start + 128, None]
        # The shell kernel increases with source radius. Intrinsic J decreases in r.
        inner_j[start : start + len(r)] = (
            np.sqrt(np.maximum(0, 1 - (r / high) ** 2)) @ coeff
        )
    light = m.intensity(radial)

    outer = np.linspace(*m.OUTER_RANGE, 65)
    central = np.linspace(0, m.CENTRAL_END, 33)
    tail = max(
        m.OUTER_UPPER**2
        * annulus(x, radial, light, True)
        / (1 - (y**2 + 2 * m.PSF**2) / m.RING**2)
        for x, y in zip(outer[:-1], outer[1:], strict=True)
    )
    limit = max(
        (annulus(x, radial, inner_j, True) + tail) / annulus(y, radial, light, False)
        for x, y in zip(central[:-1], central[1:], strict=True)
    )
    return sqrt(limit)


def added_mass_floor(central_lower, model_upper):
    """Bound the supremum of positive added spherical orbital mass inside the outer ring."""
    scale = m.lens.geometry()[0]
    radial = np.linspace(0, 6, 6001)
    brightness_lo = annulus(m.CENTRAL_END, radial, m.intensity(radial), False)
    shell = np.r_[0, np.geomspace(0.001, m.RING, 4001)]
    low, high = shell[:-1], shell[1:]
    # PSF mass inside a disk <= area times the global peak of the Gaussian.
    # After division by source radius, min(s/(2 sigma²), 1/s) bounds every kernel.
    inv_low = np.divide(1.0, low, where=low > 0, out=np.full_like(low, np.inf))
    kernel_over_s = np.minimum(high / (2 * m.PSF**2), inv_low)
    integral_hi = np.dot(m.density(low) * kernel_over_s, np.diff(shell))
    return (
        (central_lower**2 - model_upper**2)
        * brightness_lo
        * scale
        / (m.lens.G * integral_hi)
    )


if __name__ == "__main__":
    limit = upper_limit()
    print(f"Monotone shell and annulus upper sum: {limit:.6f} km/s.")
    assert limit < 262
