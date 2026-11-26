"""Test the frozen force law on independent, resolved stellar second moments."""

import argparse
import json
from math import ceil, floor, log, sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import i0e, roots_legendre

from cosmology.rar_deep_regime.double_ring_pressure import measure as lens

HERE = Path(__file__).resolve().parent
FWHM = 0.5
PSF = FWHM / sqrt(8 * log(2))
LIGHT = 10 ** (-0.4 * np.array([18.75, 17.15]))
LIGHT /= LIGHT.sum()
# The nominal width and source-redshift distance are shared with the aperture calculation.
RING = lens.THETA[1]
CENTRAL_END = 0.16
OUTER_RANGE = (0.64, 1.37)
OUTER_UPPER = 247.0


def intensity(radius):
    return sum(
        (weight * np.pi / (2 * (end - core)))
        * (1 / np.hypot(radius, core) - 1 / np.hypot(radius, end))
        for weight, core, end in zip(LIGHT, lens.CORE, lens.EDGE, strict=True)
    )


def density(radius):
    return sum(
        weight * (end + core) / ((radius**2 + core**2) * (radius**2 + end**2))
        for weight, core, end in zip(LIGHT, lens.CORE, lens.EDGE, strict=True)
    )


def circular_speed_squared(radius, log_masses, scale, acceleration):
    masses = 10 ** np.asarray(log_masses)
    enclosed = sum(
        m
        * 2
        / np.pi
        * (end * np.arctan(radius / end) - core * np.arctan(radius / core))
        / (end - core)
        for m, core, end in zip(masses, lens.CORE, lens.EDGE, strict=True)
    )
    return lens.G * enclosed / (scale * radius) + np.sqrt(
        lens.G * enclosed * acceleration * lens.KPC_KM / 1000
    )


def los_moment(sky, log_masses, scale, acceleration, cutoff, order):
    nodes, weights = roots_legendre(order)
    u, weights = (nodes + 1) / 2, weights / 2
    if cutoff:
        span = np.sqrt(np.maximum(0, cutoff**2 - sky**2))[:, None]
        z, jac = span * u, span
    else:
        # z=(u/(1-u))^2 covers the unrestricted sightline with central resolution.
        z, jac = (u / (1 - u)) ** 2, 2 * u / (1 - u) ** 3
    r = np.hypot(sky[:, None], z)
    speed = circular_speed_squared(r, log_masses, scale, acceleration)
    return (density(r) * speed * (z / r) ** 2 * jac) @ weights


def blurred(x, raw, field, seeing=PSF):
    kernel = raw * np.exp(-0.5 * ((raw - x) / seeing) ** 2) * i0e(raw * x / seeing**2)
    return simpson(field * kernel, x=raw) / seeing**2


def evaluate(log_masses, scale, acceleration, marks, resolution=2401, order=128):
    sky = np.linspace(0.0, 6.0, resolution)
    light = intensity(sky)
    full = los_moment(sky, log_masses, scale, acceleration, None, order)
    inside = los_moment(sky, log_masses, scale, acceleration, RING, order)
    predicted = np.array(
        [sqrt(blurred(r, sky, full) / blurred(r, sky, light)) for r in marks[:, 0]]
    )
    outer = np.linspace(*OUTER_RANGE, 65)
    central = np.linspace(0, CENTRAL_END, 33)
    # The positive interior force alone contributes to these moments. The discarded
    # outer-bin interior term makes the upper bound valid for all smaller interior masses.
    max_tail = max(
        OUTER_UPPER**2 * blurred(x, sky, light) / (1 - (y**2 + 2 * PSF**2) / RING**2)
        for x, y in zip(outer[:-1], outer[1:], strict=True)
    )
    limit = max(
        (blurred(x, sky, inside) + max_tail) / blurred(y, sky, light)
        for x, y in zip(central[:-1], central[1:], strict=True)
    )
    return {
        "central_full_continuation_km_s": round(float(predicted[0]), 3),
        "inner_ring_full_continuation_km_s": round(
            sqrt(
                blurred(lens.THETA[0], sky, full) / blurred(lens.THETA[0], sky, light)
            ),
            3,
        ),
        "median_predicted_minus_observed_km_s": round(
            float(np.median(predicted - marks[:, 1])), 3
        ),
        "median_absolute_residual_km_s": round(
            float(np.median(np.abs(predicted - marks[:, 1]))), 3
        ),
        "predictions_below_marker_count": int(sum(predicted < marks[:, 1])),
        "central_upper_arbitrary_positive_exterior_km_s": round(sqrt(limit), 3),
    }


def run():
    from .enclosure import added_mass_floor, upper_limit

    positive_limit = upper_limit()
    source = json.loads((HERE / "muse.json").read_text())
    points = np.array(source["radius_and_vrms"])
    assert len(points) == 41
    centre, outer = source["central_and_nearest_to_one_arcsec"]
    # Allow two image pixels in each coordinate, plus half the diagonal of a 0.2" spaxel.
    radial_pad = 0.1 * sqrt(2) + 2 / 210
    assert centre["radial_bar_arcsec"][1] + radial_pad < CENTRAL_END
    assert OUTER_RANGE[0] < outer["radial_bar_arcsec"][0] - radial_pad
    assert OUTER_RANGE[1] > outer["radial_bar_arcsec"][1] + radial_pad
    lower_central = floor(centre["moment_bar_km_s"][0] - 2 / 4.85)
    assert ceil(outer["moment_bar_km_s"][1] + 2 / 4.85) <= OUTER_UPPER
    scale, _, _ = lens.geometry()
    a0, limits = lens.acceleration()
    models = {name: (logs, a0) for name, logs in lens.LOG_MASSES.items()}
    models["upper_Chabrier"] = ((10.94, 11.32), limits[1])
    return {
        "resolved_markers": len(points),
        "monotone_shell_upper_km_s": round(positive_limit, 5),
        "outward_upper_km_s": ceil(positive_limit),
        "lower_added_interior_orbital_mass_Msun": floor(
            added_mass_floor(lower_central, ceil(positive_limit))
        ),
        "central_lower_graphical_bar_km_s": lower_central,
        "outer_upper_graphical_bar_km_s": OUTER_UPPER,
        "central_and_outer_radial_windows_arcsec": [
            [0, CENTRAL_END],
            list(OUTER_RANGE),
        ],
        "models": {
            name: evaluate(logs, scale, a, points) for name, (logs, a) in models.items()
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    path = HERE / "results.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        assert data == json.loads(path.read_text())
        print("Reproduced the resolved stellar-kinematic comparison.")
