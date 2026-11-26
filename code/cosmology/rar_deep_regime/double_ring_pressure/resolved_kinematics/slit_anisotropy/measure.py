"""Bound two independent seeing-weighted stellar measurements with stellar anisotropy."""

import argparse
import json
from math import hypot, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import chndtr

from .. import measure as muse
from .strong_radial import evaluate as strong_radial


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "keck.json").read_text())
S = muse.RING
# Every center offset of Euclidean size <= 0.03" is covered.
CENTER_ERROR = 0.03


def blur_weights(offset, radii, sigma):
    """Noncentral chi-square probability of each circular annulus after seeing."""
    cdf = chndtr((radii[None, :] / sigma) ** 2, 2, (offset[:, None] / sigma) ** 2)
    return np.column_stack((np.diff(cdf, axis=1), 1 - cdf[:, -1]))


def shell_coefficients(edges, upper):
    scale = muse.lens.geometry()[0]
    a0 = muse.lens.acceleration()[1][int(upper)]
    logs = (10.94, 11.32) if upper else (10.79, 11.19)
    mass = 10 ** np.asarray(logs)
    low, high = edges[:-1], edges[1:]
    s = high if upper else low
    enclosed = muse.lens.fractions(s[:, None])[0] @ mass
    mass_over_r = enclosed / (np.where(low > 0, low, 1) if upper else high)
    if upper:
        central_coeff = (
            2
            / (3 * pi)
            * sum(
                mass
                * (muse.lens.EDGE + muse.lens.CORE)
                / (muse.lens.EDGE * muse.lens.CORE) ** 2
            )
        )
        enclosed[0] = central_coeff * high[0] ** 3
        mass_over_r[0] = central_coeff * high[0] ** 2
    vc2 = muse.lens.G * mass_over_r / scale + np.sqrt(
        muse.lens.G * a0 * muse.lens.KPC_KM / 1000 * enclosed
    )
    return np.diff(edges) * vc2 * muse.density(low if upper else high)


class Survey:
    """Positive sums over source shells, image annuli and rectangular slit cells."""

    def __init__(self):
        self.sky = np.linspace(0, 6, 3001)
        shell = np.r_[0, np.geomspace(0.001, S, 2001)]
        self.high = shell[1:]
        interior = self.sky < S
        self.inner_count = sum(interior)
        x = self.sky[interior, None]
        self.upper_kernel = np.sqrt(
            np.maximum(0, 1 - (x / self.high) ** 2)
        ) * shell_coefficients(shell, True)
        low = shell[:-1].copy()
        low[0] = 1.0  # The first lower shell has zero mass.
        self.lower_kernel = np.sqrt(
            np.maximum(0, 1 - (x / low) ** 2)
        ) * shell_coefficients(shell, False)
        self.outer = np.linspace(*muse.OUTER_RANGE, 65)
        self.m_weights = blur_weights(self.outer, self.sky, muse.PSF)
        self.light = muse.intensity(self.sky)
        self.outer_i_hi = self.m_weights[:-1] @ self.light
        self.offsets = np.linspace(0, 1.40, 281)
        self.k_weights = blur_weights(
            self.offsets, self.sky, DATA["seeing_fwhm_arcsec"] / sqrt(8 * log(2))
        )
        self.keck_i_lo = self.k_weights @ self.lower_values(self.light)

    @staticmethod
    def lower_values(values):
        return np.r_[values[1:], 0.0]

    def positive_image(self, interior):
        field = np.zeros(len(self.sky))
        field[: self.inner_count] = interior
        return field

    def aperture(self, offset, numerator, tail):
        # Reflection of a centered circular field lets us integrate just half the slit width.
        width = np.linspace(0, DATA["slit_width_arcsec"] / 2, 33)
        height = offset + np.linspace(-0.5, 0.5, 17) * DATA["extraction_height_arcsec"]
        y_low = np.where(
            height[:-1] * height[1:] <= 0,
            0,
            np.minimum(abs(height[:-1]), abs(height[1:])),
        )
        y_hi = np.maximum(abs(height[:-1]), abs(height[1:]))
        near = np.maximum(0, np.hypot(width[:-1, None], y_low) - CENTER_ERROR)
        far = np.hypot(width[1:, None], y_hi) + CENTER_ERROR
        left = np.searchsorted(self.offsets, near, side="right") - 1
        right = np.searchsorted(self.offsets, far)
        return sqrt((np.mean(numerator[left]) + tail) / np.mean(self.keck_i_lo[right]))

    def predictions(self, family, param, centroid=None):
        if family == "symmetric":
            b = param
            upper = (1 + b) / (1 - 2 * b)
            lower = (1 - b) / (1 + 2 * b)
            # B is the exterior isotropic shell strength after scaling by the uniform upper factor.
            outer_field = (
                np.sqrt(np.maximum(0, 1 - (self.sky / S) ** 2)) * lower / upper
            )
        else:
            a = param
            upper = 1 + (self.high / a) ** 2
            lower = 1 / (1 + (self.high / a) ** 2)
            outer_field = (
                np.sqrt(np.maximum(0, 1 - (self.sky / S) ** 2))
                / (1 + (self.sky / a) ** 2) ** 2
            )
        j_upper = self.positive_image(
            self.upper_kernel @ np.broadcast_to(upper, len(self.high))
        )
        j_lower = self.positive_image(
            self.lower_kernel @ np.broadcast_to(lower, len(self.high))
        )
        outer_j_lo = self.m_weights[1:] @ self.lower_values(j_lower)
        outer_q_lo = self.m_weights[1:] @ self.lower_values(outer_field)
        a = outer_j_lo / self.outer_i_hi
        k = outer_q_lo / self.outer_i_hi
        tail = self.exterior_bound(a, k, centroid)
        assert tail > 0
        keck_j_hi = self.k_weights @ j_upper
        return np.asarray(
            [
                self.aperture(abs(row[0]), keck_j_hi, tail)
                for row in DATA["measurements"]
            ]
        )

    def exterior_bound(self, a, k, centroid):
        y = self.outer[1:]
        upper = muse.OUTER_UPPER**2
        if centroid is None:
            return np.max((upper - a) / k)
        # Extreme probability measures with one lower mean constraint have at most two atoms.
        above = y >= centroid
        single = np.max(((upper - a) / k)[above])
        below = ~above
        lo, hi = y[below, None], y[None, above]
        weight = (hi - centroid) / (hi - lo)
        av = weight * a[below, None] + (1 - weight) * a[None, above]
        kv = weight * k[below, None] + (1 - weight) * k[None, above]
        return max(single, np.max((upper - av) / kv))

    def run(self):
        rows = DATA["measurements"]
        observed = np.array(
            [
                hypot(max(0, abs(v) - dv), sigma - dsigma)
                for _, v, dv, sigma, dsigma in rows
            ]
        )
        # Register the comparison using the central slit aperture, independent of its speed.
        centre = next(i for i, row in enumerate(rows) if row[0] == 0)

        def threshold(family, bracket, centroid=None):
            return brentq(
                lambda x: self.predictions(family, x, centroid)[centre]
                - observed[centre],
                *bracket,
            )

        nominal = self.predictions("symmetric", 0)
        return {
            "central_lower_keck_km_s": float(observed[centre]),
            "strong_radial_jeans_control": strong_radial(),
            "muse_outer_upper_km_s": muse.OUTER_UPPER,
            "rectangular_center_tolerance_arcsec": CENTER_ERROR,
            "same_weight_outer_centroid": {
                "minimum_mean_radius_arcsec": 0.85,
                "isotropic_upper_km_s": round(
                    self.predictions("symmetric", 0, 0.85)[centre], 3
                ),
                "threshold_abs_beta": round(threshold("symmetric", (0, 0.25), 0.85), 5),
                "threshold_radial_radius_arcsec": round(
                    threshold("radial", (0.25, 20), 0.85), 5
                ),
                "central_upper_beta_003": round(
                    self.predictions("symmetric", 0.03, 0.85)[centre], 3
                ),
                "central_upper_radial_25_arcsec": round(
                    self.predictions("radial", 2.5, 0.85)[centre], 3
                ),
            },
            "slit_offsets_arcsec": [r[0] for r in rows],
            "keck_lower_rms_endpoint_by_offset": np.round(observed, 2).tolist(),
            "isotropic_upper_by_offset_km_s": np.round(nominal, 3).tolist(),
            "uniform_anisotropy_envelope": {
                "central_threshold_abs_beta": round(
                    threshold("symmetric", (0, 0.25)), 5
                ),
                "central_upper_km_s": {
                    str(b): round(self.predictions("symmetric", b)[centre], 3)
                    for b in (0, 0.025, 0.05, 0.075, 0.10)
                },
            },
            "radial_anisotropy_envelope": {
                "central_threshold_radius_arcsec": round(
                    threshold("radial", (0.25, 20)), 5
                ),
                "central_upper_km_s": {
                    str(a): round(self.predictions("radial", a)[centre], 3)
                    for a in (1, 1.5, 2, 3, 5, 10)
                },
            },
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = Survey().run()
    path = HERE / "results.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        assert out == json.loads(path.read_text())
        print("Reproduced the two-instrument orbital anisotropy bounds.")
