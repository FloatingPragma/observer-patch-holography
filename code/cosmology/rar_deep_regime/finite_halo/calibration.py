"""Invert the shear reading for a scale without using a rotation calibration."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from . import measure
from .. import lensing_profile
from ..udg_transfer.measure import GM_SUN, KPC


def shear_scale(mean, b, weights, cap=1.0):
    core = 4 / math.pi * (weights @ b) / lensing_profile.FOUR_G
    factor = cap * (weights @ np.sqrt(b)) / lensing_profile.FOUR_G
    return (max(0.0, mean - core) / factor) ** 2


def run():
    raw = json.loads(lensing_profile.RECEIPT.read_text())
    records = json.loads((measure.HERE / "results.json").read_text())["surveys"]
    inversions = {}
    for name, survey in records.items():
        profile = raw["profiles"][name]
        lo = np.asarray(profile["lower_bin_edge_m_s2"])
        keep = lo >= 1e-13
        b = np.asarray(profile["g_bar_m_s2"])[keep]
        high = 2 * b - lo[keep]
        w = np.asarray(survey["positive_weights"])
        observed, sigma = survey["weighted_signal"], survey["statistical_se"]
        inversions[name] = dict(
            centre_a0=shear_scale(observed, b, w),
            central_inputs_two_error_range=[
                shear_scale(observed + z * sigma, b, w) for z in (-2.0, 2.0)
            ],
            generous_shear_floor=shear_scale(
                observed - 2 * sigma, 10**0.2 * high, w, measure.MAX_SHEAR_FACTOR
            ),
        )
    ring = json.loads(
        (measure.HERE.parent / "udg_transfer" / "rings.json").read_text()
    )["panels"]["figure8"]["rings"][-1]
    v = ring["circular_kms"] * 1000
    baryon = ring["baryonic_kms"] * 1000
    r = ring["radius_kpc"] * KPC
    vmax = (
        (v + 2000 * ring["plotted_error_kms"])
        * math.sin(math.radians(31))
        / math.sin(math.radians(27))
    )
    return dict(
        shear_inversions=inversions,
        dwarf=dict(
            spherical_untrimmed_upper=v**4 / (GM_SUN * 1.47e9),
            spherical_generous_upper=vmax**4 / (GM_SUN * 1.19e9),
            measured_disk_force_inversion=(v * v - baryon * baryon) ** 2
            / (r * baryon**2),
            radius_kpc=ring["radius_kpc"],
            untrimmed_dark_speed_from_lens_kms=[
                (name, (GM_SUN * 1.47e9 * obs["centre_a0"]) ** 0.25 / 1000)
                for name, obs in inversions.items()
            ],
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as out:
        json.dump(run(), out, indent=2, allow_nan=False)
        out.write("\n")
