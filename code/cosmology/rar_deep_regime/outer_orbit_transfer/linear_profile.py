"""Find a constructive shared linear mass profile at catalogue masses and inclinations."""

import json
import math
from pathlib import Path

import numpy as np

from .. import lensing_profile
from ..udg_transfer.measure import GM_SUN


def orbit_interval(row, k, vary_speed=True):
    distance_error = k * row["distance_error_mpc"] / row["distance_mpc"]
    dv = k * row["speed_error_kms"] if vary_speed else 0.0
    stellar = 0.8 * (row["vdisk_unit_kms"] ** 2 + row["vbul_unit_kms"] ** 2)
    mass = 0.8 * row["stellar_unit_msun"] + row["helium_hi_msun"]
    factor = 1e6 / math.sqrt(GM_SUN * mass)
    # The fitted coordinate here is sqrt(a0); a negative upper endpoint is infeasible.
    low = max(
        0.0,
        factor
        * (max(0.0, row["speed_kms"] - dv) ** 2 / (1 + distance_error) - stellar),
    )
    high = (
        math.inf
        if distance_error >= 1
        else factor * ((row["speed_kms"] + dv) ** 2 / (1 - distance_error) - stellar)
    )
    return low, high


def lens_intervals(profiles, surveys):
    intervals = {}
    for name, survey in surveys.items():
        profile = profiles[name]
        mask = np.array(profile["lower_bin_edge_m_s2"]) >= 1e-13
        g = np.array(profile["g_bar_m_s2"])[mask]
        w = np.array(survey["positive_weights"])
        offset = 4 / math.pi * (w @ g)
        slope = w @ np.sqrt(g)
        mean, se = survey["weighted_signal"], survey["statistical_se"]

        def amplitude(value):
            return max(0.0, (value * lensing_profile.FOUR_G - offset) / slope)

        intervals[name] = (amplitude(mean - 2 * se), amplitude(mean + 2 * se))
    return intervals


def compatibility(rows, lens, k, vary_speed):
    intervals = {r["name"]: orbit_interval(r, k, vary_speed) for r in rows}
    intervals.update(lens)
    lower, left = max((lo, name) for name, (lo, _) in intervals.items())
    upper, right = min((hi, name) for name, (_, hi) in intervals.items())
    return lower, upper, left, right


def solve(rows, lens, vary_speed):
    lo, hi = 0.0, 1.0
    while True:
        a, b, _, _ = compatibility(rows, lens, hi, vary_speed)
        if a <= b:
            break
        hi *= 2
        if hi > 128:
            return dict(feasible=False)
    a, b, _, _ = compatibility(rows, lens, 0.0, vary_speed)
    if a <= b:
        hi = 0.0
    for _ in range(60):
        midpoint = (lo + hi) / 2
        a, b, _, _ = compatibility(rows, lens, midpoint, vary_speed)
        if a <= b:
            hi = midpoint
        else:
            lo = midpoint
    a, b, left, right = compatibility(rows, lens, hi, vary_speed)
    amplitude = (a + b) / 2
    galaxies = {}
    for r in rows:
        m = 0.8 * r["stellar_unit_msun"] + r["helium_hi_msun"]
        predicted_sq = (
            0.8 * (r["vdisk_unit_kms"] ** 2 + r["vbul_unit_kms"] ** 2)
            + math.sqrt(GM_SUN * m) * amplitude / 1e6
        )
        best_distance = r["distance_mpc"] * r["speed_kms"] ** 2 / predicted_sq
        delta = hi * r["distance_error_mpc"]
        distance = min(
            max(best_distance, r["distance_mpc"] - delta), r["distance_mpc"] + delta
        )
        galaxies[r["name"]] = dict(
            distance_mpc=distance,
            predicted_speed_kms=math.sqrt(predicted_sq * distance / r["distance_mpc"]),
        )
    return dict(
        feasible=True,
        error_multiplier=hi,
        a0_m_s2=amplitude**2,
        lower=left,
        upper=right,
        witness=galaxies,
    )


def run():
    here = Path(__file__).resolve().parent
    profiles = json.loads(lensing_profile.RECEIPT.read_text())["profiles"]
    surveys = json.loads((here.parent / "finite_halo/results.json").read_text())[
        "surveys"
    ]
    lens = lens_intervals(profiles, surveys)
    selected = json.loads((here / "results.json").read_text())["panels"]
    data = json.loads((here / "inputs.json").read_text())["galaxies"]
    panels = {}
    for name, panel in selected.items():
        rows = [r for r in data if r["name"] in panel["names"]]
        if rows:
            panels[name] = dict(
                distance_only=solve(rows, lens, False),
                distance_and_speed=solve(rows, lens, True),
            )
    common_lower = max(low for low, _ in lens.values())
    common_upper = min(high for _, high in lens.values())
    forecasts = {}
    if common_lower <= common_upper:
        for r in data:
            stellar_sq = 0.8 * (r["vdisk_unit_kms"] ** 2 + r["vbul_unit_kms"] ** 2)
            gm = GM_SUN * (0.8 * r["stellar_unit_msun"] + r["helium_hi_msun"])
            forecasts[r["name"]] = [
                math.sqrt(stellar_sq + math.sqrt(gm) * amplitude / 1e6)
                for amplitude in (common_lower, common_upper)
            ]
    return dict(
        lens_only_forecast_at_central_distance_kms=forecasts,
        lens_two_statistical_error_intervals_a0_m_s2={
            n: [a * a, b * b] for n, (a, b) in lens.items()
        },
        panels=panels,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
