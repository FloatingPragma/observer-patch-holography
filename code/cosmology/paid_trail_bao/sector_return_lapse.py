"""Forecast differential home-clock rates directly from the mapped galaxy fillings."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import sector_mass
from . import sector_return_redshift as sky

HERE = sky.HERE


def clock_rates(p, table):
    ends, ind, src = table
    j = np.searchsorted(ends, p)
    use = 6 * p**2 * (1 - p) ** 2 * 448 / 2110
    norm = 1 + use[:, None] * (ind[j] - 1)
    shift = (p**2 * use)[:, None] * src[j]
    return norm, shift


def ring_rates(r, rpix, p, table):
    mid = (r[1:] + r[:-1]) / 2
    edges = np.r_[max(0, 2 * r[0] - mid[0]), mid, 2 * r[-1] - mid[-1]]
    ids = np.searchsorted(edges, rpix, side="right") - 1
    norm, shift = clock_rates(p, table)
    rows = []
    for i in range(len(r)):
        take = ids == i
        if not np.any(take):
            rows.append(None)
            continue
        ordinary = norm[take].mean(0)
        cross = shift[take].mean(0)
        ratios = np.divide.outer(cross, ordinary)
        rows.append(
            {
                "mapped_pixels": int(take.sum()),
                "ordinary_log_frequency": np.log(ordinary).tolist(),
                "extra_log_frequency": np.log1p([ratios.min(), ratios.max()]).tolist(),
            }
        )
    return rows


def gain_interval(model, target):
    # Necessary gain >=0 for the two closed intervals to intersect.
    l, h = model
    a, b = target
    lower, upper = 0.0, np.inf
    for coeff, rhs in ((h, a), (-l, -b)):
        if coeff > 0:
            lower = max(lower, rhs / coeff)
        elif coeff < 0:
            upper = min(upper, rhs / coeff)
        elif rhs > 0:
            return 1.0, 0.0
    return lower, upper


def compare(t, rows):
    r, v, e, b = (
        np.array(t[k])
        for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
    )
    phi = np.r_[
        0, np.cumsum((b[:-1] ** 2 + b[1:] ** 2) / 2 * np.log(r[1:] / r[:-1]) / sky.C**2)
    ]
    # Necessary separate bounds valid on each annulus. Even if all mass
    # increases sit between tabulated radii, the endpoint limits bound
    # the integral. Upper inconsistencies are reported separately.
    dr = 1 / r[:-1] - 1 / r[1:]
    floor = np.maximum.accumulate(
        r * np.maximum(0, np.maximum(0, v - 2 * e) ** 2 - b**2)
    )
    ceiling = np.minimum.accumulate((r * ((v + 2 * e) ** 2 - b**2))[::-1])[::-1]
    low = np.r_[0, np.cumsum(floor[:-1] * dr / sky.C**2)]
    # Drop an inconsistent upper endpoint from the upper-side test only.
    # The lower-side exclusion does not require any upper endpoint.
    high = np.r_[0, np.cumsum(np.maximum(0, ceiling[1:]) * dr / sky.C**2)]
    counts = {
        "covered_pairs": 0,
        "below_extra_floor": 0,
        "above_consistent_extra_ceiling": 0,
        "ordinary_unit_misses": 0,
        "positive_floor_with_wrong_sign": 0,
        "spherical_boxes_inconsistent": 0,
    }
    gains = {"ordinary": [0.0, np.inf], "extra": [0.0, np.inf]}
    worst = None
    for i, first in enumerate(rows):
        if first is None:
            continue
        for j in range(i + 1, len(r)):
            if rows[j] is None:
                continue
            o, q = (
                np.array(rows[j][k]) - first[k][::-1]
                for k in ("ordinary_log_frequency", "extra_log_frequency")
            )
            a, h = low[j] - low[i], high[j] - high[i]
            consistent = bool(np.all(floor[i : j + 1] <= ceiling[i : j + 1] + 1e-9))
            counts["covered_pairs"] += 1
            counts["spherical_boxes_inconsistent"] += int(not consistent)
            counts["below_extra_floor"] += int(q[1] < a)
            counts["above_consistent_extra_ceiling"] += int(consistent and q[0] > h)
            counts["positive_floor_with_wrong_sign"] += int(a > 0 and q[1] <= 0)
            counts["ordinary_unit_misses"] += int(not o[0] <= phi[j] - phi[i] <= o[1])
            for tag, pred, target in (
                ("ordinary", o, [phi[j] - phi[i]] * 2),
                ("extra", q, [a, h if consistent else np.inf]),
            ):
                l, u = gain_interval(pred, target)
                gains[tag] = [max(gains[tag][0], l), min(gains[tag][1], u)]
            gap = a - q[1]
            if worst is None or gap > worst["floor_minus_prediction"]:
                worst = {
                    "radii_kpc": [float(r[i]), float(r[j])],
                    "required_lower": a,
                    "predicted_extra_log_ratio": q.tolist(),
                    "floor_minus_prediction": gap,
                }
    feasible = {k: bool(lo <= hi and hi > 0) for k, (lo, hi) in gains.items()}
    return {
        "counts": counts,
        "positive_scale_possible": feasible,
        "one_positive_scale_for_both": bool(
            max(v[0] for v in gains.values()) <= min(v[1] for v in gains.values())
            and all(feasible.values())
        ),
        "largest_lower_bound_gap": worst,
    }


def run():
    lab = sector_mass.inventory.LAB / "sector_return_cadence.json"
    source = HERE / "source_profile.json"
    data = json.loads(source.read_text())["resolved_maps"]["inputs"]
    rates = json.loads(lab.read_text())["bins"]
    table = [
        np.array([r[k] for r in rates])
        for k in (
            "rate",
            "independent_frequency_over_blank",
            "distinct_minus_independent_per_squared_rate",
        )
    ]
    table[0] = table[0][:, 1]
    all_cases = {}
    profiles = {}
    with np.load(sector_mass.pooled.sky.MAPS) as charts:
        for tracer in ("mapped_baryons", "atomic"):
            for capacity in sector_mass.inventory.rp.CAPACITY:
                examples = {}
                for name, t in data.items():
                    p = charts[f"{name}_{tracer}"] / capacity
                    r = np.array(t["radius_kpc"])
                    rows = ring_rates(
                        r, charts[f"{name}_{tracer}_radius_kpc"], p, table
                    )
                    examples[name] = compare(t, rows)
                    if tracer == "mapped_baryons" and capacity == 1:
                        profiles[name] = {"radius_kpc": r.tolist(), "clocks": rows}
                all_cases[f"{tracer}/{capacity}"] = {
                    "galaxies": examples,
                    "with_lower_bound_miss": sum(
                        x["counts"]["below_extra_floor"] > 0 for x in examples.values()
                    ),
                    "with_wrong_sign": sum(
                        x["counts"]["positive_floor_with_wrong_sign"] > 0
                        for x in examples.values()
                    ),
                    "joint_positive_gain": sum(
                        x["one_positive_scale_for_both"] for x in examples.values()
                    ),
                }
    return {
        "input_sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (lab, source, sector_mass.pooled.sky.MAPS)
        },
        "profiles": profiles,
        "comparisons": all_cases,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    res = run()
    dst = HERE / "sector_return_lapse.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(res, indent=2, allow_nan=False) + "\n")
    else:
        sector_mass.inventory.rp.galaxy.check(json.loads(dst.read_text()), res)
    print(
        [
            (
                k,
                v["with_lower_bound_miss"],
                v["with_wrong_sign"],
                v["joint_positive_gain"],
            )
            for k, v in res["comparisons"].items()
        ]
    )
