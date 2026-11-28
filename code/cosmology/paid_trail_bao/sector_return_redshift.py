"""Infer minimum static clock differences from rotation data without fitting a halo profile."""

import argparse
import json
from hashlib import sha256
from itertools import product
from statistics import NormalDist

import numpy as np
from cosmology.rar_deep_regime.individual_trail_budget import measure as catalogue
from scipy.constants import c

from . import sector_home_mass

HERE = sector_home_mass.HERE
C = c / 1000


def shell_bound(radius, low, high):
    radius = np.asarray(radius)
    assert len(radius) and np.all(radius > 0) and np.all(np.diff(radius) > 0)
    # r*v_extra^2 is G*M_enclosed. Its running lower envelope is the
    # smallest nondecreasing enclosed source compatible with the lower data.
    mass = np.maximum.accumulate(radius * np.maximum(low, 0))
    upper = radius * high
    drops = (mass - upper) > 1e-7
    pieces = mass[:-1] * (1 / radius[:-1] - 1 / radius[1:])
    return {
        "radii_kpc": radius.tolist(),
        "minimum_log_rate_difference": float(pieces.sum() / C**2),
        "largest_single_orbit_excess_kms2": float(max(max(low), 0)),
        "minimum_cumulative_log_rate_difference": np.r_[
            0, np.cumsum(pieces) / C**2
        ].tolist(),
        "failed_upper_constraints": int(drops.sum()),
    }


def speed_residuals(meta, rows, k, wide):
    v = np.array([r["vobs_kms"] for r in rows])
    err = np.array([r["e_vobs_kms"] for r in rows])
    gas, star = (
        np.array([catalogue.signed_square(r[key]) for r in rows])
        for key in ("vgas_kms", "vdisk_kms")
    )
    if wide:
        distance = meta["dist_mpc"]
        delta = 2 * meta["e_dist_mpc"] / distance
        ds = [max(0, 1 - delta), 1 + delta]
        inc = meta["inclination_deg"]
        ie = 2 * meta["e_inclination_deg"]
        floor = np.sin(np.deg2rad(inc)) / np.sin(np.deg2rad(min(90, inc + ie)))
        # If the lower inclination reaches zero there is no finite upper
        # speed constraint from deprojection within the sensitivity range.
        hi_inc = max(0, inc - ie)
        ceiling = (
            np.inf
            if hi_inc == 0
            else np.sin(np.deg2rad(inc)) / np.sin(np.deg2rad(hi_inc))
        )
        lum = 2 * meta["e_L36_gsun"] / meta["L36_gsun"]
        ml = [0.3 * max(0, 1 - lum), 0.8 * (1 + lum)]
        gas_scale = [0.7, 1.5]
    else:
        ds, ml, gas_scale = [1], [0.5], [1]
        floor = ceiling = 1
    baryons = np.array(
        [d * (a * star + g * gas) for d, a, g in product(ds, ml, gas_scale)]
    )
    lo = (floor * np.maximum(0, v - k * err)) ** 2 - baryons.max(0)
    hi = (ceiling * (v + k * err)) ** 2 - baryons.min(0)
    return lo, hi


def summarize(rows):
    lens = [r for r in rows.values() if len(r["radii_kpc"]) >= 2]
    pos = [
        r["minimum_log_rate_difference"]
        for r in lens
        if r["minimum_log_rate_difference"] > 0
    ]
    ok = [r for r in lens if r["failed_upper_constraints"] == 0]
    return {
        "galaxies": len(rows),
        "two_distinct_radii": len(lens),
        "positive_required": len(pos),
        "consistent_spherical_and_positive": sum(
            r["minimum_log_rate_difference"] > 0 for r in ok
        ),
        "inconsistent_spherical": len(lens) - len(ok),
        "range_positive_inferred_log_rate": [min(pos), max(pos)] if pos else [],
        "median_all_two_radius_bounds": float(
            np.median([r["minimum_log_rate_difference"] for r in lens])
        ),
    }


def run():
    meta, all_rows = catalogue.catalogue()
    selected, cuts = catalogue.select(meta, all_rows)
    # Use the full catalogue count even though the fixed cuts retain fewer
    # points. This multiplicity allowance addresses the formal speed errors.
    total = sum(map(len, all_rows.values()))
    simultaneous = NormalDist().inv_cdf(1 - 0.05 / (2 * total))
    sets = {}
    for tag, k, wide in (
        ("fixed_two_errors", 2, False),
        ("fixed_catalogue_adjusted", simultaneous, False),
        ("systematic_box", simultaneous, True),
    ):
        result = {}
        for name, rows in selected.items():
            r = np.array([x["rad_kpc"] for x in rows])
            lower, upper = speed_residuals(meta[name], rows, k, wide)
            result[name] = shell_bound(r, lower, upper)
        sets[tag] = {"summary": summarize(result), "galaxies": result}
    source = HERE / "source_profile.json"
    imaged = json.loads(source.read_text())["resolved_maps"]["inputs"]
    result = {}
    for name, t in imaged.items():
        v, e, b = (np.array(t[k]) for k in ("observed_kms", "error_kms", "baryons_kms"))
        result[name] = shell_bound(
            t["radius_kpc"],
            np.maximum(0, v - 2 * e) ** 2 - b**2,
            (v + 2 * e) ** 2 - b**2,
        )
    return {
        "sha256": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (
                catalogue.sparc.DATA / "table1.dat",
                catalogue.sparc.DATA / "table2.dat",
                source,
            )
        },
        "catalogue_points": total,
        "selection": cuts,
        "statistical_error_multiplier": simultaneous,
        "bulgeless": sets,
        "ten_imaged": {"summary": summarize(result), "galaxies": result},
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    dst = HERE / "sector_return_redshift.json"
    result = run()
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        sector_home_mass.inventory.rp.galaxy.check(json.loads(dst.read_text()), result)
    print(
        [v["summary"] for v in result["bulgeless"].values()],
        result["ten_imaged"]["summary"],
    )
