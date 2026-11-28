"""SPARC tests of an exterior potential inferred from mean first-passage clocks."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from math import radians, sin

import numpy as np

from . import occupied_work as orbit
from . import source_profile as profiles

HERE = orbit.HERE
galaxies = profiles.galaxies


def mass_test(tab):
    r = np.asarray(tab["radius_kpc"])
    assert np.all(np.diff(r) > 0)
    idx = np.flatnonzero(2 * r >= r[-1])
    if len(idx) < 3:
        return dict(rings=len(idx), included=False)
    lower, upper = orbit.orbits.budget.bounds(tab)
    lo, hi = r * lower / galaxies.G, r * upper / galaxies.G
    i, j = max(combinations(idx, 2), key=lambda ij: lo[ij[1]] - hi[ij[0]])
    # The least nonincreasing majorant of all lower endpoints is a witness.
    witness = np.maximum.accumulate(np.maximum(0, lo[idx])[::-1])[::-1]
    feasible = bool(np.all(witness <= hi[idx]))
    deficit = float(lo[j] - hi[i])
    negative = bool(np.min(hi[idx]) < 0)
    assert feasible == (not negative and deficit <= 0)
    return dict(
        included=True,
        rings=len(idx),
        nonincreasing_possible=feasible,
        required_positive=bool(np.max(lo[idx]) > 0),
        negative_upper=negative,
        growth_conflict=deficit > 0,
        witness_radii_kpc=[float(r[i]), float(r[j])],
        inner_upper_msun=float(hi[i]),
        outer_lower_msun=float(lo[j]),
        minimum_growth_msun=max(0, deficit),
    )


def shared_growth(meta, rows):
    r = np.array([p["rad_kpc"] for p in rows])
    idx = np.flatnonzero(2 * r >= r[-1])
    i0 = meta["inclination_deg"]
    ilo, ihi = (
        max(0, i0 - 2 * meta["e_inclination_deg"]),
        min(90, i0 + 2 * meta["e_inclination_deg"]),
    )
    hmin = (sin(radians(i0)) / sin(radians(ihi))) ** 2
    hmax = (sin(radians(i0)) / sin(radians(ilo))) ** 2 if ilo > 0 else np.inf
    width = 2 * meta["e_dist_mpc"] / meta["dist_mpc"]
    dmin, dmax = max(0, 1 - width), 1 + width
    best, pair = 0, None
    for i, j in combinations(idx, 2):
        ri, rj = r[i], r[j]
        a, b = rows[i], rows[j]
        A = (
            rj * max(0, b["vobs_kms"] - 2 * b["e_vobs_kms"]) ** 2
            - ri * (a["vobs_kms"] + 2 * a["e_vobs_kms"]) ** 2
        )
        stellar, gas = (
            rj * galaxies.signed_square(b[key]) - ri * galaxies.signed_square(a[key])
            for key in ("vdisk_kms", "vgas_kms")
        )
        cap = max(0.3 * stellar, 0.8 * stellar) + max(0.7 * gas, 1.5 * gas)
        floor = A * (hmin if A >= 0 else hmax) - max(dmin * cap, dmax * cap)
        if floor > best:
            best, pair = floor, [float(ri), float(rj)]
    return dict(
        positive_growth_for_every_setting=pair is not None,
        radii_kpc=pair,
        mass_growth_lower_msun=float(dmin * best / galaxies.G),
    )


def summarize(rows):
    good = [r for r in rows.values() if r["included"]]
    return dict(
        catalog=len(rows),
        eligible=len(good),
        rings=sum(r["rings"] for r in good),
        feasible=sum(r["nonincreasing_possible"] for r in good),
        growth_conflicts=sum(r["growth_conflict"] for r in good),
        negative_upper=sum(r["negative_upper"] for r in good),
        positive_required=sum(r["required_positive"] for r in good),
        shared_box_growth=sum(
            r.get("shared_box", {}).get("positive_growth_for_every_setting", False)
            for r in good
        ),
    )


def run():
    clock = orbit.LAB / "closed_passage.json"
    info = json.loads(clock.read_bytes())
    levels = list(map(Q, info["closed_minus_open_mean_partner_time_up_to_constant"]))
    radii, work = orbit.empirical()
    steps = [levels[i] - levels[i + 1] for i in range(2, 6)]
    assert min(steps) > 0
    intervals = [
        [x / float(f) for x in row] for row, f in zip(work, steps, strict=True)
    ]
    imaged, *_ = orbit.orbits.budget.load_inputs()
    images = {k: mass_test(v) for k, v in imaged.items()}
    meta, groups = galaxies.catalogue()
    selected, _ = galaxies.select(meta, groups)
    bulgeless = {}
    for name, rows in selected.items():
        tab = dict(
            radius_kpc=[p["rad_kpc"] for p in rows],
            observed_kms=[p["vobs_kms"] for p in rows],
            error_kms=[p["e_vobs_kms"] for p in rows],
            baryons_kms=np.sqrt(profiles.bare_speed2(rows)),
        )
        result = mass_test(tab)
        if result["included"]:
            result["shared_box"] = shared_growth(meta[name], rows)
        bulgeless[name] = result
    files = [
        clock,
        HERE / "PASSAGE_MASS.md",
        HERE / "source_profile.json",
        galaxies.sparc.DATA / "table1.dat",
        galaxies.sparc.DATA / "table2.dat",
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        ngc2841=dict(
            radii_kpc=radii,
            work_km2_s2=work,
            fitted_gain_intervals=intervals,
            common_gain=max(v[0] for v in intervals) <= min(v[1] for v in intervals),
            clock_outer_per_inner=float(
                (levels[4] - levels[6]) / (levels[2] - levels[4])
            ),
            observed_ratio_lower=sum(x[0] for x in work[2:])
            / sum(x[1] for x in work[:2]),
        ),
        imaged=dict(summary=summarize(images), galaxies=images),
        bulgeless=dict(summary=summarize(bulgeless), galaxies=bulgeless),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "passage_mass.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    print(result["imaged"]["summary"], result["bulgeless"]["summary"])
