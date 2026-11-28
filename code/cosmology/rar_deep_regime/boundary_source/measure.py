"""Assess preregistered common-gain source laws on disjoint galaxies."""

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from cosmology.rar_deep_regime.rar_deep_regime_diagnostic import DATA, read_table2

from .source import HERE, assignments, enroll, profile

G = 4.30091727003628e-6  # kpc (km/s)² / solar mass


def endpoint(row, source, h2=True):
    baryons = sum(
        factor * row[col] * abs(row[col])
        for col, factor in (("vgas_kms", 1), ("vdisk_kms", 0.5), ("vbul_kms", 0.7))
    )
    molecular = G * source["molecular_with_he_msun"] / row["rad_kpc"]
    vb2 = baryons + (molecular if h2 else 0)
    return {
        "observed_kms": row["vobs_kms"],
        "e_kms": row["e_vobs_kms"],
        "radius_kpc": row["rad_kpc"],
        "sparc_baryon_v2": baryons,
        "additional_molecular_v2": molecular,
        "baryon_v2": vb2,
        "excess_v2": row["vobs_kms"] ** 2 - vb2,
    }


def predictors(s):
    m = s["mass_msun"]
    xs = {
        "boundary": s["S_msun_kpc"],
        "sqrt_mass": np.sqrt(m),
        "exp_proxy": m / s["rd_kpc"],
    }
    return xs | {f"{k}_finite": val * s["F"] for k, val in xs.items()}


def scores(sources, targets, train, held):
    x = {n: predictors(src) for n, src in sources.items()}
    result = {}
    for law in x[train[0]]:
        gain = max(
            0.0, float(np.median([targets[n]["excess_v2"] / x[n][law] for n in train]))
        )
        rows = {}
        for n in held:
            pred = x[n][law] * gain
            t = targets[n]
            pred_v = np.sqrt(t["baryon_v2"] + pred)
            rows[n] = {
                "predicted_excess_v2": pred,
                "predicted_total_kms": pred_v,
                "total_residual_kms": pred_v - t["observed_kms"],
                "excess_residual_v2": pred - t["excess_v2"],
            }
        errors = np.array([r["total_residual_kms"] for r in rows.values()])
        result[law] = {
            "gain": gain,
            "median_abs_kms": np.median(abs(errors)),
            "rms_kms": np.sqrt(np.mean(errors**2)),
            "held": rows,
        }
    return result


def residuals(sources, targets):
    ids = sorted(k for k, t in targets.items() if t["excess_v2"] > 0)
    mass = np.log([sources[n]["mass_msun"] for n in ids])
    x = np.log([sources[n]["S_msun_kpc"] for n in ids]) - mass / 2
    y = np.log([targets[n]["excess_v2"] for n in ids]) - mass / 2
    slope = np.polyfit(x, y, 1)[0]
    return {
        "positive_galaxies": ids,
        "fixed_sqrt_mass_residual_slope": slope,
        "fixed_sqrt_mass_residual_spearman": spearmanr(x, y).statistic,
        "fixed_sqrt_mass_slope_leave_one_out": [
            np.polyfit(np.delete(x, i), np.delete(y, i), 1)[0] for i in range(len(ids))
        ],
        "free_mass_exponent_and_surface_exponent": np.linalg.lstsq(
            np.column_stack((np.ones(len(ids)), mass, x)), y + mass / 2, rcond=None
        )[0][1:].tolist(),
    }


def pair_control(kept, targets, pair):
    """Endpoint errors with fixed catalogue forces, rectangular nuisance sensitivities."""
    per = {}
    for n in pair:
        meta, lrow, rows, r, src = kept[n]
        t = targets[n]
        cases = {}
        for label, vary_geometry in [
            ("speed_only", False),
            ("speed_distance_inclination", True),
        ]:
            record = []
            for e_v, e_d, e_i in itertools.product(
                (-1, 1),
                (-1, 1) if vary_geometry else (0,),
                (-1, 1) if vary_geometry else (0,),
            ):
                d = meta["dist_mpc"] + e_d * meta["e_dist_mpc"]
                ratio = d / meta["dist_mpc"]
                inc = meta["inclination_deg"] + e_i * meta["e_inclination_deg"]
                speed = (
                    (t["observed_kms"] + e_v * t["e_kms"])
                    * np.sin(np.deg2rad(meta["inclination_deg"]))
                    / np.sin(np.deg2rad(inc))
                )
                y = speed**2 - ratio * t["baryon_v2"]
                alt = profile(meta, lrow, rows, r * ratio, distance=d)
                record.append((y, y / alt["S_msun_kpc"]))
            cases[label] = {
                k: [float(np.min(a)), float(np.max(a))]
                for k, a in zip(
                    ("excess_v2", "required_boundary_gain"),
                    np.transpose(record),
                    strict=True,
                )
            }
        per[n] = cases
    a, b = pair
    ratio = {}
    for label in per[a]:
        ratio[label] = {}
        for quantity in per[a][label]:
            lo, hi = per[a][label][quantity]
            blo, bhi = per[b][label][quantity]
            ratio[label][quantity] = (
                [lo / bhi, hi / blo] if blo > 0 and lo >= 0 else None
            )
    return {"endpoints": per, "ratio_bounds": ratio}


def summary_for(sources, targets, train, held, pair):
    a, b = pair
    score = scores(sources, targets, train, held)
    xs = {k: predictors(sources[k]) for k in pair}
    return {
        "held_median_abs_kms": {k: r["median_abs_kms"] for k, r in score.items()},
        "pair_source_ratios": {k: v / xs[b][k] for k, v in xs[a].items()},
        "pair_observed_excess_ratio": targets[a]["excess_v2"] / targets[b]["excess_v2"],
    }


def run():
    kept, dropped = enroll()
    sources = {k: v[-1] for k, v in kept.items()}
    train, held, pair = assignments(sources)
    orbits = {
        r["name"]: r
        for r in read_table2()
        if r["name"] in kept and r["rad_kpc"] == kept[r["name"]][3]
    }
    if set(kept) != set(orbits):
        raise ValueError("Missing selected orbital endpoint")
    targets = {k: endpoint(orbits[k], s) for k, s in sources.items()}
    variants = {
        "coarse_1p6_kpc": {"width": 1.6},
        "half_CO": {"molecular": 0.5},
        "double_CO": {"molecular": 2.0},
        "no_CO_source_or_target": {"molecular": 0.0},
        "low_stellar_source": {"star": 0.8},
        "high_stellar_source": {"star": 1.2},
        "published_source_common_catalogue_inclination": {},
    }
    controls = {}
    for key, params in variants.items():
        alt = {}
        for n, (meta, lrow, rows, r, _) in kept.items():
            if key == "published_source_common_catalogue_inclination":
                factor = np.cos(np.deg2rad(meta["inclination_deg"])) / np.cos(
                    np.deg2rad(float(lrow["i"]))
                )
                params = dict(star=factor, atomic=factor, molecular=factor)
            alt[n] = profile(meta, lrow, rows, r, **params)
        tar = {n: endpoint(orbits[n], s) for n, s in alt.items()}
        controls[key] = summary_for(alt, tar, train, held, pair)
    controls["published_target_without_CO"] = summary_for(
        sources,
        {n: endpoint(orbits[n], s, h2=False) for n, s in sources.items()},
        train,
        held,
        pair,
    )
    files = [
        HERE / "leroy_table7.dat",
        HERE / "leroy_sample.dat",
        DATA / "table1.dat",
        DATA / "table2.dat",
        HERE / "PROTOCOL.md",
    ]
    return {
        "units": "radii kpc, masses solar, speeds km/s, v2 (km/s)^2, S solar masses/kpc",
        "sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files
        },
        "excluded": dropped,
        "train": train,
        "held": held,
        "pair_high_low": pair,
        "profiles": sources,
        "targets": targets,
        "pair_and_primary_summary": summary_for(sources, targets, train, held, pair),
        "pair_geometry_controls": pair_control(kept, targets, pair),
        "held_predictions": scores(sources, targets, train, held),
        "surface_residual": residuals(sources, targets),
        "controls": controls,
    }


if __name__ == "__main__":
    Path(HERE / "results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
