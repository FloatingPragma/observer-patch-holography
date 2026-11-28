"""Test the most favorable spherical enclosed trail mass at measured outer orbits."""

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import radians, sin, sqrt
from pathlib import Path

import numpy as np

from cosmology.cold_acoustic_polarization.trails import RECEIPT
from cosmology.rar_deep_regime import rar_deep_regime_diagnostic as sparc

HERE = Path(__file__).resolve().parent
# kpc (km/s)^2 per solar mass
G = sparc.M_SUN_GM / sparc.KPC_M / 1e6
BOXES = {
    "nominal": (0, 0.5, 0.5, 1.0, 1.0),
    "one_error": (1, 0.3, 0.8, 0.7, 1.5),
    "two_errors": (2, 0.3, 0.8, 0.7, 1.5),
    "extreme": (3, 0.3, 1.0, 0.7, 2.0),
}


def native_caps():
    source = json.loads(RECEIPT.read_text())
    out = {}
    for board in source["boards"]:
        n = board["live_seams"]
        e = Fraction(board["initial_energy"])
        out[str(n)] = {
            "owners": board["owners"],
            "cap": float(Fraction(17 * n, 2) / e),
            "cap_exact": str(Fraction(17 * n, 2) / e),
            "net": float(Fraction(board["net_work"]) / e),
        }
    return out


def signed_square(v):
    return v * abs(v)


def catalogue():
    extra = {}
    for line in (sparc.DATA / "table1.dat").read_text().splitlines():
        extra[line[:11].strip()] = {
            "e_L36_gsun": float(line[48:55]),
            "distance_method": int(line[28]),
            "reference_codes": line[116:130].strip(),
        }
    meta = {r["name"]: r | extra[r["name"]] for r in sparc.read_table1()}
    groups = defaultdict(list)
    for row in sparc.read_table2():
        groups[row["name"]].append(row)
    return meta, groups


def select(meta, groups):
    kept = {}
    tally = defaultdict(int)
    for name, m in meta.items():
        rows = groups[name]
        if m["quality"] not in (1, 2) or m["inclination_deg"] < 30:
            reason = "quality_or_inclination"
        elif m["L36_gsun"] <= 0 or m["MHI_gsun"] <= 0:
            reason = "missing_mass_proxy"
        elif any(r["vbul_kms"] != 0 for r in rows):
            reason = "resolved_bulge"
        else:
            good = [
                r
                for r in rows
                if r["rad_kpc"] > 0
                and r["vobs_kms"] > 0
                and 0 < r["e_vobs_kms"] <= 0.10 * r["vobs_kms"]
                and 0.5 * signed_square(r["vdisk_kms"]) + signed_square(r["vgas_kms"])
                > 0
            ]
            if good:
                kept[name] = sorted(good, key=lambda r: r["rad_kpc"])
                reason = "included"
            else:
                reason = "no_eligible_point"
        tally[reason] += 1
    return kept, dict(tally)


def bound(r, m, box, caps):
    k, ml_min, ml_max, gas_min, gas_max = box
    D = m["dist_mpc"]
    df = 1 + k * m["e_dist_mpc"] / D
    i = m["inclination_deg"]
    incl = min(90.0, i + k * m["e_inclination_deg"])
    v_floor = (
        max(0.0, r["vobs_kms"] - k * r["e_vobs_kms"])
        * sin(radians(i))
        / sin(radians(incl))
    )
    L = m["L36_gsun"] * 1e9
    gas = sparc.GAS_HELIUM * m["MHI_gsun"] * 1e9
    g_lum = k * m["e_L36_gsun"] / m["L36_gsun"]
    star_lo = ml_min * max(0.0, 1 - g_lum)
    star_hi = ml_max * (1 + g_lum)
    R = r["rad_kpc"] * df
    states = []
    for star, gas_factor in product((star_lo, star_hi), (gas_min, gas_max)):
        speed2 = df * (
            star * signed_square(r["vdisk_kms"])
            + gas_factor * signed_square(r["vgas_kms"])
        )
        mass = (star * L + gas_factor * gas) * df**2
        unit = G * mass / R
        states.append(
            dict(
                vbar2=speed2,
                unit=unit,
                mass=mass,
                star=star,
                gas_factor=gas_factor,
                q=(v_floor**2 - speed2) / unit,
            )
        )
    best = min(states, key=lambda x: x["q"])
    caps_result = {}
    for name, stock in caps.items():
        q = stock["cap"]
        state = max(states, key=lambda x: x["vbar2"] + q * x["unit"])
        max2 = state["vbar2"] + q * state["unit"]
        missing = max(0.0, v_floor**2 - max2) * R / G / (1 + q)
        caps_result[name] = {
            "maximum_speed_kms": sqrt(max(0, max2)),
            "deficit_kms": v_floor - sqrt(max(0, max2)),
            "additional_point_mass_seed_Msun": missing,
        }
    return {
        "q_required_min": best["q"],
        "v_floor_kms": v_floor,
        "radius_kpc": R,
        "distance_mpc": D * df,
        "inclination_deg": incl,
        "stars_multiplier_at_min": best["star"],
        "gas_multiplier_at_min": best["gas_factor"],
        "mass_at_min_Msun": best["mass"],
        "bar2_at_min_kms2": best["vbar2"],
        "caps": caps_result,
    }


def fold(name):
    # Match the pre-existing held_galaxies.fold_for split, including its NUL.
    h = hashlib.sha256(b"OPH-deep-held-v1\0" + name.encode()).digest()
    return int.from_bytes(h[:8], "big") % 5


def summarize(rows, caps):
    out = {"galaxies": len(rows)}
    for box in BOXES:
        vals = np.asarray([r["bounds"][box]["q_required_min"] for r in rows])
        out[box] = {
            "median_q_required": float(np.median(vals)),
            "exceed_cap": {k: int(sum(vals > v["cap"])) for k, v in caps.items()},
            "exceed_net": {k: int(sum(vals > v["net"])) for k, v in caps.items()},
            "largest_required": sorted(
                [(r["galaxy"], r["bounds"][box]["q_required_min"]) for r in rows],
                key=lambda x: -x[1],
            )[:5],
        }
    return out


def ngc3741_controls(meta, rows, caps):
    # Gentile et al. section 3.6 (2007): final pressure-support corrections.
    sensitivity = []
    for ring, correction in zip(rows[-2:], (2.7, 4.8), strict=True):
        shifted = ring | {"vobs_kms": ring["vobs_kms"] - correction}
        sensitivity.append(
            {
                "catalogue_radius_kpc": ring["rad_kpc"],
                "removed_speed_kms": correction,
                "q_required": {
                    key: bound(shifted, meta, box, caps)["q_required_min"]
                    for key, box in BOXES.items()
                },
            }
        )
    result = {}
    r = rows[-1]
    q = caps["32"]["cap"]
    L = meta["L36_gsun"] * 1e9
    gas = sparc.GAS_HELIUM * meta["MHI_gsun"] * 1e9
    for label in ("nominal", "two_errors"):
        fav = bound(r, meta, BOXES[label], caps)
        df = fav["distance_mpc"] / meta["dist_mpc"]
        star = fav["stars_multiplier_at_min"]
        grav = G * q / r["rad_kpc"]
        gas_required = (
            fav["v_floor_kms"] ** 2 / df
            - star * (signed_square(r["vdisk_kms"]) + grav * L)
        ) / (signed_square(r["vgas_kms"]) + grav * gas)
        result[label] = {
            "gas_multiplier_for_32_seams": gas_required,
            "extra_gas_over_box_Msun_at_given_distance": (
                gas_required - fav["gas_multiplier_at_min"]
            )
            * gas
            * df**2,
            "distance_mpc_for_32_seams": fav["distance_mpc"]
            * (fav["v_floor_kms"] / fav["caps"]["32"]["maximum_speed_kms"]) ** 2,
        }
    return {"removed_pressure_correction": sensitivity, "thresholds": result}


def run():
    meta, groups = catalogue()
    kept, tally = select(meta, groups)
    caps = native_caps()
    rows = []
    for name, good in kept.items():
        r = good[-1]
        m = meta[name]
        bar = 0.5 * signed_square(r["vdisk_kms"]) + signed_square(r["vgas_kms"])
        rows.append(
            {
                "galaxy": name,
                "historical_fold": fold(name),
                "eligible_rings": len(good),
                "catalogued_rings": len(groups[name]),
                "last_catalogue_radius_kpc": max(i["rad_kpc"] for i in groups[name]),
                "meta": m,
                "ring": r,
                "gbar_m_s2": bar * 1e6 / (r["rad_kpc"] * sparc.KPC_M),
                "bounds": {key: bound(r, m, box, caps) for key, box in BOXES.items()},
            }
        )
    main = summarize(rows, caps)
    worst = main["two_errors"]["largest_required"][0][0]
    return {
        "source_native_sha256": hashlib.sha256(RECEIPT.read_bytes()).hexdigest(),
        "source_sparc_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (sparc.DATA / "table1.dat", sparc.DATA / "table2.dat")
        },
        "G_kpc_kms2_per_Msun": G,
        "source_caps": caps,
        "sample_flow": tally,
        "all": main,
        "deep": summarize([r for r in rows if r["gbar_m_s2"] < 1.2e-11], caps),
        "legacy_folds": {
            str(i): summarize([r for r in rows if r["historical_fold"] == i], caps)
            for i in range(5)
        },
        "galaxies": rows,
        "worst_two_error_example": worst,
        "NGC3741_controls": ngc3741_controls(meta["NGC3741"], kept["NGC3741"], caps),
        "worst_example_rotation_curve": [
            {
                "r_kpc": r["rad_kpc"],
                "v_obs": r["vobs_kms"],
                "sigma": r["e_vobs_kms"],
                "q_required": {
                    key: bound(r, meta[worst], box, caps)["q_required_min"]
                    for key, box in BOXES.items()
                },
            }
            for r in kept[worst]
        ],
    }


if __name__ == "__main__":
    (HERE / "results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
