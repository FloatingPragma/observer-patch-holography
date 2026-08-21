#!/usr/bin/env python3
"""Deep-regime comparison of the repair-charge condensate galaxy law with SPARC.

The dark-matter paper derives, on static spherical condensed-phase premises
with a cubic link energy, the deep-gradient law

    a_R = sqrt(a_b a_0),        v^4 = G M_b a_0,

with one dimensional constant ``a_0 = beta_b^3 / (4 pi G kappa_R)`` that the
source does not fix.  This module compares that law with the SPARC rotation
curve sample in the regime where its premises apply (``a_b << a_0``) and
reports

* the deep radial-acceleration exponent and acceleration constant fitted on
  points with ``g_bar < f a_0`` (self-consistent in ``a_0``);
* the baryonic Tully--Fisher exponent and the acceleration constant implied
  by its normalisation;
* the agreement between the two constants, which the condensate law forces
  to coincide;
* the residual of the additive all-gradient extension
  ``g = g_bar + sqrt(g_bar a_0)`` against the full sample and its Solar-System
  anomaly, which the paper's high-gradient branch is required to remove.

The acceleration constant is a fitted comparison value.  The exponents are the
parameter-free content of the law.  Nothing here is a source-derived OPH
output, and the receipt says so.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "oph.cosmology.rar_deep_regime_diagnostic.v1"
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

KPC_M = 3.0856775814913673e19
KM = 1.0e3
G_SI = 6.67430e-11
M_SUN = 1.98841e30
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7
GAS_HELIUM = 1.33  # already inside Vgas in table2; used for M_gas from M_HI
AU_M = 1.495978707e11
M_SUN_GM = G_SI * M_SUN


def _f(s: str) -> float:
    s = s.strip()
    return float(s) if s else math.nan


def read_table1(path: Path = DATA / "table1.dat") -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "hubble_type": int(line[12:14]),
                "dist_mpc": _f(line[15:21]),
                "inclination_deg": _f(line[30:34]),
                "e_inclination_deg": _f(line[35:39]),
                "L36_gsun": _f(line[40:47]),
                "MHI_gsun": _f(line[86:93]),
                "vflat_kms": _f(line[100:105]),
                "e_vflat_kms": _f(line[106:111]),
                "quality": int(line[112:115]),
            }
        )
    return rows


def read_table2(path: Path = DATA / "table2.dat") -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "rad_kpc": _f(line[19:25]),
                "vobs_kms": _f(line[26:32]),
                "e_vobs_kms": _f(line[33:38]),
                "vgas_kms": _f(line[39:45]),
                "vdisk_kms": _f(line[46:52]),
                "vbul_kms": _f(line[53:59]),
            }
        )
    return rows


def accelerations(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    max_quality: int = 2,
    min_inclination: float = 30.0,
    max_rel_error: float = 0.10,
) -> dict[str, np.ndarray]:
    """Point-wise g_obs and g_bar in m/s^2 under the standard SPARC cuts."""
    meta = {r["name"]: r for r in t1}
    g_obs, g_bar, names = [], [], []
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > max_quality or m["inclination_deg"] < min_inclination:
            continue
        if r["vobs_kms"] <= 0 or r["e_vobs_kms"] / r["vobs_kms"] > max_rel_error:
            continue
        R = r["rad_kpc"] * KPC_M
        vo = r["vobs_kms"] * KM
        vg, vd, vb = (r[k] * KM for k in ("vgas_kms", "vdisk_kms", "vbul_kms"))
        gb = (vg * abs(vg) + UPSILON_DISK * vd * abs(vd) + UPSILON_BULGE * vb * abs(vb)) / R
        if gb <= 0:
            continue
        g_obs.append(vo * vo / R)
        g_bar.append(gb)
        names.append(r["name"])
    return {
        "g_obs": np.array(g_obs),
        "g_bar": np.array(g_bar),
        "names": np.array(names),
    }


def deep_fit(g_obs: np.ndarray, g_bar: np.ndarray, fraction: float = 0.1,
             iterations: int = 50) -> dict[str, Any]:
    """Free-exponent and fixed-exponent fits on the self-consistent deep subset."""
    a0 = 1.2e-10
    for _ in range(iterations):
        mask = g_bar < fraction * a0
        x = np.log10(g_bar[mask])
        y = np.log10(g_obs[mask])
        # fixed exponent 1/2: log g_obs = 0.5 log g_bar + 0.5 log a0
        a0_new = 10 ** (2 * np.mean(y - 0.5 * x))
        if abs(a0_new / a0 - 1) < 1e-12:
            a0 = a0_new
            break
        a0 = a0_new
    mask = g_bar < fraction * a0
    x = np.log10(g_bar[mask])
    y = np.log10(g_obs[mask])
    slope, intercept = np.polyfit(x, y, 1)
    resid_fixed = y - (0.5 * x + 0.5 * math.log10(a0))
    resid_free = y - (slope * x + intercept)
    # bootstrap on galaxies would be better; point bootstrap gives a floor
    rng = np.random.default_rng(20260821)
    slopes = []
    for _ in range(400):
        idx = rng.integers(0, len(x), len(x))
        slopes.append(np.polyfit(x[idx], y[idx], 1)[0])
    return {
        "deep_fraction_of_a0": fraction,
        "n_points": int(mask.sum()),
        "a0_fixed_exponent_m_s2": float(a0),
        "rms_scatter_dex_fixed_exponent": float(np.sqrt(np.mean(resid_fixed**2))),
        "free_exponent": float(slope),
        "free_exponent_bootstrap_sd": float(np.std(slopes)),
        "rms_scatter_dex_free_exponent": float(np.sqrt(np.mean(resid_free**2))),
        "law_exponent": 0.5,
    }


def btfr_fit(t1: list[dict[str, Any]], max_quality: int = 2,
             min_inclination: float = 30.0) -> dict[str, Any]:
    Mb, Vf, names = [], [], []
    for r in t1:
        if r["quality"] > max_quality or r["inclination_deg"] < min_inclination:
            continue
        if not (r["vflat_kms"] > 0):
            continue
        mstar = UPSILON_DISK * r["L36_gsun"] * 1e9
        mgas = GAS_HELIUM * r["MHI_gsun"] * 1e9
        Mb.append((mstar + mgas) * M_SUN)
        Vf.append(r["vflat_kms"] * KM)
        names.append(r["name"])
    Mb = np.array(Mb)
    Vf = np.array(Vf)
    x = np.log10(Mb)
    y = np.log10(Vf)
    slope_v_on_m, _ = np.polyfit(x, y, 1)
    exponent = 1.0 / slope_v_on_m  # M_b ∝ V^exponent
    # fixed exponent 4: a0 = V^4 / (G M_b)
    a0_points = Vf**4 / (G_SI * Mb)
    a0 = float(10 ** np.mean(np.log10(a0_points)))
    resid = np.log10(a0_points) - math.log10(a0)
    rng = np.random.default_rng(20260821)
    exps = []
    for _ in range(400):
        idx = rng.integers(0, len(x), len(x))
        s, _ = np.polyfit(x[idx], y[idx], 1)
        exps.append(1.0 / s)
    return {
        "n_galaxies": int(len(Mb)),
        "free_exponent_M_of_V": float(exponent),
        "free_exponent_bootstrap_sd": float(np.std(exps)),
        "a0_fixed_exponent_m_s2": a0,
        "rms_scatter_dex_in_a0": float(np.sqrt(np.mean(resid**2))),
        "law_exponent": 4.0,
        "upsilon_disk": UPSILON_DISK,
    }


def additive_extension(g_obs: np.ndarray, g_bar: np.ndarray, a0: float) -> dict[str, Any]:
    """The cubic link law applied at all gradients: g = g_bar + sqrt(g_bar a0)."""
    model = g_bar + np.sqrt(g_bar * a0)
    resid = np.log10(g_obs) - np.log10(model)
    # McGaugh-Lelli-Schombert 2016 form for comparison, same a0 and own fit
    mls = g_bar / (1 - np.exp(-np.sqrt(g_bar / a0)))
    resid_mls = np.log10(g_obs) - np.log10(mls)
    grid = a0 * 10 ** np.linspace(-0.5, 0.5, 401)
    scat = [
        np.sqrt(np.mean((np.log10(g_obs) - np.log10(g_bar / (1 - np.exp(-np.sqrt(g_bar / a))))) ** 2))
        for a in grid
    ]
    best = int(np.argmin(scat))
    # Solar-System anomaly of the additive form at 1 au
    g_sun = M_SUN_GM / AU_M**2
    anomaly = math.sqrt(a0 / g_sun)
    return {
        "n_points": int(len(g_obs)),
        "rms_scatter_dex_additive_all_points": float(np.sqrt(np.mean(resid**2))),
        "mean_residual_dex_additive_all_points": float(np.mean(resid)),
        "rms_scatter_dex_mls16_same_a0": float(np.sqrt(np.mean(resid_mls**2))),
        "mls16_own_fit_a0_m_s2": float(grid[best]),
        "rms_scatter_dex_mls16_own_fit": float(scat[best]),
        "solar_system_relative_anomaly_at_1au": anomaly,
        "solar_system_bound_cassini_relative": 1.0e-8,
        "additive_form_passes_solar_system": bool(anomaly < 1.0e-8),
    }


def run() -> dict[str, Any]:
    t1 = read_table1()
    t2 = read_table2()
    acc = accelerations(t1, t2)
    deep = deep_fit(acc["g_obs"], acc["g_bar"])
    deep_sweep = [deep_fit(acc["g_obs"], acc["g_bar"], fraction=f) for f in (0.3, 0.1, 0.03)]
    btfr = btfr_fit(t1)
    a0_deep = deep["a0_fixed_exponent_m_s2"]
    a0_btfr = btfr["a0_fixed_exponent_m_s2"]
    additive = additive_extension(acc["g_obs"], acc["g_bar"], a0_deep)
    consistency = {
        "a0_deep_rar_m_s2": a0_deep,
        "a0_btfr_m_s2": a0_btfr,
        "log10_ratio_btfr_over_deep": math.log10(a0_btfr / a0_deep),
        "within_combined_scatter": bool(
            abs(math.log10(a0_btfr / a0_deep))
            < math.hypot(deep["rms_scatter_dex_fixed_exponent"], btfr["rms_scatter_dex_in_a0"])
        ),
    }
    data_sha = hashlib.sha256(
        (DATA / "table1.dat").read_bytes() + (DATA / "table2.dat").read_bytes()
    ).hexdigest()
    return {
        "schema": SCHEMA,
        "scope": "diagnostic_postdiction_one_fitted_constant",
        "physical_claim": False,
        "source_derived_output": False,
        "sample": {
            "catalogue": "SPARC, CDS J/AJ/152/157 (Lelli, McGaugh, Schombert 2016)",
            "data_sha256": data_sha,
            "cuts": "quality <= 2, inclination >= 30 deg, e_Vobs/Vobs <= 0.10, g_bar > 0",
            "n_points_after_cuts": int(len(acc["g_obs"])),
            "n_galaxies_after_cuts": int(len(set(acc["names"].tolist()))),
            "upsilon_disk": UPSILON_DISK,
            "upsilon_bulge": UPSILON_BULGE,
        },
        "law": {
            "deep_radial_acceleration": "a_R = sqrt(a_b a_0)",
            "baryonic_tully_fisher": "v^4 = G M_b a_0",
            "parameter_free_content": "exponents 1/2 and 4; a_0 is a fitted comparison value",
        },
        "deep_radial_acceleration": deep,
        "deep_radial_acceleration_fraction_sweep": deep_sweep,
        "baryonic_tully_fisher": btfr,
        "constant_consistency": consistency,
        "additive_all_gradient_extension": additive,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
