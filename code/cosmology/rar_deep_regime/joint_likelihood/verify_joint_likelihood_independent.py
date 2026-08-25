#!/usr/bin/env python3
"""Independent verifier for the joint penalized-profile objective receipt.

This script imports nothing from ``joint_rar_likelihood.py``.  It re-parses
the committed SPARC snapshot with its own fixed-width parser built from the
CDS ReadMe byte ranges, replays the declared computation (per-point Gaussian
velocity data likelihood, per-galaxy nuisance-grid penalties and profiling,
equal-correlation block covariance scan, grid-anchored delta-objective
contours, matched channel-B
estimator, paired bootstrap with the declared seed), serializes the result
as canonical JSON (sorted keys, two-space indent, trailing newline), and
byte-compares it against ``runtime/joint_likelihood_receipt.json``.

Byte equality of floating-point output requires the same arithmetic
evaluation order as the producer contract; the declared formulas below fix
that order.  The verification certifies that the committed receipt is
exactly the declared computation on the committed snapshot bytes, replayed
by code written separately from the producer.

What is not proved here.  Replay equality does not authenticate provenance
or custody of the snapshot, does not calibrate the covariance family, does
not make the postdiction a preregistered comparison, and does not turn any
number in the receipt into evidence for or against OPH.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RECEIPT_PATH = HERE / "runtime" / "joint_likelihood_receipt.json"

# Declared constants, restated independently of the producer.
KPC_M = 3.0856775814913673e19
KM = 1.0e3
G_SI = 6.67430e-11
M_SUN = 1.98841e30
UPSILON_DISK_CENTER = 0.5
UPSILON_BULGE = 0.7
GAS_HELIUM = 1.33
REFERENCE_A0 = 1.2e-10
ML_PRIOR_WIDTH_DEX = 0.1
DIST_FACTOR_MIN = 0.05
INCLINATION_CLIP = (10.0, 89.9)
THRESHOLD_68 = 1.0
THRESHOLD_95 = 3.84
SEED = 20260825
SCHEMA = "oph.cosmology.joint_rar_btfr_penalized_profile_objective.v2"

LOG10_A0_MIN = -11.0
LOG10_A0_MAX = -9.5
N_A0 = 301
N_UPSILON = 9
N_DISTANCE = 7
N_INCLINATION = 7
HALFWIDTH_SIGMA = 2.5
RHO_GRID = (0.0, 0.2, 0.4, 0.6)
FRACTIONS = (0.3, 0.1)
REPLICATES = 1000
MAX_QUALITY = 2
MIN_INCLINATION = 30.0
MAX_REL_VOBS_ERROR = 0.10


def _num(line: str, lo: int, hi: int) -> float:
    """Numeric field at CDS one-based byte range [lo, hi], blank as NaN."""
    text = line[lo - 1 : hi].strip()
    return float(text) if text else math.nan


def parse_table1(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "dist_mpc": _num(line, 16, 21),
                "e_dist_mpc": _num(line, 23, 27),
                "inclination_deg": _num(line, 31, 34),
                "e_inclination_deg": _num(line, 36, 39),
                "L36_gsun": _num(line, 41, 47),
                "MHI_gsun": _num(line, 87, 93),
                "quality": int(line[112:115]),
            }
        )
    return rows


def parse_table2(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "rad_kpc": _num(line, 20, 25),
                "vobs_kms": _num(line, 27, 32),
                "e_vobs_kms": _num(line, 34, 38),
                "vgas_kms": _num(line, 40, 45),
                "vdisk_kms": _num(line, 47, 52),
                "vbul_kms": _num(line, 54, 59),
            }
        )
    return rows


def assemble_galaxies(
    t1: list[dict[str, Any]], t2: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], int]:
    meta = {r["name"]: r for r in t1}
    collected: dict[str, dict[str, list[float]]] = {}
    n_bad_ev = 0
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > MAX_QUALITY:
            continue
        if m["inclination_deg"] < MIN_INCLINATION:
            continue
        if r["vobs_kms"] <= 0:
            continue
        if r["e_vobs_kms"] <= 0:
            n_bad_ev += 1
            continue
        if r["e_vobs_kms"] / r["vobs_kms"] > MAX_REL_VOBS_ERROR:
            continue
        vg = r["vgas_kms"] * KM
        vd = r["vdisk_kms"] * KM
        vb = r["vbul_kms"] * KM
        gas_term = vg * abs(vg)
        disk_term = vd * abs(vd)
        bul_term = UPSILON_BULGE * vb * abs(vb)
        r_m = r["rad_kpc"] * KPC_M
        vbar2_cat = gas_term + UPSILON_DISK_CENTER * disk_term + bul_term
        g_bar_cat = vbar2_cat / r_m
        if g_bar_cat <= 0:
            continue
        rec = collected.setdefault(
            r["name"],
            {
                "rad_kpc": [],
                "vobs_kms": [],
                "e_vobs_kms": [],
                "gas_term": [],
                "disk_term": [],
                "bul_term": [],
                "g_bar_cat": [],
            },
        )
        rec["rad_kpc"].append(r["rad_kpc"])
        rec["vobs_kms"].append(r["vobs_kms"])
        rec["e_vobs_kms"].append(r["e_vobs_kms"])
        rec["gas_term"].append(gas_term)
        rec["disk_term"].append(disk_term)
        rec["bul_term"].append(bul_term)
        rec["g_bar_cat"].append(g_bar_cat)
    galaxies: list[dict[str, Any]] = []
    for r in t1:
        rec = collected.get(r["name"])
        if rec is None:
            continue
        gal: dict[str, Any] = {"name": r["name"], "meta": r}
        for key, values in rec.items():
            gal[key] = np.array(values)
        gal["idx_out"] = int(np.argmax(gal["rad_kpc"]))
        galaxies.append(gal)
    return galaxies, n_bad_ev


def make_grids(meta: dict[str, Any]) -> dict[str, Any]:
    offsets_u = np.linspace(-HALFWIDTH_SIGMA, HALFWIDTH_SIGMA, N_UPSILON)
    log10_u = math.log10(UPSILON_DISK_CENTER) + ML_PRIOR_WIDTH_DEX * offsets_u
    upsilon = 10.0 ** log10_u
    pen_u = ((log10_u - math.log10(UPSILON_DISK_CENTER)) / ML_PRIOR_WIDTH_DEX) ** 2

    sd_d = (
        meta["e_dist_mpc"] / meta["dist_mpc"]
        if meta["dist_mpc"] > 0
        and not math.isnan(meta["e_dist_mpc"])
        and meta["e_dist_mpc"] > 0
        else None
    )
    if sd_d is None:
        d_grid = np.array([1.0])
        pen_d = np.array([0.0])
        distance_fixed = True
    else:
        offsets_d = np.linspace(-HALFWIDTH_SIGMA, HALFWIDTH_SIGMA, N_DISTANCE)
        d_grid = np.maximum(1.0 + sd_d * offsets_d, DIST_FACTOR_MIN)
        pen_d = ((d_grid - 1.0) / sd_d) ** 2
        distance_fixed = False

    i_cat = meta["inclination_deg"]
    e_i = meta["e_inclination_deg"]
    if math.isnan(e_i) or e_i <= 0:
        i_grid = np.array([i_cat])
        pen_i = np.array([0.0])
        inclination_fixed = True
    else:
        offsets_i = np.linspace(-HALFWIDTH_SIGMA, HALFWIDTH_SIGMA, N_INCLINATION)
        i_grid = np.clip(
            i_cat + e_i * offsets_i, INCLINATION_CLIP[0], INCLINATION_CLIP[1]
        )
        pen_i = ((i_grid - i_cat) / e_i) ** 2
        inclination_fixed = False

    penalty = pen_u[:, None, None] + pen_d[None, :, None] + pen_i[None, None, :]
    return {
        "upsilon": upsilon,
        "d_grid": d_grid,
        "i_grid": i_grid,
        "penalty": penalty,
        "distance_fixed": distance_fixed,
        "inclination_fixed": inclination_fixed,
    }


def block_stats(
    gal: dict[str, Any],
    mask: np.ndarray,
    log10_a0_grid: np.ndarray,
    grids: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    r_m = gal["rad_kpc"][mask] * KPC_M
    gas = gal["gas_term"][mask]
    disk = gal["disk_term"][mask]
    bul = gal["bul_term"][mask]
    vobs = gal["vobs_kms"][mask] * KM
    ev = gal["e_vobs_kms"][mask] * KM
    npts = int(mask.sum())
    upsilon = grids["upsilon"]
    d_grid = grids["d_grid"]
    i_grid = grids["i_grid"]
    a0 = 10.0 ** log10_a0_grid

    vbar2 = gas[None, :] + upsilon[:, None] * disk[None, :] + bul[None, :]
    valid_u = np.all(vbar2 > 0.0, axis=1)
    vb_safe = np.where(vbar2 > 0.0, vbar2, 1.0)
    m0 = vb_safe[None, :, :] + np.sqrt(
        vb_safe[None, :, :] * a0[:, None, None] * r_m[None, None, :]
    )

    S1 = np.empty((len(a0), len(upsilon), len(d_grid), len(i_grid)))
    S2 = np.empty((len(a0), len(upsilon), len(d_grid), len(i_grid)))
    sin_cat = math.sin(math.radians(gal["meta"]["inclination_deg"]))
    varlog = np.empty(len(i_grid))
    for j, d in enumerate(d_grid):
        v_model = np.sqrt(d * m0)
        for l, i_val in enumerate(i_grid):
            s = sin_cat / math.sin(math.radians(i_val))
            u = (vobs * s - v_model) / (ev * s)
            S1[:, :, j, l] = (u * u).sum(axis=2)
            S2[:, :, j, l] = u.sum(axis=2) ** 2
            varlog[l] = 2.0 * npts * math.log(s)
    S1[:, ~valid_u, :, :] = np.inf
    return S1, S2, varlog, npts


def profile_over_nuisances(
    S1: np.ndarray,
    S2: np.ndarray,
    varlog: np.ndarray,
    penalty: np.ndarray,
    npts: int,
    rho: float,
) -> dict[str, np.ndarray]:
    if not (0.0 <= rho < 1.0):
        raise ValueError(f"correlation parameter out of range [0, 1): {rho}")
    w = rho / (1.0 + (npts - 1) * rho)
    chi2 = (S1 - w * S2) / (1.0 - rho)
    total = chi2 + penalty[None, :, :, :] + varlog[None, None, None, :]
    n_a0 = total.shape[0]
    flat = total.reshape(n_a0, -1)
    argmin_flat = np.argmin(flat, axis=1)
    rows = np.arange(n_a0)
    return {
        "curve": flat[rows, argmin_flat],
        "chi2_data": chi2.reshape(n_a0, -1)[rows, argmin_flat],
        "argmin_flat": argmin_flat,
    }


def vertex(x: np.ndarray, y: np.ndarray, m: int) -> float:
    if 0 < m < len(x) - 1:
        y0, y1, y2 = y[m - 1], y[m], y[m + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom > 0.0:
            step = x[1] - x[0]
            return float(x[m] + 0.5 * step * (y0 - y2) / denom)
    return float(x[m])


def curve_contour(log10_a0_grid: np.ndarray, curve: np.ndarray) -> dict[str, Any]:
    m = int(np.argmin(curve))
    x_min = float(log10_a0_grid[m])
    step = float(log10_a0_grid[1] - log10_a0_grid[0])
    delta = curve - curve[m]

    def components(threshold: float) -> list[dict[str, Any]]:
        inside = delta <= threshold
        result: list[dict[str, Any]] = []
        start: int | None = None
        for idx in range(len(inside) + 1):
            is_inside = bool(inside[idx]) if idx < len(inside) else False
            if is_inside and start is None:
                start = idx
            if not is_inside and start is not None:
                end = idx - 1
                if start == 0:
                    lo, lo_pinned = float(log10_a0_grid[0]), True
                else:
                    x0, x1 = log10_a0_grid[start - 1], log10_a0_grid[start]
                    d0, d1 = delta[start - 1], delta[start]
                    lo = float(x0 + (x1 - x0) * (threshold - d0) / (d1 - d0))
                    lo_pinned = False
                if end == len(delta) - 1:
                    hi, hi_pinned = float(log10_a0_grid[-1]), True
                else:
                    x0, x1 = log10_a0_grid[end], log10_a0_grid[end + 1]
                    d0, d1 = delta[end], delta[end + 1]
                    hi = float(x0 + (x1 - x0) * (threshold - d0) / (d1 - d0))
                    hi_pinned = False
                result.append(
                    {
                        "log10": [lo, hi],
                        "a0_m_s2": [10.0 ** lo, 10.0 ** hi],
                        "pinned_at_grid_boundary": [lo_pinned, hi_pinned],
                        "grid_resolution_limited": bool(hi - lo < step),
                        "contains_objective_min": start <= m <= end,
                    }
                )
                start = None
        return result

    out: dict[str, Any] = {
        "log10_a0_objective_min": x_min,
        "a0_objective_min_m_s2": 10.0 ** x_min,
        "a0_objective_min_parabola_refined_m_s2": 10.0 ** vertex(log10_a0_grid, curve, m),
        "objective_min_grid_index": m,
        "min_penalized_objective": float(curve[m]),
        "reference_sublevel_sets": {},
    }
    for label, threshold in (("68", THRESHOLD_68), ("95", THRESHOLD_95)):
        all_components = components(threshold)
        containing_min = [c for c in all_components if c["contains_objective_min"]]
        if len(containing_min) != 1:
            raise AssertionError("objective minimum must lie in exactly one contour component")
        lo, hi = containing_min[0]["log10"]
        lo_pinned, hi_pinned = containing_min[0]["pinned_at_grid_boundary"]
        out[f"a0_delta_objective_contour{label}_m_s2"] = [10.0 ** lo, 10.0 ** hi]
        out[f"delta_objective_contour{label}_log10"] = [lo, hi]
        out[f"contour{label}_pinned_at_grid_boundary"] = [lo_pinned, hi_pinned]
        out[f"contour{label}_grid_resolution_limited"] = bool(hi - lo < step)
        out[f"contour{label}_all_components_log10"] = [
            component["log10"] for component in all_components
        ]
        out[f"contour{label}_n_components"] = len(all_components)
        out[f"contour{label}_has_disconnected_components"] = len(all_components) > 1
        out["reference_sublevel_sets"][f"delta_{'1' if label == '68' else '3p84'}"] = {
            "threshold": threshold,
            "components": all_components,
            "n_components": len(all_components),
            "global_min_component_index": next(
                idx
                for idx, component in enumerate(all_components)
                if component["contains_objective_min"]
            ),
            "disconnected": len(all_components) > 1,
        }
    return out


def channel_b_at(
    gal: dict[str, Any], grids: dict[str, Any], flat_index: int
) -> dict[str, Any]:
    shape = (len(grids["upsilon"]), len(grids["d_grid"]), len(grids["i_grid"]))
    u_idx, j_idx, l_idx = np.unravel_index(flat_index, shape)
    upsilon = float(grids["upsilon"][u_idx])
    d = float(grids["d_grid"][j_idx])
    i_val = float(grids["i_grid"][l_idx])
    meta = gal["meta"]
    k = gal["idx_out"]
    s = math.sin(math.radians(meta["inclination_deg"])) / math.sin(
        math.radians(i_val)
    )
    vobs = gal["vobs_kms"][k] * KM * s
    vbar2 = d * (
        gal["gas_term"][k] + upsilon * gal["disk_term"][k] + gal["bul_term"][k]
    )
    v_a2 = vobs * vobs - vbar2
    mass = (
        (upsilon * meta["L36_gsun"] + GAS_HELIUM * meta["MHI_gsun"])
        * 1e9
        * M_SUN
        * d
        * d
    )
    return {
        "v_a2": float(v_a2),
        "mass": float(mass),
        "a0_b": float(v_a2 * v_a2 / (G_SI * mass)) if v_a2 > 0 and mass > 0 else None,
    }


def verdict_label(interval_low: float, interval_high: float) -> str:
    if interval_low > 0.0:
        return "TENSION_JOINT_CHANNEL_B_ABOVE_CHANNEL_A"
    if interval_high < 0.0:
        return "TENSION_JOINT_CHANNEL_B_BELOW_CHANNEL_A"
    return "CONSISTENT_INTERVAL_CONTAINS_ZERO"


VERDICT_RULE = (
    "The 95 percent paired bootstrap interval for log10(a0_B/a0_A) decides "
    "the label: entirely above zero gives "
    "TENSION_JOINT_CHANNEL_B_ABOVE_CHANNEL_A, entirely below zero gives "
    "TENSION_JOINT_CHANNEL_B_BELOW_CHANNEL_A, otherwise "
    "CONSISTENT_INTERVAL_CONTAINS_ZERO. The rule is symmetric in direction. "
    "An interval containing zero is compatible with every model that "
    "reproduces both relations, including the standard null, and is not "
    "evidence for OPH. A tension counts only against this declared submodel "
    "with its declared conventions, never against the framework beyond that "
    "scope."
)


def bootstrap_block(
    log10_a0_grid: np.ndarray,
    common_curves: np.ndarray,
    log_a0_b: np.ndarray,
) -> dict[str, Any]:
    n = common_curves.shape[0]
    if n < 2:
        raise ValueError("paired bootstrap needs at least two common galaxies")
    total = common_curves.sum(axis=0)
    x_ml = float(log10_a0_grid[int(np.argmin(total))])
    a0_a = 10.0 ** x_ml
    a0_b = float(10.0 ** np.mean(log_a0_b))
    point_ratio = math.log10(a0_b / a0_a)
    rng = np.random.default_rng(SEED)
    ratios = np.empty(REPLICATES)
    for t in range(REPLICATES):
        idx = rng.integers(0, n, n)
        x_rep = float(log10_a0_grid[int(np.argmin(common_curves[idx].sum(axis=0)))])
        b_rep = float(np.mean(log_a0_b[idx]))
        ratios[t] = b_rep - x_rep
    interval = np.percentile(ratios, [2.5, 97.5])
    n_le = int(np.count_nonzero(ratios <= 0.0))
    n_ge = int(np.count_nonzero(ratios >= 0.0))
    return {
        "n_common_galaxies": n,
        "channel_a_common_a0_m_s2": float(a0_a),
        "channel_b_common_a0_m_s2": a0_b,
        "point_log10_ratio_b_over_a": point_ratio,
        "bootstrap_unit": "common-set galaxy, both channels resampled together",
        "bootstrap_replicates": REPLICATES,
        "bootstrap_seed": SEED,
        "log10_ratio_95pct": [float(interval[0]), float(interval[1])],
        "count_log10_ratio_le_zero": n_le,
        "le_zero_plus_one_fraction": float((n_le + 1) / (REPLICATES + 1)),
        "count_log10_ratio_ge_zero": n_ge,
        "ge_zero_plus_one_fraction": float((n_ge + 1) / (REPLICATES + 1)),
        "verdict": verdict_label(float(interval[0]), float(interval[1])),
    }


def paired_block(
    per_gal: list[dict[str, Any]],
    gal_by_name: dict[str, dict[str, Any]],
    grids_by_name: dict[str, dict[str, Any]],
    curves: np.ndarray,
    log10_a0_grid: np.ndarray,
    ml_index: int,
    rho: float,
    fraction: float,
) -> dict[str, Any]:
    n_not_deep = 0
    n_nonpositive_va2 = 0
    n_nonpositive_mb = 0
    excluded_bulge_luminosity_ambiguous: list[str] = []
    used_names: list[str] = []
    used_logs: list[float] = []
    used_rows: list[int] = []
    n_candidates = 0
    for row_idx, g in enumerate(per_gal):
        gal = gal_by_name[g["name"]]
        n_candidates += 1
        k = gal["idx_out"]
        if not (gal["g_bar_cat"][k] < fraction * REFERENCE_A0):
            n_not_deep += 1
            continue
        if bool(np.any(gal["bul_term"] != 0.0)):
            excluded_bulge_luminosity_ambiguous.append(g["name"])
            continue
        grids = grids_by_name[g["name"]]
        flat_index = int(g["rho"][rho]["argmin_flat"][ml_index])
        est = channel_b_at(gal, grids, flat_index)
        if est["v_a2"] <= 0:
            n_nonpositive_va2 += 1
            continue
        if est["mass"] <= 0:
            n_nonpositive_mb += 1
            continue
        used_names.append(g["name"])
        used_logs.append(math.log10(est["a0_b"]))
        used_rows.append(row_idx)
    common_curves = curves[np.array(used_rows, dtype=int)]
    log_a0_b = np.array(used_logs)
    paired = bootstrap_block(log10_a0_grid, common_curves, log_a0_b)
    return {
        "n_candidate_galaxies": n_candidates,
        "n_outermost_not_deep": n_not_deep,
        "n_nonpositive_anomalous_speed2": n_nonpositive_va2,
        "n_nonpositive_baryonic_mass": n_nonpositive_mb,
        "n_excluded_bulge_luminosity_ambiguous": len(
            excluded_bulge_luminosity_ambiguous
        ),
        "excluded_bulge_luminosity_ambiguous": sorted(
            excluded_bulge_luminosity_ambiguous
        ),
        "n_used": len(used_names),
        "a0_b_unweighted_log_mean_m_s2": float(10.0 ** np.mean(log_a0_b)),
        **paired,
    }


def compute() -> dict[str, Any]:
    t1 = parse_table1(DATA / "table1.dat")
    t2 = parse_table2(DATA / "table2.dat")
    galaxies, n_bad_ev = assemble_galaxies(t1, t2)
    log10_a0_grid = np.linspace(LOG10_A0_MIN, LOG10_A0_MAX, N_A0)
    grids_by_name = {gal["name"]: make_grids(gal["meta"]) for gal in galaxies}
    n_fallback_distance = sum(
        1 for g in grids_by_name.values() if g["distance_fixed"]
    )
    n_fallback_inclination = sum(
        1 for g in grids_by_name.values() if g["inclination_fixed"]
    )

    subsets: list[tuple[str, float | None]] = [("full", None)]
    for f in FRACTIONS:
        subsets.append((f"deep_f_{str(f).replace('.', 'p')}", f))

    subset_records: dict[str, dict[str, Any]] = {}
    for subset_name, fraction in subsets:
        per_gal: list[dict[str, Any]] = []
        n_points = 0
        for gal in galaxies:
            if fraction is None:
                mask = np.ones(len(gal["rad_kpc"]), dtype=bool)
            else:
                mask = gal["g_bar_cat"] < fraction * REFERENCE_A0
            npts = int(mask.sum())
            if npts == 0:
                continue
            grids = grids_by_name[gal["name"]]
            S1, S2, varlog, _ = block_stats(gal, mask, log10_a0_grid, grids)
            rec: dict[str, Any] = {"name": gal["name"], "npts": npts, "rho": {}}
            for rho in RHO_GRID:
                rec["rho"][rho] = profile_over_nuisances(
                    S1, S2, varlog, grids["penalty"], npts, rho
                )
            per_gal.append(rec)
            n_points += npts
        subset_records[subset_name] = {
            "fraction": fraction,
            "per_gal": per_gal,
            "n_points": n_points,
            "n_galaxies": len(per_gal),
        }

    n_free_by_subset = {}
    for subset_name, rec in subset_records.items():
        total = 0
        for g in rec["per_gal"]:
            grids = grids_by_name[g["name"]]
            total += 1
            total += 0 if grids["distance_fixed"] else 1
            total += 0 if grids["inclination_fixed"] else 1
        n_free_by_subset[subset_name] = total

    gal_by_name = {gal["name"]: gal for gal in galaxies}
    subset_results: dict[str, Any] = {}
    for subset_name, rec in subset_records.items():
        fraction = rec["fraction"]
        per_gal = rec["per_gal"]
        per_rho_rows: list[dict[str, Any]] = []
        for rho in RHO_GRID:
            curves = np.stack([g["rho"][rho]["curve"] for g in per_gal])
            chi2_data = np.stack([g["rho"][rho]["chi2_data"] for g in per_gal])
            total_curve = curves.sum(axis=0)
            contour = curve_contour(log10_a0_grid, total_curve)
            m = contour["objective_min_grid_index"]
            chi2_data_at_min = float(chi2_data[:, m].sum())
            n_edge = 0
            for g in per_gal:
                grids = grids_by_name[g["name"]]
                shape = (
                    len(grids["upsilon"]),
                    len(grids["d_grid"]),
                    len(grids["i_grid"]),
                )
                u_idx, j_idx, l_idx = np.unravel_index(
                    int(g["rho"][rho]["argmin_flat"][m]), shape
                )
                at_edge = u_idx in (0, shape[0] - 1)
                if shape[1] > 1:
                    at_edge = at_edge or j_idx in (0, shape[1] - 1)
                if shape[2] > 1:
                    at_edge = at_edge or l_idx in (0, shape[2] - 1)
                if at_edge:
                    n_edge += 1
            dof_upper = rec["n_points"] - 1
            dof_lower = rec["n_points"] - 1 - n_free_by_subset[subset_name]
            row: dict[str, Any] = {
                "rho": rho,
                **contour,
                "chi2_data_at_objective_min": chi2_data_at_min,
                "dof_no_nuisance_count": dof_upper,
                "dof_full_nuisance_count": dof_lower,
                "reduced_chi2_dof_no_nuisance": chi2_data_at_min / dof_upper,
                "reduced_chi2_dof_full_nuisance": (
                    chi2_data_at_min / dof_lower if dof_lower > 0 else None
                ),
                "dof_convention": (
                    "two nominal sensitivity denominators, with the profiled "
                    "nuisances counted either fixed or free; penalties and "
                    "grid boundaries mean neither denominator is a calibrated "
                    "effective dof or a proved bound"
                ),
                "n_galaxies_profiled_nuisance_at_grid_edge_at_objective_min": n_edge,
            }
            if fraction is not None:
                row["paired_btfr"] = paired_block(
                    per_gal,
                    gal_by_name,
                    grids_by_name,
                    curves,
                    log10_a0_grid,
                    m,
                    rho,
                    fraction,
                )
            per_rho_rows.append(row)
        subset_results[subset_name] = {
            "fraction_of_reference": fraction,
            "deep_cut_g_bar_max_m_s2": (
                None if fraction is None else fraction * REFERENCE_A0
            ),
            "n_points": rec["n_points"],
            "n_galaxies": rec["n_galaxies"],
            "n_free_nuisance_parameters": n_free_by_subset[subset_name],
            "per_rho": per_rho_rows,
        }

    observed = {}
    for subset_name, block in subset_results.items():
        rows = []
        for row in block["per_rho"]:
            rows.append(
                {
                    "rho": row["rho"],
                    "reduced_chi2_dof_no_nuisance": row[
                        "reduced_chi2_dof_no_nuisance"
                    ],
                    "reduced_chi2_dof_full_nuisance": row[
                        "reduced_chi2_dof_full_nuisance"
                    ],
                    "contour95_pinned_at_grid_boundary": row[
                        "contour95_pinned_at_grid_boundary"
                    ],
                    "n_galaxies_profiled_nuisance_at_grid_edge_at_objective_min": row[
                        "n_galaxies_profiled_nuisance_at_grid_edge_at_objective_min"
                    ],
                }
            )
        observed[subset_name] = rows

    return {
        "schema": SCHEMA,
        "scope": (
            "labeled_postdiction_penalized_profile_objective_declared_conventions_"
            "fixed_absolute_cuts"
        ),
        "physical_claim": False,
        "source_derived_output": False,
        "seen_data_postdiction": True,
        "sample": {
            "catalogue": "SPARC, CDS J/AJ/152/157 (Lelli, McGaugh, Schombert 2016)",
            "cuts": (
                "quality <= 2, catalogue inclination >= 30 deg, v_obs > 0, "
                "e_Vobs > 0, e_Vobs/v_obs <= 0.10, g_bar_cat > 0 at the "
                "committed mass-to-light values"
            ),
            "n_retained_points": subset_records["full"]["n_points"],
            "n_retained_galaxies": subset_records["full"]["n_galaxies"],
            "n_points_excluded_nonpositive_e_vobs": n_bad_ev,
            "outermost_point": (
                "largest radius among the retained points of each galaxy"
            ),
            "data_sha256": hashlib.sha256(
                (DATA / "table1.dat").read_bytes()
                + (DATA / "table2.dat").read_bytes()
            ).hexdigest(),
        },
        "model": {
            "acceleration_law": "g_model = g_bar + sqrt(g_bar a0)",
            "velocity_model": (
                "v_model^2(r) = d [v_bar_cat^2(r) + "
                "sqrt(v_bar_cat^2(r) a0 r_cat)] at r = d r_cat"
            ),
            "baryonic_speed2": (
                "v_bar_cat^2 = v_gas|v_gas| + Upsilon_d v_disk|v_disk| + "
                "0.7 v_bul|v_bul|, signed gas convention, catalogue "
                "components at unit mass-to-light"
            ),
            "distance_scaling": (
                "r = d r_cat from fixed angular sizes; M proportional to "
                "D^2 from photometry, so v_bar^2 = d v_bar_cat^2 and the "
                "anomalous term sqrt(v_bar^2 a0 r) also scales as d; the "
                "table-1 baryonic mass scales as d^2"
            ),
            "inclination_scaling": (
                "v_obs(i) = v_obs_cat sin(i_cat)/sin(i) and e_Vobs scales by "
                "the same factor; the differential form dv/v = -cot(i) di "
                "matches the committed cot(i) error factor"
            ),
        },
        "error_model": {
            "per_point": (
                "Gaussian in velocity with sigma = e_Vobs sin(i_cat)/sin(i); "
                "the 2 n ln s(i) variance-normalization term is kept in the "
                "penalized objective and a0-independent constants are dropped"
            ),
            "intra_galaxy_covariance": (
                "equal-correlation block R = (1-rho) I + rho J per galaxy, "
                "closed-form quadratic form; rho is a declared family "
                "parameter scanned over a fixed grid; a full covariance "
                "calibration for SPARC rotation curves is open, so contours "
                "are reported at every rho and no single rho is selected"
            ),
            "rho_grid": list(RHO_GRID),
            "not_propagated": (
                "luminosity and HI-mass measurement errors, "
                "stellar-population systematics beyond the declared "
                "mass-to-light prior, baryonic-model decomposition errors, "
                "and inter-galaxy correlations"
            ),
        },
        "nuisances": {
            "disk_mass_to_light": {
                "prior": (
                    "log-normal centered at 0.5 with width 0.1 dex, declared"
                ),
                "grid": (
                    f"{N_UPSILON} points over +-{HALFWIDTH_SIGMA} prior sigma"
                ),
                "status": "declared",
            },
            "bulge_mass_to_light": {
                "value": UPSILON_BULGE,
                "status": (
                    "declared fixed fallback; the snapshot carries no "
                    "per-galaxy bulge datum"
                ),
            },
            "distance": {
                "prior": (
                    "Gaussian factor prior at 1 with the catalogue relative "
                    "error e_Dist/Dist per galaxy"
                ),
                "grid": (
                    f"{N_DISTANCE} points over +-{HALFWIDTH_SIGMA} sigma, "
                    f"clipped below at {DIST_FACTOR_MIN}"
                ),
                "status": "catalogue-supplied",
                "n_galaxies_fallback_fixed": n_fallback_distance,
            },
            "inclination": {
                "prior": (
                    "Gaussian at the catalogue inclination with the "
                    "catalogue error e_i per galaxy"
                ),
                "grid": (
                    f"{N_INCLINATION} points over +-{HALFWIDTH_SIGMA} sigma, "
                    f"clipped to [{INCLINATION_CLIP[0]}, "
                    f"{INCLINATION_CLIP[1]}] deg"
                ),
                "status": "catalogue-supplied",
                "n_galaxies_fallback_fixed": n_fallback_inclination,
            },
            "subset_membership": (
                "deep subsets fixed before fitting by g_bar_cat < f * "
                "1.2e-10 m/s^2 at catalogue nuisance values; nuisance shifts "
                "do not re-select points"
            ),
        },
        "grids": {
            "log10_a0": {
                "min": LOG10_A0_MIN,
                "max": LOG10_A0_MAX,
                "n": N_A0,
                "step_dex": (LOG10_A0_MAX - LOG10_A0_MIN) / (N_A0 - 1),
            },
        },
        "calibration": {
            "reference_delta_objective_thresholds": {
                "delta_1": THRESHOLD_68,
                "delta_3p84": THRESHOLD_95,
            },
            "status": (
                "uncalibrated reference contours on a penalized profile "
                "objective; the mass-to-light term is an astrophysical prior, "
                "not an identified auxiliary-data likelihood, so Wilks "
                "coverage and maximum-likelihood labels are not licensed; "
                "delta is measured from the grid minimum, the parabola vertex "
                "is a numerical refinement row, and contours narrower than "
                "one grid step carry a grid-resolution flag"
            ),
            "confidence_interval_claim": False,
            "coverage_calibration_open": True,
        },
        "seed": SEED,
        "subset_results": subset_results,
        "btfr_channel": {
            "estimator": (
                "a0_B = (v_A^2)^2 / (G M_b d^2), "
                "v_A^2 = (v_obs(r_out) s(i))^2 - d v_bar_cat^2(r_out, U), "
                "M_b = (U L_[3.6] + 1.33 M_HI) 1e9 M_sun at catalogue "
                "distance"
            ),
            "bulge_exclusion": (
                "channel B excludes every galaxy with any nonzero retained "
                "V_bulge because table 1 L_[3.6] is total luminosity and the "
                "snapshot supplies no disk/bulge luminosity split; retaining "
                "such galaxies would mix a 0.7 bulge subtraction with a "
                "disk-Upsilon total-luminosity denominator"
            ),
            "nuisance_application": (
                "per galaxy the profiled (U, d, i) at the subset "
                "penalized-objective-minimizing a0 for the same rho; "
                "channel-B nuisance "
                "values are held at these point estimates inside the "
                "bootstrap, a declared simplification, while channel A "
                "re-profiles a0 exactly per replicate"
            ),
            "deep_requirement": (
                "outermost retained point passes the same fixed absolute "
                "cut g_bar_cat < f * 1.2e-10 m/s^2 as the channel-A subset"
            ),
        },
        "misfit_statement": {
            "rule": (
                "Under this model and error family, misfit would appear as "
                "reduced chi-square far above one at every rho on the "
                "declared grid, as contour endpoints pinned at the declared "
                "grid boundary, or as profiled nuisances stacked at their "
                "grid edges across the sample. These are objective/residual "
                "diagnostics; no confidence or posterior probability is "
                "computed or implied."
            ),
            "observed": observed,
        },
        "verdict_rule": VERDICT_RULE,
        "inference_boundary": {
            "status": (
                "labeled postdiction penalized objective on the seen committed "
                "snapshot; no new data, no arming, no discharge, no scored "
                "comparison, no frozen contract"
            ),
            "a0_source_status": (
                "the OPH source does not fix a0; nothing here derives a "
                "source value, so a0 is a fitted comparison parameter and "
                "the source-value question is open"
            ),
            "covariance_status": (
                "the equal-correlation family is a declared stand-in; a "
                "calibrated covariance model for SPARC rotation curves is "
                "open, and contour widths at different rho values diagnose "
                "the sensitivity without furnishing confidence coverage"
            ),
            "full_sample_scope": (
                "the full-sample row extrapolates the deep-regime profile "
                "to all gradients; its misfit bears on that additive "
                "extension submodel under the declared error family, never "
                "on the framework beyond that scope, and the deep subsets "
                "are the derived regime of the candidate law"
            ),
            "neutrality": (
                "consistency of the channels or agreement of the contour "
                "with any external value is shared with every model that "
                "reproduces the same relations, including the standard "
                "null, and is not evidence for OPH"
            ),
            "not_done_here": (
                "no relativistic lensing, no cosmological observables, no "
                "cluster or CMB likelihood, no covariance calibration, no "
                "source derivation of a0, no preregistration"
            ),
        },
    }


def first_difference(a: Any, b: Any, path: str = "$") -> str | None:
    if type(a) is not type(b):
        return f"{path}: type {type(a).__name__} vs {type(b).__name__}"
    if isinstance(a, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a:
                return f"{path}.{key}: only in receipt"
            if key not in b:
                return f"{path}.{key}: only in recomputation"
            found = first_difference(a[key], b[key], f"{path}.{key}")
            if found:
                return found
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            found = first_difference(x, y, f"{path}[{i}]")
            if found:
                return found
        return None
    if a != b:
        return f"{path}: {a!r} vs {b!r}"
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    args = ap.parse_args()
    committed = args.receipt.read_bytes()
    fresh_obj = compute()
    fresh = (json.dumps(fresh_obj, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if fresh == committed:
        digest = hashlib.sha256(committed).hexdigest()
        print(f"PASS: byte-identical recomputation, sha256 {digest}")
        return 0
    diff = first_difference(json.loads(committed), fresh_obj)
    print(f"FAIL: recomputation differs from receipt at {diff}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
