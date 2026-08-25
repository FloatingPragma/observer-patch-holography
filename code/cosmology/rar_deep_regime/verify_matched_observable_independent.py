#!/usr/bin/env python3
"""Independent verifier for the matched-observable dark-sector receipt.

This script imports neither ``matched_observable_diagnostic.py`` nor the
committed ``rar_deep_regime_diagnostic.py``.  It re-parses the snapshot with
its own fixed-width parser built from the CDS ReadMe byte ranges, recomputes
every matched-observable number in the receipt from the raw tables (the
three quoted mixed-proxy values of the superseded block are carried over
verbatim as the committed record they quote), serializes the result as
canonical JSON (sorted keys, two-space indent, trailing newline), and
byte-compares it against ``runtime/matched_observable_receipt.json``.

Byte equality of floating-point output requires the same arithmetic
evaluation order as the producer contract; the formulas below fix that
order.  The verification therefore certifies that the committed receipt is
exactly the declared computation on the committed snapshot bytes, replayed
by code written separately from the producer.

What is not proved here.  Replay equality does not authenticate provenance
or custody of the snapshot, does not make the diagnostic a likelihood, and
does not turn the postdiction into evidence for or against OPH.
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
from scipy import optimize

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RECEIPT_PATH = HERE / "runtime" / "matched_observable_receipt.json"
MIXED_PROXY_RECEIPT = HERE / "receipts" / "sparc_deep_regime_diagnostic.json"

# Physical and convention constants, restated independently of the producer.
KPC_M = 3.0856775814913673e19
KM = 1.0e3
G_SI = 6.67430e-11
M_SUN = 1.98841e30
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7
GAS_HELIUM = 1.33
REFERENCE_A0 = 1.2e-10
LOG10_A0_BOUNDS = (-12.0, -9.0)
PRIMARY_FRACTION = 0.1
SENSITIVITY_FRACTION = 0.3
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260825
SCHEMA = "oph.cosmology.matched_observable_dark_sector_diagnostic.v1"


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


def baryonic_squared_speed_m2s2(row: dict[str, Any]) -> float:
    vg = row["vgas_kms"] * KM
    vd = row["vdisk_kms"] * KM
    vb = row["vbul_kms"] * KM
    return vg * abs(vg) + UPSILON_DISK * vd * abs(vd) + UPSILON_BULGE * vb * abs(vb)


def point_accelerations(
    t1: list[dict[str, Any]], t2: list[dict[str, Any]]
) -> dict[str, np.ndarray]:
    meta = {r["name"]: r for r in t1}
    g_obs, g_bar, names = [], [], []
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > 2 or m["inclination_deg"] < 30.0:
            continue
        if r["vobs_kms"] <= 0 or r["e_vobs_kms"] / r["vobs_kms"] > 0.10:
            continue
        radius = r["rad_kpc"] * KPC_M
        vo = r["vobs_kms"] * KM
        vg = r["vgas_kms"] * KM
        vd = r["vdisk_kms"] * KM
        vb = r["vbul_kms"] * KM
        gb = (
            vg * abs(vg) + UPSILON_DISK * vd * abs(vd) + UPSILON_BULGE * vb * abs(vb)
        ) / radius
        if gb <= 0:
            continue
        g_obs.append(vo * vo / radius)
        g_bar.append(gb)
        names.append(r["name"])
    return {
        "g_obs": np.array(g_obs),
        "g_bar": np.array(g_bar),
        "names": np.array(names),
    }


def fit_total_a0(g_obs: np.ndarray, g_bar: np.ndarray) -> float:
    if len(g_obs) < 2:
        raise ValueError("at least two points are required")

    def objective(log10_a0: float) -> float:
        a0 = 10.0 ** log10_a0
        model = g_bar + np.sqrt(g_bar * a0)
        residual = np.log10(g_obs) - np.log10(model)
        return float(np.mean(residual**2))

    result = optimize.minimize_scalar(
        objective, bounds=LOG10_A0_BOUNDS, method="bounded",
        options={"xatol": 1.0e-14},
    )
    if not result.success:
        raise RuntimeError(f"a0 minimization failed: {result.message}")
    return float(10.0 ** result.x)


def outermost_rows(
    t1: list[dict[str, Any]], t2: list[dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    meta = {r["name"]: r for r in t1}
    rows: dict[str, dict[str, Any]] = {}
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > 2 or m["inclination_deg"] < 30.0:
            continue
        if r["vobs_kms"] <= 0 or r["e_vobs_kms"] / r["vobs_kms"] > 0.10:
            continue
        vo = r["vobs_kms"] * KM
        vbar2 = baryonic_squared_speed_m2s2(r)
        v_a2 = vo * vo - vbar2
        radius_m = r["rad_kpc"] * KPC_M
        g_bar = vbar2 / radius_m
        if g_bar <= 0:
            continue
        cur = rows.get(r["name"])
        if cur is None or r["rad_kpc"] > cur["rad_kpc"]:
            rows[r["name"]] = dict(
                r, v_a2=v_a2, vbar2=vbar2, g_bar=g_bar, radius_m=radius_m
            )
    return rows


def channel_b(
    t1: list[dict[str, Any]], t2: list[dict[str, Any]], fraction: float
) -> dict[str, Any]:
    meta = {r["name"]: r for r in t1}
    rows = outermost_rows(t1, t2)
    names: list[str] = []
    log_a0: list[float] = []
    sigma_log: list[float] = []
    log_a0_same_radius: list[float] = []
    n_not_deep = 0
    n_nonpositive_va2 = 0
    n_nonpositive_mb = 0
    for name, r in rows.items():
        if not (r["g_bar"] < fraction * REFERENCE_A0):
            n_not_deep += 1
            continue
        if r["v_a2"] <= 0:
            n_nonpositive_va2 += 1
            continue
        m = meta[name]
        mstar = UPSILON_DISK * m["L36_gsun"] * 1e9
        mgas = GAS_HELIUM * m["MHI_gsun"] * 1e9
        mass = (mstar + mgas) * M_SUN
        if mass <= 0:
            n_nonpositive_mb += 1
            continue
        a0_b = r["v_a2"] * r["v_a2"] / (G_SI * mass)
        vo = r["vobs_kms"] * KM
        ev = r["e_vobs_kms"] * KM
        inc_rad = math.radians(m["inclination_deg"])
        e_inc_rad = math.radians(m["e_inclination_deg"])
        sigma_va2 = math.hypot(
            2.0 * vo * ev, 2.0 * vo * vo * e_inc_rad / math.tan(inc_rad)
        )
        names.append(name)
        log_a0.append(math.log10(a0_b))
        sigma_log.append((2.0 / math.log(10.0)) * sigma_va2 / r["v_a2"])
        log_a0_same_radius.append(
            math.log10(r["v_a2"] * r["v_a2"] / (r["vbar2"] * r["radius_m"]))
        )
    names_arr = np.array(names)
    log_arr = np.array(log_a0)
    sigma_arr = np.array(sigma_log)
    same_radius_arr = np.array(log_a0_same_radius)
    if len(names) == 0:
        empty_summary = {
            "n_candidate_galaxies": int(len(rows)),
            "n_outermost_not_deep": n_not_deep,
            "n_nonpositive_anomalous_speed2": n_nonpositive_va2,
            "n_nonpositive_baryonic_mass": n_nonpositive_mb,
            "n_used": 0,
            "a0_unweighted_log_mean_m_s2": None,
            "empty_selection": (
                "no galaxy passes the channel-B requirements at this cut; "
                "no estimate is reported"
            ),
        }
        return {"names": names_arr, "log_a0": log_arr, "summary": empty_summary}
    weights = 1.0 / sigma_arr**2
    scatter = float(np.std(log_arr))
    median_sigma = float(np.median(sigma_arr))
    return {
        "names": names_arr,
        "log_a0": log_arr,
        "summary": {
            "estimator": "a0_B = (v_obs(r_out)^2 - v_bar(r_out)^2)^2 / (G M_b)",
            "baryonic_mass": (
                "M_b = 0.5 L_[3.6] + 1.33 M_HI from table 1, committed "
                "declared mass-to-light and helium choices"
            ),
            "deep_outermost_requirement_g_bar_max_m_s2": fraction * REFERENCE_A0,
            "n_candidate_galaxies": int(len(rows)),
            "n_outermost_not_deep": n_not_deep,
            "n_nonpositive_anomalous_speed2": n_nonpositive_va2,
            "n_nonpositive_baryonic_mass": n_nonpositive_mb,
            "n_used": int(len(names)),
            "a0_unweighted_log_mean_m_s2": float(10.0 ** np.mean(log_arr)),
            "scatter_dex": scatter,
            "median_propagated_sigma_log10_dex": median_sigma,
            "scatter_dominated_by_unpropagated_terms": bool(
                scatter > 2.0 * median_sigma
            ),
            "error_model": (
                "sigma(v_A^2) = hypot(2 v_obs e_Vobs, 2 v_obs^2 cot(i) e_i); "
                "distance, mass-to-light, and baryonic-model errors are not "
                "propagated"
            ),
            "a0_inverse_variance_weighted_m_s2": float(
                10.0 ** (np.sum(weights * log_arr) / np.sum(weights))
            ),
            "inverse_variance_weighting_caution": (
                "sensitivity row only: the observed scatter is dominated by "
                "unpropagated terms, and measurement-error weights correlate "
                "with the estimate, so the primary combination is unweighted"
            ),
            "a0_same_radius_mass_variant_m_s2": float(
                10.0 ** np.mean(same_radius_arr)
            ),
            "same_radius_mass_variant": (
                "G M_b replaced by v_bar^2(r_out) r_out; exact in the "
                "exterior point-mass limit; sensitivity row only"
            ),
        },
    }


def verdict_label(interval_low: float, interval_high: float) -> str:
    if interval_low > 0.0:
        return "TENSION_MATCHED_CHANNEL_B_ABOVE_CHANNEL_A"
    if interval_high < 0.0:
        return "TENSION_MATCHED_CHANNEL_B_BELOW_CHANNEL_A"
    return "CONSISTENT_INTERVAL_CONTAINS_ZERO"


VERDICT_RULE = (
    "The 95 percent paired interval for log10(a0_B/a0_A) decides the label: "
    "entirely above zero gives TENSION_MATCHED_CHANNEL_B_ABOVE_CHANNEL_A, "
    "entirely below zero gives TENSION_MATCHED_CHANNEL_B_BELOW_CHANNEL_A, "
    "otherwise CONSISTENT_INTERVAL_CONTAINS_ZERO. The rule is symmetric in "
    "direction. An interval containing zero is compatible with every model "
    "that reproduces both relations, including the standard null, and is "
    "not evidence for OPH."
)


def paired_common_bootstrap(
    acc: dict[str, np.ndarray], b_channel: dict[str, Any], fraction: float
) -> dict[str, Any]:
    mask = acc["g_bar"] < fraction * REFERENCE_A0
    deep_obs = acc["g_obs"][mask]
    deep_bar = acc["g_bar"][mask]
    deep_names = acc["names"][mask]
    common = np.intersect1d(np.unique(deep_names), b_channel["names"])
    common_mask = np.isin(deep_names, common)
    d_obs = deep_obs[common_mask]
    d_bar = deep_bar[common_mask]
    d_names = deep_names[common_mask]
    b_keep = np.isin(b_channel["names"], common)
    b_names = b_channel["names"][b_keep]
    b_log = b_channel["log_a0"][b_keep]

    a0_a_common = fit_total_a0(d_obs, d_bar)
    a0_b_common = float(10.0 ** np.mean(b_log))
    point_ratio = math.log10(a0_b_common / a0_a_common)

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    ratios: list[float] = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sampled = rng.choice(common, size=len(common), replace=True)
        point_indices: list[np.ndarray] = []
        sampled_logs: list[float] = []
        for name in sampled:
            point_indices.append(np.flatnonzero(d_names == name))
            sampled_logs.append(float(b_log[np.flatnonzero(b_names == name)[0]]))
        joined = np.concatenate(point_indices)
        a0_a = fit_total_a0(d_obs[joined], d_bar[joined])
        a0_b = float(10.0 ** np.mean(sampled_logs))
        ratios.append(math.log10(a0_b / a0_a))
    ratios_arr = np.asarray(ratios)
    interval = np.percentile(ratios_arr, [2.5, 97.5])
    n_le_zero = int(np.count_nonzero(ratios_arr <= 0.0))
    n_ge_zero = int(np.count_nonzero(ratios_arr >= 0.0))
    return {
        "n_common_galaxies": int(len(common)),
        "n_channel_a_points_common": int(len(d_obs)),
        "channel_a_common_a0_m_s2": float(a0_a_common),
        "channel_b_common_a0_m_s2": a0_b_common,
        "point_log10_ratio_b_over_a": point_ratio,
        "bootstrap_unit": "common-set galaxy, both channels resampled together",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "log10_ratio_95pct": [float(interval[0]), float(interval[1])],
        "count_log10_ratio_le_zero": n_le_zero,
        "le_zero_plus_one_fraction": float(
            (n_le_zero + 1) / (BOOTSTRAP_REPLICATES + 1)
        ),
        "count_log10_ratio_ge_zero": n_ge_zero,
        "ge_zero_plus_one_fraction": float(
            (n_ge_zero + 1) / (BOOTSTRAP_REPLICATES + 1)
        ),
        "verdict": verdict_label(float(interval[0]), float(interval[1])),
    }


def fraction_block(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    acc: dict[str, np.ndarray],
    fraction: float,
) -> dict[str, Any]:
    mask = acc["g_bar"] < fraction * REFERENCE_A0
    a0_a = fit_total_a0(acc["g_obs"][mask], acc["g_bar"][mask])
    b_res = channel_b(t1, t2, fraction)
    paired = paired_common_bootstrap(acc, b_res, fraction)
    return {
        "fraction_of_reference": fraction,
        "deep_cut_g_bar_max_m_s2": fraction * REFERENCE_A0,
        "channel_a_rar": {
            "model": "g_model = g_bar + sqrt(g_bar a0), unweighted log residuals",
            "a0_m_s2": float(a0_a),
            "n_points": int(mask.sum()),
            "n_galaxies": int(len(np.unique(acc["names"][mask]))),
        },
        "channel_b_matched": b_res["summary"],
        "paired_common_set": paired,
    }


def resolution_statement(
    primary_verdict: str, mixed_interval: list[float], matched_interval: list[float]
) -> str:
    mixed = (
        f"the superseded mixed-proxy paired interval was "
        f"[{mixed_interval[0]:.3f}, {mixed_interval[1]:.3f}] dex"
    )
    matched = (
        f"the matched-observable paired interval is "
        f"[{matched_interval[0]:.3f}, {matched_interval[1]:.3f}] dex"
    )
    if primary_verdict == "CONSISTENT_INTERVAL_CONTAINS_ZERO":
        return (
            f"{matched} and contains zero, while {mixed} and excluded zero. "
            "Under this diagnostic the mixed-proxy displacement is "
            "attributable to the finite-radius total-speed proxy, not to a "
            "normalization disagreement between the two relations. This "
            "consistency is shared with the standard null and is not "
            "evidence for OPH."
        )
    if primary_verdict == "TENSION_MATCHED_CHANNEL_B_ABOVE_CHANNEL_A":
        return (
            f"{matched} and stays above zero; {mixed}. The proxy mismatch "
            "does not account for the full displacement, and a positive "
            "tension persists under this diagnostic."
        )
    return (
        f"{matched} and stays below zero; {mixed}. The displacement changes "
        "sign under the matched observable, and a negative tension is "
        "present under this diagnostic."
    )


def compute() -> dict[str, Any]:
    t1 = parse_table1(DATA / "table1.dat")
    t2 = parse_table2(DATA / "table2.dat")
    acc = point_accelerations(t1, t2)
    primary = fraction_block(t1, t2, acc, PRIMARY_FRACTION)
    sensitivity = fraction_block(t1, t2, acc, SENSITIVITY_FRACTION)
    data_sha = hashlib.sha256(
        (DATA / "table1.dat").read_bytes() + (DATA / "table2.dat").read_bytes()
    ).hexdigest()
    mixed_bytes = MIXED_PROXY_RECEIPT.read_bytes()
    mixed = json.loads(mixed_bytes)
    mixed_paired = mixed["constant_consistency"]["paired_parent_galaxy_bootstrap"]
    mixed_interval = mixed_paired["log10_ratio_btfr_over_rar_95pct"]
    matched_interval = primary["paired_common_set"]["log10_ratio_95pct"]
    primary_verdict = primary["paired_common_set"]["verdict"]
    mixed_point = mixed["constant_consistency"]["log10_ratio_btfr_over_deep"]
    matched_point = primary["paired_common_set"]["point_log10_ratio_b_over_a"]
    return {
        "schema": SCHEMA,
        "scope": (
            "labeled_postdiction_diagnostic_matched_anomalous_observable_"
            "fixed_absolute_cuts"
        ),
        "physical_claim": False,
        "source_derived_output": False,
        "seen_data_postdiction": True,
        "sample": {
            "catalogue": "SPARC, CDS J/AJ/152/157 (Lelli, McGaugh, Schombert 2016)",
            "data_sha256": data_sha,
            "cuts": (
                "quality <= 2, inclination >= 30 deg, e_Vobs/Vobs <= 0.10, "
                "g_bar > 0; identical for both channels"
            ),
            "upsilon_disk": UPSILON_DISK,
            "upsilon_bulge": UPSILON_BULGE,
            "gas_helium_factor": GAS_HELIUM,
            "mass_to_light_status": (
                "committed declared comparison choices, not source-derived"
            ),
            "outermost_point": (
                "largest radius among the retained points of each galaxy, so "
                "the same point-level cuts govern both channels"
            ),
        },
        "matched_observable": {
            "definition": (
                "v_A^2(r_out) = v_obs(r_out)^2 - v_bar(r_out)^2 with "
                "v_bar^2 = v_gas|v_gas| + 0.5 v_disk|v_disk| + "
                "0.7 v_bul|v_bul|"
            ),
            "channel_a": (
                "the same anomalous component fitted pointwise inside the "
                "total model on the fixed absolute deep subset"
            ),
            "channel_b": (
                "a0_B = (v_A^2(r_out))^2 / (G M_b) per galaxy, restricted to "
                "outermost points passing the same fixed deep cut"
            ),
            "derivation": (
                "v_A^2(r) = r sqrt(a_0 G M_b / r^2) = sqrt(G M_b a_0) in the "
                "deep exterior regime; inversion a_0 = (v_A^2)^2/(G M_b); "
                "asymptotic agreement with v^4 = G M_b a_0; dimension and "
                "sign checks in the module docstring"
            ),
        },
        "deep_cut": {
            "reference_a0_m_s2": REFERENCE_A0,
            "convention": (
                "g_bar < f * reference, fixed before fitting, identical for "
                "the channel-A subset and the channel-B outermost-point "
                "requirement; the reference is a comparison-only imported "
                "scale, not a fitted OPH output"
            ),
        },
        "primary_fraction_0p1": primary,
        "sensitivity_fraction_0p3": sensitivity,
        "verdict_rule": VERDICT_RULE,
        "diagnostic_verdict": primary_verdict,
        "superseded_mixed_proxy_diagnostic": {
            "receipt_path": "receipts/sparc_deep_regime_diagnostic.json",
            "receipt_sha256": hashlib.sha256(mixed_bytes).hexdigest(),
            "proxy_defect": (
                "catalogue Vflat is a finite-radius total speed used as an "
                "asymptotic anomalous-speed proxy without baryonic "
                "subtraction"
            ),
            "mixed_proxy_point_log10_ratio": float(mixed_point),
            "mixed_proxy_log10_ratio_95pct": [
                float(mixed_interval[0]),
                float(mixed_interval[1]),
            ],
            "matched_point_log10_ratio": float(matched_point),
            "displacement_removed_dex": float(mixed_point - matched_point),
            "resolution_statement": resolution_statement(
                primary_verdict, mixed_interval, matched_interval
            ),
        },
        "inference_boundary": {
            "status": (
                "labeled postdiction diagnostic on the seen committed "
                "snapshot; no new data, no arming, no discharge, no scored "
                "comparison"
            ),
            "not_a_preregistered_likelihood": (
                "the joint likelihood needs a measurement-error and nuisance "
                "covariance model, an asymptotic-velocity treatment, and a "
                "frozen contract before exposure to any new data; that work "
                "is open"
            ),
            "error_model_limits": (
                "only e_Vobs and inclination error through the observed term "
                "are propagated; distance, mass-to-light, "
                "stellar-population, correlated rotation-curve, and "
                "baryonic-model uncertainties are not"
            ),
            "enclosed_mass": (
                "the exterior point-mass limit at r_out is not separately "
                "enforced; enclosed mass below the total biases a0_B "
                "downward; the same-radius mass variant bounds this "
                "sensitivity"
            ),
            "bulge_mass_convention": (
                "table-1 masses apply the disk mass-to-light value to the "
                "total luminosity while the table-2 subtraction separates "
                "bulges; inherited from the committed comparison"
            ),
            "neutrality": (
                "consistency with zero is shared with the standard null and "
                "does not support OPH; a tension in either direction would "
                "count against the shared normalization only within this "
                "declared diagnostic"
            ),
            "verdict_scope": (
                "applies to this declared diagnostic, not to the recovered "
                "OPH core or the generic dark-sector theorem"
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
