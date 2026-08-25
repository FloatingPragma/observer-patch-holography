#!/usr/bin/env python3
"""Matched-observable dark-sector diagnostic on the committed SPARC snapshot.

The committed mixed-proxy diagnostic (``rar_deep_regime_diagnostic.py``,
receipt ``receipts/sparc_deep_regime_diagnostic.json``) compared a deep-regime
radial-acceleration fit against a Tully--Fisher normalization built from the
catalogue ``Vflat``.  Catalogue ``Vflat`` is a finite-radius TOTAL speed used
there as an asymptotic ANOMALOUS speed without baryonic subtraction, so the
two channels measured different observables.  This module builds both channels
from one observable, the anomalous squared circular speed

    v_A^2(r) = v_obs(r)^2 - v_bar(r)^2,
    v_bar^2  = v_gas|v_gas| + 0.5 v_disk|v_disk| + 0.7 v_bul|v_bul|,

with the committed declared mass-to-light choices (0.5 disk, 0.7 bulge) and
the committed helium factor 1.33 that the catalogue gas term contains.

Channel-B estimator derivation (independent re-derivation, recorded per the
sign/frame discipline).  The dark-matter paper's deep-regime relations are

    a_A = sqrt(a_b a_0)          (anomalous acceleration),
    v^2(r) = a(r) r              (circular motion, per additive component).

Outside the bulk of the baryonic mass, a_b(r) = G M_b / r^2, so

    v_A^2(r) = a_A(r) r = r sqrt(a_0 G M_b / r^2) = sqrt(G M_b a_0),

independent of radius.  Inverting,

    a_0 = (v_A^2)^2 / (G M_b).

Cross-check against the committed baryonic Tully--Fisher law: as r grows the
baryonic term v_bar^2 = G M_b / r vanishes, the total squared speed tends to
v_A^2, and (v_A^2)^2 = v_inf^4 = G M_b a_0, the committed normalization.
Dimension check: (m^2 s^-2)^2 / (m^3 kg^-1 s^-2 * kg) = m s^-2.  Sign check:
the anomalous term adds to the baryonic term, so the model requires
v_obs^2 - v_bar^2 = v_A^2 > 0; a nonpositive sampled value admits no estimate
and is excluded and counted.  Bias-direction check reproducing the audited
proxy defect: at finite radius v_obs^2(r) = v_bar^2(r) + v_A^2 >= v_A^2, so
the unsubtracted proxy Vflat^4/(G M_b) is biased upward relative to
(v_A^2)^2/(G M_b) by the factor (1 + v_bar^2/v_A^2)^2 at the measurement
radius.

Finite-radius caveat and its handling.  The identity v_A^2 = sqrt(G M_b a_0)
needs two regimes at the evaluation radius: the deep regime (a_b << a_0) and
the exterior-mass regime (enclosed baryonic mass close to the total M_b).
The outermost measured radius of a rotation curve may satisfy neither.  The
deep regime is enforced with the SAME fixed absolute cut as channel A:
channel B keeps only galaxies whose outermost retained point has
g_bar(r_out) < f * 1.2e-10 m/s^2, with the cut fixed before fitting, and the
receipt reports every sample size and exclusion count.  The exterior-mass
regime is not separately enforced; where enclosed mass falls below the total,
the estimator is biased downward by the enclosed fraction under the local
reading of the profile.  A same-radius mass variant, replacing G M_b by
v_bar^2(r_out) r_out, is reported as a labeled sensitivity row; it is exact
in the exterior point-mass limit and uses the table-2 mass-to-light split at
the same radius.

Channel A is the audited corrected radial-acceleration fit, reused verbatim
from the committed module: the total model g_model = g_bar + sqrt(g_bar a_0)
fitted by unweighted log residuals on the fixed absolute deep subset.

Statistics.  Measurement errors are propagated to v_A^2 from the snapshot
columns the data carry: the per-point velocity error e_Vobs (random error
from non-circular motions and asymmetries, per the catalogue ReadMe) and the
per-galaxy inclination error e_i applied through the observed term, with
sigma(v_A^2) = hypot(2 v_obs e_Vobs, 2 v_obs^2 cot(i) e_i).  Distance,
mass-to-light, stellar-population, and correlated rotation-curve
uncertainties, and errors in the baryonic model velocities, are not
propagated.  The primary channel-B combination is the unweighted log mean,
the same combination rule as the superseded mixed-proxy channel so that only
the observable changes between the two receipts; the propagated errors are
reported per galaxy in summary, and the inverse-variance weighted mean is a
labeled sensitivity row because the observed scatter is dominated by
unpropagated terms and measurement-error weights correlate with the
estimate.  The paired comparison is a galaxy-level bootstrap of
log10(a0_B / a0_A) on the common galaxy set, with both channels recomputed
per replicate, percentile intervals, and zero tallies reported as counts with
plus-one fractions (count + 1)/(replicates + 1), never as probability zero.

The decision rule is direction-neutral: an interval entirely above zero and
an interval entirely below zero are both labeled tensions, symmetrically; an
interval containing zero is labeled consistent and is compatible with every
model that reproduces both relations, including the standard null, so it is
not evidence for OPH.

What is not proved here.  This is a labeled postdiction diagnostic on the
seen committed snapshot with fixed declared mass-to-light ratios and no full
covariance model.  It supersedes the mixed-proxy channel-B observable; it is
not the preregistered joint likelihood, which needs a measurement-error and
nuisance covariance model, an asymptotic-velocity treatment, and a frozen
contract before exposure to any new data.  Nothing here derives a_0, arms or
discharges a frozen prediction, scores a comparison, or bears on the generic
modular-charge dark-sector theorem. Table 1 gives total L_[3.6], while table 2
separates disk and bulge rotation components. Galaxies with any nonzero
retained V_bulge are therefore excluded from channel B: the snapshot does not
provide the disk/bulge luminosity split needed to combine a 0.7 bulge
subtraction with a consistent baryonic-mass denominator.
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
sys.path.insert(0, str(HERE))

import rar_deep_regime_diagnostic as base  # noqa: E402

SCHEMA = "oph.cosmology.matched_observable_dark_sector_diagnostic.v1"
RECEIPT_PATH = HERE / "runtime" / "matched_observable_receipt.json"
MIXED_PROXY_RECEIPT = HERE / "receipts" / "sparc_deep_regime_diagnostic.json"

PRIMARY_FRACTION = 0.1
SENSITIVITY_FRACTION = 0.3
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260825  # declared for this diagnostic


def anomalous_squared_speed_m2s2(
    vobs_kms: float, vgas_kms: float, vdisk_kms: float, vbul_kms: float
) -> tuple[float, float]:
    """Return (v_A^2, v_bar^2) in m^2/s^2 at one rotation-curve point.

    v_bar^2 uses the committed declared mass-to-light choices, 0.5 for the
    disk and 0.7 for the bulge, on the unit-mass-to-light catalogue
    components, with signed gas handling identical to the committed module.
    """
    vo = vobs_kms * base.KM
    vg = vgas_kms * base.KM
    vd = vdisk_kms * base.KM
    vb = vbul_kms * base.KM
    vbar2 = (
        vg * abs(vg)
        + base.UPSILON_DISK * vd * abs(vd)
        + base.UPSILON_BULGE * vb * abs(vb)
    )
    return vo * vo - vbar2, vbar2


def outermost_retained_points(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    max_quality: int = 2,
    min_inclination: float = 30.0,
    max_rel_error: float = 0.10,
) -> dict[str, dict[str, Any]]:
    """Outermost retained rotation-curve point per eligible galaxy.

    Point-level and galaxy-level cuts are the committed standard cuts, so the
    same points govern both channels.  The outermost point is taken within
    the retained sample.
    """
    meta = {r["name"]: r for r in t1}
    rows: dict[str, dict[str, Any]] = {}
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > max_quality or m["inclination_deg"] < min_inclination:
            continue
        if r["vobs_kms"] <= 0 or r["e_vobs_kms"] / r["vobs_kms"] > max_rel_error:
            continue
        v_a2, vbar2 = anomalous_squared_speed_m2s2(
            r["vobs_kms"], r["vgas_kms"], r["vdisk_kms"], r["vbul_kms"]
        )
        radius_m = r["rad_kpc"] * base.KPC_M
        g_bar = vbar2 / radius_m
        if g_bar <= 0:
            continue
        cur = rows.get(r["name"])
        has_nonzero_bulge = r["vbul_kms"] != 0.0 or bool(
            cur and cur["has_nonzero_bulge"]
        )
        if cur is None or r["rad_kpc"] > cur["rad_kpc"]:
            rows[r["name"]] = dict(
                r,
                v_a2=v_a2,
                vbar2=vbar2,
                g_bar=g_bar,
                radius_m=radius_m,
                has_nonzero_bulge=has_nonzero_bulge,
            )
        else:
            cur["has_nonzero_bulge"] = has_nonzero_bulge
    return rows


def channel_b(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    fraction: float,
    reference_a0: float = base.DEEP_CUT_REFERENCE_A0_M_S2,
) -> dict[str, Any]:
    """Per-galaxy matched flat-speed estimates a0_B = (v_A^2)^2 / (G M_b).

    Restricted to galaxies whose outermost retained point passes the same
    fixed absolute deep cut as channel A.  Also returns the same-radius mass
    variant replacing G M_b by v_bar^2(r_out) r_out.
    """
    meta = {r["name"]: r for r in t1}
    rows = outermost_retained_points(t1, t2)
    names: list[str] = []
    log_a0: list[float] = []
    sigma_log: list[float] = []
    log_a0_same_radius: list[float] = []
    n_not_deep = 0
    n_nonpositive_va2 = 0
    n_nonpositive_mb = 0
    excluded_bulge_luminosity_ambiguous: list[str] = []
    for name, r in rows.items():
        if not (r["g_bar"] < fraction * reference_a0):
            n_not_deep += 1
            continue
        if r["has_nonzero_bulge"]:
            excluded_bulge_luminosity_ambiguous.append(name)
            continue
        if r["v_a2"] <= 0:
            n_nonpositive_va2 += 1
            continue
        m = meta[name]
        mstar = base.UPSILON_DISK * m["L36_gsun"] * 1e9
        mgas = base.GAS_HELIUM * m["MHI_gsun"] * 1e9
        mass = (mstar + mgas) * base.M_SUN
        if mass <= 0:
            n_nonpositive_mb += 1
            continue
        a0_b = r["v_a2"] * r["v_a2"] / (base.G_SI * mass)
        vo = r["vobs_kms"] * base.KM
        ev = r["e_vobs_kms"] * base.KM
        inc_rad = math.radians(m["inclination_deg"])
        e_inc_rad = math.radians(m["e_inclination_deg"])
        sigma_va2 = math.hypot(2.0 * vo * ev, 2.0 * vo * vo * e_inc_rad / math.tan(inc_rad))
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
            "n_excluded_bulge_luminosity_ambiguous": len(
                excluded_bulge_luminosity_ambiguous
            ),
            "excluded_bulge_luminosity_ambiguous": sorted(
                excluded_bulge_luminosity_ambiguous
            ),
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
        "sigma_log": sigma_arr,
        "summary": {
            "estimator": "a0_B = (v_obs(r_out)^2 - v_bar(r_out)^2)^2 / (G M_b)",
            "baryonic_mass": (
                "M_b = 0.5 L_[3.6] + 1.33 M_HI from table 1, committed "
                "declared mass-to-light and helium choices"
            ),
            "deep_outermost_requirement_g_bar_max_m_s2": fraction * reference_a0,
            "n_candidate_galaxies": int(len(rows)),
            "n_outermost_not_deep": n_not_deep,
            "n_nonpositive_anomalous_speed2": n_nonpositive_va2,
            "n_nonpositive_baryonic_mass": n_nonpositive_mb,
            "n_excluded_bulge_luminosity_ambiguous": len(
                excluded_bulge_luminosity_ambiguous
            ),
            "excluded_bulge_luminosity_ambiguous": sorted(
                excluded_bulge_luminosity_ambiguous
            ),
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


def _verdict(interval_low: float, interval_high: float) -> str:
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
    acc: dict[str, np.ndarray],
    b_channel: dict[str, Any],
    fraction: float,
    reference_a0: float = base.DEEP_CUT_REFERENCE_A0_M_S2,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> dict[str, Any]:
    """Galaxy-level paired bootstrap of log10(a0_B/a0_A) on the common set.

    Each replicate resamples common-set galaxies with replacement; a sampled
    galaxy carries all its deep channel-A points and its one channel-B datum,
    so within-galaxy correlation and the channel overlap are preserved.
    """
    mask = acc["g_bar"] < fraction * reference_a0
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

    a0_a_common = base._fit_total_a0(d_obs, d_bar)
    a0_b_common = float(10.0 ** np.mean(b_log))
    point_ratio = math.log10(a0_b_common / a0_a_common)

    rng = np.random.default_rng(seed)
    ratios: list[float] = []
    for _ in range(bootstrap_replicates):
        sampled = rng.choice(common, size=len(common), replace=True)
        point_indices: list[np.ndarray] = []
        sampled_logs: list[float] = []
        for name in sampled:
            point_indices.append(np.flatnonzero(d_names == name))
            sampled_logs.append(float(b_log[np.flatnonzero(b_names == name)[0]]))
        joined = np.concatenate(point_indices)
        a0_a = base._fit_total_a0(d_obs[joined], d_bar[joined])
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
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": seed,
        "log10_ratio_95pct": [float(interval[0]), float(interval[1])],
        "count_log10_ratio_le_zero": n_le_zero,
        "le_zero_plus_one_fraction": float(
            (n_le_zero + 1) / (bootstrap_replicates + 1)
        ),
        "count_log10_ratio_ge_zero": n_ge_zero,
        "ge_zero_plus_one_fraction": float(
            (n_ge_zero + 1) / (bootstrap_replicates + 1)
        ),
        "verdict": _verdict(float(interval[0]), float(interval[1])),
    }


def fraction_block(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    acc: dict[str, np.ndarray],
    fraction: float,
) -> dict[str, Any]:
    reference_a0 = base.DEEP_CUT_REFERENCE_A0_M_S2
    mask = acc["g_bar"] < fraction * reference_a0
    a0_a = base._fit_total_a0(acc["g_obs"][mask], acc["g_bar"][mask])
    b_res = channel_b(t1, t2, fraction)
    paired = paired_common_bootstrap(acc, b_res, fraction)
    return {
        "fraction_of_reference": fraction,
        "deep_cut_g_bar_max_m_s2": fraction * reference_a0,
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
            "The combined matched-analysis changes remove the displacement "
            "and are consistent with a mixed-proxy explanation. Because the "
            "observable, radius choice, eligible sample, and combination rule "
            "change together, this comparison does not causally isolate which "
            "change is responsible. This "
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


def run() -> dict[str, Any]:
    t1 = base.read_table1()
    t2 = base.read_table2()
    acc = base.accelerations(t1, t2)
    primary = fraction_block(t1, t2, acc, PRIMARY_FRACTION)
    sensitivity = fraction_block(t1, t2, acc, SENSITIVITY_FRACTION)
    data_sha = hashlib.sha256(
        (base.DATA / "table1.dat").read_bytes()
        + (base.DATA / "table2.dat").read_bytes()
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
            "upsilon_disk": base.UPSILON_DISK,
            "upsilon_bulge": base.UPSILON_BULGE,
            "gas_helium_factor": base.GAS_HELIUM,
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
            "reference_a0_m_s2": base.DEEP_CUT_REFERENCE_A0_M_S2,
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
                "galaxies with any nonzero retained V_bulge are excluded from "
                "channel B because table 1 gives total L_[3.6] and the "
                "snapshot supplies no disk/bulge luminosity split"
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


def canonical_json(obj: dict[str, Any]) -> str:
    return json.dumps(obj, indent=2, sort_keys=True) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output", type=Path, default=RECEIPT_PATH)
    args = ap.parse_args()
    out = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(canonical_json(out), encoding="utf-8")
    paired = out["primary_fraction_0p1"]["paired_common_set"]
    print(f"receipt written: {args.output}")
    print(
        "primary paired log10(a0_B/a0_A) interval:",
        paired["log10_ratio_95pct"],
        paired["verdict"],
    )


if __name__ == "__main__":
    main()
