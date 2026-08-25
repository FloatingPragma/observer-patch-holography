#!/usr/bin/env python3
"""Deep-regime comparison of the anomalous-density galaxy law with SPARC.

The dark-matter paper derives, conditional on its profile premises, the
anomalous acceleration and baryonic Tully--Fisher laws

    a_A = sqrt(a_b a_0),        v^4 = G M_b a_0.

Because the paper's Poisson equation contains both baryonic and anomalous
density, the quantity compared with the observed circular acceleration is

    g_model = g_bar + sqrt(g_bar a_0),

not the anomalous term alone.  The source does not fix ``a_0``.  This module
therefore performs a diagnostic postdiction on fixed absolute ``g_bar`` cuts
and reports

* the anomalous-term exponent and acceleration constant fitted on points with
  ``g_bar < f a_ref``, where the comparison-only reference
  ``a_ref = 1.2e-10 m/s^2`` is fixed before fitting;
* the baryonic Tully--Fisher exponent and the acceleration constant implied
  by its normalisation;
* a galaxy-cluster bootstrap comparison of the two constants, which the
  conditional density profile requires to coincide;
* the residual of the additive all-gradient extension
  ``g = g_bar + sqrt(g_bar a_0)`` against the full sample and its Solar-System
  anomaly, which the paper's high-gradient branch is required to remove.

The fixed cut avoids the discontinuous, and for this sample sometimes cycling,
``g_bar < f a_0`` fit-and-reselect procedure used by the superseded v1
diagnostic.  The acceleration constant is a fitted comparison value.  Nothing
here is a source-derived OPH output or a full likelihood, and the receipt says
so.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import optimize

SCHEMA = "oph.cosmology.rar_deep_regime_diagnostic.v2"
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
DEEP_CUT_REFERENCE_A0_M_S2 = 1.2e-10
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260824
LOG10_A0_BOUNDS = (-12.0, -9.0)


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


def _fit_total_a0(g_obs: np.ndarray, g_bar: np.ndarray) -> float:
    """Fit the fixed-exponent total law by unweighted log-residual RMS."""
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


def _fit_total_free_exponent(
    g_obs: np.ndarray,
    g_bar: np.ndarray,
    reference_a0: float = DEEP_CUT_REFERENCE_A0_M_S2,
) -> tuple[float, float]:
    """Fit ``g_bar + sqrt(g_bar*a0)*(g_bar/a_ref)^(alpha-1/2)``.

    The total-law likelihood uses every retained point. It avoids selecting on
    the noisy sign of ``g_obs-g_bar``, which would bias a log-excess fit.
    """
    def residual(params: np.ndarray) -> np.ndarray:
        log10_a0, alpha = params
        a0 = 10.0 ** log10_a0
        anomaly = np.sqrt(g_bar * a0) * (
            g_bar / reference_a0
        ) ** (alpha - 0.5)
        return np.log10(g_obs) - np.log10(g_bar + anomaly)

    result = optimize.least_squares(
        residual,
        x0=np.array([-10.0, 0.5]),
        bounds=(
            np.array([LOG10_A0_BOUNDS[0], -1.0]),
            np.array([LOG10_A0_BOUNDS[1], 2.0]),
        ),
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=1000,
    )
    if not result.success:
        raise RuntimeError(f"free-exponent minimization failed: {result.message}")
    return float(10.0 ** result.x[0]), float(result.x[1])


def _cluster_indices(names: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    galaxies = np.unique(names)
    sampled = rng.choice(galaxies, size=len(galaxies), replace=True)
    return np.concatenate([np.flatnonzero(names == galaxy) for galaxy in sampled])


def deep_fit(
    g_obs: np.ndarray,
    g_bar: np.ndarray,
    names: np.ndarray,
    fraction: float = 0.1,
    reference_a0: float = DEEP_CUT_REFERENCE_A0_M_S2,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    """Fit the total law on a fixed absolute deep cut and cluster-bootstrap it."""
    mask = g_bar < fraction * reference_a0
    selected_obs = g_obs[mask]
    selected_bar = g_bar[mask]
    selected_names = names[mask]
    a0 = _fit_total_a0(selected_obs, selected_bar)
    free_a0, slope = _fit_total_free_exponent(
        selected_obs, selected_bar, reference_a0
    )
    model = selected_bar + np.sqrt(selected_bar * a0)
    resid_fixed = np.log10(selected_obs) - np.log10(model)
    free_model = selected_bar + np.sqrt(selected_bar * free_a0) * (
        selected_bar / reference_a0
    ) ** (slope - 0.5)
    resid_free = np.log10(selected_obs) - np.log10(free_model)

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    boot_a0: list[float] = []
    boot_slopes: list[float] = []
    for _ in range(bootstrap_replicates):
        idx = _cluster_indices(selected_names, rng)
        sample_obs = selected_obs[idx]
        sample_bar = selected_bar[idx]
        boot_a0.append(_fit_total_a0(sample_obs, sample_bar))
        _, boot_slope = _fit_total_free_exponent(
            sample_obs, sample_bar, reference_a0
        )
        boot_slopes.append(boot_slope)

    a0_interval = np.percentile(boot_a0, [2.5, 97.5])
    slope_interval = np.percentile(boot_slopes, [2.5, 97.5])
    return {
        "deep_fraction_of_reference_a0": fraction,
        "deep_cut_reference_a0_m_s2": reference_a0,
        "deep_cut_g_bar_max_m_s2": fraction * reference_a0,
        "n_points": int(mask.sum()),
        "n_galaxies": int(len(np.unique(selected_names))),
        "a0_fixed_exponent_m_s2": float(a0),
        "optimizer_log10_a0_bounds": list(LOG10_A0_BOUNDS),
        "optimizer_solution_interior": bool(
            LOG10_A0_BOUNDS[0] < math.log10(a0) < LOG10_A0_BOUNDS[1]
        ),
        "a0_galaxy_cluster_bootstrap_95pct_m_s2": [
            float(a0_interval[0]), float(a0_interval[1])
        ],
        "rms_scatter_dex_fixed_exponent": float(np.sqrt(np.mean(resid_fixed**2))),
        "free_anomalous_exponent_in_total_model": float(slope),
        "free_model_a0_m_s2": float(free_a0),
        "free_exponent_galaxy_cluster_bootstrap_sd": float(np.std(boot_slopes)),
        "free_exponent_galaxy_cluster_bootstrap_95pct": [
            float(slope_interval[0]), float(slope_interval[1])
        ],
        "rms_scatter_dex_free_total_model": float(np.sqrt(np.mean(resid_free**2))),
        "law_exponent": 0.5,
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "point_weighting": "equal retained rotation-curve points",
        "bootstrap_unit": "galaxy cluster",
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
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    exps = []
    a0_boot = []
    for _ in range(BOOTSTRAP_REPLICATES):
        idx = rng.integers(0, len(x), len(x))
        s, _ = np.polyfit(x[idx], y[idx], 1)
        exps.append(1.0 / s)
        a0_boot.append(float(10 ** np.mean(np.log10(a0_points[idx]))))
    a0_interval = np.percentile(a0_boot, [2.5, 97.5])
    return {
        "n_galaxies": int(len(Mb)),
        "free_exponent_M_of_V": float(exponent),
        "free_exponent_bootstrap_sd": float(np.std(exps)),
        "a0_fixed_exponent_m_s2": a0,
        "a0_galaxy_bootstrap_95pct_m_s2": [
            float(a0_interval[0]), float(a0_interval[1])
        ],
        "rms_scatter_dex_in_a0": float(np.sqrt(np.mean(resid**2))),
        "law_exponent": 4.0,
        "upsilon_disk": UPSILON_DISK,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
    }


def paired_consistency_bootstrap(
    t1: list[dict[str, Any]],
    acc: dict[str, np.ndarray],
    fraction: float = 0.1,
    reference_a0: float = DEEP_CUT_REFERENCE_A0_M_S2,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    """Parent-galaxy bootstrap of the shared RAR/BTFR normalization.

    All catalogue parents are resampled before either analysis cut is applied.
    A sampled galaxy carries all retained rotation-curve points and, when
    eligible, its one BTFR datum. This preserves within-galaxy correlation and
    the overlap between the two diagnostics.
    """
    parents = np.array([row["name"] for row in t1])
    deep_mask = acc["g_bar"] < fraction * reference_a0
    deep_obs = acc["g_obs"][deep_mask]
    deep_bar = acc["g_bar"][deep_mask]
    deep_names = acc["names"][deep_mask]

    btfr_a0: dict[str, float] = {}
    for row in t1:
        if row["quality"] > 2 or row["inclination_deg"] < 30.0:
            continue
        if not (row["vflat_kms"] > 0):
            continue
        mstar = UPSILON_DISK * row["L36_gsun"] * 1e9
        mgas = GAS_HELIUM * row["MHI_gsun"] * 1e9
        mass = (mstar + mgas) * M_SUN
        velocity = row["vflat_kms"] * KM
        btfr_a0[row["name"]] = float(velocity**4 / (G_SI * mass))

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rar_values: list[float] = []
    btfr_values: list[float] = []
    log_ratios: list[float] = []
    for _ in range(bootstrap_replicates):
        sampled = rng.choice(parents, size=len(parents), replace=True)
        names, multiplicities = np.unique(sampled, return_counts=True)
        rar_indices: list[np.ndarray] = []
        sampled_btfr: list[float] = []
        for name, multiplicity in zip(names, multiplicities):
            idx = np.flatnonzero(deep_names == name)
            if len(idx):
                rar_indices.extend([idx] * int(multiplicity))
            value = btfr_a0.get(str(name))
            if value is not None:
                sampled_btfr.extend([value] * int(multiplicity))
        if not rar_indices or len(sampled_btfr) < 2:
            continue
        joined = np.concatenate(rar_indices)
        a0_rar = _fit_total_a0(deep_obs[joined], deep_bar[joined])
        a0_btfr = float(10 ** np.mean(np.log10(sampled_btfr)))
        rar_values.append(a0_rar)
        btfr_values.append(a0_btfr)
        log_ratios.append(math.log10(a0_btfr / a0_rar))

    if len(log_ratios) != bootstrap_replicates:
        raise RuntimeError("paired bootstrap produced an empty diagnostic sample")
    rar_interval = np.percentile(rar_values, [2.5, 97.5])
    btfr_interval = np.percentile(btfr_values, [2.5, 97.5])
    ratio_interval = np.percentile(log_ratios, [2.5, 97.5])
    nonpositive_count = int(np.count_nonzero(np.asarray(log_ratios) <= 0.0))
    return {
        "parent_catalogue_galaxies": int(len(parents)),
        "bootstrap_unit": "parent galaxy before RAR and BTFR cuts",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "rar_a0_95pct_m_s2": [float(rar_interval[0]), float(rar_interval[1])],
        "btfr_a0_95pct_m_s2": [float(btfr_interval[0]), float(btfr_interval[1])],
        "log10_ratio_btfr_over_rar_95pct": [
            float(ratio_interval[0]), float(ratio_interval[1])
        ],
        "bootstrap_nonpositive_log10_ratio_count": nonpositive_count,
        "bootstrap_nonpositive_log10_ratio_plus_one_fraction": float(
            (nonpositive_count + 1) / (bootstrap_replicates + 1)
        ),
        "diagnostic_verdict": (
            "TENSION_UNDER_DECLARED_DIAGNOSTIC"
            if ratio_interval[0] > 0.0
            else "INCONCLUSIVE_COMPATIBLE"
        ),
    }


def additive_extension(g_obs: np.ndarray, g_bar: np.ndarray, a0: float) -> dict[str, Any]:
    """The deep density profile extrapolated to all gradients."""
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
    deep = deep_fit(acc["g_obs"], acc["g_bar"], acc["names"])
    deep_sweep = [
        deep_fit(acc["g_obs"], acc["g_bar"], acc["names"], fraction=f)
        for f in (0.3, 0.1, 0.03)
    ]
    btfr = btfr_fit(t1)
    a0_deep = deep["a0_fixed_exponent_m_s2"]
    a0_btfr = btfr["a0_fixed_exponent_m_s2"]
    additive = additive_extension(acc["g_obs"], acc["g_bar"], a0_deep)
    paired = paired_consistency_bootstrap(t1, acc)
    consistency = {
        "a0_deep_rar_m_s2": a0_deep,
        "a0_btfr_m_s2": a0_btfr,
        "log10_ratio_btfr_over_deep": math.log10(a0_btfr / a0_deep),
        "paired_parent_galaxy_bootstrap": paired,
        "diagnostic_verdict": paired["diagnostic_verdict"],
        "raw_scatter_is_not_parameter_uncertainty": True,
    }
    data_sha = hashlib.sha256(
        (DATA / "table1.dat").read_bytes() + (DATA / "table2.dat").read_bytes()
    ).hexdigest()
    return {
        "schema": SCHEMA,
        "scope": "diagnostic_postdiction_fixed_absolute_cuts_one_fitted_constant",
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
            "deep_radial_acceleration": "a_A = sqrt(a_b a_0)",
            "observed_total_acceleration": "g_model = g_bar + sqrt(g_bar a_0)",
            "baryonic_tully_fisher": "v^4 = G M_b a_0",
            "parameter_free_content": "exponents 1/2 and 4; a_0 is a fitted comparison value",
        },
        "inference_boundary": {
            "deep_cut": (
                "fixed absolute cuts g_bar < f * 1.2e-10 m/s^2; the reference "
                "is a comparison-only imported scale, not a fitted OPH output"
            ),
            "likelihood": "unweighted log-residual diagnostic, not a measurement-error likelihood",
            "systematics": (
                "fixed mass-to-light ratios; distance, inclination, correlated "
                "rotation-curve, and stellar-population uncertainties are not marginalized"
            ),
            "btfr_proxy": (
                "Vflat is used as an asymptotic-speed proxy. The receipt does not "
                "subtract the finite-radius baryonic contribution, parse a flat "
                "radius, or establish that Vflat equals V_infinity; this can bias "
                "the naive Vflat^4/(G M_b) normalization upward"
            ),
            "verdict_scope": (
                "a tension applies only to this declared density-profile diagnostic, "
                "not to the recovered OPH core or the generic dark-sector theorem"
            ),
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
