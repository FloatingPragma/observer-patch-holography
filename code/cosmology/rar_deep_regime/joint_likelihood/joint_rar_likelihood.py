#!/usr/bin/env python3
"""Joint SPARC rotation-curve and BTFR profile likelihood, declared conventions.

The word "joint" names the shared-nuisance, jointly-resampled treatment of
the two channels; the likelihood object is the per-point rotation-curve
likelihood, and no combined two-channel likelihood is formed (the BTFR
channel enters as a matched estimator compared through the paired
bootstrap).

This module replaces the mixed-proxy diagnostic combination rule with a full
per-point rotation-curve likelihood for the deep-regime candidate model

    g_model = g_bar + sqrt(g_bar a0),

equivalently, per rotation-curve point at catalogue distance,

    v_model^2(r) = v_bar^2(r) + sqrt(v_bar^2(r) a0 r),

with per-galaxy nuisance profiling (disk mass-to-light, distance,
inclination), a declared intra-galaxy equal-correlation covariance family
scanned over a grid, profile-likelihood intervals for a0 on the full retained
sample and on fixed absolute deep subsets, and a matched
baryon-subtracted BTFR channel under the same nuisance values with a
direction-neutral paired consistency statement.

Sign, frame, and unit conventions (each re-derived here, per the audit
discipline):

* Baryonic squared speed.  v_bar^2 = v_gas|v_gas| + Upsilon_d v_disk|v_disk|
  + 0.7 v_bul|v_bul| in m^2/s^2, the committed signed-gas convention: a
  negative catalogue V_gas denotes a gas contribution that subtracts from the
  centripetal balance, and the catalogue disk and bulge components are given
  at unit mass-to-light.  The helium factor 1.33 is inside the catalogue
  V_gas and multiplies M_HI in the table-1 baryonic mass.
* Distance scaling.  Angular positions are the measured frame, so the
  physical radius scales as r = d r_cat under a distance factor
  d = D / D_cat.  Photometric masses scale as M proportional to L
  proportional to D^2, so component squared speeds scale as
  v^2 = G M / r proportional to D^2 / D = D, giving
  v_bar^2 = d v_bar_cat^2 at r = d r_cat.  The anomalous term scales as
  v_A^2 = sqrt(v_bar^2 a0 r) proportional to sqrt(d * d) = d, so the whole
  model obeys v_model^2(d) = d v_model_cat^2 while the observed velocity is
  distance independent.  The table-1 baryonic mass scales as d^2.
* Inclination scaling.  The catalogue V_obs is the line-of-sight velocity
  corrected by 1/sin(i_cat).  Re-correcting to a trial inclination i gives
  v_obs(i) = v_obs_cat sin(i_cat)/sin(i), and the same factor multiplies the
  catalogue random error e_Vobs because that error is quoted on the corrected
  velocity.  The differential form dv/v = -cot(i) di reproduces the
  committed cot(i) error factor of the matched diagnostic.
* Dimension check of the anomalous term: v_A^2 = sqrt(v_bar^2 a0 r) has
  units sqrt(m^2 s^-2 * m s^-2 * m) = m^2 s^-2.  Sign check: a0 > 0 adds a
  positive anomalous term, so v_model^2 > v_bar^2 for every retained point.

Likelihood.  Per galaxy, the residual vector u_k = (v_obs_k(i) -
v_model_k(a0, U, d)) / (e_Vobs_k sin(i_cat)/sin(i)) enters the
equal-correlation Gaussian block

    -2 ln L = [sum u_k^2 - w (sum u_k)^2] / (1 - rho) + 2 n ln s(i) + const,
    w = rho / (1 + (n-1) rho),   s(i) = sin(i_cat)/sin(i),

the closed form of u^T R^-1 u for R = (1-rho) I + rho J, with the
i-dependent variance normalization retained and a0-independent constants
dropped.  The correlation rho is a declared family parameter scanned over a
fixed grid; a full covariance calibration for SPARC rotation curves is open,
so the a0 interval is reported separately at every declared rho value and no
single rho is selected.

Nuisances (profiled on declared grids, Gaussian penalties added to
-2 ln L):

* disk mass-to-light Upsilon_d, log-normal prior centered at the committed
  0.5 with declared width 0.1 dex, grid of 9 points over +-2.5 prior sigma;
* distance factor d, Gaussian prior at 1 with the catalogue relative error
  e_Dist/Dist, grid of 7 points over +-2.5 sigma, clipped below at 0.05;
* inclination i, Gaussian prior at the catalogue value with the catalogue
  error e_i, grid of 7 points over +-2.5 sigma, clipped to [10, 89.9] deg;
* bulge mass-to-light fixed at the committed 0.7 (declared fallback, the
  snapshot carries no per-galaxy bulge datum);
* when a catalogue error datum is absent or nonpositive the corresponding
  nuisance is fixed at the catalogue value and the receipt counts the
  fallback (zero galaxies on this snapshot).

Profiling and calibration.  For each a0 grid value the penalized -2 ln L is
minimized over the per-galaxy nuisance grid; the summed profiled curve gives
the maximum-likelihood a0 at the grid argmin (a three-point parabola vertex
is a labeled refinement row) and the 68.3 and 95 percent intervals from the
declared asymptotic chi-square one-parameter thresholds 1.0 and 3.84 on the
delta of the profiled curve measured from its grid minimum.  The
calibration is the Wilks asymptotic rule and is declared approximate under
grid profiling with Gaussian penalties; intervals narrower than one grid
step are flagged as grid-resolution limited; the synthetic-coverage test
exercises the construction on a known-a0 fixture.
Deep subsets are fixed before fitting by g_bar_cat < f * 1.2e-10 m/s^2 at
the catalogue nuisance values; nuisance shifts do not re-select points.

Matched BTFR channel.  Per deep-outermost galaxy the baryon-subtracted
observable v_A^2 = (v_obs(r_out) s(i*))^2 - d* v_bar_cat^2(r_out, U*) with
the profiled nuisance values (U*, d*, i*) taken at the subset
maximum-likelihood a0 for the same rho, and the estimator
a0_B = (v_A^2)^2 / (G M_b d*^2) with M_b = (U* L_[3.6] + 1.33 M_HI) at
catalogue distance.  The paired comparison is a galaxy-level bootstrap of
log10(a0_B / a0_A) on the common set with the channel-A profile re-minimized
per replicate; channel-B nuisance values are held at their point estimates
per replicate, a declared simplification.  Zero tallies are reported as
counts with plus-one fractions.  The decision rule is direction-neutral:
intervals entirely above or entirely below zero are labeled tensions
symmetrically, and an interval containing zero is labeled consistent, is
shared with the standard null, and is not evidence for OPH.

What this module does not do.  It does not derive a source value of a0 (the
OPH source does not fix a0, so a0 stays a fitted comparison parameter), does
not calibrate the intra-galaxy covariance, does not treat relativistic
lensing or cosmological observables, does not arm or discharge any frozen
prediction, and does not constitute a preregistered contract; it is a
labeled postdiction likelihood on the seen committed snapshot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "oph.cosmology.joint_rar_btfr_profile_likelihood.v1"
HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RECEIPT_PATH = HERE / "runtime" / "joint_likelihood_receipt.json"

KPC_M = 3.0856775814913673e19
KM = 1.0e3
G_SI = 6.67430e-11
M_SUN = 1.98841e30
UPSILON_DISK_CENTER = 0.5
UPSILON_BULGE = 0.7
GAS_HELIUM = 1.33
REFERENCE_A0_M_S2 = 1.2e-10
ML_PRIOR_WIDTH_DEX = 0.1
DIST_FACTOR_MIN = 0.05
INCLINATION_CLIP_DEG = (10.0, 89.9)
THRESHOLD_68 = 1.0
THRESHOLD_95 = 3.84
SEED = 20260825


@dataclass(frozen=True)
class Config:
    """Declared grid and scan resolutions; every value is part of the record."""

    log10_a0_min: float = -11.0
    log10_a0_max: float = -9.5
    n_a0: int = 301
    n_upsilon: int = 9
    n_distance: int = 7
    n_inclination: int = 7
    halfwidth_sigma: float = 2.5
    rho_grid: tuple[float, ...] = (0.0, 0.2, 0.4, 0.6)
    fractions: tuple[float, ...] = (0.3, 0.1)
    bootstrap_replicates: int = 1000
    seed: int = SEED
    max_quality: int = 2
    min_inclination_deg: float = 30.0
    max_rel_vobs_error: float = 0.10

    def log10_a0_grid(self) -> np.ndarray:
        return np.linspace(self.log10_a0_min, self.log10_a0_max, self.n_a0)


DEFAULT_CONFIG = Config()


def _f(s: str) -> float:
    s = s.strip()
    return float(s) if s else math.nan


def read_table1(path: Path = DATA / "table1.dat") -> list[dict[str, Any]]:
    """Parse table1 by the CDS ReadMe byte ranges, including error columns."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "dist_mpc": _f(line[15:21]),
                "e_dist_mpc": _f(line[22:27]),
                "inclination_deg": _f(line[30:34]),
                "e_inclination_deg": _f(line[35:39]),
                "L36_gsun": _f(line[40:47]),
                "MHI_gsun": _f(line[86:93]),
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


def build_galaxies(
    t1: list[dict[str, Any]],
    t2: list[dict[str, Any]],
    config: Config = DEFAULT_CONFIG,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Retained per-galaxy point arrays under the committed standard cuts.

    Cuts, identical to the committed diagnostics: quality <= 2, catalogue
    inclination >= 30 deg, v_obs > 0, e_Vobs/v_obs <= 0.10, g_bar_cat > 0 at
    the committed mass-to-light values.  Points with e_Vobs <= 0 admit no
    Gaussian residual, are excluded, and are counted.  Galaxies are kept in
    catalogue order.
    """
    meta = {r["name"]: r for r in t1}
    by_name: dict[str, dict[str, list[float]]] = {}
    n_nonpositive_e_vobs = 0
    for r in t2:
        m = meta[r["name"]]
        if m["quality"] > config.max_quality:
            continue
        if m["inclination_deg"] < config.min_inclination_deg:
            continue
        if r["vobs_kms"] <= 0:
            continue
        if r["e_vobs_kms"] <= 0:
            n_nonpositive_e_vobs += 1
            continue
        if r["e_vobs_kms"] / r["vobs_kms"] > config.max_rel_vobs_error:
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
        rec = by_name.setdefault(
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
        rec = by_name.get(r["name"])
        if rec is None:
            continue
        gal: dict[str, Any] = {"name": r["name"], "meta": r}
        for key, values in rec.items():
            gal[key] = np.array(values)
        gal["idx_out"] = int(np.argmax(gal["rad_kpc"]))
        galaxies.append(gal)
    return galaxies, {"n_points_excluded_nonpositive_e_vobs": n_nonpositive_e_vobs}


def nuisance_grids(meta: dict[str, Any], config: Config) -> dict[str, Any]:
    """Per-galaxy declared nuisance grids, penalties, and fallback flags."""
    offsets_u = np.linspace(
        -config.halfwidth_sigma, config.halfwidth_sigma, config.n_upsilon
    )
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
        offsets_d = np.linspace(
            -config.halfwidth_sigma, config.halfwidth_sigma, config.n_distance
        )
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
        offsets_i = np.linspace(
            -config.halfwidth_sigma, config.halfwidth_sigma, config.n_inclination
        )
        i_grid = np.clip(
            i_cat + e_i * offsets_i, INCLINATION_CLIP_DEG[0], INCLINATION_CLIP_DEG[1]
        )
        pen_i = ((i_grid - i_cat) / e_i) ** 2
        inclination_fixed = False

    penalty = (
        pen_u[:, None, None] + pen_d[None, :, None] + pen_i[None, None, :]
    )
    return {
        "upsilon": upsilon,
        "d_grid": d_grid,
        "i_grid": i_grid,
        "penalty": penalty,
        "distance_fixed": distance_fixed,
        "inclination_fixed": inclination_fixed,
    }


def galaxy_block_stats(
    gal: dict[str, Any],
    mask: np.ndarray,
    log10_a0_grid: np.ndarray,
    grids: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """S1, S2, and variance-normalization arrays on the nuisance grid.

    Returns S1[a,u,j,l] = sum_k u_k^2, S2[a,u,j,l] = (sum_k u_k)^2, and
    varlog[l] = 2 n ln s(i_l), with combos whose baryonic squared speed is
    nonpositive at any point set to S1 = +inf so the profile skips them.
    """
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

    n_a0 = len(a0)
    n_u = len(upsilon)
    n_d = len(d_grid)
    n_i = len(i_grid)
    S1 = np.empty((n_a0, n_u, n_d, n_i))
    S2 = np.empty((n_a0, n_u, n_d, n_i))
    sin_cat = math.sin(math.radians(gal["meta"]["inclination_deg"]))
    varlog = np.empty(n_i)
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


def equal_correlation_chi2(u: np.ndarray, rho: float) -> float:
    """Closed form of u^T R^-1 u for R = (1-rho) I + rho J, with validation."""
    if not (0.0 <= rho < 1.0):
        raise ValueError(f"correlation parameter out of range [0, 1): {rho}")
    n = len(u)
    if 1.0 + (n - 1) * rho <= 0.0:
        raise ValueError(f"equal-correlation matrix not positive definite: {rho}")
    w = rho / (1.0 + (n - 1) * rho)
    s1 = float(np.sum(u * u))
    s2 = float(np.sum(u)) ** 2
    return (s1 - w * s2) / (1.0 - rho)


def profiled_curve(
    S1: np.ndarray,
    S2: np.ndarray,
    varlog: np.ndarray,
    penalty: np.ndarray,
    npts: int,
    rho: float,
) -> dict[str, np.ndarray]:
    """Per-a0 profile over the nuisance grid at one correlation value."""
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


def parabola_vertex(x: np.ndarray, y: np.ndarray, m: int) -> float:
    """Three-point parabola vertex around index m, labeled refinement only."""
    if 0 < m < len(x) - 1:
        y0, y1, y2 = y[m - 1], y[m], y[m + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom > 0.0:
            step = x[1] - x[0]
            return float(x[m] + 0.5 * step * (y0 - y2) / denom)
    return float(x[m])


def interval_from_curve(
    log10_a0_grid: np.ndarray, curve: np.ndarray
) -> dict[str, Any]:
    """Profile interval by threshold crossing on the delta curve.

    The maximum-likelihood point is the grid argmin, so it lies inside its
    own intervals by construction; the parabola vertex is a labeled
    refinement row.  Delta is measured from the grid minimum.  Endpoints are
    the first linear-in-grid crossings walking outward from the minimum; an
    endpoint with no crossing inside the grid is pinned at the grid boundary
    and flagged; an interval narrower than one grid step is flagged as
    grid-resolution limited.
    """
    m = int(np.argmin(curve))
    x_ml = float(log10_a0_grid[m])
    step = float(log10_a0_grid[1] - log10_a0_grid[0])
    delta = curve - curve[m]

    def crossing(direction: int, threshold: float) -> tuple[float, bool]:
        j = m
        while 0 <= j + direction < len(delta):
            if delta[j + direction] >= threshold:
                x0 = log10_a0_grid[j]
                x1 = log10_a0_grid[j + direction]
                d0 = delta[j]
                d1 = delta[j + direction]
                return float(x0 + (x1 - x0) * (threshold - d0) / (d1 - d0)), False
            j += direction
        return float(log10_a0_grid[0 if direction < 0 else -1]), True

    out: dict[str, Any] = {
        "log10_a0_ml": x_ml,
        "a0_ml_m_s2": 10.0 ** x_ml,
        "a0_ml_parabola_refined_m_s2": 10.0
        ** parabola_vertex(log10_a0_grid, curve, m),
        "ml_grid_index": m,
        "min_minus2lnl": float(curve[m]),
    }
    for label, threshold in (("68", THRESHOLD_68), ("95", THRESHOLD_95)):
        lo, lo_pinned = crossing(-1, threshold)
        hi, hi_pinned = crossing(+1, threshold)
        out[f"a0_interval{label}_m_s2"] = [10.0 ** lo, 10.0 ** hi]
        out[f"interval{label}_log10"] = [lo, hi]
        out[f"interval{label}_pinned_at_grid_boundary"] = [lo_pinned, hi_pinned]
        out[f"interval{label}_grid_resolution_limited"] = bool(hi - lo < step)
    return out


def combo_shape(grids: dict[str, Any]) -> tuple[int, int, int]:
    return (len(grids["upsilon"]), len(grids["d_grid"]), len(grids["i_grid"]))


def channel_b_estimate(
    gal: dict[str, Any], grids: dict[str, Any], flat_index: int
) -> dict[str, Any]:
    """Baryon-subtracted outermost estimator at the profiled nuisance values.

    a0_B = (v_A^2)^2 / (G M_b d^2) with
    v_A^2 = (v_obs(r_out) s(i))^2 - d v_bar_cat^2(r_out, U) and
    M_b = (U L_[3.6] + 1.33 M_HI) 1e9 M_sun at catalogue distance; the d^2
    factor is the photometric mass scaling re-derived in the module
    docstring.
    """
    u_idx, j_idx, l_idx = np.unravel_index(flat_index, combo_shape(grids))
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
        "upsilon": upsilon,
        "distance_factor": d,
        "inclination_deg": i_val,
    }


def _verdict(interval_low: float, interval_high: float) -> str:
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


def paired_bootstrap(
    log10_a0_grid: np.ndarray,
    common_curves: np.ndarray,
    log_a0_b: np.ndarray,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    """Galaxy-level paired bootstrap on the common set.

    Each replicate resamples common galaxies with replacement via
    rng.integers(0, n, n); channel A re-minimizes the summed profiled curves
    of the sampled galaxies at grid resolution (per-galaxy nuisance
    profiling is inside the curves), channel B takes the unweighted log mean
    of the sampled fixed per-galaxy estimates.
    """
    n = common_curves.shape[0]
    if n < 2:
        raise ValueError("paired bootstrap needs at least two common galaxies")
    total = common_curves.sum(axis=0)
    x_ml = float(log10_a0_grid[int(np.argmin(total))])
    a0_a = 10.0 ** x_ml
    a0_b = float(10.0 ** np.mean(log_a0_b))
    point_ratio = math.log10(a0_b / a0_a)
    rng = np.random.default_rng(seed)
    ratios = np.empty(replicates)
    for t in range(replicates):
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
        "bootstrap_replicates": replicates,
        "bootstrap_seed": seed,
        "log10_ratio_95pct": [float(interval[0]), float(interval[1])],
        "count_log10_ratio_le_zero": n_le,
        "le_zero_plus_one_fraction": float((n_le + 1) / (replicates + 1)),
        "count_log10_ratio_ge_zero": n_ge,
        "ge_zero_plus_one_fraction": float((n_ge + 1) / (replicates + 1)),
        "verdict": _verdict(float(interval[0]), float(interval[1])),
    }


def run(
    t1: list[dict[str, Any]] | None = None,
    t2: list[dict[str, Any]] | None = None,
    config: Config = DEFAULT_CONFIG,
    include_provenance: bool = True,
) -> dict[str, Any]:
    if t1 is None:
        t1 = read_table1()
    if t2 is None:
        t2 = read_table2()
    galaxies, point_counts = build_galaxies(t1, t2, config)
    log10_a0_grid = config.log10_a0_grid()
    grids_by_name = {
        gal["name"]: nuisance_grids(gal["meta"], config) for gal in galaxies
    }
    n_fallback_distance = sum(
        1 for g in grids_by_name.values() if g["distance_fixed"]
    )
    n_fallback_inclination = sum(
        1 for g in grids_by_name.values() if g["inclination_fixed"]
    )

    subsets: list[tuple[str, float | None]] = [("full", None)]
    for f in config.fractions:
        subsets.append((f"deep_f_{str(f).replace('.', 'p')}", f))

    # Per (subset, galaxy): profiled curves for every rho, plus bookkeeping.
    subset_records: dict[str, dict[str, Any]] = {}
    for subset_name, fraction in subsets:
        per_gal: list[dict[str, Any]] = []
        n_points = 0
        for gal in galaxies:
            if fraction is None:
                mask = np.ones(len(gal["rad_kpc"]), dtype=bool)
            else:
                mask = gal["g_bar_cat"] < fraction * REFERENCE_A0_M_S2
            npts = int(mask.sum())
            if npts == 0:
                continue
            grids = grids_by_name[gal["name"]]
            S1, S2, varlog, _ = galaxy_block_stats(gal, mask, log10_a0_grid, grids)
            rec: dict[str, Any] = {"name": gal["name"], "npts": npts, "rho": {}}
            for rho in config.rho_grid:
                rec["rho"][rho] = profiled_curve(
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

    n_free_nuisances_by_subset = {}
    for subset_name, rec in subset_records.items():
        total = 0
        for g in rec["per_gal"]:
            grids = grids_by_name[g["name"]]
            total += 1  # disk mass-to-light is free for every galaxy
            total += 0 if grids["distance_fixed"] else 1
            total += 0 if grids["inclination_fixed"] else 1
        n_free_nuisances_by_subset[subset_name] = total

    gal_by_name = {gal["name"]: gal for gal in galaxies}
    subset_results: dict[str, Any] = {}
    for subset_name, rec in subset_records.items():
        fraction = rec["fraction"]
        per_gal = rec["per_gal"]
        per_rho_rows: list[dict[str, Any]] = []
        for rho in config.rho_grid:
            curves = np.stack([g["rho"][rho]["curve"] for g in per_gal])
            chi2_data = np.stack([g["rho"][rho]["chi2_data"] for g in per_gal])
            total_curve = curves.sum(axis=0)
            interval = interval_from_curve(log10_a0_grid, total_curve)
            m = interval["ml_grid_index"]
            chi2_data_ml = float(chi2_data[:, m].sum())
            n_edge = 0
            for g in per_gal:
                grids = grids_by_name[g["name"]]
                shape = combo_shape(grids)
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
            dof_lower = (
                rec["n_points"] - 1 - n_free_nuisances_by_subset[subset_name]
            )
            row: dict[str, Any] = {
                "rho": rho,
                **interval,
                "chi2_data_at_ml": chi2_data_ml,
                "dof_no_nuisance_count": dof_upper,
                "dof_full_nuisance_count": dof_lower,
                "reduced_chi2_dof_no_nuisance": chi2_data_ml / dof_upper,
                "reduced_chi2_dof_full_nuisance": (
                    chi2_data_ml / dof_lower if dof_lower > 0 else None
                ),
                "dof_convention": (
                    "bounds on the effective count: penalized profiled "
                    "nuisances sit between fully free and fully fixed"
                ),
                "n_galaxies_profiled_nuisance_at_grid_edge_at_ml": n_edge,
            }

            if fraction is not None:
                row["paired_btfr"] = _paired_btfr_block(
                    per_gal,
                    gal_by_name,
                    grids_by_name,
                    curves,
                    log10_a0_grid,
                    m,
                    rho,
                    fraction,
                    config,
                )
            per_rho_rows.append(row)
        subset_results[subset_name] = {
            "fraction_of_reference": fraction,
            "deep_cut_g_bar_max_m_s2": (
                None if fraction is None else fraction * REFERENCE_A0_M_S2
            ),
            "n_points": rec["n_points"],
            "n_galaxies": rec["n_galaxies"],
            "n_free_nuisance_parameters": n_free_nuisances_by_subset[subset_name],
            "per_rho": per_rho_rows,
        }

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": (
            "labeled_postdiction_profile_likelihood_declared_conventions_"
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
            **point_counts,
            "outermost_point": (
                "largest radius among the retained points of each galaxy"
            ),
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
                "profiled objective and a0-independent constants are dropped"
            ),
            "intra_galaxy_covariance": (
                "equal-correlation block R = (1-rho) I + rho J per galaxy, "
                "closed-form quadratic form; rho is a declared family "
                "parameter scanned over a fixed grid; a full covariance "
                "calibration for SPARC rotation curves is open, so intervals "
                "are reported at every rho and no single rho is selected"
            ),
            "rho_grid": list(config.rho_grid),
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
                    f"{config.n_upsilon} points over +-"
                    f"{config.halfwidth_sigma} prior sigma"
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
                    f"{config.n_distance} points over +-"
                    f"{config.halfwidth_sigma} sigma, clipped below at "
                    f"{DIST_FACTOR_MIN}"
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
                    f"{config.n_inclination} points over +-"
                    f"{config.halfwidth_sigma} sigma, clipped to "
                    f"[{INCLINATION_CLIP_DEG[0]}, {INCLINATION_CLIP_DEG[1]}] "
                    "deg"
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
                "min": config.log10_a0_min,
                "max": config.log10_a0_max,
                "n": config.n_a0,
                "step_dex": (
                    (config.log10_a0_max - config.log10_a0_min)
                    / (config.n_a0 - 1)
                ),
            },
        },
        "calibration": {
            "thresholds_delta_minus2lnl": {"68": THRESHOLD_68, "95": THRESHOLD_95},
            "status": (
                "asymptotic chi-square, one parameter, declared approximate "
                "under grid profiling with Gaussian penalties; delta is "
                "measured from the grid minimum, the headline "
                "maximum-likelihood point is the grid argmin so it lies "
                "inside its own intervals, the parabola vertex is a labeled "
                "refinement row, and intervals narrower than one grid step "
                "carry a grid-resolution flag"
            ),
        },
        "seed": config.seed,
        "subset_results": subset_results,
        "btfr_channel": {
            "estimator": (
                "a0_B = (v_A^2)^2 / (G M_b d^2), "
                "v_A^2 = (v_obs(r_out) s(i))^2 - d v_bar_cat^2(r_out, U), "
                "M_b = (U L_[3.6] + 1.33 M_HI) 1e9 M_sun at catalogue "
                "distance"
            ),
            "nuisance_application": (
                "per galaxy the profiled (U, d, i) at the subset "
                "maximum-likelihood a0 for the same rho; channel-B nuisance "
                "values are held at these point estimates inside the "
                "bootstrap, a declared simplification, while channel A "
                "re-profiles a0 exactly per replicate"
            ),
            "deep_requirement": (
                "outermost retained point passes the same fixed absolute "
                "cut g_bar_cat < f * 1.2e-10 m/s^2 as the channel-A subset"
            ),
        },
        "misfit_statement": _misfit_statement(subset_results, config),
        "verdict_rule": VERDICT_RULE,
        "inference_boundary": {
            "status": (
                "labeled postdiction likelihood on the seen committed "
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
                "open, and interval widths at different rho values bound "
                "the sensitivity"
            ),
            "full_sample_scope": (
                "the full-sample row extrapolates the deep-regime profile "
                "to all gradients; its misfit bears on that additive "
                "extension submodel under the declared error family, never "
                "on the framework beyond that scope, and the deep subsets "
                "are the derived regime of the candidate law"
            ),
            "neutrality": (
                "consistency of the channels or agreement of the interval "
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
    if include_provenance:
        receipt["sample"]["data_sha256"] = hashlib.sha256(
            (DATA / "table1.dat").read_bytes() + (DATA / "table2.dat").read_bytes()
        ).hexdigest()
    return receipt


def _paired_btfr_block(
    per_gal: list[dict[str, Any]],
    gal_by_name: dict[str, dict[str, Any]],
    grids_by_name: dict[str, dict[str, Any]],
    curves: np.ndarray,
    log10_a0_grid: np.ndarray,
    ml_index: int,
    rho: float,
    fraction: float,
    config: Config,
) -> dict[str, Any]:
    """Channel-B estimates and the paired bootstrap for one (subset, rho)."""
    n_not_deep = 0
    n_nonpositive_va2 = 0
    n_nonpositive_mb = 0
    used_names: list[str] = []
    used_logs: list[float] = []
    used_rows: list[int] = []
    n_candidates = 0
    for row_idx, g in enumerate(per_gal):
        gal = gal_by_name[g["name"]]
        n_candidates += 1
        k = gal["idx_out"]
        if not (gal["g_bar_cat"][k] < fraction * REFERENCE_A0_M_S2):
            n_not_deep += 1
            continue
        grids = grids_by_name[g["name"]]
        flat_index = int(g["rho"][rho]["argmin_flat"][ml_index])
        est = channel_b_estimate(gal, grids, flat_index)
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
    paired = paired_bootstrap(
        log10_a0_grid,
        common_curves,
        log_a0_b,
        config.bootstrap_replicates,
        config.seed,
    )
    return {
        "n_candidate_galaxies": n_candidates,
        "n_outermost_not_deep": n_not_deep,
        "n_nonpositive_anomalous_speed2": n_nonpositive_va2,
        "n_nonpositive_baryonic_mass": n_nonpositive_mb,
        "n_used": len(used_names),
        "a0_b_unweighted_log_mean_m_s2": float(10.0 ** np.mean(log_a0_b)),
        **paired,
    }


def _misfit_statement(
    subset_results: dict[str, Any], config: Config
) -> dict[str, Any]:
    """Likelihood-only description of what misfit looks like, with values.

    No posterior quantity appears here: the rows are reduced chi-square
    values, boundary flags, and grid-edge counts of the profiled objective.
    """
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
                    "interval95_pinned_at_grid_boundary": row[
                        "interval95_pinned_at_grid_boundary"
                    ],
                    "n_galaxies_profiled_nuisance_at_grid_edge_at_ml": row[
                        "n_galaxies_profiled_nuisance_at_grid_edge_at_ml"
                    ],
                }
            )
        observed[subset_name] = rows
    return {
        "rule": (
            "Under this model and error family, misfit would appear as "
            "reduced chi-square far above one at every rho on the declared "
            "grid, as interval endpoints pinned at the declared grid "
            "boundary, or as profiled nuisances stacked at their grid "
            "edges across the sample. These are likelihood statements; no "
            "posterior probability is computed or implied."
        ),
        "observed": observed,
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
    print(f"receipt written: {args.output}")
    for subset_name, block in out["subset_results"].items():
        for row in block["per_rho"]:
            line = (
                f"{subset_name} rho={row['rho']}: a0_ml={row['a0_ml_m_s2']:.4e} "
                f"68={row['a0_interval68_m_s2']} 95={row['a0_interval95_m_s2']} "
                f"red_chi2={row['reduced_chi2_dof_no_nuisance']:.3f}"
            )
            print(line)
            if "paired_btfr" in row:
                p = row["paired_btfr"]
                print(
                    f"  paired: ratio={p['point_log10_ratio_b_over_a']:+.4f} "
                    f"95={p['log10_ratio_95pct']} {p['verdict']}"
                )


if __name__ == "__main__":
    main()
