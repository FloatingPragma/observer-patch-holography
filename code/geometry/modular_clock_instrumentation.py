#!/usr/bin/env python3
"""Free-fermion modular-clock instrumentation of the cyclic cap-net tower.

Evaluates, on realized Gaussian MaxEnt reference states, the two receipt
families the cyclic cap-net repair run left pending at the geometry level
(compact paper, Definition `def:cap-mesh-crossratio-receipts`):

* the independently normalized geometric 2pi-KMS receipt, and
* the modular cross-ratio Cauchy receipt with identified Moebius limit.

Construction (all modular data computed, none declared):

1. Reference states: half-filled ground states of critical nearest-neighbour
   hopping (anti-periodic twist, Fermi velocity 1) on the stage-r
   cap-boundary collar rings N_r = 16, 32, 64. Correlations are closed-form
   plane-wave sums, evaluated in extended precision (mpmath, 120 digits),
   because the arc correlation spectrum reaches 1e-30 of the endpoints,
   far beyond float64.
2. Entanglement Hamiltonian of the half-ring arc: h = log((1-C_I)/C_I) by
   extended-precision eigendecomposition.
3. Modular velocity profile: the lattice EH carries odd-range hopping terms
   whose chiral weight at half filling (k_F = pi/2) is r (-1)^((r-1)/2), so
   the profile is the resummation
       beta(j) = -2 sum_{r odd} r (-1)^((r-1)/2) h_{i, i+r},
   with the range-r term centered at the bond (i = j + (1-r)/2). The
   nearest-neighbour truncation alone is off by ~9 percent and does NOT
   converge; the resummation converges to the 2pi-normalized conformal
   Bisognano-Wichmann profile
       beta_CFT(x) = 2N sin(pi(x-u)/N) sin(pi(v-x)/N) / sin(pi(v-u)/N)
   (bond j at x = j + 1/2, cuts at u = -1/2, v = m - 1/2).
4. Receipts: (KMS) median interior relative residual against beta_CFT,
   decreasing along refinement, with the wrongly normalized (x1.2) profile
   separated by a widening factor; (CR) modular transport cross-ratios
   exp(2pi sum 1/beta) over the fixed angular quadruple (1/4, 3/4 of the
   arc) converge, Cauchy along refinement, to the Moebius cross-ratio
       cr_geo = [S(b-u) S(v-a)] / [S(v-b) S(a-u)], S(x) = sin(pi x/N),
   which is stage-independent for proportional positions.

Scope, stated exactly: this is BOUNDARY-COLLAR instrumentation on the
declared Gaussian MaxEnt family (an Axiom-level input). Cap-interior
modular data, the null-net families (Cyc/NTI/weak additivity/MI), the event
families E1-E4, and the physical-identification receipts are pending.

Run:
    python3 code/geometry/modular_clock_instrumentation.py
writes code/geometry/runs/modular_clock_instrumentation_report.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
REPORT_PATH = HERE / "runs" / "modular_clock_instrumentation_report.json"

mp.mp.dps = 120
RMAX = 9  # odd-range resummation depth for the profile


# ---------------------------------------------------------------------------
# extended-precision modular data
# ---------------------------------------------------------------------------

def ring_correlation_value(n_ring: int, d: int) -> mp.mpf:
    """C(d) = (1/N) sum_{filled k} cos(k d), k = (2t+1) pi / N (anti-periodic),
    filled = negative-energy modes of -cos k (half filling)."""
    s = mp.mpf(0)
    for t in range(n_ring):
        k = (2 * t + 1) * mp.pi / n_ring
        if mp.cos(k) > 0:
            s += mp.cos(k * d)
    return s / n_ring


def arc_entanglement_hamiltonian(n_ring: int, m: int) -> np.ndarray:
    """h = log((1-C_I)/C_I) for the m-site arc, via mpmath eigendecomposition;
    returned as float64 (safe: the extreme logs are O(70), well inside
    float range, only the eigenvalues needed extended precision)."""
    vals = {d: ring_correlation_value(n_ring, d) for d in range(m)}
    c_arc = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            c_arc[i, j] = vals[abs(i - j)]
    evals, vecs = mp.eigsy(c_arc)
    diag = mp.matrix(m, m)
    for i in range(m):
        diag[i, i] = mp.log((1 - evals[i]) / evals[i])
    h_arc = vecs * diag * vecs.T
    return np.array([[float(h_arc[i, j]) for j in range(m)] for i in range(m)])


def resummed_profile(h_arc: np.ndarray, rmax: int = RMAX) -> np.ndarray:
    """beta(j) = -2 sum_{r odd <= rmax} r (-1)^((r-1)/2) h_{i,i+r}, centered
    at bond j; entries where the largest stencil does not fit use the
    largest odd range that does."""
    m = h_arc.shape[0]
    beta = np.zeros(m - 1)
    for j in range(m - 1):
        total = 0.0
        for r in range(1, rmax + 1, 2):
            i = j + (1 - r) // 2
            if i < 0 or i + r >= m:
                continue
            weight = r * (-1) ** (((r - 1) // 2) % 2)
            total += -2.0 * weight * h_arc[i, i + r]
        beta[j] = total
    return beta


def cft_profile(n_ring: int, x: np.ndarray, u: float, v: float) -> np.ndarray:
    s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
    return 2.0 * n_ring * s(x - u) * s(v - x) / s(v - u)


# ---------------------------------------------------------------------------
# receipts
# ---------------------------------------------------------------------------

def stage_data(n_ring: int) -> dict:
    m = n_ring // 2
    h_arc = arc_entanglement_hamiltonian(n_ring, m)
    beta = resummed_profile(h_arc)
    u, v = -0.5, m - 0.5
    xs = np.arange(m - 1) + 0.5
    beta_geo = cft_profile(n_ring, xs, u, v)
    return {"m": m, "beta": beta, "beta_geo": beta_geo, "u": u, "v": v}


def kms_profile_receipt(data: dict, n_ring: int,
                        interior_fraction: float = 0.5) -> dict:
    beta, beta_geo = data["beta"], data["beta_geo"]
    m1 = len(beta)
    lo = int(m1 * (1 - interior_fraction) / 2)
    hi = m1 - lo
    rel = np.abs(beta[lo:hi] / beta_geo[lo:hi] - 1.0)
    rel_wrong = np.abs(beta[lo:hi] / (1.2 * beta_geo[lo:hi]) - 1.0)
    return {
        "n_ring": n_ring,
        "arc_sites": data["m"],
        "interior_bonds": int(hi - lo),
        "median_relative_residual_2pi": float(np.median(rel)),
        "max_relative_residual_2pi": float(np.max(rel)),
        "median_relative_residual_wrong_1p2": float(np.median(rel_wrong)),
        "separation_factor": float(np.median(rel_wrong) / np.median(rel)),
    }


def crossratio_receipt(data: dict, n_ring: int) -> dict:
    """Transport cross-ratio over the fixed angular quadruple (1/4, 3/4 of
    the arc, anchored at the cuts)."""
    m, beta = data["m"], data["beta"]
    # exact quarter/three-quarter fractions of the arc: a - u = m/4,
    # b - u = 3m/4, so the Moebius target is stage-independent
    a_idx, b_idx = m // 4 - 1, (3 * m) // 4 - 1
    t = float(np.sum(1.0 / beta[a_idx:b_idx]))
    cr_latt = float(np.exp(2.0 * np.pi * t))
    s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
    u, v = data["u"], data["v"]
    a, b = a_idx + 0.5, b_idx + 0.5
    cr_geo = float((s(b - u) * s(v - a)) / (s(v - b) * s(a - u)))
    return {
        "n_ring": n_ring,
        "quadruple_bonds": [a_idx, b_idx],
        "cr_lattice": cr_latt,
        "cr_moebius": cr_geo,
        "relative_error": abs(cr_latt / cr_geo - 1.0),
    }


def instrument_tower(rings: tuple[int, ...] = (16, 32, 64)) -> dict:
    stages = {n: stage_data(n) for n in rings}
    kms = [kms_profile_receipt(stages[n], n) for n in rings]
    crs = [crossratio_receipt(stages[n], n) for n in rings]
    kms_res = [k["median_relative_residual_2pi"] for k in kms]
    cr_err = [c["relative_error"] for c in crs]
    cr_vals = [c["cr_lattice"] for c in crs]
    cauchy = [abs(cr_vals[i + 1] - cr_vals[i]) for i in range(len(cr_vals) - 1)]
    verdicts = {
        "kms_residual_decreasing": all(
            kms_res[i] > kms_res[i + 1] for i in range(len(kms_res) - 1)
        ),
        "kms_final_median_residual": kms_res[-1],
        "wrong_normalization_separated": all(
            k["separation_factor"] > 5.0 for k in kms
        ),
        "crossratio_error_decreasing": all(
            cr_err[i] > cr_err[i + 1] for i in range(len(cr_err) - 1)
        ),
        "crossratio_final_relative_error": cr_err[-1],
        "crossratio_cauchy_decreasing": all(
            cauchy[i] > cauchy[i + 1] for i in range(len(cauchy) - 1)
        ) if len(cauchy) > 1 else True,
    }
    witnessed = {
        "geometric_2pi_kms_boundary_collar": bool(
            verdicts["kms_residual_decreasing"]
            and verdicts["wrong_normalization_separated"]
        ),
        "modular_cross_ratio_boundary_collar": bool(
            verdicts["crossratio_error_decreasing"]
            and verdicts["crossratio_cauchy_decreasing"]
        ),
    }
    return {
        "artifact": "oph_modular_clock_instrumentation",
        "object_id": "ModularClockInstrumentation_Issue503",
        "issue": 503,
        "scope": (
            "boundary-collar instrumentation: declared Gaussian MaxEnt "
            "reference states (critical hopping, anti-periodic twist, half "
            "filling, v_F = 1) on the stage-r cap-boundary collar rings; "
            "modular data computed at 120-digit precision; profile extracted "
            "by odd-range chiral resummation (nearest-neighbour truncation "
            "alone is off by ~9 percent and does not converge); receipts "
            "evaluated are the boundary-circle profile-KMS and cross-ratio "
            "clauses; cap-interior modular data are pending"
        ),
        "rings": list(rings),
        "resummation_rmax": RMAX,
        "kms_profile_receipt": kms,
        "crossratio_receipt": crs,
        "verdicts": verdicts,
        "receipts_witnessed": witnessed,
        "receipts_pending": [
            "cap-interior modular data (full 2D cap algebras)",
            "Cyc/NTI/weak-additivity/MI null-net receipts",
            "E1-E4 event receipts",
            "UC/VR/scale physical-identification receipts",
        ],
    }


def main() -> None:
    report = instrument_tower()
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT_PATH, "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(report["verdicts"], indent=2))
    print(json.dumps(report["receipts_witnessed"], indent=2))


if __name__ == "__main__":
    main()
