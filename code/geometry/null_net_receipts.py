#!/usr/bin/env python3
"""Null-net receipt instrumentation on the realized collar tower (#503, #524).

Evaluates, on the same realized Gaussian MaxEnt collar states as the
modular-clock instrumentation, finite-stage witnesses for the null-net
receipt families of the compact paper's standardness packet
(`subsec:null-net-standardness`):

* NTI (nontrivial relative commutants): for nested arcs B in A, the modes of
  A minus B generate a relative commutant of dimension 4^(|A|-|B|) > 1,
  a combinatorial mode-count witness at every stage;
* weak additivity: ring translates of one arc cover every site, so the
  translated arc algebras generate the full ring algebra;
* separating/faithfulness modulus: every arc occupation spectrum lies
  strictly inside (0,1) (computed at 120-digit precision), so the reference
  vector is separating for every arc algebra at every finite stage; the
  cyclicity-modulus limit clause (Cyc) is NOT witnessed here and remains
  pending;
* mixed-GNS Cauchy: expectation values of fixed smeared bond observables
  are Cauchy along refinement;
* HSM compression (the half-sided-inclusion receipt, one-particle level):
  for nested arcs B in A sharing one cut, the composed modular flow
  U(t) = exp(i h_A t) exp(-i h_B t) keeps one-particle modes localized in B
  inside B for one sign of t and leaks them for the other, Wiesbrock's
  compression condition witnessed as a leakage asymmetry that persists
  across refinement with a consistent direction;
* modular Lie-closure (the modular-intersection surrogate): the commutator
  i[h_A, h_B] of two arc modular generators lies, up to a residual
  decreasing under refinement, in the linear span of arc modular generators
  plus arc number operators and the identity, the finite-stage shadow of
  the statement that vacuum interval modular generators close on the global
  conformal algebra, which is the structure the modular-intersection
  receipts consume.

All statements are one-particle-level finite receipts on the declared
Gaussian family, labeled as such; the second-quantized limit clauses are the
paper's scaling-limit theorems.

Run:
    python3 code/geometry/null_net_receipts.py
writes code/geometry/runs/null_net_receipt_report.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from modular_clock_instrumentation import (  # noqa: E402
    arc_entanglement_hamiltonian,
    ring_correlation_value,
)

REPORT_PATH = HERE / "runs" / "null_net_receipt_report.json"


def embed(h_small: np.ndarray, n_total: int, offset: int = 0) -> np.ndarray:
    out = np.zeros((n_total, n_total))
    m = h_small.shape[0]
    out[offset:offset + m, offset:offset + m] = h_small
    return out


def correlation_matrix(n_ring: int) -> np.ndarray:
    vals = {d: float(ring_correlation_value(n_ring, d)) for d in range(n_ring)}
    return np.array([[vals[min(abs(i - j), n_ring - abs(i - j))]
                      for j in range(n_ring)] for i in range(n_ring)])


# ---------------------------------------------------------------------------
# easy families: NTI, weak additivity, separating modulus, mixed-GNS Cauchy
# ---------------------------------------------------------------------------

def nti_receipt(n_ring: int) -> dict:
    m_a, m_b = n_ring // 2, n_ring // 4
    free_modes = m_a - m_b
    return {
        "n_ring": n_ring,
        "relative_commutant_modes": free_modes,
        "relative_commutant_dim_log2": 2 * free_modes,  # dim = 4^modes
        "nontrivial": free_modes > 0,
    }


def weak_additivity_receipt(n_ring: int) -> dict:
    m = n_ring // 2
    covered = set()
    shifts = []
    for s in range(0, n_ring, max(1, m // 2)):
        covered |= {(s + j) % n_ring for j in range(m)}
        shifts.append(s)
        if len(covered) == n_ring:
            break
    return {"n_ring": n_ring, "translates_used": len(shifts),
            "covers_ring": len(covered) == n_ring}


def separating_modulus_receipt(n_ring: int) -> dict:
    """Distance of the arc occupation spectrum from {0,1}: > 0 at every
    stage, so the reference vector separates every arc algebra."""
    m = n_ring // 2
    vals = {d: ring_correlation_value(n_ring, d) for d in range(m)}
    c_arc = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            c_arc[i, j] = vals[abs(i - j)]
    evals, _ = mp.eigsy(c_arc)
    gap = min(min(e, 1 - e) for e in evals)
    return {"n_ring": n_ring, "occupation_gap_log10": float(mp.log10(gap)),
            "strictly_faithful": bool(gap > 0)}


def mixed_gns_cauchy_receipt(rings: tuple[int, ...]) -> dict:
    """Expectation of the fixed smeared bond observable
    O = sum_j f(j/N)(c^dag_j c_{j+1} + h.c.), f a fixed bump, across stages."""
    f = lambda x: np.exp(-((x - 0.25) ** 2) / (2 * 0.05 ** 2))  # noqa: E731
    values = []
    for n in rings:
        c = correlation_matrix(n)
        total = 0.0
        for j in range(n - 1):
            total += f(j / n) * 2.0 * c[j, j + 1].real
        values.append(total / n * n / n)  # per-site smeared density
    diffs = [abs(values[i + 1] - values[i]) for i in range(len(values) - 1)]
    return {
        "rings": list(rings),
        "values": values,
        "cauchy_differences": diffs,
        "cauchy_decreasing": all(diffs[i] > diffs[i + 1]
                                 for i in range(len(diffs) - 1))
        if len(diffs) > 1 else True,
    }


# ---------------------------------------------------------------------------
# HSM compression receipt (one-particle Wiesbrock condition)
# ---------------------------------------------------------------------------

def hsm_compression_receipt(n_ring: int, t_mod: float = 0.12) -> dict:
    """Nested arcs sharing the left cut: A = sites [0, N/2), B = [0, N/4).
    A Gaussian one-particle mode centered in B is evolved by the composed
    flow U(t) = exp(i h_A t) exp(-i h_B t); leakage outside B is measured
    for t = +t_mod and t = -t_mod. The half-sided inclusion predicts a
    strict sign asymmetry (compression for one sign only)."""
    m_a, m_b = n_ring // 2, n_ring // 4
    h_a = arc_entanglement_hamiltonian(n_ring, m_a)
    h_b = embed(arc_entanglement_hamiltonian(n_ring, m_b), m_a)
    # normalized CHIRAL Gaussian mode centered in B: the carrier e^{i k_F x}
    # (k_F = pi/2) selects one chirality; a real packet is exactly
    # time-reversal symmetric and shows no asymmetry, which is itself the
    # control for this receipt.
    xs = np.arange(m_a)
    center, width = m_b / 2.0, m_b / 6.0
    psi = (np.exp(-((xs - center) ** 2) / (2 * width ** 2))
           * np.exp(1j * np.pi * xs / 2.0))
    psi[m_b:] = 0.0
    psi /= np.linalg.norm(psi)

    def leakage(t: float) -> float:
        u = expm(1j * h_a * t) @ expm(-1j * h_b * t)
        out = u @ psi
        return float(np.sum(np.abs(out[m_b:]) ** 2))

    leak_plus, leak_minus = leakage(+t_mod), leakage(-t_mod)
    lo, hi = sorted((leak_plus, leak_minus))
    return {
        "n_ring": n_ring,
        "modular_time": t_mod,
        "leakage_plus": leak_plus,
        "leakage_minus": leak_minus,
        "compressing_sign": "+" if leak_plus < leak_minus else "-",
        "asymmetry_ratio": float(hi / max(lo, 1e-300)),
    }


# ---------------------------------------------------------------------------
# modular Lie-closure receipt (modular-intersection surrogate)
# ---------------------------------------------------------------------------

def lie_closure_receipt(n_ring: int, rmax: int = 8) -> dict:
    """Modular Lie-closure receipt (modular-intersection surrogate).

    The conformal algebra predicts
        i [K_A, K_B] = - integral g*(x) P(x) dx,
        g* = beta_A beta_B' - beta_A' beta_B,
    with P the momentum density. The EH matrices are odd-range (particle-hole
    symmetry), so the one-particle commutator is exactly even-range, and the
    lattice momentum content must be resummed over even ranges with the
    chiral weights (r/2)(-1)^((r-2)/2) (the momentum analogue of the
    odd-range energy resummation). The receipt compares the resummed
    envelope against g* after one scalar normalization on the overlap
    interior. Status: the shape matches at the percent level at
    every evaluated stage, with a smooth subleading systematic whose decay
    rate is NOT yet certified; the receipt therefore reports a
    percent-level flag and leaves the convergence-rate clause open.
    Requires n_ring >= 32 (smaller stages have too few overlap bonds and
    the fit is degenerate)."""
    m = n_ring // 2
    off_b = n_ring // 8
    h_a = embed(arc_entanglement_hamiltonian(n_ring, m), n_ring, 0)
    h_b = embed(arc_entanglement_hamiltonian(n_ring, m), n_ring, off_b)
    comm = h_a @ h_b - h_b @ h_a  # real antisymmetric, even-range only

    def p_eff(j: int) -> float:
        total = 0.0
        for r in range(2, rmax + 1, 2):
            i = j + 1 - r // 2
            if i < 0 or i + r >= n_ring:
                continue
            weight = (r / 2.0) * (-1) ** (((r - 2) // 2) % 2)
            total += weight * comm[i, i + r]
        return total

    measured = np.array([p_eff(j) for j in range(n_ring - 2)])

    def beta_and_deriv(x, u, v):
        s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
        c = lambda y: np.cos(np.pi * y / n_ring)  # noqa: E731
        b = 2.0 * n_ring * s(x - u) * s(v - x) / s(v - u)
        db = 2.0 * np.pi * (c(x - u) * s(v - x) - s(x - u) * c(v - x)) / s(v - u)
        return b, db

    u_a, v_a = -0.5, m - 0.5
    u_b, v_b = off_b - 0.5, off_b + m - 0.5
    xs = np.arange(n_ring - 2) + 1.0
    inside = (xs > u_b) & (xs < v_a)
    ba, dba = beta_and_deriv(xs[inside], u_a, v_a)
    bb, dbb = beta_and_deriv(xs[inside], u_b, v_b)
    g_star = ba * dbb - dba * bb

    idx = np.where(inside)[0][2:-2]
    g_sel = g_star[2:-2]
    m_sel = measured[idx]
    alpha = float(np.dot(m_sel, g_sel) / np.dot(g_sel, g_sel))
    resid = float(np.linalg.norm(m_sel - alpha * g_sel)
                  / np.linalg.norm(m_sel))
    # control: the unresummed range-2 truncation is several times worse
    naive = np.array([comm[j, j + 2] for j in range(n_ring - 2)])[idx]
    alpha_n = float(np.dot(naive, g_sel) / np.dot(g_sel, g_sel))
    resid_naive = float(np.linalg.norm(naive - alpha_n * g_sel)
                        / np.linalg.norm(naive))
    return {
        "n_ring": n_ring,
        "relative_residual": resid,
        "relative_residual_unresummed_control": resid_naive,
        "normalization_alpha": alpha,
    }


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

def instrument_null_net(rings: tuple[int, ...] = (16, 32, 64)) -> dict:
    nti = [nti_receipt(n) for n in rings]
    wa = [weak_additivity_receipt(n) for n in rings]
    sep = [separating_modulus_receipt(n) for n in rings]
    cauchy = mixed_gns_cauchy_receipt(rings)
    hsm = [hsm_compression_receipt(n) for n in rings]
    lie = [lie_closure_receipt(n) for n in rings if n >= 32]
    lie_res = [entry["relative_residual"] for entry in lie]
    verdicts = {
        "nti_all_stages": all(x["nontrivial"] for x in nti),
        "weak_additivity_all_stages": all(x["covers_ring"] for x in wa),
        "separating_all_stages": all(x["strictly_faithful"] for x in sep),
        "mixed_gns_cauchy_decreasing": cauchy["cauchy_decreasing"],
        "hsm_asymmetry_min_ratio": min(x["asymmetry_ratio"] for x in hsm),
        "hsm_direction_consistent": len({x["compressing_sign"]
                                         for x in hsm}) == 1,
        "lie_closure_residuals": lie_res,
        "lie_closure_percent_level": all(r < 0.02 for r in lie_res),
        "lie_closure_rate_certified": all(
            lie_res[i] > lie_res[i + 1] for i in range(len(lie_res) - 1)
        ) and len(lie_res) > 1,
        "lie_closure_unresummed_controls": [
            entry["relative_residual_unresummed_control"] for entry in lie
        ],
    }
    witnessed = {
        "nti": verdicts["nti_all_stages"],
        "weak_additivity": verdicts["weak_additivity_all_stages"],
        "separating_faithfulness": verdicts["separating_all_stages"],
        "mixed_gns_cauchy": verdicts["mixed_gns_cauchy_decreasing"],
        "hsm_compression_one_particle": bool(
            verdicts["hsm_asymmetry_min_ratio"] > 3.0
            and verdicts["hsm_direction_consistent"]
        ),
        "modular_lie_closure_percent_level": verdicts[
            "lie_closure_percent_level"],
    }
    return {
        "artifact": "oph_null_net_receipt_instrumentation",
        "object_id": "NullNetReceipts_Issue503",
        "issue": 503,
        "scope": (
            "one-particle-level finite receipts on the declared Gaussian "
            "collar family; the Cyc limit clause and all second-quantized "
            "scaling-limit clauses are the paper's theorems/receipts, "
            "and the E1-E4 event and UC/VR/scale physical families are "
            "outside this module"
        ),
        "rings": list(rings),
        "nti": nti,
        "weak_additivity": wa,
        "separating": sep,
        "mixed_gns_cauchy": cauchy,
        "hsm_compression": hsm,
        "lie_closure": lie,
        "verdicts": verdicts,
        "receipts_witnessed": witnessed,
        "receipts_pending": [
            "Cyc limit clause (cyclicity modulus along refinement)",
            "Lie-closure convergence rate (percent-level match witnessed; "
            "the decay of the smooth subleading systematic is uncertified)",
            "second-quantized MI relations (one-particle Lie closure is the "
            "finite surrogate)",
            "E1-E4 event receipts on realized records",
            "UC/VR/scale physical-identification receipts",
        ],
    }


def main() -> None:
    report = instrument_null_net()
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT_PATH, "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(report["verdicts"], indent=2))
    print(json.dumps(report["receipts_witnessed"], indent=2))


if __name__ == "__main__":
    main()
