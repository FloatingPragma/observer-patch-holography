#!/usr/bin/env python3
"""Check the Cs-133 Feshbach scalarization theorems on a synthetic fixture.

A finite rotationally invariant self-adjoint Hamiltonian carries an isolated
ground manifold P = P3 (+) P4 with ranks 7 and 9, multiplicity-free, coupled
to a Q complement.  The Feshbach-Schur operator

    F_P(z) = P(H-z)P - P H Q [Q(H-z)Q]^{-1} Q H P

is scalar on each channel, F_P(z) = f3(z) P3 + f4(z) P4, each channel scalar
has derivative at or below -1 on the sampled resolvent interval, and one
sign bracket per channel isolates the unique root.  The two roots reproduce
the exact P-channel eigenvalues of H, and the clock normal form gives
epsilon = E4 - E3 = 4*a_Cs with E0 = (9*E4 + 7*E3)/16 on the K eigenvalues
-9/4 and 7/4.

This module is a synthetic theorem checker.  It consumes no cesium
measurement and no measured particle mass.  Physical promotion requires the
open source packets listed in the artifact.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "clock"
    / "cs133_feshbach_scalarization.json"
)

DIM_F3 = 7
DIM_F4 = 9
K_F3 = -9.0 / 4.0
K_F4 = 7.0 / 4.0

# Real z samples inside the complementary resolvent set of the fixture.
RESOLVENT_SAMPLE_GRID = (-0.4, -0.2, 0.0, 0.2, 0.4)

F3_BRACKET = (-0.2, 0.2)
F4_BRACKET = (-0.1, 0.3)

OPEN_SOURCE_PACKETS = (
    "R_ALPHA_ATOMIC_SCHEME_PACKET",
    "R_ELECTRON_ABSOLUTE_RATIO_PACKET",
    "R_QCD_NUCLEAR_133CS_PACKET",
    "R_ATOMIC_FESHBACH_SCALARS_133CS",
    "R_CLOCK_REFINEMENT_LIMIT",
    "R_NO_G_CLOCK_DAG",
    "R_CLOCK_PROSPECTIVE_FREEZE",
)


@dataclass(frozen=True)
class RootResult:
    lo: float
    hi: float
    mid: float
    f_lo: float
    f_hi: float
    iterations: int


def theorem_registry() -> dict[str, str]:
    """One-sentence statements for the ported Cs-133 clock theorems."""

    return {
        "CS133_FESHBACH_SPECTRAL_EQUIVALENCE_THEOREM": (
            "For z with Q(H-z)Q invertible, z is an eigenvalue of H exactly "
            "when 0 is an eigenvalue of F_P(z), with equal geometric "
            "multiplicity and an explicit eigenvector bijection."
        ),
        "CS133_FESHBACH_CHANNEL_SCALARIZATION_THEOREM": (
            "Under rotational invariance with multiplicity-free "
            "P = P3 (+) P4, F_P(z) = f3(z) P3 + f4(z) P4 for unique real "
            "scalar functions given by normalized channel traces."
        ),
        "CS133_FESHBACH_STRICT_MONOTONICITY_THEOREM": (
            "On every real complementary-resolvent interval, "
            "dF_P/dz <= -P, so f3' <= -1 and f4' <= -1 and each channel "
            "scalar is strictly decreasing."
        ),
        "CS133_SCALAR_SIGN_BRACKET_ROOT_THEOREM": (
            "A sign bracket f_F(L) > 0 > f_F(U) on one "
            "complementary-resolvent interval proves existence and "
            "uniqueness of the channel root in [L, U]."
        ),
        "CS133_CLOCK_SCALAR_NORMAL_FORM_THEOREM": (
            "With K eigenvalues -9/4 and 7/4, E0 = (9*E4 + 7*E3)/16 and "
            "a_Cs = (E4 - E3)/4 give H_clock = E0 P + a_Cs K, so "
            "epsilon_Cs = E4 - E3 = 4*a_Cs."
        ),
        "CS133_FULL_OPERATOR_TO_GAP_INTERVAL_THEOREM": (
            "An operator-norm enclosure ||H - H_tilde|| <= delta with "
            "stable Feshbach channel labels encloses each root within delta "
            "and the gap within 2*delta."
        ),
        "CURRENT_OPH_CS133_GAP_NONENTAILMENT_THEOREM": (
            "A source that leaves alpha_at, mu_e, the nuclear "
            "current/Compton packet, or the atomic remainder functional "
            "unfixed admits extensions with different Feshbach scalars and "
            "different roots, so no unique epsilon_Cs is entailed."
        ),
        "SOURCE_GEV_REQUIRES_OPERATIONAL_CLOCK_RATIO_THEOREM": (
            "E/GeV factors as E/(h*nu_Cs) times the exact unit conversion "
            "h*nu_Cs/GeV, so a source-only numerical GeV prediction is "
            "equivalent to a source-only prediction of an operational clock "
            "ratio."
        ),
    }


def make_synthetic_clock_model() -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Build the rotationally invariant fixture with an isolated P manifold.

    Basis order: P3(7), P4(9), Q3(7), Q4(9).  Every m-substate in one F
    channel carries the same 2x2 block, so the Hamiltonian commutes with the
    diagonal rotation representation and scalarization is exact.  The exact
    P-channel eigenvalues are the lower roots of the two 2x2 blocks.
    """

    d3, d4 = DIM_F3, DIM_F4
    n = 2 * (d3 + d4)
    h = np.zeros((n, n), dtype=float)
    s_p3 = slice(0, d3)
    s_p4 = slice(d3, d3 + d4)
    s_q3 = slice(d3 + d4, 2 * d3 + d4)
    s_q4 = slice(2 * d3 + d4, n)

    e3, e4 = 0.10, 0.20
    q3, q4 = 2.00, 2.50
    g3, g4 = 0.30, 0.40

    h[s_p3, s_p3] = e3 * np.eye(d3)
    h[s_p4, s_p4] = e4 * np.eye(d4)
    h[s_q3, s_q3] = q3 * np.eye(d3)
    h[s_q4, s_q4] = q4 * np.eye(d4)
    h[s_p3, s_q3] = g3 * np.eye(d3)
    h[s_q3, s_p3] = g3 * np.eye(d3)
    h[s_p4, s_q4] = g4 * np.eye(d4)
    h[s_q4, s_p4] = g4 * np.eye(d4)

    p3 = np.zeros((n, n))
    p4 = np.zeros((n, n))
    p3[s_p3, s_p3] = np.eye(d3)
    p4[s_p4, s_p4] = np.eye(d4)

    def lower(a: float, b: float, g: float) -> float:
        return (a + b - math.sqrt((b - a) ** 2 + 4 * g * g)) / 2

    exact = {
        "E3": lower(e3, q3, g3),
        "E4": lower(e4, q4, g4),
    }
    exact["epsilon"] = exact["E4"] - exact["E3"]
    exact["a_Cs"] = exact["epsilon"] / 4
    return h, p3, p4, exact


def orthonormal_basis_of_projector(p: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    vals, vecs = np.linalg.eigh((p + p.conj().T) / 2)
    keep = vals > 0.5
    if np.any((vals > tol) & (vals < 1 - tol)):
        raise ValueError("matrix is not a numerical orthogonal projector")
    return vecs[:, keep]


def projector_checks(p: np.ndarray, tol: float = 1e-10) -> dict[str, Any]:
    return {
        "hermitian_defect": float(np.linalg.norm(p - p.conj().T, 2)),
        "idempotence_defect": float(np.linalg.norm(p @ p - p, 2)),
        "rank": int(round(float(np.real(np.trace(p))))),
        "pass": bool(
            np.linalg.norm(p - p.conj().T, 2) <= tol
            and np.linalg.norm(p @ p - p, 2) <= tol
        ),
    }


def rotation_commutant_check(h: np.ndarray, seed: int = 0) -> dict[str, Any]:
    """Commutator of H with a random channel unitary acting jointly on P and Q."""

    rng = np.random.default_rng(seed)
    o3 = np.linalg.qr(rng.standard_normal((DIM_F3, DIM_F3)))[0]
    o4 = np.linalg.qr(rng.standard_normal((DIM_F4, DIM_F4)))[0]
    n = h.shape[0]
    u = np.zeros((n, n))
    u[0:7, 0:7] = o3
    u[7:16, 7:16] = o4
    u[16:23, 16:23] = o3
    u[23:32, 23:32] = o4
    return {
        "unitary_defect": float(np.linalg.norm(u @ u.T - np.eye(n), 2)),
        "commutator_norm": float(np.linalg.norm(h @ u - u @ h, 2)),
    }


def feshbach_operator(
    h: np.ndarray,
    u3: np.ndarray,
    u4: np.ndarray,
    uq: np.ndarray,
    z: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return F_P(z) and dF_P/dz in the orthonormal P-channel coordinates."""

    up = np.concatenate([u3, u4], axis=1)
    hpp = up.conj().T @ h @ up
    hpq = up.conj().T @ h @ uq
    hqp = hpq.conj().T
    hqq = uq.conj().T @ h @ uq
    aq = hqq - z * np.eye(hqq.shape[0])
    inv = np.linalg.inv(aq)
    f = hpp - z * np.eye(hpp.shape[0]) - hpq @ inv @ hqp
    derivative = -np.eye(hpp.shape[0]) - hpq @ (inv @ inv) @ hqp
    return (f + f.conj().T) / 2, (derivative + derivative.conj().T) / 2


def channel_scalar(block: np.ndarray) -> tuple[float, float]:
    scalar = float(np.real(np.trace(block)) / block.shape[0])
    defect = float(np.linalg.norm(block - scalar * np.eye(block.shape[0]), 2))
    return scalar, defect


def channel_functions(
    h: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
    z: float,
) -> dict[str, Any]:
    """Evaluate f3, f4, their derivatives, and the scalarization defects at z."""

    n = h.shape[0]
    p = p3 + p4
    q = np.eye(n) - p
    u3 = orthonormal_basis_of_projector(p3)
    u4 = orthonormal_basis_of_projector(p4)
    uq = orthonormal_basis_of_projector(q)
    f, df = feshbach_operator(h, u3, u4, uq, z)
    n3 = u3.shape[1]
    b3 = f[:n3, :n3]
    b4 = f[n3:, n3:]
    d3 = df[:n3, :n3]
    d4 = df[n3:, n3:]
    f3, defect3 = channel_scalar(b3)
    f4, defect4 = channel_scalar(b4)
    df3, ddefect3 = channel_scalar(d3)
    df4, ddefect4 = channel_scalar(d4)
    off = float(np.linalg.norm(f[:n3, n3:], 2))
    return {
        "z": z,
        "f3": f3,
        "f4": f4,
        "df3": df3,
        "df4": df4,
        "scalar_defect_F3": defect3,
        "scalar_defect_F4": defect4,
        "derivative_scalar_defect_F3": ddefect3,
        "derivative_scalar_defect_F4": ddefect4,
        "cross_channel_defect": off,
        "derivative_max_eigenvalue": float(np.max(np.linalg.eigvalsh(df))),
    }


def bisect_channel(
    h: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
    channel: str,
    lo: float,
    hi: float,
    tol: float = 1e-13,
    max_iter: int = 200,
) -> RootResult:
    """Bisect one strictly decreasing channel scalar on a sign bracket."""

    key = "f3" if channel == "F3" else "f4"

    def f(x: float) -> float:
        return float(channel_functions(h, p3, p4, x)[key])

    flo, fhi = f(lo), f(hi)
    if not (flo > 0 and fhi < 0):
        raise ValueError(f"expected decreasing sign bracket, got {flo}, {fhi}")
    it = 0
    while hi - lo > tol and it < max_iter:
        mid = (lo + hi) / 2
        fm = f(mid)
        if fm > 0:
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
        it += 1
    return RootResult(lo, hi, (lo + hi) / 2, flo, fhi, it)


def gap_normal_form_checks(e3: float, e4: float) -> dict[str, Any]:
    """Exact clock normal form H_clock = E0 P + a_Cs K on the P manifold."""

    diag3 = [1.0] * DIM_F3 + [0.0] * DIM_F4
    diag4 = [0.0] * DIM_F3 + [1.0] * DIM_F4
    p3 = np.diag(diag3)
    p4 = np.diag(diag4)
    p = p3 + p4
    k = K_F3 * p3 + K_F4 * p4
    e0 = (9.0 * e4 + 7.0 * e3) / 16.0
    a = (e4 - e3) / 4.0
    h_clock = e3 * p3 + e4 * p4
    return {
        "E0": e0,
        "a_Cs": a,
        "epsilon": e4 - e3,
        "K_eigenvalues": [K_F3, K_F4],
        "epsilon_minus_four_a": (e4 - e3) - 4.0 * a,
        "channel_identity_residual_F3": abs(e0 + K_F3 * a - e3),
        "channel_identity_residual_F4": abs(e0 + K_F4 * a - e4),
        "normal_form_residual": float(np.linalg.norm(h_clock - e0 * p - a * k, 2)),
    }


def build_artifact() -> dict[str, Any]:
    h, p3, p4, exact = make_synthetic_clock_model()
    projector_p3 = projector_checks(p3)
    projector_p4 = projector_checks(p4)
    orthogonality_defect = float(np.linalg.norm(p3 @ p4, 2))
    selfadjoint_defect = float(np.linalg.norm(h - h.T, 2))
    rotation = rotation_commutant_check(h)

    samples = [channel_functions(h, p3, p4, z) for z in RESOLVENT_SAMPLE_GRID]
    max_scalar_defect = max(
        max(s["scalar_defect_F3"], s["scalar_defect_F4"]) for s in samples
    )
    max_cross_defect = max(s["cross_channel_defect"] for s in samples)
    max_derivative = max(
        max(s["df3"], s["df4"], s["derivative_max_eigenvalue"]) for s in samples
    )

    r3 = bisect_channel(h, p3, p4, "F3", *F3_BRACKET)
    r4 = bisect_channel(h, p3, p4, "F4", *F4_BRACKET)
    values = np.linalg.eigvalsh(h)
    e3_multiplicity = int(np.sum(np.isclose(values, exact["E3"], atol=1e-11)))
    e4_multiplicity = int(np.sum(np.isclose(values, exact["E4"], atol=1e-11)))
    normal_form = gap_normal_form_checks(r3.mid, r4.mid)

    checks = {
        "H_selfadjoint": selfadjoint_defect < 1e-12,
        "H_commutes_with_channel_rotations": rotation["commutator_norm"] < 1e-12,
        "P3_is_rank_7_projector": bool(projector_p3["pass"]) and projector_p3["rank"] == 7,
        "P4_is_rank_9_projector": bool(projector_p4["pass"]) and projector_p4["rank"] == 9,
        "P3_P4_orthogonal": orthogonality_defect < 1e-12,
        "channel_scalarization_residual_below_1e_12": max_scalar_defect < 1e-12,
        "cross_channel_blocks_vanish": max_cross_defect < 1e-12,
        "derivatives_at_or_below_minus_one": max_derivative <= -1.0 + 1e-12,
        "F3_sign_bracket_isolates_one_root": r3.f_lo > 0.0 > r3.f_hi,
        "F4_sign_bracket_isolates_one_root": r4.f_lo > 0.0 > r4.f_hi,
        "F3_root_matches_exact_eigenvalue": abs(r3.mid - exact["E3"]) < 1e-11,
        "F4_root_matches_exact_eigenvalue": abs(r4.mid - exact["E4"]) < 1e-11,
        "E3_multiplicity_7_in_full_spectrum": e3_multiplicity == 7,
        "E4_multiplicity_9_in_full_spectrum": e4_multiplicity == 9,
        "gap_positive": r4.lo - r3.hi > 0.0,
        "epsilon_equals_four_a": abs(normal_form["epsilon_minus_four_a"]) < 1e-15,
        "normal_form_residual_below_1e_12": normal_form["normal_form_residual"] < 1e-12,
    }
    checks_pass = all(bool(value) for value in checks.values())

    return {
        "artifact": "oph_cs133_feshbach_scalarization",
        "status": "SYNTHETIC_FIXTURE_THEOREM_CHECKS_ONLY",
        "synthetic_fixture": True,
        "measured_clock_data_consumed": False,
        "measured_particle_data_consumed": False,
        "physical_source_prediction_ready": False,
        "public_clock_gap_promotion_allowed": False,
        "theorem_registry": theorem_registry(),
        "synthetic_model": {
            "basis_order": "P3(7), P4(9), Q3(7), Q4(9)",
            "P_manifold": "multiplicity-free V3 (+) V4 with ranks 7 and 9",
            "isolation": (
                "The Q block energies 2.0 and 2.5 sit above the P block "
                "energies 0.1 and 0.2, so the dressed ground manifold is "
                "isolated below the Q branch."
            ),
            "rotational_invariance": (
                "Every m-substate in one F channel carries the same 2x2 "
                "block, so H commutes with the diagonal rotation "
                "representation."
            ),
            "P_energies": {"e3": 0.10, "e4": 0.20},
            "Q_energies": {"q3": 2.00, "q4": 2.50},
            "couplings": {"g3": 0.30, "g4": 0.40},
            "exact_channel_eigenvalues": {"E3": exact["E3"], "E4": exact["E4"]},
            "exact_epsilon": exact["epsilon"],
            "exact_a_Cs": exact["a_Cs"],
        },
        "projector_checks": {
            "P3": projector_p3,
            "P4": projector_p4,
            "P3_P4_orthogonality_defect": orthogonality_defect,
            "H_selfadjoint_defect": selfadjoint_defect,
            "rotation_commutant": rotation,
        },
        "scalarization_checks": {
            "resolvent_sample_grid": list(RESOLVENT_SAMPLE_GRID),
            "samples": samples,
            "max_scalar_defect": max_scalar_defect,
            "max_cross_channel_defect": max_cross_defect,
            "max_derivative_value": max_derivative,
        },
        "root_checks": {
            "F3_bracket": list(F3_BRACKET),
            "F4_bracket": list(F4_BRACKET),
            "F3_root": dataclasses.asdict(r3),
            "F4_root": dataclasses.asdict(r4),
            "epsilon_interval": [r4.lo - r3.hi, r4.hi - r3.lo],
            "a_interval": [(r4.lo - r3.hi) / 4, (r4.hi - r3.lo) / 4],
            "E3_multiplicity_in_full_spectrum": e3_multiplicity,
            "E4_multiplicity_in_full_spectrum": e4_multiplicity,
        },
        "gap_normal_form": normal_form,
        "reduction": {
            "pipeline": (
                "full cesium Hamiltonian -> two scalar interval evaluators "
                "f3, f4 -> root enclosures E3, E4 -> eps_Cs = E4 - E3"
            ),
            "public_packet_form": (
                "The public atomic packet may export two scalar interval "
                "evaluators with resolvent certificates and sign brackets "
                "in place of a 55-electron matrix."
            ),
        },
        "open_source_packets": list(OPEN_SOURCE_PACKETS),
        "checks": checks,
        "checks_pass": checks_pass,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"], "checks_pass": artifact["checks_pass"]}, indent=2))
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
