#!/usr/bin/env python3
"""Analyze the strongest current D12 quark mass branch and CKM residual.

Chain role: keep the honest D12 quark continuation program explicit without
promoting the light-quark selector value to recovered-core status.

Mathematics: apply a one-scalar D12 light-quark overlap selector candidate to
the current forward Yukawas, evaluate the resulting mass branch, and compute the
matrix-valued left-unitary residual generator needed to reach the target CKM
surface while preserving that mass branch.

OPH-derived inputs: the current forward Yukawa artifact, the quark exactness
audit, and the already-emitted spread package.

Output: a diagnostic continuation artifact carrying the best current D12
mass-side candidate and the explicit remaining CKM/CP generator residual.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FORWARD_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
AUDIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _left_diag(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = matrix @ matrix.conjugate().T
    eig_vals, eig_vecs = np.linalg.eigh(hermitian)
    order = np.argsort(eig_vals)
    eig_vals = eig_vals[order]
    eig_vecs = eig_vecs[:, order]
    return np.sqrt(np.clip(eig_vals, 0.0, None)), eig_vecs


def _jarlskog(v_ckm: np.ndarray) -> float:
    return float(np.imag(v_ckm[0, 0] * v_ckm[1, 1] * np.conjugate(v_ckm[0, 1]) * np.conjugate(v_ckm[1, 0])))


def _principal_unitary_log(matrix: np.ndarray) -> np.ndarray:
    eig_vals, eig_vecs = np.linalg.eig(matrix)
    phases = np.angle(eig_vals)
    log_diag = np.diag(1j * phases)
    return eig_vecs @ log_diag @ np.linalg.inv(eig_vecs)


def _apply_delta(delta_value: float, y_u: np.ndarray, y_d: np.ndarray, sigma_u: float, sigma_d: float) -> dict[str, Any]:
    b_ord = np.asarray([-1.0, 0.0, 1.0], dtype=float)
    tau_u = 0.5 * delta_value * sigma_d / (sigma_u + sigma_d)
    tau_d = 0.5 * delta_value * sigma_u / (sigma_u + sigma_d)
    d_u = np.diag(np.exp(tau_u * b_ord))
    d_d = np.diag(np.exp(tau_d * b_ord))
    y_u_trial = d_u @ y_u @ d_u
    y_d_trial = d_d @ y_d @ d_d
    m_u, u_left = _left_diag(y_u_trial)
    m_d, d_left = _left_diag(y_d_trial)
    v_ckm = u_left.conjugate().T @ d_left
    return {
        "Delta_ud_overlap": float(delta_value),
        "Lambda_ud_B_transport": float((sigma_u * sigma_d / (2.0 * (sigma_u + sigma_d))) * delta_value),
        "tau_u_log_per_side": float(tau_u),
        "tau_d_log_per_side": float(tau_d),
        "m_u": [float(value) for value in m_u.tolist()],
        "m_d": [float(value) for value in m_d.tolist()],
        "abs_V_CKM": np.abs(v_ckm).tolist(),
        "jarlskog": _jarlskog(v_ckm),
        "U_u_left": u_left,
        "U_d_left": d_left,
    }


def _rms_log_error(m_u: list[float], m_d: list[float], target_u: np.ndarray, target_d: np.ndarray) -> float:
    residual = np.concatenate([np.log(np.asarray(m_u) / target_u), np.log(np.asarray(m_d) / target_d)])
    return float(np.sqrt(np.mean(residual * residual)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 quark mass-branch and CKM residual artifact.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--audit", default=str(AUDIT_JSON))
    parser.add_argument("--delta-overlap-candidate", type=float, default=0.6695617711471163 / 5.0)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward = _load_json(Path(args.forward))
    audit = _load_json(Path(args.audit))
    y_u = _complex_matrix(forward["Y_u"])
    y_d = _complex_matrix(forward["Y_d"])
    target_u = np.asarray(audit["reference_targets"]["singular_values_u"], dtype=float)
    target_d = np.asarray(audit["reference_targets"]["singular_values_d"], dtype=float)
    sigma_u = float(audit["spread_emitter_audit"]["sigma_u_total_log_per_side"])
    sigma_d = float(audit["spread_emitter_audit"]["sigma_d_total_log_per_side"])

    candidate = _apply_delta(args.delta_overlap_candidate, y_u, y_d, sigma_u, sigma_d)
    candidate["rms_log_error_vs_reference_targets"] = _rms_log_error(candidate["m_u"], candidate["m_d"], target_u, target_d)

    best: dict[str, Any] | None = None
    for delta_value in np.linspace(-0.4, 0.4, 4001):
        payload = _apply_delta(float(delta_value), y_u, y_d, sigma_u, sigma_d)
        err = _rms_log_error(payload["m_u"], payload["m_d"], target_u, target_d)
        if best is None or err < best["rms_log_error_vs_reference_targets"]:
            best = {key: value for key, value in payload.items() if key not in {"U_u_left", "U_d_left"}}
            best["rms_log_error_vs_reference_targets"] = err

    s12, s13, s23, delta_phase = 0.22501, 0.003732, 0.04183, 1.147
    c12 = math.sqrt(1.0 - s12 * s12)
    c13 = math.sqrt(1.0 - s13 * s13)
    c23 = math.sqrt(1.0 - s23 * s23)
    v_target = np.asarray(
        [
            [c12 * c13, s12 * c13, s13 * np.exp(-1j * delta_phase)],
            [
                -s12 * c23 - c12 * s23 * s13 * np.exp(1j * delta_phase),
                c12 * c23 - s12 * s23 * s13 * np.exp(1j * delta_phase),
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13 * np.exp(1j * delta_phase),
                -c12 * s23 - s12 * c23 * s13 * np.exp(1j * delta_phase),
                c23 * c13,
            ],
        ],
        dtype=complex,
    )

    u_left = candidate["U_u_left"]
    d_left = candidate["U_d_left"]
    generator_input = u_left @ v_target @ d_left.conjugate().T
    k_ckm = _principal_unitary_log(generator_input)

    result = {
        "artifact": "oph_quark_d12_mass_branch_and_ckm_residual",
        "generated_utc": _timestamp(),
        "status": "diagnostic_only_d12_continuation",
        "public_promotion_allowed": False,
        "theorem_tier": "D12_continuation_only",
        "candidate_selector_value_source": "uploaded_lane5_d12_continuation_t1_over_5_candidate",
        "candidate_mass_branch_from_t1_over_5": {
            key: value for key, value in candidate.items() if key not in {"U_u_left", "U_d_left"}
        },
        "best_honest_one_scalar_mass_point_on_same_family": best,
        "ckm_cp_residual_generator": {
            "definition": "K_CKM = log(U_u_left(candidate) @ V_CKM_target @ U_d_left(candidate)^dagger)",
            "real": np.real(k_ckm).tolist(),
            "imag": np.imag(k_ckm).tolist(),
            "abs": np.abs(k_ckm).tolist(),
            "off_diagonal_abs": {
                "12": float(abs(k_ckm[0, 1])),
                "13": float(abs(k_ckm[0, 2])),
                "23": float(abs(k_ckm[1, 2])),
            },
            "off_diagonal_phase_radians": {
                "12": float(np.angle(k_ckm[0, 1])),
                "13": float(np.angle(k_ckm[0, 2])),
                "23": float(np.angle(k_ckm[1, 2])),
            },
        },
        "notes": [
            "This artifact records the strongest current D12 continuation candidate for the light-quark split without overriding the recovered-core no-go.",
            "The mass-side candidate materially improves the light sector, but CKM and CP remain badly under-produced on this one-scalar family.",
            "The strongest smaller constructive object for CKM/CP is therefore a matrix-valued left-unitary transport generator rather than another scalar selector guess.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
            default=lambda value: value.tolist() if isinstance(value, np.ndarray) else value,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
