#!/usr/bin/env python3
"""Analyze the strongest current D12 quark mass branch and CKM closure.

Chain role: keep the honest D12 quark continuation program explicit without
promoting the light-quark selector value to recovered-core status.

Mathematics: apply a one-scalar D12 light-quark overlap selector candidate to
the current forward Yukawas, evaluate the resulting mass branch, and read off
the same-label left-transport unitary `V_CKM^fwd = U_u^dagger U_d` together
with its principal anti-Hermitian logarithm.

OPH-derived inputs: the current forward Yukawa artifact, the quark exactness
audit, and the already-emitted spread package.

Output: a D12 continuation artifact carrying the strongest current mass-side
sample point on the emitted D12 ray, the honest forward CKM transport unitary,
its principal generator, and the sharper branch-repair boundary beneath the
physical CKM shell.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from sigma_ud_orbit_provider import load_sigma_ud_singleton_uniqueness_witness


ROOT = Path(__file__).resolve().parents[2]
FORWARD_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
AUDIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
TARGET_THETA_12 = 0.2256
TARGET_THETA_23 = 0.0438
TARGET_THETA_13 = 0.00347


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _left_diag(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = matrix @ matrix.conjugate().T
    eig_vals, eig_vecs = np.linalg.eigh(hermitian)
    order = np.argsort(eig_vals)
    eig_vals = eig_vals[order]
    eig_vecs = eig_vecs[:, order]
    return np.sqrt(np.clip(eig_vals, 0.0, None)), eig_vecs


def _jarlskog(v_ckm: np.ndarray) -> float:
    return float(np.imag(v_ckm[0, 0] * v_ckm[1, 1] * np.conjugate(v_ckm[0, 1]) * np.conjugate(v_ckm[1, 0])))


def _matrix_exp(matrix: np.ndarray) -> np.ndarray:
    eig_vals, eig_vecs = np.linalg.eig(matrix)
    return eig_vecs @ np.diag(np.exp(eig_vals)) @ np.linalg.inv(eig_vecs)


def _principal_unitary_log(matrix: np.ndarray) -> np.ndarray:
    eig_vals, eig_vecs = np.linalg.eig(matrix)
    phases = np.angle(eig_vals)
    log_diag = np.diag(1j * phases)
    raw = eig_vecs @ log_diag @ np.linalg.inv(eig_vecs)
    return 0.5 * (raw - raw.conjugate().T)


def _standard_ckm_gauge_fix(v_ckm: np.ndarray) -> dict[str, Any]:
    """Choose a standard representative with V_ud, V_us, V_cs, V_cb, V_tb real-positive."""
    row_phases = np.zeros(3, dtype=float)
    col_phases = np.zeros(3, dtype=float)
    row_phases[0] = float(np.angle(v_ckm[0, 0]))
    col_phases[1] = float(row_phases[0] - np.angle(v_ckm[0, 1]))
    col_phases[2] = float(np.angle(v_ckm[1, 1]) + col_phases[1] - np.angle(v_ckm[1, 2]))
    row_phases[1] = float(np.angle(v_ckm[1, 2]) + col_phases[2])
    row_phases[2] = float(np.angle(v_ckm[2, 2]) + col_phases[2])
    left = np.diag(np.exp(-1j * row_phases))
    right = np.diag(np.exp(1j * col_phases))
    return {
        "matrix": left @ v_ckm @ right,
        "row_phase_radians": row_phases.tolist(),
        "col_phase_radians": col_phases.tolist(),
        "conditions": [
            "V_ud real-positive",
            "V_us real-positive",
            "V_cs real-positive",
            "V_cb real-positive",
            "V_tb real-positive",
        ],
    }


def _generator_gauge_fix(generator: np.ndarray, matrix: np.ndarray) -> dict[str, Any]:
    """Diagonal-conjugate so K_12 and K_23 are real-positive."""
    phase_shifts = np.zeros(3, dtype=float)
    phase_shifts[1] = float(-np.angle(generator[0, 1]))
    stage = np.diag(np.exp(-1j * phase_shifts)) @ generator @ np.diag(np.exp(1j * phase_shifts))
    phase_shifts[2] = float(phase_shifts[1] - np.angle(stage[1, 2]))
    left = np.diag(np.exp(-1j * phase_shifts))
    right = np.diag(np.exp(1j * phase_shifts))
    return {
        "matrix": left @ generator @ right,
        "matrix_surface": left @ matrix @ right,
        "phase_shifts_radians": phase_shifts.tolist(),
    }


def _standard_ckm_parameters(v_ckm: np.ndarray) -> dict[str, float]:
    s13 = min(1.0, max(0.0, float(abs(v_ckm[0, 2]))))
    c13 = math.sqrt(max(0.0, 1.0 - s13 * s13))
    s12 = min(1.0, max(0.0, float(abs(v_ckm[0, 1]) / c13)))
    s23 = min(1.0, max(0.0, float(abs(v_ckm[1, 2]) / c13)))
    theta_12 = math.asin(s12)
    theta_23 = math.asin(s23)
    theta_13 = math.asin(s13)
    c12 = math.sqrt(max(0.0, 1.0 - s12 * s12))
    c23 = math.sqrt(max(0.0, 1.0 - s23 * s23))
    numerator = (s12 * s23) ** 2 + (c12 * c23 * s13) ** 2 - float(abs(v_ckm[2, 0])) ** 2
    denominator = 2.0 * s12 * s23 * c12 * c23 * s13
    cos_delta = max(-1.0, min(1.0, numerator / denominator))
    delta = math.acos(cos_delta)
    jarlskog = _jarlskog(v_ckm)
    if jarlskog < 0.0:
        delta = 2.0 * math.pi - delta
    return {
        "theta_12": theta_12,
        "theta_23": theta_23,
        "theta_13": theta_13,
        "delta_ckm": delta,
        "jarlskog": jarlskog,
    }


def _j_max(theta_12: float, theta_23: float, theta_13: float) -> float:
    s12 = math.sin(theta_12)
    s23 = math.sin(theta_23)
    s13 = math.sin(theta_13)
    c12 = math.cos(theta_12)
    c23 = math.cos(theta_23)
    c13 = math.cos(theta_13)
    return c12 * c23 * (c13**2) * s12 * s23 * s13


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
        "V_CKM_forward": v_ckm,
        "U_u_left": u_left,
        "U_d_left": d_left,
    }


def _rms_log_error(m_u: list[float], m_d: list[float], target_u: np.ndarray, target_d: np.ndarray) -> float:
    residual = np.concatenate([np.log(np.asarray(m_u) / target_u), np.log(np.asarray(m_d) / target_d)])
    return float(np.sqrt(np.mean(residual * residual)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 quark mass-branch and CKM closure artifact.")
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
            best = {
                key: value
                for key, value in payload.items()
                if key not in {"U_u_left", "U_d_left", "V_CKM_forward"}
            }
            best["rms_log_error_vs_reference_targets"] = err

    v_forward = candidate["V_CKM_forward"]
    standard_gauge = _standard_ckm_gauge_fix(v_forward)
    v_standard = standard_gauge["matrix"]
    eigenphase_radians = [float(value) for value in np.angle(np.linalg.eigvals(v_standard)).tolist()]
    k_principal = _principal_unitary_log(v_standard)
    generator_gauge = _generator_gauge_fix(k_principal, v_standard)
    k_ckm = generator_gauge["matrix"]
    v_generator_surface = generator_gauge["matrix_surface"]
    generator_invariants = {
        "theta_12_K": float(abs(k_ckm[0, 1])),
        "theta_23_K": float(abs(k_ckm[1, 2])),
        "theta_13_K": float(abs(k_ckm[0, 2])),
        "phi_K": float(np.angle(k_ckm[0, 2])),
        "chi_diagonal_imag": [float(value) for value in np.imag(np.diag(k_ckm)).tolist()],
    }
    standard_parameters = _standard_ckm_parameters(v_standard)
    closure_residual = float(np.linalg.norm(_matrix_exp(k_ckm) - v_generator_surface, ord="fro"))
    uniqueness = load_sigma_ud_singleton_uniqueness_witness()
    selector_value = uniqueness.get("selected_sigma") if bool(uniqueness.get("theorem_grade_select")) else None
    theta_12 = standard_parameters["theta_12"]
    theta_23 = standard_parameters["theta_23"]
    theta_13 = standard_parameters["theta_13"]
    j_max = _j_max(theta_12, theta_23, theta_13)

    result = {
        "artifact": "oph_quark_d12_mass_branch_and_ckm_closure",
        "generated_utc": _timestamp(),
        "status": "d12_current_sheet_ckm_cp_transport_closed_wrong_branch_no_go",
        "public_promotion_allowed": False,
        "theorem_tier": "D12_continuation_only",
        "branch_key": ["D12", None],
        "quark_relative_sheet_selector": selector_value,
        "current_sheet_status": "single_local_reference_sheet_only",
        "physical_branch_status": "current_d12_sheet_is_strict_no_go_for_physical_ckm_shell",
        "sample_selector_value_source": "sample_same_family_point_on_D12_ud_mass_ray",
        "sample_same_family_point": {
            key: value
            for key, value in candidate.items()
            if key not in {"U_u_left", "U_d_left", "V_CKM_forward"}
        },
        "same_sheet_mass_comparison_scan": {
            "status": "comparison_only_not_a_relative_sheet_scan",
            "grid_kind": "np.linspace(-0.4, 0.4, 4001)",
            "scan_coordinate": "Delta_ud_overlap",
            "score": "RMS log error against reference_targets",
            "same_sheet_only": True,
            "uses_reference_targets": True,
            "disqualified_for_sigma_ud_selection": True,
            "reason": "This finite scan varies only the same-sheet mass-side overlap coordinate against target masses, so it cannot enumerate Sigma_ud or emit quark_relative_sheet_selector.",
        },
        "comparison_only_best_same_family_point": dict(best, status="comparison_only_not_promotable"),
        "forward_same_label_transport": {
            "definition": "V_CKM^fwd = U_u_left(candidate)^dagger @ U_d_left(candidate)",
            "standard_rephasing_gauge": {
                "conditions": standard_gauge["conditions"],
                "row_phase_radians": standard_gauge["row_phase_radians"],
                "col_phase_radians": standard_gauge["col_phase_radians"],
            },
            "V_CKM_forward_standard_gauge": _encode_complex_matrix(v_standard),
            "abs_V_CKM": np.abs(v_standard).tolist(),
            "eigenphase_principal_strip_radians": eigenphase_radians,
            "principal_log_exists_uniquely": all(abs(value) < math.pi for value in eigenphase_radians),
            "jarlskog": standard_parameters["jarlskog"],
        },
        "same_label_transport_generator": {
            "definition": "K_CKM = Log_pr(V_CKM^fwd)",
            "generator_gauge_convention": "diagonal conjugation chosen so that K_12 and K_23 are real-positive",
            "generator_gauge_phase_radians": generator_gauge["phase_shifts_radians"],
            "V_CKM_generator_gauge": _encode_complex_matrix(v_generator_surface),
            "real": np.real(k_ckm).tolist(),
            "imag": np.imag(k_ckm).tolist(),
            "abs": np.abs(k_ckm).tolist(),
            "generator_invariants": generator_invariants,
        },
        "standard_ckm_parameters": standard_parameters,
        "physical_ckm_comparison_shell": {
            "theta_12": TARGET_THETA_12,
            "theta_23": TARGET_THETA_23,
            "theta_13": TARGET_THETA_13,
            "absolute_misses": {
                "theta_12": TARGET_THETA_12 - theta_12,
                "theta_23": TARGET_THETA_23 - theta_23,
                "theta_13": TARGET_THETA_13 - theta_13,
            },
            "undershoot_factors": {
                "theta_12": TARGET_THETA_12 / theta_12 if theta_12 > 0.0 else None,
                "theta_23": TARGET_THETA_23 / theta_23 if theta_23 > 0.0 else None,
                "theta_13": TARGET_THETA_13 / theta_13 if theta_13 > 0.0 else None,
            },
            "jarlskog_fraction_of_max_allowed_by_current_angles": (
                abs(standard_parameters["jarlskog"]) / j_max if j_max > 0.0 else None
            ),
            "selection_loss": (
                (theta_12 - TARGET_THETA_12) ** 2
                + (theta_23 - TARGET_THETA_23) ** 2
                + (theta_13 - TARGET_THETA_13) ** 2
            ),
        },
        "closure_residual": {
            "definition": "||exp(K_CKM) - V_CKM^fwd||_F on the generator-gauge surface",
            "fro_norm": closure_residual,
        },
        "remaining_open_objects": (
            [
                "intrinsic_scale_law_D12",
                "quark_exact_mean_split_value_law_or_carrier_repair",
            ]
            if selector_value is not None
            else [
                "quark_relative_sheet_selector",
                "intrinsic_scale_law_D12",
                "quark_exact_mean_split_value_law_or_carrier_repair",
            ]
        ),
        "debug_only_target_seeded_generator": {
            "status": "retired",
            "reason": "the honest same-label transport unitary is emitted directly by the forward Yukawa step and no target CKM seed is needed",
        },
        "notes": [
            "This artifact records the strongest current D12 continuation sample point for the light-quark split without overriding the recovered-core no-go.",
            "On the D12 continuation branch the CKM/CP lane closes honestly once the forward Yukawa step is reached, because the same-label transport unitary is already V_CKM^fwd = U_u^dagger U_d.",
            "But the current D12 sheet is not the physical quark branch: the selector is already fixed to sigma_ref on the emitted same-label left-handed solver surface, while the current theorem-grade rigidity result is only same-sheet rephasing rigidity and does not prove CKM rigidity under ray_modulus motion on D12_ud_mass_ray. The emitted angles on this sheet still undershoot the comparison shell substantially.",
            "The only finite local scan on disk is a same-sheet Delta_ud_overlap scan against reference targets; it is comparison-only and cannot be repurposed as a Sigma_ud orbit scan.",
            (
                "The solver-side same-label left-handed orbit now closes to the singleton sigma_ref, so the exact next object is the intrinsic scale law on D12_ud_mass_ray; that selected branch still inherits the current CKM-shell no-go."
                if selector_value is not None
                else "The exact next object is therefore one discrete quark_relative_sheet_selector; mass-side scale fixing on the selected branch remains a separate issue after that branch shift."
            ),
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
