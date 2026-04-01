#!/usr/bin/env python3
"""Concrete helpers for the same-label left-handed quark orbit frontier.

The current local corpus still does not emit the full finite sigma_ud orbit.
But it does emit one real same-label left-handed reference-sheet evaluation on
the D12 branch. This module exposes that singleton honestly so downstream
artifacts can distinguish "no provider output at all" from the sharper remaining
blocker: no distinct same-label left-handed relative-sheet class beyond the
reference representative.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from sigma_ud_orbit_provider_interface import (
    CKMTuple,
    CanonicalToken,
    OrbitElement,
    SigmaUDOrbitProvider,
    build_sigma_ud_orbit,
)


ROOT = Path(__file__).resolve().parents[2]
FORWARD_YUKAWAS_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
D12_AUDIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"
D12_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
LOCAL_BASIS_ORBIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_local_basis_orbit_diagnostic.json"
TRANSPORT_FRAME_DIAGNOSTIC_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_transport_frame_diagnostic_orbit.json"
REFERENCE_SHEET_TOKEN = "D12::same_label_left::reference_sheet"
REFERENCE_SHEET_SIGMA_ID = "sigma_ref"
REFERENCE_FAMILY_BRANCH_KEY = ["D12", None]
REFERENCE_SELECTED_BRANCH_KEY = ["D12", REFERENCE_SHEET_SIGMA_ID]
STANDARD_GAUGE_ANCHORS = {
    "V_ud": (0, 0),
    "V_us": (0, 1),
    "V_cs": (1, 1),
    "V_cb": (1, 2),
    "V_tb": (2, 2),
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, list[list[float]]]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _jarlskog(v_ckm: np.ndarray) -> float:
    return float(np.imag(v_ckm[0, 0] * v_ckm[1, 1] * np.conjugate(v_ckm[0, 1]) * np.conjugate(v_ckm[1, 0])))


def _left_diag(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = matrix @ matrix.conjugate().T
    eig_vals, eig_vecs = np.linalg.eigh(hermitian)
    order = np.argsort(eig_vals)
    eig_vals = eig_vals[order]
    eig_vecs = eig_vecs[:, order]
    return np.sqrt(np.clip(eig_vals, 0.0, None)), eig_vecs


def _apply_delta(delta_value: float, y_u: np.ndarray, y_d: np.ndarray, sigma_u: float, sigma_d: float) -> dict[str, Any]:
    b_ord = np.asarray([-1.0, 0.0, 1.0], dtype=float)
    tau_u = 0.5 * delta_value * sigma_d / (sigma_u + sigma_d)
    tau_d = 0.5 * delta_value * sigma_u / (sigma_u + sigma_d)
    d_u = np.diag(np.exp(tau_u * b_ord))
    d_d = np.diag(np.exp(tau_d * b_ord))
    y_u_trial = d_u @ y_u @ d_u
    y_d_trial = d_d @ y_d @ d_d
    _, u_left = _left_diag(y_u_trial)
    _, d_left = _left_diag(y_d_trial)
    return {
        "U_u_left": u_left,
        "U_d_left": d_left,
        "V_CKM": u_left.conjugate().T @ d_left,
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
    cos_delta = 1.0 if denominator == 0.0 else max(-1.0, min(1.0, numerator / denominator))
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


def load_emitted_reference_sheet_evaluation(
    forward_path: Path = FORWARD_YUKAWAS_JSON,
    audit_path: Path = D12_AUDIT_JSON,
    branch_path: Path = D12_BRANCH_JSON,
) -> dict[str, Any]:
    if not forward_path.exists() or not audit_path.exists() or not branch_path.exists():
        return {
            "available": False,
            "reason_missing": "forward Yukawas, quark audit, or D12 mass-branch transport artifact has not been emitted yet",
            "required_inputs": [
                "code/particles/runs/flavor/forward_yukawas.json",
                "code/particles/runs/flavor/quark_current_family_exactness_audit.json",
                "code/particles/runs/flavor/quark_d12_mass_branch_and_ckm_residual.json",
            ],
        }

    forward = _load_json(forward_path)
    audit = _load_json(audit_path)
    branch = _load_json(branch_path)
    delta_value = float(branch["sample_same_family_point"]["Delta_ud_overlap"])
    sigma_u = float(audit["spread_emitter_audit"]["sigma_u_total_log_per_side"])
    sigma_d = float(audit["spread_emitter_audit"]["sigma_d_total_log_per_side"])
    reconstructed = _apply_delta(
        delta_value,
        _decode_complex_matrix(forward["Y_u"]),
        _decode_complex_matrix(forward["Y_d"]),
        sigma_u,
        sigma_d,
    )

    u_u_left = reconstructed["U_u_left"]
    u_d_left = reconstructed["U_d_left"]
    gauge = branch["forward_same_label_transport"]["standard_rephasing_gauge"]
    row_phases = np.asarray(gauge["row_phase_radians"], dtype=float)
    col_phases = np.asarray(gauge["col_phase_radians"], dtype=float)

    # Apply the branch's published diagonal rephasing so the emitted pair and
    # V_CKM live on the same standard-gauge representative.
    u_u_left = u_u_left @ np.diag(np.exp(1j * row_phases))
    u_d_left = u_d_left @ np.diag(np.exp(1j * col_phases))
    v_ckm = u_u_left.conjugate().T @ u_d_left
    published_v_ckm = _decode_complex_matrix(branch["forward_same_label_transport"]["V_CKM_forward_standard_gauge"])
    consistency_residual = float(np.linalg.norm(v_ckm - published_v_ckm, ord="fro"))
    if consistency_residual > 1.0e-10:
        raise ValueError("reference-sheet gauge reconstruction drifted from published D12 standard CKM representative")

    ckm = branch.get("standard_ckm_parameters") or _standard_ckm_parameters(v_ckm)
    return {
        "available": True,
        "sigma_id": REFERENCE_SHEET_SIGMA_ID,
        "canonical_token": REFERENCE_SHEET_TOKEN,
        "branch_key": list(REFERENCE_FAMILY_BRANCH_KEY),
        "provider_status": "singleton_reference_sheet_only",
        "selection_proof": {
            "theorem_grade_select": False,
            "reason": (
                "This is the emitted D12 reference-sheet representative only. No distinct same-label left-handed "
                "relative-sheet class or intrinsic uniqueness theorem is exposed yet."
            ),
        },
        "coverage_status": "reference_sheet_representative_only",
        "U_u_left": _encode_complex_matrix(u_u_left),
        "U_d_left": _encode_complex_matrix(u_d_left),
        "V_CKM": _encode_complex_matrix(v_ckm),
        "ckm_invariants": {
            "theta_12": float(ckm["theta_12"]),
            "theta_23": float(ckm["theta_23"]),
            "theta_13": float(ckm["theta_13"]),
            "delta_ckm": float(ckm["delta_ckm"]),
            "jarlskog": float(ckm["jarlskog"]),
        },
        "gauge_fix": {
            "kind": "published_d12_standard_rephasing_gauge",
            "row_phase_radians": gauge["row_phase_radians"],
            "col_phase_radians": gauge["col_phase_radians"],
            "consistency_residual_fro_norm": consistency_residual,
        },
        "source_artifacts": {
            "forward_yukawas": "code/particles/runs/flavor/forward_yukawas.json",
            "quark_current_family_exactness_audit": "code/particles/runs/flavor/quark_current_family_exactness_audit.json",
            "d12_mass_branch_and_ckm_residual": "code/particles/runs/flavor/quark_d12_mass_branch_and_ckm_residual.json",
        },
    }


def build_emitted_reference_sheet_orbit_elements() -> list[dict[str, Any]]:
    payload = load_emitted_reference_sheet_evaluation()
    if not payload.get("available"):
        return []

    uniqueness = load_sigma_ud_singleton_uniqueness_witness()
    provider = CurrentReferenceSheetSigmaUDProvider(payload)
    return [
        {
            "sigma_id": item.sigma_id,
            "canonical_token": item.canonical_token.token,
            "branch_key": (
                uniqueness["selected_sigma"]["branch_key"]
                if bool(uniqueness.get("theorem_grade_select"))
                else payload["branch_key"]
            ),
            "family_branch_key": (
                payload["branch_key"] if bool(uniqueness.get("theorem_grade_select")) else None
            ),
            "coverage_status": payload["coverage_status"],
            "selection_proof": {
                "theorem_grade_select": bool(uniqueness.get("theorem_grade_select")),
                "reason": (
                    str(uniqueness["selection_reason"])
                    if bool(uniqueness.get("theorem_grade_select"))
                    else payload["selection_proof"]["reason"]
                ),
                "theorem_artifact": "code/particles/runs/flavor/quark_sigma_ud_singleton_uniqueness_theorem.json",
            },
            "provider_status": payload["provider_status"],
            "U_u_left": payload["U_u_left"],
            "U_d_left": payload["U_d_left"],
            "V_CKM": payload["V_CKM"],
            "ckm_invariants": {
                "theta_12": item.ckm_invariants.theta_12,
                "theta_23": item.ckm_invariants.theta_23,
                "theta_13": item.ckm_invariants.theta_13,
                "delta_ckm": item.ckm_invariants.delta_ckm,
                "jarlskog": item.ckm_invariants.jarlskog,
            },
            "source_artifacts": payload["source_artifacts"],
            "gauge_fix": payload["gauge_fix"],
        }
        for item in build_sigma_ud_orbit(provider)
    ]


def load_already_local_diagnostic_orbit(path: Path = LOCAL_BASIS_ORBIT_JSON) -> dict[str, Any]:
    if not path.exists():
        return {
            "available": False,
            "artifact": "code/particles/runs/flavor/quark_local_basis_orbit_diagnostic.json",
            "reason_missing": "already-local chirality-basis diagnostic has not been emitted yet",
        }

    payload = _load_json(path)
    physical_reference = next(
        (
            item
            for item in payload.get("elements", [])
            if item.get("basis_u") == "L" and item.get("basis_d") == "L"
        ),
        None,
    )
    return {
        "available": True,
        "artifact": "code/particles/runs/flavor/quark_local_basis_orbit_diagnostic.json",
        "scope": payload.get("scope"),
        "theorem_use": payload.get("theorem_use"),
        "physical_reference_element": physical_reference,
        "best_nonphysical_candidate": payload.get("best_nonphysical_candidate"),
        "notes": payload.get("notes", []),
    }


def load_transport_frame_diagnostic_orbit(path: Path = TRANSPORT_FRAME_DIAGNOSTIC_JSON) -> dict[str, Any]:
    if not path.exists():
        return {
            "available": False,
            "artifact": "code/particles/runs/flavor/quark_transport_frame_diagnostic_orbit.json",
            "reason_missing": "transport-frame diagnostic orbit has not been emitted yet",
        }

    payload = _load_json(path)
    return {
        "available": True,
        "artifact": "code/particles/runs/flavor/quark_transport_frame_diagnostic_orbit.json",
        "self_overlap_symbol": ((payload.get("self_overlap") or {}).get("symbol")),
        "ckm_invariants": payload.get("ckm_invariants"),
        "debug_log_shell_loss": payload.get("debug_log_shell_loss"),
        "missing_sector_attachment": payload.get("missing_sector_attachment"),
        "why_not_promotable": payload.get("why_not_promotable"),
    }


def load_sigma_ud_singleton_uniqueness_witness(
    local_basis_path: Path = LOCAL_BASIS_ORBIT_JSON,
    branch_path: Path = D12_BRANCH_JSON,
) -> dict[str, Any]:
    if not local_basis_path.exists() or not branch_path.exists():
        return {
            "available": False,
            "proof_status": "supporting_artifacts_missing",
            "reason_missing": "local basis diagnostic or D12 branch artifact has not been emitted yet",
        }

    local_basis = _load_json(local_basis_path)
    branch = _load_json(branch_path)
    physical = [item for item in local_basis.get("elements", []) if bool(item.get("physical_admissible"))]
    only_left_left_survives = (
        len(physical) == 1
        and physical[0].get("basis_u") == "L"
        and physical[0].get("basis_d") == "L"
    )

    v_standard = _decode_complex_matrix(branch["forward_same_label_transport"]["V_CKM_forward_standard_gauge"])
    anchor_checks = []
    equations = []
    for name, (row_idx, col_idx) in STANDARD_GAUGE_ANCHORS.items():
        value = complex(v_standard[row_idx, col_idx])
        anchor_checks.append(
            {
                "name": name,
                "position": [row_idx, col_idx],
                "abs": float(abs(value)),
                "real": float(np.real(value)),
                "imag": float(np.imag(value)),
                "is_real_positive": abs(np.imag(value)) < 1.0e-12 and float(np.real(value)) > 0.0,
            }
        )
        row = np.zeros(6, dtype=float)
        row[row_idx] = -1.0
        row[3 + col_idx] = 1.0
        equations.append(row)

    constraint_matrix = np.vstack(equations)
    phase_constraint_rank = int(np.linalg.matrix_rank(constraint_matrix))
    phase_nullity = int(constraint_matrix.shape[1] - phase_constraint_rank)
    global_phase_vector = np.ones(6, dtype=float)
    global_phase_is_kernel = bool(np.allclose(constraint_matrix @ global_phase_vector, 0.0, atol=1.0e-12))
    standard_gauge_representative_unique = (
        all(bool(item["is_real_positive"]) for item in anchor_checks)
        and phase_nullity == 1
        and global_phase_is_kernel
    )

    theorem_grade_select = bool(only_left_left_survives and standard_gauge_representative_unique)
    return {
        "available": True,
        "proof_status": (
            "same_label_left_handed_local_orbit_singleton_closed"
            if theorem_grade_select
            else "same_label_left_handed_local_orbit_singleton_not_closed"
        ),
        "theorem_grade_select": theorem_grade_select,
        "scope": "emitted_local_solver_surface_only",
        "selected_sigma": {
            "sigma_id": REFERENCE_SHEET_SIGMA_ID,
            "canonical_token": REFERENCE_SHEET_TOKEN,
            "branch_key": list(REFERENCE_SELECTED_BRANCH_KEY),
        },
        "theorem_statement": (
            "On the emitted local solver surface, the only currently exposed same-label left-handed orbit representative is "
            "sigma_ref. The local chirality-basis diagnostic leaves only the ordered L/L choice physically admissible, and "
            "the published standard CKM gauge with V_ud, V_us, V_cs, V_cb, V_tb real-positive fixes the remaining diagonal "
            "rephasing to the trivial global phase. A separate ordered-basis-to-sheet identification witness is not yet "
            "surfaced beyond that emitted reference representative."
        ),
        "selection_reason": (
            "The current emitted local same-label left-handed solver surface exposes only sigma_ref, the chirality-basis "
            "diagnostic leaves only the ordered L/L side physically admissible, and the published five-anchor standard CKM "
            "gauge removes residual diagonal rephasing except the trivial global phase."
            if theorem_grade_select
            else "Either a non-left chirality basis remains admissible or the published five-anchor standard CKM gauge does not uniquely fix the reference representative."
        ),
        "local_basis_admissibility": {
            "physical_admissible_count": len(physical),
            "only_left_left_survives": only_left_left_survives,
            "physical_reference": physical[0] if physical else None,
            "artifact": "code/particles/runs/flavor/quark_local_basis_orbit_diagnostic.json",
        },
        "standard_rephasing_gauge_uniqueness": {
            "anchors": anchor_checks,
            "phase_constraint_rank": phase_constraint_rank,
            "phase_nullity": phase_nullity,
            "global_phase_is_kernel": global_phase_is_kernel,
            "standard_gauge_representative_unique": standard_gauge_representative_unique,
            "artifact": "code/particles/runs/flavor/quark_d12_mass_branch_and_ckm_residual.json",
        },
        "limitations": [
            "This closes the emitted local same-label left-handed orbit only.",
            "It does not manufacture any nonlocal sheet provider beyond the current solver surface.",
            "It does not yet expose a separate ordered-basis-to-sheet identification witness between the L/L diagnostic basis and the emitted sigma_ref representative.",
            "The selected singleton remains the current D12 reference sheet, so the physical CKM-shell no-go still applies.",
        ],
    }


def build_sigma_ud_provider_frontier(path: Path = LOCAL_BASIS_ORBIT_JSON) -> dict[str, Any]:
    diagnostic = load_already_local_diagnostic_orbit(path)
    transport_frame = load_transport_frame_diagnostic_orbit()
    reference_sheet = load_emitted_reference_sheet_evaluation()
    uniqueness = load_sigma_ud_singleton_uniqueness_witness(path, D12_BRANCH_JSON)
    if bool(uniqueness.get("theorem_grade_select")) and reference_sheet.get("available"):
        return {
            "status": "same_label_left_handed_local_orbit_singleton_closed",
            "required_provider_methods": [
                "enumerate_relative_sheets_d12()",
                "evaluate_relative_sheet(token)",
            ],
            "blocking_reason": None,
            "smallest_missing_runtime_object": None,
            "orbit_completion_status": "singleton_closed_by_uniqueness_theorem",
            "emitted_reference_sheet": {
                "sigma_id": reference_sheet["sigma_id"],
                "canonical_token": reference_sheet["canonical_token"],
                "branch_key": list(uniqueness["selected_sigma"]["branch_key"]),
                "family_branch_key": reference_sheet["branch_key"],
                "coverage_status": reference_sheet["coverage_status"],
                "ckm_invariants": reference_sheet["ckm_invariants"],
            },
            "uniqueness_theorem": {
                "proof_status": uniqueness["proof_status"],
                "selection_reason": uniqueness["selection_reason"],
                "scope": uniqueness["scope"],
            },
            "next_exact_object_after_orbit_closure": "quark_d12_t1_value_law",
            "already_local_diagnostic_orbit_available": bool(diagnostic.get("available")),
            "already_local_diagnostic_orbit_artifact": diagnostic.get("artifact"),
            "transport_frame_diagnostic_orbit_available": bool(transport_frame.get("available")),
            "transport_frame_diagnostic_orbit_artifact": transport_frame.get("artifact"),
            "transport_frame_diagnostic_orbit_loss": (
                (transport_frame.get("debug_log_shell_loss") or {}).get("transport_frame_self_overlap")
            ),
            "why_transport_frame_diagnostic_orbit_is_insufficient": (
                transport_frame.get("why_not_promotable")
                or "The transport-frame witness is compare-only and not sector-attached."
            ),
        }
    if reference_sheet.get("available"):
        return {
            "status": "reference_sheet_singleton_emitted_missing_distinct_same_label_sheet",
            "required_provider_methods": [
                "enumerate_relative_sheets_d12()",
                "evaluate_relative_sheet(token)",
            ],
            "blocking_reason": (
                "The current local corpus emits only one same-label left-handed D12 reference-sheet evaluation. "
                "No distinct relative-sheet token or non-reference sigma -> CKM evaluation is exposed yet."
            ),
            "smallest_missing_runtime_object": (
                "one additional same-label left-handed relative-sheet evaluation not gauge-equivalent to the emitted "
                "D12 reference representative, or an intrinsic uniqueness theorem collapsing the orbit to that singleton"
            ),
            "orbit_completion_status": "singleton_reference_sheet_only",
            "emitted_reference_sheet": {
                "sigma_id": reference_sheet["sigma_id"],
                "canonical_token": reference_sheet["canonical_token"],
                "branch_key": reference_sheet["branch_key"],
                "coverage_status": reference_sheet["coverage_status"],
                "ckm_invariants": reference_sheet["ckm_invariants"],
            },
            "already_local_diagnostic_orbit_available": bool(diagnostic.get("available")),
            "already_local_diagnostic_orbit_artifact": diagnostic.get("artifact"),
            "transport_frame_diagnostic_orbit_available": bool(transport_frame.get("available")),
            "transport_frame_diagnostic_orbit_artifact": transport_frame.get("artifact"),
            "transport_frame_diagnostic_orbit_loss": (
                (transport_frame.get("debug_log_shell_loss") or {}).get("transport_frame_self_overlap")
            ),
            "why_diagnostic_orbit_is_insufficient": (
                "The already-local chirality-basis orbit can move CKM moduli compare-only, but it is not the same-label "
                "left-handed relative-sheet orbit and therefore cannot emit sigma_ud."
            ),
            "why_transport_frame_diagnostic_orbit_is_insufficient": (
                transport_frame.get("why_not_promotable")
                or "The transport-frame witness is compare-only and not sector-attached."
            ),
        }

    return {
        "status": "same_label_left_handed_sigma_ud_orbit_unemitted",
        "required_provider_methods": [
            "enumerate_relative_sheets_d12()",
            "evaluate_relative_sheet(token)",
        ],
        "blocking_reason": (
            "No emitted same-label left-handed relative-sheet enumerator or sigma-to-CKM evaluator exists on the current local corpus."
        ),
        "smallest_missing_runtime_object": (
            "first non-empty provider output Sigma_ud_orbit.elements = "
            "[{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM, ckm_invariants}]"
        ),
        "already_local_diagnostic_orbit_available": bool(diagnostic.get("available")),
        "already_local_diagnostic_orbit_artifact": diagnostic.get("artifact"),
        "transport_frame_diagnostic_orbit_available": bool(transport_frame.get("available")),
        "transport_frame_diagnostic_orbit_artifact": transport_frame.get("artifact"),
        "transport_frame_diagnostic_orbit_loss": (
            (transport_frame.get("debug_log_shell_loss") or {}).get("transport_frame_self_overlap")
        ),
        "why_diagnostic_orbit_is_insufficient": (
            "The already-local chirality-basis orbit can move CKM moduli compare-only, but it is not the same-label left-handed relative-sheet orbit and therefore cannot emit sigma_ud."
        ),
        "why_transport_frame_diagnostic_orbit_is_insufficient": (
            transport_frame.get("why_not_promotable")
            or "The transport-frame witness is compare-only and not sector-attached."
        ),
    }


class MissingSigmaUDOrbitProvider(SigmaUDOrbitProvider):
    """Explicit placeholder until the live solver can emit same-label left-handed orbit elements."""

    def enumerate_relative_sheets_d12(self) -> Sequence[CanonicalToken]:
        raise NotImplementedError(
            "Current local corpus still lacks an emitted finite same-label left-handed Sigma_ud orbit."
        )

    def evaluate_relative_sheet(self, token: CanonicalToken) -> OrbitElement:
        raise NotImplementedError(
            "Current local corpus still lacks a concrete same-label left-handed sigma -> CKM evaluator."
        )


class CurrentReferenceSheetSigmaUDProvider(SigmaUDOrbitProvider):
    """Emit the one already-local same-label left-handed D12 reference-sheet evaluation."""

    def __init__(self, payload: dict[str, Any] | None = None) -> None:
        self._payload = payload or load_emitted_reference_sheet_evaluation()
        if not self._payload.get("available"):
            raise ValueError("reference-sheet provider requires emitted forward Yukawa and D12 transport artifacts")

    def enumerate_relative_sheets_d12(self) -> Sequence[CanonicalToken]:
        return [CanonicalToken(token=self._payload["canonical_token"])]

    def evaluate_relative_sheet(self, token: CanonicalToken) -> OrbitElement:
        if token.token != self._payload["canonical_token"]:
            raise KeyError(f"unknown D12 reference-sheet token: {token.token}")
        ckm = self._payload["ckm_invariants"]
        return OrbitElement(
            sigma_id=self._payload["sigma_id"],
            canonical_token=CanonicalToken(token=self._payload["canonical_token"]),
            U_u_left=self._payload["U_u_left"],
            U_d_left=self._payload["U_d_left"],
            V_CKM=self._payload["V_CKM"],
            ckm_invariants=CKMTuple(
                theta_12=float(ckm["theta_12"]),
                theta_23=float(ckm["theta_23"]),
                theta_13=float(ckm["theta_13"]),
                delta_ckm=float(ckm["delta_ckm"]),
                jarlskog=float(ckm["jarlskog"]),
            ),
        )
