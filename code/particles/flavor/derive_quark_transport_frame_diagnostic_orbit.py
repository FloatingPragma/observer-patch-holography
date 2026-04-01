#!/usr/bin/env python3
"""Emit the compare-only transport-frame diagnostic orbit for the quark lane.

Chain role: expose the strongest already-local non-sector-attached CKM-like
signal downstream of the common-refinement line-lift, without promoting it to a
same-label left-handed relative-sheet orbit element.

Mathematics: extract deterministic source/target transport frames from the
rank-one same-label projectors on the line-lift artifact, phase-lock each
target line by the corresponding same-label transport map, then evaluate the
gauge-invariant self-overlap ``F0^dagger F1`` as a compare-only CKM-like
matrix.

OPH-derived inputs: the projective same-label common-refinement line-lift and
the current D12 quark transport closure artifact.

Output: a compare-only transport-frame diagnostic orbit witness together with
its exact non-promotion reason on the sigma_ud frontier.
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
DEFAULT_LINE_LIFT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
DEFAULT_TRANSPORT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_transport_frame_diagnostic_orbit.json"
TARGET = {
    "theta_12": 0.2256,
    "theta_23": 0.0438,
    "theta_13": 0.00347,
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _principal_vector_from_projector(projector: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(projector)
    vector = eigenvectors[:, int(np.argmax(eigenvalues))]
    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 1.0e-15:
        vector = vector * np.exp(-1j * np.angle(vector[pivot]))
    return vector


def _align_target_phase(source_vec: np.ndarray, target_vec: np.ndarray, transport_map: np.ndarray) -> np.ndarray:
    image = transport_map @ source_vec
    amplitude = np.vdot(target_vec, image)
    if abs(amplitude) > 1.0e-15:
        target_vec = target_vec * np.exp(1j * np.angle(amplitude))
    return target_vec


def _ckm_tuple(matrix: np.ndarray) -> dict[str, float]:
    s13 = min(1.0, max(0.0, float(abs(matrix[0, 2]))))
    c13 = math.sqrt(max(0.0, 1.0 - s13 * s13))
    s12 = min(1.0, max(0.0, float(abs(matrix[0, 1]) / c13)))
    s23 = min(1.0, max(0.0, float(abs(matrix[1, 2]) / c13)))

    theta12 = math.asin(s12)
    theta23 = math.asin(s23)
    theta13 = math.asin(s13)

    c12 = math.sqrt(max(0.0, 1.0 - s12 * s12))
    c23 = math.sqrt(max(0.0, 1.0 - s23 * s23))
    numerator = (s12 * s23) ** 2 + (c12 * c23 * s13) ** 2 - float(abs(matrix[2, 0])) ** 2
    denominator = 2.0 * s12 * s23 * c12 * c23 * s13
    cos_delta = 1.0 if denominator == 0.0 else max(-1.0, min(1.0, numerator / denominator))
    delta = math.acos(cos_delta)
    jarlskog = float(np.imag(matrix[0, 0] * matrix[1, 1] * np.conj(matrix[0, 1]) * np.conj(matrix[1, 0])))
    if jarlskog < 0.0:
        delta = 2.0 * math.pi - delta

    return {
        "theta_12": theta12,
        "theta_23": theta23,
        "theta_13": theta13,
        "delta_ckm": delta,
        "jarlskog": jarlskog,
    }


def _debug_log_shell_loss(ckm: dict[str, Any]) -> float:
    return sum(math.log(float(ckm[key]) / TARGET[key]) ** 2 for key in TARGET)


def build_artifact(line_lift: dict[str, Any], transport: dict[str, Any]) -> dict[str, Any]:
    items = list(line_lift.get("transport_partial_isometry_by_label_and_refinement_pair") or [])
    if len(items) != 3:
        raise ValueError("line-lift artifact must expose exactly three same-label transport partial isometries")

    source_frame = []
    target_frame = []
    per_label = []
    for item in items:
        source_projector = _decode_complex_matrix(item["source_projector"])
        target_projector = _decode_complex_matrix(item["target_projector"])
        transport_map = _decode_complex_matrix(item["transport_map"])

        source_vector = _principal_vector_from_projector(source_projector)
        target_vector = _principal_vector_from_projector(target_projector)
        target_vector = _align_target_phase(source_vector, target_vector, transport_map)

        source_frame.append(source_vector)
        target_frame.append(target_vector)
        per_label.append(
            {
                "label": item["label"],
                "source_overlap_norm": float(np.linalg.norm(source_projector @ source_vector)),
                "target_overlap_norm": float(np.linalg.norm(target_projector @ target_vector)),
            }
        )

    frame_source = np.column_stack(source_frame)
    frame_target = np.column_stack(target_frame)
    self_overlap = frame_source.conj().T @ frame_target
    ckm = _ckm_tuple(self_overlap)

    current_same_sheet = dict(transport["standard_ckm_parameters"])
    current_same_sheet_loss = _debug_log_shell_loss(current_same_sheet)
    transport_frame_loss = _debug_log_shell_loss(ckm)

    return {
        "artifact": "oph_quark_transport_frame_diagnostic_orbit",
        "generated_utc": _timestamp(),
        "proof_status": "compare_only_not_sector_attached",
        "public_promotion_allowed": False,
        "scope": "common_refinement_transport_frame_only",
        "input_artifacts": {
            "line_lift": "code/particles/runs/flavor/overlap_edge_line_lift.json",
            "d12_transport": "code/particles/runs/flavor/quark_d12_mass_branch_and_ckm_residual.json",
        },
        "line_lift_labels": [item["label"] for item in items],
        "frame_construction": {
            "source_frame": "principal rank-one source-projector eigenvectors ordered by same-label line-lift labels",
            "target_frame": "principal rank-one target-projector eigenvectors phase-aligned by same-label transport_map action",
            "phase_lock": "same_label_transport_map_phase_alignment",
            "per_label_norm_checks": per_label,
        },
        "self_overlap": {
            "symbol": "F0^dagger F1",
            "matrix": _encode_complex_matrix(self_overlap),
            "abs_matrix": np.abs(self_overlap).tolist(),
        },
        "ckm_invariants": ckm,
        "comparison_shell": TARGET,
        "debug_log_shell_loss": {
            "transport_frame_self_overlap": transport_frame_loss,
            "current_same_sheet": current_same_sheet_loss,
            "improvement_factor_vs_current_same_sheet": (
                current_same_sheet_loss / transport_frame_loss if transport_frame_loss > 0.0 else None
            ),
        },
        "current_same_sheet_reference": current_same_sheet,
        "residual_gauge_quotient": {
            "transport_group": line_lift.get("transport_group"),
            "gauge_quotient": line_lift.get("gauge_quotient"),
            "vertex_rephasing_gauge_class": line_lift.get("vertex_rephasing_gauge_class"),
            "presentation_independence_status": line_lift.get("presentation_independence_status"),
            "action_on_self_overlap": "M -> D_u^dagger M D_d with D_u, D_d in U(1)^3",
            "entrywise_moduli_invariant": True,
            "ckm_moduli_invariant": True,
        },
        "smallest_exact_blocker": {
            "id": "sector_attached_transport_frame_to_sigma_ud_lift_equivalence",
            "statement": (
                "No sector-attached lift or equivalence currently identifies the common-refinement self-overlap "
                "matrix with an emitted same-label left-handed Sigma_ud element carrying "
                "{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM}."
            ),
        },
        "missing_sector_attachment": {
            "cannot_emit": [
                "sigma_id",
                "canonical_token",
                "U_u_left",
                "U_d_left",
                "V_CKM",
            ],
            "reason": (
                "The construction uses only common-refinement labels f1,f2,f3 and same-label line-lift transport maps. "
                "It does not emit a finite relative-sheet token or sector-attached left-handed up/down evaluations."
            ),
        },
        "why_not_promotable": (
            "The self-overlap F0^dagger F1 is a genuine already-local compare-only witness, but it lives on common-refinement "
            "transport frames rather than on emitted left-handed quark sheet evaluations. The residual objectwise-U(1) quotient "
            "acts only by diagonal rephasings M -> D_u^dagger M D_d and therefore preserves all entrywise moduli |M_ij|, so the "
            "remaining obstruction is not moving CKM moduli but the lack of a sector-attached lift/equivalence to an emitted "
            "same-label left-handed Sigma_ud element. Without that lift the construction cannot populate "
            "Sigma_ud_orbit.elements = [{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM, ckm_invariants}] or emit sigma_ud."
        ),
        "notes": [
            "This is the strongest already-local CKM-like signal downstream of the common-refinement line-lift on the current corpus.",
            "Its self-overlap is derived, not hardcoded: rerunning this emitter recomputes the compare-only matrix from the line-lift transport maps.",
            "The witness sharpens the sigma_ud frontier by separating a real local diagnostic from the still-missing sector-attached relative-sheet provider.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark transport-frame diagnostic orbit artifact.")
    parser.add_argument("--line-lift", default=str(DEFAULT_LINE_LIFT))
    parser.add_argument("--transport", default=str(DEFAULT_TRANSPORT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.line_lift)), _load_json(Path(args.transport)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
