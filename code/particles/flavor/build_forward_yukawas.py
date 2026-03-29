#!/usr/bin/env python3
"""Turn local quark tensor candidates into masses and mixing diagnostics.

Chain role: take the active `/particles` quark tensors and expose their
singular-value masses, left bases, CKM matrix, and invariants on one surface.

Mathematics: matrix assembly, SVD, left-diagonalization, commutator
invariants, and Jarlskog extraction.

OPH-derived inputs: explicit local Yukawa tensors from the flavor descent.

Output: the live quark forward candidate plus certification metadata used by
the public table and exactness audits.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_complex_matrix(name: str, payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        real = np.asarray(payload["real"], dtype=float)
        imag = np.asarray(payload["imag"], dtype=float)
        matrix = real + 1j * imag
    else:
        matrix = np.asarray(payload, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix")
    return matrix


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _build_from_tensors(payload: dict[str, Any], prefix: str) -> np.ndarray:
    S = np.asarray(payload[f"S_{prefix}"], dtype=float)
    phi = np.asarray(payload[f"Phi_{prefix}"], dtype=float)
    A = np.asarray(payload.get(f"A_{prefix}", np.ones((3, 3))), dtype=float)
    if S.shape != (3, 3) or phi.shape != (3, 3) or A.shape != (3, 3):
        raise ValueError(f"S_{prefix}, Phi_{prefix}, and A_{prefix} must all be 3x3")
    return A * np.exp(-S) * np.exp(1j * phi)


def _build_from_factorized(payload: dict[str, Any], prefix: str) -> np.ndarray:
    suppression = np.asarray(payload[f"S_{prefix}"], dtype=float)
    phase = np.asarray(payload[f"Phi_{prefix}"], dtype=float)
    log_left = np.asarray(payload[f"D_{prefix}_left_log"], dtype=float)
    log_right = np.asarray(payload[f"D_{prefix}_right_log"], dtype=float)
    even_log = np.asarray(payload.get(f"E_{prefix}_log", [0.0, 0.0, 0.0]), dtype=float)
    g_scalar = float(payload[f"g_{prefix}"])
    if suppression.shape != (3, 3) or phase.shape != (3, 3):
        raise ValueError(f"S_{prefix} and Phi_{prefix} must both be 3x3")
    if log_left.shape != (3,) or log_right.shape != (3,) or even_log.shape != (3,):
        raise ValueError(f"D_{prefix}_left_log, D_{prefix}_right_log, and E_{prefix}_log must all be length-3")
    return (
        g_scalar
        * np.diag(np.exp(log_left + even_log))
        @ (np.exp(-suppression) * np.exp(1j * phase))
        @ np.diag(np.exp(log_right + even_log))
    )


def _sorted_left_diagonalizer(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = matrix @ matrix.conj().T
    evals, evecs = np.linalg.eigh(hermitian)
    order = np.argsort(evals)
    evals = np.real_if_close(evals[order])
    evecs = evecs[:, order]
    singular_values = np.sqrt(np.clip(evals, 0.0, None))
    return singular_values, evecs


def _jarlskog(ckm: np.ndarray) -> float:
    return float(np.imag(ckm[0, 0] * ckm[1, 1] * np.conj(ckm[0, 1]) * np.conj(ckm[1, 0])))


def _commutator_invariant(y_u: np.ndarray, y_d: np.ndarray) -> dict[str, float]:
    h_u = y_u @ y_u.conj().T
    h_d = y_d @ y_d.conj().T
    commutator = h_u @ h_d - h_d @ h_u
    det_value = np.linalg.det(commutator)
    return {
        "fro_norm": float(np.linalg.norm(commutator, ord="fro")),
        "det_abs": float(abs(det_value)),
        "det_imag": float(np.imag(det_value)),
    }


def _quark_certification(payload: dict[str, Any]) -> tuple[bool, bool, bool, list[str]]:
    blockers: list[str] = []
    factorized_mode = all(
        key in payload
        for key in (
            "D_u_left_log",
            "D_u_right_log",
            "D_d_left_log",
            "D_d_right_log",
            "g_u",
            "g_d",
        )
    )
    dense_entrywise_amplitude_used = not factorized_mode and ("A_u" in payload or "A_d" in payload)
    template_amplitude_fallback_used = not factorized_mode and ("A_u" not in payload or "A_d" not in payload)
    up_down_sector_distinct = False
    if all(key in payload for key in ("S_u", "S_d", "Phi_u", "Phi_d")):
        s_u = np.asarray(payload["S_u"], dtype=float)
        s_d = np.asarray(payload["S_d"], dtype=float)
        phi_u = np.asarray(payload["Phi_u"], dtype=float)
        phi_d = np.asarray(payload["Phi_d"], dtype=float)
        up_down_sector_distinct = not (np.allclose(s_u, s_d) and np.allclose(phi_u, phi_d))
    if payload.get("sector_distinctness_witness", {}).get("value") is not None:
        up_down_sector_distinct = up_down_sector_distinct or float(payload["sector_distinctness_witness"]["value"]) > 1.0e-12
    if template_amplitude_fallback_used:
        blockers.append("template_amplitude_fallback")
    if dense_entrywise_amplitude_used:
        blockers.append("dense_entrywise_amplitude")
    if not up_down_sector_distinct:
        blockers.append("u_d_sector_clone")
    if not bool(payload.get("uses_full_projector_algebra", False)):
        blockers.append("missing_projector_resolved_descent")
    if bool(payload.get("degenerate_placeholder_fallback_used", False)):
        blockers.append("odd_splitter_degenerate_fallback")
    if str(payload.get("quark_descent_proof_status", "")) != "closed":
        blockers.append("quark_descent_candidate_only")
    noncentrality_witness = dict(payload.get("noncentrality_witness", {}))
    if float(noncentrality_witness.get("fro_norm", 0.0)) <= 1.0e-18:
        blockers.append("noncentrality_missing")
    if dict(payload.get("budget_neutrality_certificate", {})).get("status") != "closed":
        blockers.append("budget_neutrality_open")
    even_excitation_present = payload.get("E_u_log") is not None and payload.get("E_d_log") is not None
    if not even_excitation_present:
        if payload.get("exact_missing_object") == "oph_family_excitation_gap_map":
            blockers.append("quark_gap_pair_map_missing")
        else:
            blockers.append("quark_even_excitation_evaluator_missing")
    elif str(payload.get("even_excitation_proof_status", "")) not in {"closed", "exact"}:
        blockers.append("quark_even_excitation_candidate_only")
    if payload.get("exact_missing_object") not in {None, ""}:
        blockers.append("quark_source_mode_unresolved")
    if (
        payload.get("b_odd_source_scalar_evaluator_artifact") == "oph_quark_diagonal_B_odd_source_scalar_evaluator"
        and (
            payload.get("J_B_source_u") is None
            or payload.get("J_B_source_d") is None
        )
    ):
        blockers.append("quark_b_odd_source_values_open")
    if payload.get("u_raw_channel_norm_candidate") in {None, 1.0} or payload.get("d_raw_channel_norm_candidate") in {None, 1.0}:
        blockers.append("quark_shared_norm_unavailable")
    return template_amplitude_fallback_used, up_down_sector_distinct, not blockers, blockers


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    tensor_ready = all(key in payload for key in ("S_u", "S_d", "Phi_u", "Phi_d"))
    if "Y_u" in payload and "Y_d" in payload and not tensor_ready:
        y_u = _decode_complex_matrix("Y_u", payload["Y_u"])
        y_d = _decode_complex_matrix("Y_d", payload["Y_d"])
        source_mode = "explicit_matrices"
        template_amplitude_fallback_used = False
        up_down_sector_distinct = not np.allclose(y_u, y_d)
        forward_certified = True
        promotion_blockers: list[str] = []
    else:
        factorized_mode = all(
            key in payload
            for key in (
                "D_u_left_log",
                "D_u_right_log",
                "D_d_left_log",
                "D_d_right_log",
                "g_u",
                "g_d",
            )
        )
        if factorized_mode:
            y_u = _build_from_factorized(payload, "u")
            y_d = _build_from_factorized(payload, "d")
            source_mode = "factorized_descent"
        else:
            y_u = _build_from_tensors(payload, "u")
            y_d = _build_from_tensors(payload, "d")
            source_mode = "tensor_derived"
        (
            template_amplitude_fallback_used,
            up_down_sector_distinct,
            forward_certified,
            promotion_blockers,
        ) = _quark_certification(payload)

    masses_u, u_left = _sorted_left_diagonalizer(y_u)
    masses_d, d_left = _sorted_left_diagonalizer(y_d)
    ckm = u_left.conj().T @ d_left
    certification_status = "forward_matrix_certified" if forward_certified else "placeholder_unpromotable"
    public_surface_candidate_allowed = (
        source_mode == "factorized_descent"
        and not template_amplitude_fallback_used
        and up_down_sector_distinct
        and not bool(payload.get("dense_entrywise_amplitude_used", False))
    )

    return {
        "artifact": "oph_forward_yukawas",
        "generated_utc": _timestamp(),
        "source_mode": source_mode,
        "template_amplitude_fallback_used": template_amplitude_fallback_used,
        "dense_entrywise_amplitude_used": "A_u" in payload or "A_d" in payload,
        "up_down_sector_distinct": up_down_sector_distinct,
        "uses_full_projector_algebra": payload.get("uses_full_projector_algebra"),
        "quark_descent_proof_status": payload.get("quark_descent_proof_status"),
        "shared_scalarization_law_status": payload.get("shared_scalarization_law_status"),
        "delta_logg_q_status": payload.get("delta_logg_q_status"),
        "sector_mean_split_artifact": payload.get("sector_mean_split_artifact"),
        "sector_mean_split_status": payload.get("sector_mean_split_status"),
        "degenerate_placeholder_fallback_used": payload.get("degenerate_placeholder_fallback_used"),
        "predictive_promotion_allowed": forward_certified,
        "public_surface_candidate_allowed": public_surface_candidate_allowed,
        "public_surface_candidate_scope": ["m_u", "m_d", "m_s", "m_c", "m_b"],
        "public_surface_candidate_status": (
            "best_available_reference_free_current_family_candidate"
            if public_surface_candidate_allowed
            else "not_allowed"
        ),
        "forward_certified": forward_certified,
        "certification_status": certification_status,
        "promotion_blockers": promotion_blockers,
        "Y_u": _encode_complex_matrix(y_u),
        "Y_d": _encode_complex_matrix(y_d),
        "singular_values_u": [float(x) for x in masses_u.tolist()],
        "singular_values_d": [float(x) for x in masses_d.tolist()],
        "U_u_left": _encode_complex_matrix(u_left),
        "U_d_left": _encode_complex_matrix(d_left),
        "V_CKM": _encode_complex_matrix(ckm),
        "jarlskog": _jarlskog(ckm),
        "commutator_invariant": _commutator_invariant(y_u, y_d),
        "noncentrality_witness": payload.get("noncentrality_witness"),
        "metadata": {
            **dict(payload.get("metadata", {})),
            "note": "This artifact stays a sandbox placeholder until the quark descent law closes. The exact reduced even family is already fixed on the current spectrum class; if a current-family spread pair is populated, the builder now emits a real hierarchy candidate, but it remains non-promotable until the spread emitter itself is theorem-fed. Factorized-only mode is allowed, and dense entrywise amplitudes are not promotable in certified quark mode.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build forward Yukawa artifacts from matrices or tensors.")
    parser.add_argument("--input", required=True, help="Input JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
