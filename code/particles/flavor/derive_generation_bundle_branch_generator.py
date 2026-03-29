#!/usr/bin/env python3
"""Derive the intrinsic generation-bundle branch-generator artifact."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_complex_matrix(name: str, payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        matrix = np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
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


def _spectral_projectors(matrix: np.ndarray) -> tuple[list[np.ndarray], list[float], list[float]]:
    evals, evecs = np.linalg.eigh(matrix)
    order = np.argsort(evals)[::-1]
    evals = np.real_if_close(evals[order])
    evecs = evecs[:, order]
    projectors: list[np.ndarray] = []
    for idx in range(3):
        vec = evecs[:, idx : idx + 1]
        projectors.append(vec @ vec.conj().T)
    gaps = [float(abs(evals[idx] - evals[idx + 1])) for idx in range(2)]
    return projectors, [float(item) for item in evals.tolist()], gaps


def _same_label_overlap_amplitudes(projectors_by_refinement: list[list[np.ndarray]]) -> list[float]:
    if len(projectors_by_refinement) < 2:
        return []
    left = projectors_by_refinement[-2]
    right = projectors_by_refinement[-1]
    overlaps: list[float] = []
    for left_proj, right_proj in zip(left, right, strict=True):
        overlaps.append(float(abs(np.trace(right_proj @ left_proj))))
    return overlaps


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    refinements = list(payload.get("refinements", []))
    if not refinements:
        raise ValueError("family transport kernel must include refinements")

    decoded_refinements: list[dict[str, Any]] = []
    projectors_by_refinement: list[list[np.ndarray]] = []
    for refinement in refinements:
        if "hermitian_descendant" in refinement:
            compressed = _decode_complex_matrix("hermitian_descendant", refinement["hermitian_descendant"])
            generator_kind = "compressed_hermitian_descendant_proxy"
        else:
            transport = _decode_complex_matrix("transport_operator", refinement["transport_operator"])
            compressed = 0.5 * (transport + transport.conj().T)
            generator_kind = "compressed_transport_hermitian_part_proxy"
        trace_scalar = np.trace(compressed) / 3.0
        centered = compressed - trace_scalar * np.eye(3, dtype=complex)
        projectors, eigenvalues, spectral_gaps = _spectral_projectors(centered)
        decoded_refinements.append(
            {
                "level": int(refinement["level"]),
                "compressed": compressed,
                "centered": centered,
                "projectors": projectors,
                "eigenvalues": eigenvalues,
                "spectral_gaps": spectral_gaps,
                "generator_kind": generator_kind,
            }
        )
        projectors_by_refinement.append(projectors)

    latest = decoded_refinements[-1]
    level = int(latest["level"])
    compressed = np.asarray(latest["compressed"], dtype=complex)
    centered = np.asarray(latest["centered"], dtype=complex)
    projectors = list(latest["projectors"])
    eigenvalues = list(latest["eigenvalues"])
    spectral_gaps = list(latest["spectral_gaps"])
    generator_kind = str(latest["generator_kind"])
    min_gap = min(spectral_gaps) if spectral_gaps else 0.0
    refinement_stability = dict(payload.get("refinement_stability", {}))
    defect_sup = float(refinement_stability.get("conjugacy_defect_sup", 0.0))
    same_label_overlap_amplitudes = _same_label_overlap_amplitudes(projectors_by_refinement)
    trace_square = float(np.trace(centered @ centered).real)
    spectral_span = float(max(eigenvalues) - min(eigenvalues))
    conservative_lower_bound = float(min_gap - 2.0 * defect_sup)
    defect_gap_ratio = float(defect_sup / min_gap) if min_gap > 0.0 else None

    return {
        "artifact": "oph_intrinsic_generation_bundle_branch_generator",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only",
        "operator_theorem_candidate": "oph_generation_bundle_branch_generator_splitting",
        "carrier_kind": "persistent_charged_generation_multiplicity_bundle",
        "carrier_dimension": 3,
        "realized_generation_count": 3,
        "source_chain": [
            "ax:maxent",
            "internal_refinement_directed_system",
            "prop:regulatedrealization",
            "thm:markovcollar",
            "thm:master_ng_equals_3",
        ],
        "boundary_action_kind": "compact_unitary_groupoid_or_central_extension",
        "charged_block_projector_kind": "same_label_sector_block",
        "generation_bundle_projector_kind": "realized_ng_equals_3_bundle",
        "branch_generator_kind": generator_kind,
        "current_proxy_kind": generator_kind,
        "compressed_branch_generator": _encode_complex_matrix(compressed),
        "centered_compressed_branch_generator": _encode_complex_matrix(centered),
        "centering_rule": "K_gen - tr(K_gen)/3 * I",
        "selfadjoint": bool(np.allclose(centered, centered.conj().T)),
        "canonicality_certificate": {
            "carrier_kind": "persistent_charged_generation_multiplicity_bundle",
            "realized_generation_count": 3,
            "same_label_only": True,
            "exact_markov_required": False,
            "canonicality_mode": "up_to_refinement_conjugacy_and_objectwise_u1_projectivization",
            "line_lift_is_readout_of": "oph_intrinsic_generation_bundle_branch_generator",
        },
        "noncentrality_certificate": {
            "selfadjoint": bool(np.allclose(centered, centered.conj().T)),
            "trace_free": True,
            "is_scalar_multiple_identity": bool(np.allclose(centered, np.zeros((3, 3), dtype=complex))),
            "eigenvalues": eigenvalues,
            "trace_square": trace_square,
            "spectral_span": spectral_span,
        },
        "simple_spectrum_certificate": {
            "eigenvalues": eigenvalues,
            "spectral_gaps": spectral_gaps,
            "min_gap": float(min_gap),
            "simple_spectrum": bool(min_gap > 1.0e-12),
            "uniform_gap_candidate": bool(refinement_stability.get("uniform_three_cluster_gap", False)),
            "projector_order_preserved": bool(refinement_stability.get("projector_order_preserved", False)),
        },
        "persistent_gap_certificate": {
            "base_min_gap": float(min_gap),
            "conjugacy_defect_sup": defect_sup,
            "defect_gap_ratio": defect_gap_ratio,
            "projector_order_preserved": bool(refinement_stability.get("projector_order_preserved", False)),
            "riesz_bound_passes": bool(conservative_lower_bound > 0.0),
            "conservative_gap_lower_bound": conservative_lower_bound,
        },
        "schur_diagonal_pairing_certificate": {
            "same_label_only": True,
            "exact_markov_required": False,
        },
        "family_projectors_by_label": [
            {"label": label, "projector": _encode_complex_matrix(projector)}
            for label, projector in zip(("f1", "f2", "f3"), projectors, strict=True)
        ],
        "projective_readout_certificate": {
            "family_labels": ["f1", "f2", "f3"],
            "same_label_overlap_amplitudes": same_label_overlap_amplitudes,
        },
        "actual_generator_transfer_candidate": {
            "candidate_id": "centered_branch_generator_commutator_transfer",
            "central_split_rule": "K_actual = K_centered + c*I with the central term removed before comparison",
            "actual_generator_centering_rule": "K_actual_centered = K_actual - tr(K_actual)/3 * I",
            "bridge_strategy": "control the compression-descendant commutator after the central split",
            "actual_proxy_centered_residual_kind": "compression_descendant_commutator",
            "first_order_residual_after_central_split": "vanishes_if_commutator_zero",
            "descended_commutator_control_mode": "exact_zero_or_uniform_quadratic",
            "desired_outcome": "commutator vanishes exactly or is uniformly quadratic-small on realized refinement arrows",
            "quadratic_residual_bound_candidate": "||[C_actual_centered,C_proxy_centered]|| <= O(defect^2)",
            "quadratic_factorization_claim": "all surviving centered P->P corrections factor through P->Q->P",
            "smaller_exact_missing_clause": "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split",
            "transfer_if_closed_effect": "proxy_defect_vs_gap_estimate_lifts_to_actual_generator",
            "actual_generator_promotion_margin_if_bound_closed": "reuse_proxy_simple_spectrum_and_riesz_margin",
            "status": "candidate_only",
        },
        "readout_status": "line_lift_is_downstream_serialization",
        "remaining_missing_theorem": "oph_generation_bundle_branch_generator_splitting",
        "promotion_gate": {
            "proof_status": "candidate_only",
            "exact_missing_ingredient": "OPH proof that the actual compressed branch generator satisfies the same defect-vs-gap estimate as the current proxy",
            "exact_transfer_bridge": "central split plus compression-descendant commutator control",
            "smaller_exact_missing_clause": "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split",
        },
        "latest_refinement_level": level,
        "metadata": {
            "transport_artifact": payload.get("artifact"),
            "note": "This artifact promotes the flavor origin from a vague missing transport functor to a concrete centered compressed branch-generator candidate on the realized three-generation charged bundle. The standard Kato/Riesz persistence shell is now explicit on the current proxy; the remaining OPH-only burden is proving that the actual compressed branch generator satisfies the same defect-vs-gap estimate, with the current constructive bridge being a central split plus compression-descendant commutator control where the first surviving residual vanishes when the descended commutator vanishes and otherwise is treated as uniformly quadratic through P->Q->P factorization.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the generation-bundle branch-generator artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
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
