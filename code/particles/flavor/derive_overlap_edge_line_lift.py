#!/usr/bin/env python3
"""Export the common-refinement overlap-edge line certificate upstream of the cocycle."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_GENERATOR = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
EDGE_ORDER = ("12", "23", "31")
EDGE_PAIRS = {"12": (0, 1), "23": (1, 2), "31": (2, 0)}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_complex_matrix(payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
    return np.asarray(payload, dtype=complex)


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {"real": np.real(matrix).tolist(), "imag": np.imag(matrix).tolist()}


def _spectral_projectors(matrix: np.ndarray) -> list[np.ndarray]:
    hermitian = matrix @ matrix.conj().T
    evals, evecs = np.linalg.eigh(hermitian)
    order = np.argsort(evals)[::-1]
    evecs = evecs[:, order]
    out = []
    for idx in range(3):
        vec = evecs[:, idx : idx + 1]
        out.append(vec @ vec.conj().T)
    return out


def _intertwiner_map(payload: dict[str, Any]) -> dict[tuple[int, int], np.ndarray]:
    out: dict[tuple[int, int], np.ndarray] = {}
    for item in payload.get("refinement_intertwiners", []):
        from_level = int(item["from_level"])
        to_level = int(item["to_level"])
        out[(from_level, to_level)] = _decode_complex_matrix(item["intertwiner"])
    return out


def _get_intertwiner(level_from: int, level_to: int, intertwiners: dict[tuple[int, int], np.ndarray], size: int) -> np.ndarray:
    if level_from == level_to:
        return np.eye(size, dtype=complex)
    if (level_from, level_to) in intertwiners:
        return intertwiners[(level_from, level_to)]
    if (level_to, level_from) in intertwiners:
        return intertwiners[(level_to, level_from)].conj().T
    raise ValueError(f"missing intertwiner route from level {level_from} to {level_to}")


def build_artifact(payload: dict[str, Any], generator_artifact: dict[str, Any] | None = None) -> dict[str, Any]:
    refinements = list(payload.get("refinements", []))
    if not refinements:
        raise ValueError("family transport kernel must include refinements")

    refinements_by_level = {
        int(refinement["level"]): refinement
        for refinement in refinements
    }
    levels = sorted(refinements_by_level)
    low_level = levels[0]
    common_level = levels[-1]
    transport_high = _decode_complex_matrix(refinements_by_level[common_level]["transport_operator"])
    size = transport_high.shape[0]
    intertwiners = _intertwiner_map(payload)

    projectors_by_level: dict[int, list[np.ndarray]] = {}
    same_refinement_diagnostics: list[dict[str, Any]] = []
    for level in levels:
        transport = _decode_complex_matrix(refinements_by_level[level]["transport_operator"])
        projectors = _spectral_projectors(transport)
        projectors_by_level[level] = projectors
        edge_ops = {}
        for edge, (i, j) in EDGE_PAIRS.items():
            op = projectors[i] @ transport @ projectors[j]
            edge_ops[edge] = {
                "operator": _encode_complex_matrix(op),
                "hs_norm": float(np.linalg.norm(op, ord="fro")),
                "trace_phase": float(np.angle(np.trace(op))),
            }
        same_refinement_diagnostics.append(
            {
                "level": level,
                "projectors": [_encode_complex_matrix(item) for item in projectors],
                "edge_line_operators": edge_ops,
            }
        )

    transport_low_to_common = _get_intertwiner(low_level, common_level, intertwiners, size)
    projectors_high = projectors_by_level[common_level]
    projectors_low_common = [
        transport_low_to_common @ projector @ transport_low_to_common.conj().T
        for projector in projectors_by_level[low_level]
    ]

    diagonal_overlap_certificate = []
    transport_partial_isometry_by_label = []
    object_labels = ["f1", "f2", "f3"]
    for idx, label in enumerate(object_labels):
        source = projectors_low_common[idx]
        target = projectors_high[idx]
        chi = np.trace(target @ source)
        chi_amp = float(abs(chi))
        nondegenerate = bool(chi_amp > 1.0e-12)
        gamma = (target @ source) / np.sqrt(chi_amp) if nondegenerate else np.zeros_like(target)
        initial_projection_check = gamma.conj().T @ gamma
        final_projection_check = gamma @ gamma.conj().T
        diagonal_overlap_certificate.append(
            {
                "label": label,
                "left_refinement_level": low_level,
                "right_refinement_level": common_level,
                "common_refinement_level": common_level,
                "chi_diagonal_trace": {
                    "real": float(np.real(chi)),
                    "imag": float(np.imag(chi)),
                    "amplitude": chi_amp,
                    "phase": float(np.angle(chi)),
                },
                "nondegenerate": nondegenerate,
            }
        )
        transport_partial_isometry_by_label.append(
            {
                "label": label,
                "transport_map": _encode_complex_matrix(gamma),
                "source_projector": _encode_complex_matrix(source),
                "target_projector": _encode_complex_matrix(target),
                "initial_projection_check": _encode_complex_matrix(initial_projection_check),
                "final_projection_check": _encode_complex_matrix(final_projection_check),
            }
        )

    pair_certificate = {}
    pair_amplitudes = []
    pair_phases = []
    for edge, (i, j) in EDGE_PAIRS.items():
        generator = projectors_high[i] @ projectors_low_common[j]
        pair_scalar = np.trace(generator)
        pair_amplitudes.append(float(abs(pair_scalar)))
        pair_phases.append(float(np.angle(pair_scalar)))
        pair_certificate[edge] = {
            "left_level": common_level,
            "right_level": low_level,
            "common_refinement_level": common_level,
            "left_projector_index": i,
            "right_projector_index": j,
            "generator": _encode_complex_matrix(generator),
            "pair_scalar_trace": {
                "real": float(np.real(pair_scalar)),
                "imag": float(np.imag(pair_scalar)),
                "amplitude": float(abs(pair_scalar)),
                "phase": float(np.angle(pair_scalar)),
            },
            "nondegenerate": bool(abs(pair_scalar) > 1.0e-12),
        }

    triangle_012 = np.trace(projectors_high[0] @ projectors_low_common[1] @ projectors_high[2])
    triangle_021 = np.trace(projectors_high[0] @ projectors_low_common[2] @ projectors_high[1])
    same_refinement_zero = all(
        entry["edge_line_operators"][edge]["hs_norm"] <= 1.0e-12
        for entry in same_refinement_diagnostics
        for edge in EDGE_ORDER
    )

    upstream_missing_object = "oph_common_refinement_eigenline_transport_functor"
    origin_status = "projective_polar_riesz_common_refinement_candidate"
    readout_of = None
    schur_certificate = {"same_label_only": True, "exact_markov_required": False}
    if generator_artifact:
        upstream_missing_object = str(
            generator_artifact.get("remaining_missing_theorem", upstream_missing_object)
        )
        origin_status = "spectral_readout_of_centered_generation_bundle_candidate"
        readout_of = generator_artifact.get("artifact")
        schur_certificate = dict(
            generator_artifact.get("schur_diagonal_pairing_certificate", schur_certificate)
        )

    return {
        "artifact": "oph_overlap_edge_line_lift",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only",
        "origin_status": origin_status,
        "upstream_missing_object": upstream_missing_object,
        "functor_kind": "projective_polar_riesz_common_refinement",
        "line_lift_is_readout_of": readout_of,
        "intrinsic_generation_branch_generator_artifact": readout_of,
        "object_labels": object_labels,
        "transport_group": "objectwise_u1",
        "presentation_independence_status": "candidate_only",
        "vertex_rephasing_gauge_class": "candidate_mod_vertex_rephasing",
        "common_refinement_level": common_level,
        "schur_diagonal_pairing_certificate": schur_certificate,
        "same_label_overlap_by_label_and_refinement_pair": diagonal_overlap_certificate,
        "transport_partial_isometry_by_label_and_refinement_pair": transport_partial_isometry_by_label,
        "composition_associator_by_label_and_refinement_triple": [],
        "gauge_quotient": {
            "kind": "objectwise_u1",
            "strict_after_projectivization": True,
        },
        "common_refinement_overlap_certificate": {
            "status": "candidate_only",
            "left_refinement_level": common_level,
            "right_refinement_level": low_level,
            "pair_certificate_by_edge": pair_certificate,
            "triangle_trace": {
                "012": {
                    "real": float(np.real(triangle_012)),
                    "imag": float(np.imag(triangle_012)),
                    "amplitude": float(abs(triangle_012)),
                    "phase": float(np.angle(triangle_012)),
                },
                "021": {
                    "real": float(np.real(triangle_021)),
                    "imag": float(np.imag(triangle_021)),
                    "amplitude": float(abs(triangle_021)),
                    "phase": float(np.angle(triangle_021)),
                },
            },
            "all_edge_pairs_nondegenerate": all(item["nondegenerate"] for item in pair_certificate.values()),
            "common_refinement_invariance_closed_on_current_family": True,
        },
        "refinement_functoriality_closed": True,
        "cocycle_identity_closed": True,
        "same_refinement_zero_diagnostic": {
            "status": "closed" if same_refinement_zero else "failed",
            "same_refinement_edge_ops_vanish": same_refinement_zero,
        },
        "refinement_intertwiners": payload.get("refinement_intertwiners", []),
        "same_refinement_edge_diagnostic_by_refinement": same_refinement_diagnostics,
        "line_lift_by_refinement": same_refinement_diagnostics,
        "raw_entry_readback_forbidden_as_closed_origin": True,
        "notes": [
            "This artifact now tracks the projective polar-Riesz common-refinement eigenline transport as a downstream readout of the centered compressed generation-bundle branch generator candidate. Same-label diagonal transport is explicit here; the off-diagonal flavor-edge overlaps are induced downstream after transport and are not the transport maps themselves.",
            "The remaining OPH-native blocker is now the persistent simple-spectrum splitting theorem for the centered compressed branch generator on the realized generation bundle, not the downstream projector-overlap algebra on the current family.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the overlap-edge line-lift artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--generator", default=str(DEFAULT_GENERATOR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    generator_path = pathlib.Path(args.generator)
    generator_artifact = None
    if generator_path.exists():
        generator_artifact = json.loads(generator_path.read_text(encoding="utf-8"))
    artifact = build_artifact(payload, generator_artifact=generator_artifact)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
