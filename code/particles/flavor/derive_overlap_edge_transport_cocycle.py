#!/usr/bin/env python3
"""Derive an overlap-edge transport cocycle candidate from projector overlaps on a common refinement."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_LINE_LIFT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
EDGE_ORDER = ("12", "23", "31")
EDGE_MATRIX_INDEX = {"12": (0, 1), "23": (1, 2), "31": (2, 0)}


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


def _wrap_phase(value: float) -> float:
    return float(math.atan2(math.sin(value), math.cos(value)))


def _pairwise_suppression_from_edges(edge_payload: dict[str, dict[str, float]]) -> list[list[float]]:
    eps = 1.0e-12
    lookup = {
        (0, 1): edge_payload["12"]["amplitude"],
        (1, 0): edge_payload["12"]["amplitude"],
        (1, 2): edge_payload["23"]["amplitude"],
        (2, 1): edge_payload["23"]["amplitude"],
        (2, 0): edge_payload["31"]["amplitude"],
        (0, 2): edge_payload["31"]["amplitude"],
    }
    out: list[list[float]] = []
    for i in range(3):
        row: list[float] = []
        for j in range(3):
            if i == j:
                row.append(0.0)
                continue
            row.append(float(-0.5 * np.log(max(float(lookup[(i, j)]), eps))))
        out.append(row)
    return out


def build_artifact(payload: dict[str, Any], line_lift: dict[str, Any] | None = None) -> dict[str, Any]:
    refinements = list(payload.get("refinements", []))
    if not refinements:
        raise ValueError("family transport kernel must include refinements")
    if line_lift is None:
        raise ValueError("line-lift artifact is required for projector-overlap cocycle derivation")

    descendant_transport_by_refinement: list[dict[str, Any]] = []
    transport_norms: list[dict[str, float]] = []
    hermitian_descendants: list[np.ndarray] = []
    for refinement in refinements:
        level = int(refinement.get("level"))
        matrix = _decode_complex_matrix("transport_operator", refinement.get("transport_operator"))
        hermitian = matrix @ matrix.conj().T
        descendant_transport_by_refinement.append(
            {
                "level": level,
                "transport_operator": {
                    "real": np.real(matrix).tolist(),
                    "imag": np.imag(matrix).tolist(),
                },
                "hermitian_descendant": {
                    "real": np.real(hermitian).tolist(),
                    "imag": np.imag(hermitian).tolist(),
                },
            }
        )
        hermitian_descendants.append(hermitian)
        transport_norms.append({"level": level, "operator_norm": float(np.linalg.norm(matrix, ord=2))})

    line_certificate = dict(line_lift.get("common_refinement_overlap_certificate", {}))
    pair_by_edge = dict(line_certificate.get("pair_certificate_by_edge", {}))
    latest_edge_payload: dict[str, dict[str, float]] = {}
    for edge in EDGE_ORDER:
        pair_scalar = dict(pair_by_edge[edge]["pair_scalar_trace"])
        latest_edge_payload[edge] = {
            "real": float(pair_scalar["real"]),
            "imag": float(pair_scalar["imag"]),
            "amplitude": float(pair_scalar["amplitude"]),
            "phase": float(pair_scalar["phase"]),
        }

    triangle_trace = dict(line_certificate.get("triangle_trace", {}))
    cycle_holonomy = {
        "012": float(dict(triangle_trace.get("012", {})).get("phase", 0.0)),
        "021": -float(dict(triangle_trace.get("012", {})).get("phase", 0.0)),
    }
    pairwise_suppression = _pairwise_suppression_from_edges(latest_edge_payload)

    stability = dict(payload.get("refinement_stability", {}))
    theorem_gap_gamma = float(stability.get("min_spectral_gap", 0.0))
    conjugacy_defect_sup = float(stability.get("conjugacy_defect_sup", 0.0))
    defect_gap_ratio = None
    if theorem_gap_gamma > 0.0:
        defect_gap_ratio = conjugacy_defect_sup / theorem_gap_gamma

    h_margin = None
    if len(refinements) >= 2 and theorem_gap_gamma > 0.0:
        t0 = _decode_complex_matrix("transport_operator", refinements[0].get("transport_operator"))
        t1 = _decode_complex_matrix("transport_operator", refinements[-1].get("transport_operator"))
        h0 = hermitian_descendants[0]
        h1 = hermitian_descendants[-1]
        epsilon = float(np.linalg.norm(t1 - t0, ord=2))
        direct_h_gap = float(np.linalg.norm(h1 - h0, ord=2))
        upper_bound = float((np.linalg.norm(t0, ord=2) + np.linalg.norm(t1, ord=2)) * epsilon)
        gamma_half = theorem_gap_gamma / 2.0
        h_margin = {
            "t0_operator_norm": float(np.linalg.norm(t0, ord=2)),
            "t1_operator_norm": float(np.linalg.norm(t1, ord=2)),
            "epsilon_transport_gap": epsilon,
            "gamma": theorem_gap_gamma,
            "gamma_half": gamma_half,
            "hermitian_descendant_norm_upper_bound": upper_bound,
            "hermitian_descendant_norm_direct": direct_h_gap,
            "passes": direct_h_gap < gamma_half,
        }

    amplitudes = [latest_edge_payload[edge]["amplitude"] for edge in EDGE_ORDER]
    amplitude_floor = 1.0e-12
    non_floor_amplitude_closed = bool(amplitudes) and all(value > amplitude_floor for value in amplitudes)
    all_equal_off_diagonal = len({round(float(value), 12) for value in amplitudes}) <= 1

    cycle_identity_closed = bool(line_lift.get("cocycle_identity_closed", False))
    refinement_functoriality_closed = bool(line_lift.get("refinement_functoriality_closed", False))

    return {
        "artifact": "oph_overlap_edge_transport_cocycle",
        "generated_utc": _timestamp(),
        "proof_status": "common_refinement_projector_overlap_candidate",
        "persistence_proof_status": "conditional_standard_math_closed" if h_margin and h_margin["passes"] else "candidate_only",
        "cocycle_origin_status": "induced_from_projective_eigenline_transport",
        "upstream_missing_object": line_lift.get("upstream_missing_object", "oph_common_refinement_eigenline_transport_functor"),
        "overlap_edge_line_lift_artifact": line_lift.get("artifact"),
        "projective_transport_functor_kind": line_lift.get("functor_kind"),
        "transport_artifact": payload.get("artifact"),
        "labels": ["f1", "f2", "f3"],
        "edge_basis": ["12", "23", "31"],
        "edge_transport_by_refinement": [
            {
                "level": int(line_lift.get("common_refinement_level", refinements[-1]["level"])),
                "edges": latest_edge_payload,
            }
        ],
        "cycle_holonomy_by_refinement": [
            {
                "level": int(line_lift.get("common_refinement_level", refinements[-1]["level"])),
                "omega_012": cycle_holonomy["012"],
                "omega_021": cycle_holonomy["021"],
            }
        ],
        "projector_overlap_line_certificate": line_certificate,
        "descendant_transport_operator_by_refinement": descendant_transport_by_refinement,
        "transport_operator_norms_by_refinement": transport_norms,
        "latest_level": int(line_lift.get("common_refinement_level", refinements[-1]["level"])),
        "derived_pairwise_suppression": pairwise_suppression,
        "derived_cycle_holonomy": cycle_holonomy,
        "conjugacy_defect_sup": conjugacy_defect_sup,
        "theorem_gap_gamma": theorem_gap_gamma,
        "defect_gap_ratio": defect_gap_ratio,
        "riesz_threshold": 0.5,
        "riesz_bound_passes": defect_gap_ratio is not None and defect_gap_ratio < 0.5,
        "hermitian_descendant_riesz_margin": h_margin,
        "cocycle_identity_closed": cycle_identity_closed,
        "refinement_functoriality_closed": refinement_functoriality_closed,
        "vertex_rephasing_gauge_class": line_lift.get("vertex_rephasing_gauge_class", "candidate_mod_vertex_rephasing"),
        "common_refinement_invariance_closed": bool(line_certificate.get("common_refinement_invariance_closed_on_current_family", False)),
        "non_floor_amplitude_certificate": {
            "status": "closed" if non_floor_amplitude_closed else "failed",
            "amplitude_floor": amplitude_floor,
            "all_equal_off_diagonal": all_equal_off_diagonal,
        },
        "metadata": {
            "note": "Overlap-edge transport cocycle candidate is now treated as an induced off-diagonal edge object downstream of the projective same-label eigenline transport readout, which itself is now a serialization of the centered compressed generation-bundle branch-generator candidate rather than a primitive transport law. The remaining OPH-native blocker is the persistent simple-spectrum splitting theorem behind that generator.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Derive an overlap-edge transport cocycle candidate.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input family transport-kernel JSON path.")
    parser.add_argument("--line-lift", default=str(DEFAULT_LINE_LIFT), help="Overlap-edge line-lift JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output cocycle JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    line_lift = json.loads(pathlib.Path(args.line_lift).read_text(encoding="utf-8"))
    artifact = build_artifact(payload, line_lift=line_lift)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
