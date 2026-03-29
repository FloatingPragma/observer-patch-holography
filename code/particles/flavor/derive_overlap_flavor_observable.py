#!/usr/bin/env python3
"""Reduce a candidate transport-kernel artifact into a flavor observable."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _default_template() -> dict[str, Any]:
    return {
        "artifact": "oph_flavor_observable",
        "generated_utc": _timestamp(),
        "status": "template",
        "observable_kind": "persistent_spectral_triple",
        "proof_status": "candidate_only",
        "labels": ["f1", "f2", "f3"],
        "intrinsic_label_order": ["f1", "f2", "f3"],
        "family_projectors": [],
        "family_eigenvalues": [0.0, 0.0, 0.0],
        "spectral_gaps": [0.0, 0.0],
        "pairwise_suppression": [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        "cycle_phases": {
            "012": 0.0,
            "021": 0.0,
        },
        "persistent_projector_certificate": {
            "status": "unavailable",
            "stable_projector_count": 0,
            "uniform_three_cluster_gap": False,
            "conjugation_class_stable": False,
        },
        "persistent_spectral_triple": {},
        "refinement_stability": {},
        "gauge_invariance_checks": {},
        "metadata": {
            "note": "Replace this template with a derived overlap/defect flavor observable.",
        },
    }


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


def _require_square_3x3(name: str, value: Any) -> list[list[float]]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{name} must be a 3x3 nested list")
    out: list[list[float]] = []
    for row in value:
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError(f"{name} must be a 3x3 nested list")
        out.append([float(item) for item in row])
    return out


def _projectors_from_transport(matrix: np.ndarray) -> tuple[list[dict[str, Any]], list[float], list[float]]:
    hermitian = matrix @ matrix.conj().T
    evals, evecs = np.linalg.eigh(hermitian)
    order = np.argsort(evals)[::-1]
    evals = np.real_if_close(evals[order])
    evecs = evecs[:, order]
    projectors: list[dict[str, Any]] = []
    for idx in range(3):
        vec = evecs[:, idx : idx + 1]
        projector = vec @ vec.conj().T
        projectors.append(_encode_complex_matrix(projector))
    gaps = [float(abs(evals[idx] - evals[idx + 1])) for idx in range(2)]
    return projectors, [float(x) for x in evals.tolist()], gaps


def _pairwise_suppression(matrix: np.ndarray, projectors: list[np.ndarray]) -> list[list[float]]:
    out: list[list[float]] = []
    eps = 1e-12
    for i in range(3):
        row: list[float] = []
        for j in range(3):
            if i == j:
                row.append(0.0)
                continue
            amplitude = abs(np.trace(projectors[i] @ matrix @ projectors[j]))
            row.append(float(-np.log(max(amplitude, eps))))
        out.append(row)
    return out


def _cycle_phases(matrix: np.ndarray, projectors: list[np.ndarray]) -> dict[str, float]:
    i, j, k = (0, 1, 2)
    value = np.trace(projectors[i] @ matrix @ projectors[j] @ matrix @ projectors[k] @ matrix)
    phase = float(np.angle(value))
    return {
        "012": phase,
        "021": -phase,
    }


def _from_transport_payload(payload: dict[str, Any], cocycle: dict[str, Any] | None = None) -> dict[str, Any]:
    refinements = payload.get("refinements", [])
    if not isinstance(refinements, list) or not refinements:
        raise ValueError("transport-kernel payload must include at least one refinement")
    latest = refinements[-1]
    matrix = _decode_complex_matrix("transport_operator", latest.get("transport_operator"))
    projector_payloads, eigenvalues, gaps = _projectors_from_transport(matrix)
    projectors = [_decode_complex_matrix("family_projector", item) for item in projector_payloads]
    refinement_stability = dict(payload.get("refinement_stability", {}))
    gauge_checks = dict(payload.get("gauge_invariance_checks", {}))
    persistent_projector_certificate = {
        "status": str(payload.get("proof_status", "candidate_only")),
        "stable_projector_count": int(refinement_stability.get("stable_projector_count", 0)),
        "uniform_three_cluster_gap": bool(refinement_stability.get("uniform_three_cluster_gap", False)),
        "projector_order_preserved": bool(refinement_stability.get("projector_order_preserved", False)),
        "conjugation_class_stable": bool(gauge_checks.get("conjugation_class_stable", False)),
        "conjugacy_defect_sup": refinement_stability.get("conjugacy_defect_sup"),
    }
    if cocycle:
        persistent_projector_certificate["theorem_gap_gamma"] = cocycle.get("theorem_gap_gamma")
        persistent_projector_certificate["defect_gap_ratio"] = cocycle.get("defect_gap_ratio")
        persistent_projector_certificate["riesz_bound_passes"] = cocycle.get("riesz_bound_passes")
        persistent_projector_certificate["hermitian_descendant_riesz_margin"] = cocycle.get("hermitian_descendant_riesz_margin")
        persistent_projector_certificate["persistence_proof_status"] = cocycle.get("persistence_proof_status")
        pairwise_suppression = cocycle.get("derived_pairwise_suppression")
        cycle_phases = cocycle.get("derived_cycle_holonomy")
    else:
        pairwise_suppression = _pairwise_suppression(matrix, projectors)
        cycle_phases = _cycle_phases(matrix, projectors)
    return {
        "artifact": "oph_flavor_observable",
        "generated_utc": _timestamp(),
        "status": str(payload.get("status", "derived_from_transport")),
        "observable_kind": "persistent_spectral_triple",
        "proof_status": str(payload.get("proof_status", "conjugacy_riesz_candidate")),
        "labels": ["f1", "f2", "f3"],
        "intrinsic_label_order": ["f1", "f2", "f3"],
        "family_projectors": projector_payloads,
        "family_eigenvalues": eigenvalues,
        "spectral_gaps": gaps,
        "pairwise_suppression": pairwise_suppression,
        "cycle_phases": cycle_phases,
        "persistent_projector_certificate": persistent_projector_certificate,
        "persistent_spectral_triple": {
            "projectors": projector_payloads,
            "pairwise_suppression": pairwise_suppression,
            "cycle_holonomy": cycle_phases,
            "edge_transport_cocycle": None if cocycle is None else cocycle.get("artifact"),
        },
        "refinement_stability": refinement_stability,
        "gauge_invariance_checks": gauge_checks,
        "metadata": {
            "transport_artifact": payload.get("artifact", "unknown"),
            "cocycle_artifact": None if cocycle is None else cocycle.get("artifact"),
            **dict(payload.get("metadata", {})),
        },
    }


def normalize_payload(payload: dict[str, Any], cocycle: dict[str, Any] | None = None) -> dict[str, Any]:
    if payload.get("artifact") == "oph_family_transport_kernel":
        return _from_transport_payload(payload, cocycle=cocycle)
    labels = payload.get("labels", ["f1", "f2", "f3"])
    if not isinstance(labels, list) or len(labels) != 3:
        raise ValueError("labels must be a length-3 list")
    cycle_phases = payload.get("cycle_phases", {})
    if not isinstance(cycle_phases, dict):
        raise ValueError("cycle_phases must be a mapping")
    return {
        "artifact": "oph_flavor_observable",
        "generated_utc": _timestamp(),
        "status": str(payload.get("status", "normalized_input")),
        "observable_kind": str(payload.get("observable_kind", "persistent_spectral_triple")),
        "proof_status": str(payload.get("proof_status", "candidate_only")),
        "labels": [str(item) for item in labels],
        "intrinsic_label_order": [str(item) for item in payload.get("intrinsic_label_order", labels)],
        "family_projectors": list(payload.get("family_projectors", [])),
        "family_eigenvalues": [float(item) for item in payload.get("family_eigenvalues", [0.0, 0.0, 0.0])],
        "spectral_gaps": [float(item) for item in payload.get("spectral_gaps", [0.0, 0.0])],
        "pairwise_suppression": _require_square_3x3(
            "pairwise_suppression",
            payload.get("pairwise_suppression", _default_template()["pairwise_suppression"]),
        ),
        "cycle_phases": {str(key): float(value) for key, value in cycle_phases.items()},
        "persistent_projector_certificate": dict(payload.get("persistent_projector_certificate", {})),
        "persistent_spectral_triple": dict(payload.get("persistent_spectral_triple", {})),
        "refinement_stability": dict(payload.get("refinement_stability", {})),
        "gauge_invariance_checks": dict(payload.get("gauge_invariance_checks", {})),
        "metadata": dict(payload.get("metadata", {})),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Reduce a candidate transport kernel into a flavor observable.")
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
        help="Input JSON path. Defaults to the family transport-kernel artifact; omit with --template to emit a template.",
    )
    parser.add_argument(
        "--cocycle",
        default=str(DEFAULT_COCYCLE),
        help="Optional overlap-edge transport cocycle JSON path.",
    )
    parser.add_argument("--template", action="store_true", help="Emit a template observable instead of reading input.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    if args.template:
        normalized = _default_template()
    elif args.input:
        payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
        cocycle_path = pathlib.Path(args.cocycle)
        cocycle = None
        if cocycle_path.exists():
            cocycle = json.loads(cocycle_path.read_text(encoding="utf-8"))
        normalized = normalize_payload(payload, cocycle=cocycle)
    else:
        normalized = _default_template()

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
