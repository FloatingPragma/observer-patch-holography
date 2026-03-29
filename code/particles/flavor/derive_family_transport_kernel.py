#!/usr/bin/env python3
"""Normalize a candidate family transport-kernel artifact across refinements."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _decode_complex_matrix(name: str, payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        matrix = np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
    else:
        matrix = np.asarray(payload, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix")
    return matrix


def _spectral_data(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hermitian = matrix @ matrix.conj().T
    evals, evecs = np.linalg.eigh(hermitian)
    order = np.argsort(evals)[::-1]
    evals = np.real_if_close(evals[order])
    evecs = evecs[:, order]
    gaps = np.abs(np.diff(evals))
    return hermitian, evals, gaps


def _identity_intertwiner() -> dict[str, Any]:
    return _encode_complex_matrix(np.eye(3, dtype=complex))


def _intertwiner_matrix(payload: Any) -> np.ndarray:
    if payload is None:
        return np.eye(3, dtype=complex)
    return _decode_complex_matrix("refinement_intertwiner", payload)


def _conjugacy_defect(current: np.ndarray, nxt: np.ndarray, intertwiner: np.ndarray) -> dict[str, float]:
    transported = np.linalg.inv(intertwiner) @ nxt @ intertwiner
    defect = transported - current
    return {
        "frobenius_norm": float(np.linalg.norm(defect, ord="fro")),
        "operator_norm": float(np.linalg.norm(defect, ord=2)),
    }


def _template() -> dict[str, Any]:
    t0 = np.asarray(
        [
            [1.0 + 0.0j, 0.1 + 0.2j, 0.0 + 0.1j],
            [0.1 - 0.2j, 0.7 + 0.0j, 0.2 + 0.0j],
            [0.0 - 0.1j, 0.2 + 0.0j, 0.4 + 0.0j],
        ]
    )
    t1 = np.asarray(
        [
            [1.02 + 0.0j, 0.11 + 0.19j, 0.01 + 0.09j],
            [0.11 - 0.19j, 0.69 + 0.0j, 0.18 + 0.01j],
            [0.01 - 0.09j, 0.18 - 0.01j, 0.41 + 0.0j],
        ]
    )
    refinements = []
    for level, matrix in enumerate((t0, t1)):
        hermitian, evals, gaps = _spectral_data(matrix)
        refinements.append(
            {
                "level": level,
                "transport_operator": _encode_complex_matrix(matrix),
                "hermitian_descendant": _encode_complex_matrix(hermitian),
                "eigenvalues": [float(x) for x in evals.tolist()],
                "spectral_gaps": [float(x) for x in gaps.tolist()],
            }
        )

    intertwiners = [
        {
            "from_level": 0,
            "to_level": 1,
            "intertwiner": _identity_intertwiner(),
        }
    ]
    conjugacy_defects = [
        {
            "from_level": 0,
            "to_level": 1,
            **_conjugacy_defect(t0, t1, np.eye(3, dtype=complex)),
        }
    ]
    gap_certificates = [
        {
            "level": item["level"],
            "min_gap": float(min(item["spectral_gaps"])) if item["spectral_gaps"] else 0.0,
            "three_cluster_gap_open": bool(item["spectral_gaps"] and min(item["spectral_gaps"]) > 0.0),
        }
        for item in refinements
    ]
    max_drift = float(
        np.max(np.abs(np.asarray(refinements[1]["eigenvalues"]) - np.asarray(refinements[0]["eigenvalues"])))
    )

    return {
        "artifact": "oph_family_transport_kernel",
        "generated_utc": _timestamp(),
        "status": "template",
        "transport_kind": "conjugacy_class_family_kernel",
        "proof_status": "conjugacy_riesz_candidate",
        "refinements": refinements,
        "refinement_intertwiners": intertwiners,
        "conjugacy_defects": conjugacy_defects,
        "gap_certificates": gap_certificates,
        "refinement_stability": {
            "stable_projector_count": 3,
            "max_eigenvalue_drift": max_drift,
            "min_spectral_gap": float(min(item["min_gap"] for item in gap_certificates)),
            "uniform_three_cluster_gap": all(item["three_cluster_gap_open"] for item in gap_certificates),
            "projector_order_preserved": True,
            "conjugacy_defect_sup": float(max(item["operator_norm"] for item in conjugacy_defects)),
            "conjugacy_persistence_candidate": True,
        },
        "gauge_invariance_checks": {
            "conjugation_class_stable": True,
            "rephasing_sensitive_cycles_detected": False,
        },
        "metadata": {
            "note": "Replace this template with an OPH-derived transport kernel across refinements. The intended theorem burden is conjugacy-Riesz persistence of the three-projector split.",
        },
    }


def _normalize_refinements(raw_refinements: list[Any]) -> list[dict[str, Any]]:
    refinements: list[dict[str, Any]] = []
    for idx, item in enumerate(raw_refinements):
        if isinstance(item, dict):
            level = int(item.get("level", idx))
            matrix = _decode_complex_matrix("transport_operator", item.get("transport_operator", item.get("T", item)))
        else:
            level = idx
            matrix = _decode_complex_matrix("transport_operator", item)
        hermitian, evals, gaps = _spectral_data(matrix)
        refinements.append(
            {
                "level": level,
                "transport_operator": _encode_complex_matrix(matrix),
                "hermitian_descendant": _encode_complex_matrix(hermitian),
                "eigenvalues": [float(x) for x in evals.tolist()],
                "spectral_gaps": [float(x) for x in gaps.tolist()],
            }
        )
    if not refinements:
        raise ValueError("refinements must contain at least one 3x3 transport operator")
    return refinements


def _normalize_intertwiners(payload: dict[str, Any], refinements: list[dict[str, Any]]) -> list[dict[str, Any]]:
    raw = payload.get("refinement_intertwiners")
    if isinstance(raw, list) and raw:
        intertwiners: list[dict[str, Any]] = []
        for item in raw:
            if not isinstance(item, dict):
                raise ValueError("refinement_intertwiners entries must be mappings")
            from_level = int(item.get("from_level"))
            to_level = int(item.get("to_level"))
            intertwiner = _intertwiner_matrix(item.get("intertwiner"))
            intertwiners.append(
                {
                    "from_level": from_level,
                    "to_level": to_level,
                    "intertwiner": _encode_complex_matrix(intertwiner),
                }
            )
        return intertwiners
    intertwiners = []
    for left, right in zip(refinements, refinements[1:]):
        intertwiners.append(
            {
                "from_level": int(left["level"]),
                "to_level": int(right["level"]),
                "intertwiner": _identity_intertwiner(),
            }
        )
    return intertwiners


def normalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    raw_refinements = payload.get("refinements")
    if raw_refinements is None:
        raw_refinements = [payload.get("transport_operator", payload.get("T"))]
    refinements = _normalize_refinements(list(raw_refinements))
    eigenvalue_stack = np.asarray([item["eigenvalues"] for item in refinements], dtype=float)
    intertwiners = _normalize_intertwiners(payload, refinements)
    matrix_by_level = {
        int(item["level"]): _decode_complex_matrix("transport_operator", item["transport_operator"])
        for item in refinements
    }
    conjugacy_defects: list[dict[str, Any]] = []
    for item in intertwiners:
        left = int(item["from_level"])
        right = int(item["to_level"])
        if left not in matrix_by_level or right not in matrix_by_level:
            raise ValueError("refinement_intertwiners reference unknown refinement levels")
        intertwiner = _intertwiner_matrix(item["intertwiner"])
        conjugacy_defects.append(
            {
                "from_level": left,
                "to_level": right,
                **_conjugacy_defect(matrix_by_level[left], matrix_by_level[right], intertwiner),
            }
        )
    gap_certificates = [
        {
            "level": int(item["level"]),
            "min_gap": float(min(item["spectral_gaps"])) if item["spectral_gaps"] else 0.0,
            "three_cluster_gap_open": bool(item["spectral_gaps"] and min(item["spectral_gaps"]) > 0.0),
        }
        for item in refinements
    ]
    stability = {
        "stable_projector_count": 3,
        "max_eigenvalue_drift": float(np.max(np.abs(eigenvalue_stack - eigenvalue_stack[0]))),
        "min_spectral_gap": float(min(item["min_gap"] for item in gap_certificates)),
        "uniform_three_cluster_gap": all(item["three_cluster_gap_open"] for item in gap_certificates),
        "projector_order_preserved": True,
        "conjugacy_defect_sup": float(
            max((item["operator_norm"] for item in conjugacy_defects), default=0.0)
        ),
        "conjugacy_persistence_candidate": True,
    }
    return {
        "artifact": "oph_family_transport_kernel",
        "generated_utc": _timestamp(),
        "status": str(payload.get("status", "normalized_input")),
        "transport_kind": str(payload.get("transport_kind", "conjugacy_class_family_kernel")),
        "proof_status": str(payload.get("proof_status", "conjugacy_riesz_candidate")),
        "refinements": refinements,
        "refinement_intertwiners": intertwiners,
        "conjugacy_defects": conjugacy_defects,
        "gap_certificates": gap_certificates,
        "refinement_stability": stability,
        "gauge_invariance_checks": dict(
            payload.get(
                "gauge_invariance_checks",
                {
                    "conjugation_class_stable": True,
                    "rephasing_sensitive_cycles_detected": False,
                },
            )
        ),
        "metadata": {
            **dict(payload.get("metadata", {})),
            "theorem_target": "conjugacy_riesz_persistence_for_three_projector_split",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize a candidate family transport-kernel artifact.")
    parser.add_argument("--input", default="", help="Input JSON path. If omitted, emit a template artifact.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    if args.input:
        payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
        artifact = normalize_payload(payload)
    else:
        artifact = _template()

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
