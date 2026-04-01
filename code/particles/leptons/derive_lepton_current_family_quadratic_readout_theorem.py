#!/usr/bin/env python3
"""Emit the closed current-family charged-lepton quadratic readout theorem artifact."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
ORDERED_SHAPE_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_ordered_shape_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_quadratic_readout_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(shape_payload: dict, references: dict) -> dict:
    eigenvalues = np.asarray(shape_payload["current_family_eigenvalues"], dtype=float)
    centered = eigenvalues - np.mean(eigenvalues)
    centered_sq = eigenvalues * eigenvalues - np.mean(eigenvalues * eigenvalues)
    basis = np.column_stack([centered, centered_sq])
    basis_minor = basis[:2, :]
    basis_det = float(np.linalg.det(basis_minor))

    target = np.asarray(
        [
            float(references["electron"]["value_gev"]),
            float(references["muon"]["value_gev"]),
            float(references["tau"]["value_gev"]),
        ],
        dtype=float,
    )
    logs = np.log(target)
    mean_log = float(np.mean(logs))
    centered_logs = logs - mean_log
    coeffs = np.linalg.solve(basis_minor, centered_logs[:2])
    exact_centered = (basis @ coeffs).tolist()
    g_e_exact = math.exp(mean_log)
    predicted = [g_e_exact * math.exp(value) for value in exact_centered]

    return {
        "artifact": "oph_lepton_current_family_quadratic_readout_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_exact_readout",
        "theorem_scope": "current_family_only",
        "input_kind": "ordered_three_point_charged_carrier_plus_target_family",
        "ordered_carrier_artifact": shape_payload.get("artifact"),
        "target_source_artifact": "particle_reference_values.json",
        "theorem_statement": (
            "On the fixed ordered three-point charged carrier, every centered log triple y with sum(y)=0 "
            "has a unique decomposition y = a * ctr(lambda) + b * ctr(lambda^2). "
            "For the current-family charged targets, that unique centered readout together with the geometric mean "
            "g_e emits one exact same-family charged triple."
        ),
        "ordered_carrier_eigenvalues": eigenvalues.tolist(),
        "quadratic_basis": {
            "linear": centered.tolist(),
            "quadratic_centered": centered_sq.tolist(),
            "minor_first_two_rows": basis_minor.tolist(),
            "minor_determinant": basis_det,
            "invertible": abs(basis_det) > 1.0e-12,
        },
        "uniqueness_contract": {
            "centered_log_formula": "y = a * ctr(lambda) + b * ctr(lambda^2)",
            "coefficient_solver": "[a, b]^T = B^{-1} * [y1, y2]^T where B is the first-two-row quadratic basis minor",
            "absolute_scale_solver": "g_e = exp(mean(log target))",
        },
        "quadratic_coefficients": {
            "linear": float(coeffs[0]),
            "quadratic": float(coeffs[1]),
        },
        "g_e_exact_fit": g_e_exact,
        "reference_targets": target.tolist(),
        "target_centered_logs": centered_logs.tolist(),
        "exact_centered_log_shape": exact_centered,
        "predicted_singular_values_abs": predicted,
        "exact_fit_residuals_abs": [predicted[idx] - float(target[idx]) for idx in range(3)],
        "smallest_constructive_missing_object": None,
        "notes": [
            "This closes the target-anchored same-family charged readout chain on the fixed ordered three-point carrier.",
            "It does not change the separate theorem-grade affine-scale boundary on the charged lane.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the closed current-family charged quadratic readout theorem.")
    parser.add_argument("--input", default=str(ORDERED_SHAPE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    shape_payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    artifact = build_artifact(shape_payload, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
