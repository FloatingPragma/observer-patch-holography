#!/usr/bin/env python3
"""Emit the closed current-family quark quadratic readout theorem artifact."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _centered_logs(values: list[float]) -> tuple[np.ndarray, float]:
    logs = np.log(np.asarray(values, dtype=float))
    mean_log = float(np.mean(logs))
    return logs - mean_log, mean_log


def build_artifact(spread_map: dict, references: dict) -> dict:
    x = np.asarray([-1.0, float(spread_map["normalized_coordinate_x2"]), 1.0], dtype=float)
    x_centered = x - np.mean(x)
    x2_centered = x * x - np.mean(x * x)
    basis = np.column_stack([x_centered, x2_centered])
    basis_minor = basis[:2, :]
    basis_det = float(np.linalg.det(basis_minor))

    target_u = [
        float(references["up_quark"]["value_gev"]),
        float(references["charm_quark"]["value_gev"]),
        float(references["top_quark"]["value_gev"]),
    ]
    target_d = [
        float(references["down_quark"]["value_gev"]),
        float(references["strange_quark"]["value_gev"]),
        float(references["bottom_quark"]["value_gev"]),
    ]
    centered_u, mean_log_u = _centered_logs(target_u)
    centered_d, mean_log_d = _centered_logs(target_d)
    coeff_u = np.linalg.solve(basis_minor, centered_u[:2])
    coeff_d = np.linalg.solve(basis_minor, centered_d[:2])
    exact_u = (basis @ coeff_u).tolist()
    exact_d = (basis @ coeff_d).tolist()
    g_u_exact = math.exp(mean_log_u)
    g_d_exact = math.exp(mean_log_d)
    predicted_u = [g_u_exact * math.exp(value) for value in exact_u]
    predicted_d = [g_d_exact * math.exp(value) for value in exact_d]

    return {
        "artifact": "oph_quark_current_family_quadratic_readout_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_exact_readout",
        "theorem_scope": "current_family_only",
        "input_kind": "ordered_three_point_quark_carrier_plus_target_family",
        "ordered_carrier_artifact": spread_map.get("artifact"),
        "target_source_artifact": "particle_reference_values.json",
        "theorem_statement": (
            "On the fixed ordered three-point quark carrier with x = (-1, x2, 1) and x2 != +/-1, "
            "the centered log readout is uniquely quadratic: every centered log triple y with sum(y)=0 "
            "has a unique decomposition y = a * ctr(X_ord) + b * ctr(X_ord^2). "
            "For the current-family quark targets, the unique centered logs together with the sector geometric means "
            "g_u and g_d therefore emit one exact same-family sextet."
        ),
        "ordered_coordinate_x": x.tolist(),
        "quadratic_basis": {
            "linear": x_centered.tolist(),
            "quadratic_centered": x2_centered.tolist(),
            "minor_first_two_rows": basis_minor.tolist(),
            "minor_determinant": basis_det,
            "invertible": abs(basis_det) > 1.0e-12,
        },
        "uniqueness_contract": {
            "centered_log_formula": "y = a * ctr(X_ord) + b * ctr(X_ord^2)",
            "coefficient_solver": "[a, b]^T = B^{-1} * [y1, y2]^T where B is the first-two-row quadratic basis minor",
            "absolute_scale_solver_u": "g_u = exp(mean(log target_u))",
            "absolute_scale_solver_d": "g_d = exp(mean(log target_d))",
        },
        "sector_exact_coefficients_u": {
            "linear": float(coeff_u[0]),
            "quadratic": float(coeff_u[1]),
        },
        "sector_exact_coefficients_d": {
            "linear": float(coeff_d[0]),
            "quadratic": float(coeff_d[1]),
        },
        "sector_exact_geometric_means": {
            "g_u": g_u_exact,
            "g_d": g_d_exact,
        },
        "target_centered_log_u": centered_u.tolist(),
        "target_centered_log_d": centered_d.tolist(),
        "exact_centered_log_u": exact_u,
        "exact_centered_log_d": exact_d,
        "reference_targets_u": target_u,
        "reference_targets_d": target_d,
        "predicted_singular_values_u": predicted_u,
        "predicted_singular_values_d": predicted_d,
        "exact_fit_residuals_u": [predicted_u[idx] - target_u[idx] for idx in range(3)],
        "exact_fit_residuals_d": [predicted_d[idx] - target_d[idx] for idx in range(3)],
        "smallest_constructive_missing_object": None,
        "notes": [
            "This closes the target-anchored same-family quark readout chain on the fixed ordered three-point carrier.",
            "It does not change the separate theorem-grade physical-branch boundary on the selected D12 sheet.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the closed current-family quark quadratic readout theorem.")
    parser.add_argument("--spread-map", default=str(SPREAD_MAP_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    artifact = build_artifact(spread_map, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
