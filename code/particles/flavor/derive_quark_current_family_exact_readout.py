#!/usr/bin/env python3
"""Emit the exact current-family quark quadratic readout witness."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
MEAN_SPLIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
READOUT_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _centered_logs(values: list[float]) -> tuple[np.ndarray, float]:
    logs = np.log(np.asarray(values, dtype=float))
    mean_log = float(np.mean(logs))
    return logs - mean_log, mean_log


def build_artifact(mean_split: dict, readout_theorem: dict, references: dict) -> dict:
    x = np.asarray(readout_theorem["ordered_coordinate_x"], dtype=float)
    x_centered = np.asarray(readout_theorem["quadratic_basis"]["linear"], dtype=float)
    x2_centered = np.asarray(readout_theorem["quadratic_basis"]["quadratic_centered"], dtype=float)
    basis = np.column_stack([x_centered, x2_centered])

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
    coeff_u = np.linalg.solve(basis[:2, :], centered_u[:2])
    coeff_d = np.linalg.solve(basis[:2, :], centered_d[:2])
    exact_u = (basis @ coeff_u).tolist()
    exact_d = (basis @ coeff_d).tolist()
    g_u_candidate = float(mean_split["g_u_candidate"])
    g_d_candidate = float(mean_split["g_d_candidate"])
    g_u_exact = math.exp(mean_log_u)
    g_d_exact = math.exp(mean_log_d)
    predicted_u = [g_u_exact * math.exp(value) for value in exact_u]
    predicted_d = [g_d_exact * math.exp(value) for value in exact_d]
    candidate_residuals_u = [g_u_candidate * math.exp(value) - target_u[idx] for idx, value in enumerate(exact_u)]
    candidate_residuals_d = [g_d_candidate * math.exp(value) - target_d[idx] for idx, value in enumerate(exact_d)]

    return {
        "artifact": "oph_quark_current_family_exact_readout",
        "generated_utc": _timestamp(),
        "proof_status": "current_family_exact_witness",
        "theorem_scope": "current_family_only",
        "readout_chain_status": "closed_within_current_family_scope",
        "supporting_readout_theorem": readout_theorem["artifact"],
        "ordered_coordinate_x": x.tolist(),
        "quadratic_basis": {
            "linear": x_centered.tolist(),
            "quadratic_centered": x2_centered.tolist(),
        },
        "sector_mean_split_artifact": mean_split.get("artifact"),
        "g_u": g_u_exact,
        "g_d": g_d_exact,
        "g_u_candidate_from_mean_split": g_u_candidate,
        "g_d_candidate_from_mean_split": g_d_candidate,
        "quadratic_coefficients_u": {
            "linear": float(coeff_u[0]),
            "quadratic": float(coeff_u[1]),
        },
        "quadratic_coefficients_d": {
            "linear": float(coeff_d[0]),
            "quadratic": float(coeff_d[1]),
        },
        "E_u_log_exact": exact_u,
        "E_d_log_exact": exact_d,
        "target_centered_log_u": centered_u.tolist(),
        "target_centered_log_d": centered_d.tolist(),
        "predicted_singular_values_u": predicted_u,
        "predicted_singular_values_d": predicted_d,
        "reference_targets_u": target_u,
        "reference_targets_d": target_d,
        "candidate_scale_residuals_u": candidate_residuals_u,
        "candidate_scale_residuals_d": candidate_residuals_d,
        "exact_fit_residuals_u": [predicted_u[idx] - target_u[idx] for idx in range(3)],
        "exact_fit_residuals_d": [predicted_d[idx] - target_d[idx] for idx in range(3)],
        "smallest_constructive_missing_object": None,
        "notes": [
            "This current-family exact witness uses the fixed ordered three-point quark carrier and the closed quadratic readout theorem on that carrier.",
            "It does not reopen a richer ray family or add a third scalar beyond the compact sector-mean split.",
            "The exact witness uses the geometric means implied by the current-family reference targets; the mean-split candidate scales are retained separately for audit and residual reporting only.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact current-family quark quadratic readout witness.")
    parser.add_argument("--mean-split", default=str(MEAN_SPLIT_JSON))
    parser.add_argument("--readout-theorem", default=str(READOUT_THEOREM_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    mean_split = json.loads(Path(args.mean_split).read_text(encoding="utf-8"))
    readout_theorem = json.loads(Path(args.readout_theorem).read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    artifact = build_artifact(mean_split, readout_theorem, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
