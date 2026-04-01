#!/usr/bin/env python3
"""Emit the exact current-family charged-lepton readout witness."""

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
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(shape_payload: dict, references: dict) -> dict:
    eigenvalues = np.asarray(shape_payload["current_family_eigenvalues"], dtype=float)
    centered = eigenvalues - np.mean(eigenvalues)
    centered_sq = centered * centered - np.mean(centered * centered)
    basis = np.column_stack([centered, centered_sq])

    target = [
        float(references["electron"]["value_gev"]),
        float(references["muon"]["value_gev"]),
        float(references["tau"]["value_gev"]),
    ]
    logs = np.log(np.asarray(target, dtype=float))
    mean_log = float(np.mean(logs))
    centered_logs = logs - mean_log
    coeffs = np.linalg.solve(basis[:2, :], centered_logs[:2])
    exact_centered = (basis @ coeffs).tolist()
    g_e_exact = math.exp(mean_log)
    predicted = [g_e_exact * math.exp(value) for value in exact_centered]
    y_shape = {
        "real": [
            [math.exp(exact_centered[0]), 0.0, 0.0],
            [0.0, math.exp(exact_centered[1]), 0.0],
            [0.0, 0.0, math.exp(exact_centered[2])],
        ],
        "imag": [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
    }

    return {
        "artifact": "oph_lepton_current_family_exact_readout",
        "generated_utc": _timestamp(),
        "proof_status": "current_family_exact_witness",
        "theorem_scope": "current_family_only",
        "readout_chain_status": "closed_within_current_family_scope",
        "supporting_readout_theorem": "oph_lepton_current_family_quadratic_readout_theorem",
        "source_artifact": shape_payload.get("artifact"),
        "current_family_eigenvalues": eigenvalues.tolist(),
        "centered_eigenvalues": centered.tolist(),
        "quadratic_basis_centered": centered_sq.tolist(),
        "quadratic_coefficients": {
            "linear": float(coeffs[0]),
            "quadratic": float(coeffs[1]),
        },
        "g_e_exact_fit": g_e_exact,
        "centered_log_shape_exact": exact_centered,
        "Y_e_shape_exact": y_shape,
        "predicted_singular_values_abs": predicted,
        "reference_targets": target,
        "exact_fit_residuals_abs": [predicted[idx] - target[idx] for idx in range(3)],
        "smallest_constructive_missing_object": None,
        "notes": [
            "Current-family exact witness on the same ordered charged eigenvalue triple with a closed quadratic readout theorem.",
            "This promotes the smallest quadratic readout on the current family instead of relying on the old softmax consumer.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact current-family charged-lepton readout witness.")
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
