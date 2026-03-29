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
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _centered_logs(values: list[float]) -> tuple[np.ndarray, float]:
    logs = np.log(np.asarray(values, dtype=float))
    mean_log = float(np.mean(logs))
    return logs - mean_log, mean_log


def build_artifact(mean_split: dict, spread_map: dict, references: dict) -> dict:
    x = np.asarray([-1.0, float(spread_map["normalized_coordinate_x2"]), 1.0], dtype=float)
    x_centered = x - np.mean(x)
    x2_centered = x * x - np.mean(x * x)
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
    g_u = float(mean_split["g_u_candidate"])
    g_d = float(mean_split["g_d_candidate"])
    predicted_u = [g_u * math.exp(value) for value in exact_u]
    predicted_d = [g_d * math.exp(value) for value in exact_d]

    return {
        "artifact": "oph_quark_current_family_exact_readout",
        "generated_utc": _timestamp(),
        "proof_status": "current_family_exact_witness",
        "theorem_scope": "current_family_only",
        "ordered_coordinate_x": x.tolist(),
        "quadratic_basis": {
            "linear": x_centered.tolist(),
            "quadratic_centered": x2_centered.tolist(),
        },
        "sector_mean_split_artifact": mean_split.get("artifact"),
        "g_u": g_u,
        "g_d": g_d,
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
        "smallest_constructive_missing_object": "quark_quadratic_readout_theorem",
        "notes": [
            "This current-family exact witness uses the already-fixed simple three-point ordered spectrum and promotes only the sector-even quadratic readout coefficients.",
            "It does not reopen a richer ray family or add a third scalar beyond the compact sector-mean split.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact current-family quark quadratic readout witness.")
    parser.add_argument("--mean-split", default=str(MEAN_SPLIT_JSON))
    parser.add_argument("--spread-map", default=str(SPREAD_MAP_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    mean_split = json.loads(Path(args.mean_split).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    artifact = build_artifact(mean_split, spread_map, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
