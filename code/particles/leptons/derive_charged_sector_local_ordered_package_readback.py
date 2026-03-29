#!/usr/bin/env python3
"""Emit the charged sector-local ordered-package readback artifact."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_EMISSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_source_emission.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_readback.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _diag_payload(diagonal_values: list[float]) -> dict[str, list[list[float]]]:
    real = [[0.0, 0.0, 0.0] for _ in range(3)]
    imag = [[0.0, 0.0, 0.0] for _ in range(3)]
    for idx, value in enumerate(diagonal_values):
        real[idx][idx] = float(value)
    return {"real": real, "imag": imag}


def build_artifact(source_emission: dict[str, object]) -> dict[str, object]:
    ordered_package = [float(value) for value in source_emission["source_side_ordered_package_log_per_side_emitted"]]
    x2 = float(source_emission["ordered_family_coordinate"][1])
    mean_log = sum(ordered_package) / 3.0
    centered = [value - mean_log for value in ordered_package]
    r1, r2, r3 = ordered_package
    quadratic_midpoint_defect = (((1.0 - x2) * r1) + ((1.0 + x2) * r3)) / 2.0 - r2
    denom = 1.0 - x2 * x2
    b_e = quadratic_midpoint_defect / denom if abs(denom) > 1.0e-15 else 0.0
    a_e = (r3 - r1) / 2.0
    sigma_e = r3 - r1
    eta_e = 2.0 * quadratic_midpoint_defect
    shape_singular_values = [math.exp(value) for value in centered]
    reconstruction = [
        a_e * (-(3.0 + x2) / 3.0) + b_e * ((1.0 - x2 * x2) / 3.0),
        a_e * ((2.0 * x2) / 3.0) + b_e * (-(2.0 * (1.0 - x2 * x2)) / 3.0),
        a_e * ((3.0 - x2) / 3.0) + b_e * ((1.0 - x2 * x2) / 3.0),
    ]
    residual = [centered[idx] - reconstruction[idx] for idx in range(3)]

    return {
        "artifact": "oph_charged_sector_local_ordered_package_readback",
        "generated_utc": _timestamp(),
        "proof_status": "current_package_readback_closed",
        "predictive_promotion_allowed": False,
        "input_artifact": source_emission.get("artifact"),
        "input_kind": "sector_local_ordered_package_log_per_side",
        "output_kind": "charged_affine_quadratic_source_readback",
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "source_side_ordered_package_log_per_side_emitted": ordered_package,
        "source_side_ordered_package_mean_log_per_side_emitted": mean_log,
        "source_side_ordered_package_centered_log_emitted": centered,
        "source_side_quadratic_midpoint_defect_log_per_side_emitted": quadratic_midpoint_defect,
        "a_e_log_coeff": a_e,
        "b_e_log_coeff": b_e,
        "sigma_e_total_log_per_side_emitted": sigma_e,
        "eta_e_split_log_per_side_emitted": eta_e,
        "shape_singular_values_emitted": shape_singular_values,
        "Y_e_shape_emitted": _diag_payload(shape_singular_values),
        "a_e_log_coeff_formula": "(r3 - r1) / 2",
        "b_e_log_coeff_formula": "(((1 - x2) * r1 + (1 + x2) * r3) / 2 - r2) / (1 - x2^2)",
        "sigma_e_total_log_per_side_emitted_formula": "r3 - r1",
        "eta_e_split_log_per_side_emitted_formula": "2 * (((1 - x2) * r1 + (1 + x2) * r3) / 2 - r2)",
        "centered_source_formula": "[r1 - r_bar, r2 - r_bar, r3 - r_bar]",
        "shape_singular_values_emitted_formula": "exp(source_side_ordered_package_centered_log_emitted)",
        "affine_quadratic_reconstruction_residual": residual,
        "current_package_linear_subray_only": abs(quadratic_midpoint_defect) < 1.0e-12,
        "same_support_exhausted": abs(quadratic_midpoint_defect) < 1.0e-12,
        "smallest_constructive_missing_object": "oph_charged_sector_local_current_support_obstruction_certificate",
        "notes": [
            "This artifact reads back the affine-quadratic charged coefficients from the current sector-local ordered package.",
            "The current package still lands on the linear subray when the midpoint defect vanishes.",
            "The next remaining charged gap is now the obstruction certificate that shows the current support itself is too small, not another same-support midpoint-defect scalar.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged sector-local ordered-package readback artifact.")
    parser.add_argument("--source-emission", default=str(DEFAULT_SOURCE_EMISSION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_emission = json.loads(Path(args.source_emission).read_text(encoding="utf-8"))
    artifact = build_artifact(source_emission)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
