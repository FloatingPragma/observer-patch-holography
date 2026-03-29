#!/usr/bin/env python3
"""Emit the charged source ordered-package artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEAN_CERT = ROOT / "particles" / "runs" / "flavor" / "charged_mean_eigenvalue_certificate.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "flavor" / "family_excitation_evaluator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_source_emission.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(mean_certificate: dict, family: dict) -> dict:
    unordered = [float(value) for value in mean_certificate["current_family_eigenvalues"]]
    ordering = sorted(range(len(unordered)), key=lambda idx: unordered[idx])
    ordered = [unordered[idx] for idx in ordering]
    mean_log = sum(ordered) / 3.0
    x2 = float(family["family_coordinate_x"][1])
    return {
        "artifact": "oph_charged_sector_local_ordered_package_source_emission",
        "generated_utc": _timestamp(),
        "proof_status": "current_family_source_emission_closed",
        "predictive_promotion_allowed": False,
        "source_artifact": mean_certificate.get("artifact"),
        "source_package_origin": "common_refinement_eigenvalue_readout",
        "left_common_family_eigenvalues": mean_certificate.get("left_common_family_eigenvalues"),
        "right_common_family_eigenvalues": mean_certificate.get("right_common_family_eigenvalues"),
        "current_family_eigenvalues": unordered,
        "ordering_permutation": ordering,
        "source_side_unordered_package_log_per_side_emitted": unordered,
        "source_side_ordered_package_log_per_side_emitted": ordered,
        "source_side_ordered_package_mean_log_per_side_emitted": mean_log,
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "rho_ord": float(family["rho_ord"]),
        "readout_formulas": {
            "unordered_package": "current_family_eigenvalues",
            "ordered_package": "sort(current_family_eigenvalues)",
            "mean_log": "sum(source_side_ordered_package_log_per_side_emitted) / 3",
        },
        "smallest_constructive_missing_object": "oph_charged_sector_local_ordered_package_value_law",
        "notes": [
            "This artifact exposes the current charged ordered package directly from the common-refinement eigenvalue package.",
            "It provides the upstream ordered package without adding a larger charged carrier family.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged sector-local ordered-package source emission artifact.")
    parser.add_argument("--mean-certificate", default=str(DEFAULT_MEAN_CERT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    mean_certificate = json.loads(Path(args.mean_certificate).read_text(encoding="utf-8"))
    family = json.loads(Path(args.family).read_text(encoding="utf-8"))
    artifact = build_artifact(mean_certificate, family)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
