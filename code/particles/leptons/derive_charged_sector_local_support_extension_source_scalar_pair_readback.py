#!/usr/bin/env python3
"""Emit the smaller charged support-extension source-scalar pair readback.

Chain role: collect the first two same-carrier support-extension scalar
readbacks into one machine-readable primitive beneath the larger completion
shell.

Mathematics: invariant readback on the fixed charged carrier. `eta` is the
weighted midpoint-defect invariant and `sigma` is the endpoint-ratio
invariant.

OPH-derived inputs: the charged support-extension eta source-readback, the
charged endpoint-ratio breaker, and the full support-extension completion law.

Output: a smaller source-scalar pair readback artifact that exposes the ordered
eta-then-sigma frontier without pretending those values are emitted.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ETA_READBACK = (
    ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_eta_source_readback.json"
)
DEFAULT_ENDPOINT_RATIO_BREAKER = (
    ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_endpoint_ratio_breaker.json"
)
DEFAULT_COMPLETION_LAW = (
    ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_source_scalar_pair_readback.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(eta_readback: dict, endpoint_ratio_breaker: dict, completion_law: dict) -> dict:
    return {
        "artifact": "oph_charged_sector_local_support_extension_source_scalar_pair_readback",
        "generated_utc": _timestamp(),
        "status": "closed_smaller_primitive",
        "proof_status": "source_scalar_pair_readback_formula_closed_values_open",
        "smaller_than_full_completion_shell": True,
        "carrier_not_enlarged": True,
        "scalar_order": [
            "eta_source_support_extension_log_per_side",
            "sigma_source_support_extension_total_log_per_side",
        ],
        "next_single_residual_object": "eta_source_support_extension_log_per_side",
        "next_single_residual_object_after_eta": "sigma_source_support_extension_total_log_per_side",
        "eta_source_support_extension_log_per_side": None,
        "sigma_source_support_extension_total_log_per_side": None,
        "eta_readback_invariant_name": eta_readback.get("eta_readback_invariant_name"),
        "eta_equivalence_formula": eta_readback.get("eta_equivalence_formula"),
        "sigma_readback_invariant_name": endpoint_ratio_breaker.get("endpoint_ratio_breaker_invariant_name"),
        "sigma_equivalence_formula": endpoint_ratio_breaker.get("sigma_equivalence_formula"),
        "endpoint_ratio_formula": endpoint_ratio_breaker.get("endpoint_ratio_formula"),
        "endpoint_ratio_breaker_invariant_formula": endpoint_ratio_breaker.get("endpoint_ratio_breaker_invariant_formula"),
        "weighted_midpoint_formula": eta_readback.get("weighted_midpoint_formula"),
        "midpoint_defect_log_per_side_formula": eta_readback.get("midpoint_defect_log_per_side_formula"),
        "source_side_ordered_package_log_per_side_ext_formula": completion_law.get(
            "source_side_ordered_package_log_per_side_ext_formula"
        ),
        "centered_extension_formula": completion_law.get("centered_extension_formula"),
        "e_log_centered_ext_formula": completion_law.get("e_log_centered_ext_formula"),
        "mu_log_centered_ext_formula": completion_law.get("mu_log_centered_ext_formula"),
        "tau_log_centered_ext_formula": completion_law.get("tau_log_centered_ext_formula"),
        "e_ext_formula": completion_law.get("e_ext_formula"),
        "mu_ext_formula": completion_law.get("mu_ext_formula"),
        "tau_ext_formula": completion_law.get("tau_ext_formula"),
        "ordered_family_coordinate": completion_law.get("ordered_family_coordinate"),
        "linear_basis_vector_centered": completion_law.get("linear_basis_vector_centered"),
        "extension_basis_vector_centered": completion_law.get("extension_basis_vector_centered"),
        "g_active_candidate": completion_law.get("g_active_candidate"),
        "mu_source_log_per_side": completion_law.get("mu_source_log_per_side"),
        "notes": [
            "This is the smallest same-carrier readback primitive that keeps both charged support-extension source scalars together.",
            "It closes the invariant formulas for the ordered eta-then-sigma pair without promoting any numerical source values.",
            "The first residual is still eta_source_support_extension_log_per_side; once eta is emitted, sigma_source_support_extension_total_log_per_side remains as the unique endpoint-ratio breaker on the same carrier.",
        ],
        "parent_artifacts": {
            "eta_source_readback": eta_readback.get("artifact"),
            "endpoint_ratio_breaker": endpoint_ratio_breaker.get("artifact"),
            "support_extension_completion_law": completion_law.get("artifact"),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged support-extension source-scalar pair readback artifact.")
    parser.add_argument("--eta-readback", default=str(DEFAULT_ETA_READBACK))
    parser.add_argument("--endpoint-ratio-breaker", default=str(DEFAULT_ENDPOINT_RATIO_BREAKER))
    parser.add_argument("--completion-law", default=str(DEFAULT_COMPLETION_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    eta_readback = json.loads(Path(args.eta_readback).read_text(encoding="utf-8"))
    endpoint_ratio_breaker = json.loads(Path(args.endpoint_ratio_breaker).read_text(encoding="utf-8"))
    completion_law = json.loads(Path(args.completion_law).read_text(encoding="utf-8"))
    artifact = build_artifact(eta_readback, endpoint_ratio_breaker, completion_law)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
