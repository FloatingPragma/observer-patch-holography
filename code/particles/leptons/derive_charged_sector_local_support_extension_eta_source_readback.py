#!/usr/bin/env python3
"""Emit the eta source-readback primitive beneath the charged support shell.

Chain role: isolate the first support-extension scalar on the active charged
carrier without enlarging the family.

Mathematics: weighted midpoint-defect / source-readback invariant on the
ordered three-point support-extension family.

OPH-derived inputs: the support-extension completion law, the minimal
support-extension emitter, the current support value law, and the forward
charged artifact.

Output: the eta source-readback primitive, leaving
`sigma_source_support_extension_total_log_per_side` as the next and only
same-carrier residual object.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_COMPLETION_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"
DEFAULT_MINIMAL_EXTENSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"
DEFAULT_CURRENT_VALUE_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"
DEFAULT_FORWARD = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
DEFAULT_ENDPOINT_RATIO_BREAKER = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_endpoint_ratio_breaker.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_eta_source_readback.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    completion_law: dict,
    minimal_extension: dict,
    current_value_law: dict | None = None,
    forward: dict | None = None,
    endpoint_ratio_breaker: dict | None = None,
) -> dict:
    current_value_law = current_value_law or {}
    forward = forward or {}
    endpoint_ratio_breaker = endpoint_ratio_breaker or {}

    x2 = float(completion_law["ordered_family_coordinate"][1])
    mu_source = float(completion_law["mu_source_log_per_side"])
    g_active = completion_law.get("g_active_candidate")

    zero_check: dict[str, object] = {}
    ordered = current_value_law.get("source_side_ordered_package_log_per_side_emitted")
    if ordered is not None and len(ordered) == 3:
        r1_ext, r2_ext, r3_ext = [float(value) for value in ordered]
        weighted_midpoint = (((1.0 - x2) * r1_ext) + ((1.0 + x2) * r3_ext)) / 2.0
        midpoint_defect = weighted_midpoint - r2_ext
        eta_from_ordered = ((1.0 - x2) * r1_ext) + ((1.0 + x2) * r3_ext) - (2.0 * r2_ext)
        zero_check.update(
            {
                "current_support_weighted_midpoint_log_per_side": weighted_midpoint,
                "current_support_midpoint_defect_log_per_side": midpoint_defect,
                "current_support_eta_from_ordered_package": eta_from_ordered,
                "current_support_source_side_ordered_package_log_per_side": [r1_ext, r2_ext, r3_ext],
            }
        )

    singular_values_abs = forward.get("singular_values_abs")
    if singular_values_abs is not None and len(singular_values_abs) == 3:
        masses = [float(value) for value in singular_values_abs]
        eta_from_mass_ratio = math.log((masses[0] ** (1.0 - x2) * masses[2] ** (1.0 + x2)) / (masses[1] ** 2))
        logs = [math.log(value) for value in masses]
        mean_log = sum(logs) / 3.0
        centered = [value - mean_log for value in logs]
        eta_from_centered = ((1.0 - x2) * centered[0]) + ((1.0 + x2) * centered[2]) - (2.0 * centered[1])
        zero_check.update(
            {
                "current_support_eta_from_centered_logs": eta_from_centered,
                "current_support_eta_from_mass_ratio": eta_from_mass_ratio,
                "current_support_singular_values_abs": masses,
            }
        )

    return {
        "artifact": "oph_charged_sector_local_support_extension_eta_source_readback",
        "generated_utc": _timestamp(),
        "status": "closed_smaller_primitive",
        "proof_status": "eta_source_readback_formula_closed_value_open",
        "predictive_promotion_allowed": False,
        "parent_artifacts": {
            "support_extension_completion_law": completion_law.get("artifact"),
            "minimal_source_support_extension_emitter": minimal_extension.get("artifact"),
            "endpoint_ratio_breaker": endpoint_ratio_breaker.get("artifact"),
            "current_support_value_law": current_value_law.get("artifact"),
            "forward_charged_leptons": forward.get("artifact"),
        },
        "ordered_family_coordinate": completion_law.get("ordered_family_coordinate"),
        "mu_source_log_per_side": mu_source,
        "g_active_candidate": g_active,
        "linear_basis_vector_centered": completion_law.get("linear_basis_vector_centered"),
        "extension_basis_vector_centered": completion_law.get("extension_basis_vector_centered"),
        "support_extension_denominator": completion_law.get("support_extension_denominator"),
        "eta_source_support_extension_log_per_side": None,
        "eta_readback_invariant_name": "J_mid_e_mu_tau_support_extension",
        "eta_readback_invariant_formula": "(1 - x2) * r1_ext + (1 + x2) * r3_ext - 2 * r2_ext",
        "eta_equivalence_formula": "eta_source_support_extension_log_per_side = J_mid_e_mu_tau_support_extension",
        "weighted_midpoint_formula": "((1 - x2) * r1_ext + (1 + x2) * r3_ext) / 2",
        "midpoint_defect_log_per_side_formula": "weighted_midpoint_ext - r2_ext",
        "eta_from_midpoint_defect_formula": "2 * (weighted_midpoint_ext - r2_ext)",
        "centered_log_readback_formula": "(1 - x2) * e_log_centered_ext + (1 + x2) * tau_log_centered_ext - 2 * mu_log_centered_ext",
        "mass_ratio_readback_formula": "log((e_ext^(1 - x2) * tau_ext^(1 + x2)) / mu_ext^2)",
        "shape_ratio_readback_formula": "log((shape_e_ext^(1 - x2) * shape_tau_ext^(1 + x2)) / shape_mu_ext^2)",
        "common_shift_invariant": True,
        "common_scale_invariant": True,
        "carrier_not_enlarged": True,
        "smaller_than_full_completion_shell": True,
        "source_side_ordered_package_log_per_side_ext_formula": {
            "r1_ext": completion_law["source_side_ordered_package_log_per_side_ext_formula"]["r1"],
            "r2_ext": completion_law["source_side_ordered_package_log_per_side_ext_formula"]["r2"],
            "r3_ext": completion_law["source_side_ordered_package_log_per_side_ext_formula"]["r3"],
        },
        "e_log_centered_ext_formula": completion_law.get("e_log_centered_ext_formula"),
        "mu_log_centered_ext_formula": completion_law.get("mu_log_centered_ext_formula"),
        "tau_log_centered_ext_formula": completion_law.get("tau_log_centered_ext_formula"),
        "eta_zero_on_current_support_certificate": zero_check,
        "next_single_residual_object_after_eta": "sigma_source_support_extension_total_log_per_side",
        "downstream_sigma_artifact": endpoint_ratio_breaker.get("artifact"),
        "notes": [
            "This isolates the first support-extension scalar itself as a weighted midpoint-defect invariant on the current charged carrier.",
            "It does not reuse the quarantined rigid eta candidate.",
            "It is common-shift and common-scale invariant, so it does not depend on the unresolved absolute writeback.",
            "Once eta_source_support_extension_log_per_side is emitted through this primitive, the only remaining same-carrier residual is sigma_source_support_extension_total_log_per_side via the existing endpoint-ratio breaker.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged eta source-readback primitive.")
    parser.add_argument("--completion-law", default=str(DEFAULT_COMPLETION_LAW))
    parser.add_argument("--minimal-extension", default=str(DEFAULT_MINIMAL_EXTENSION))
    parser.add_argument("--current-value-law", default=str(DEFAULT_CURRENT_VALUE_LAW))
    parser.add_argument("--forward", default=str(DEFAULT_FORWARD))
    parser.add_argument("--endpoint-ratio-breaker", default=str(DEFAULT_ENDPOINT_RATIO_BREAKER))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    completion_law = json.loads(Path(args.completion_law).read_text(encoding="utf-8"))
    minimal_extension = json.loads(Path(args.minimal_extension).read_text(encoding="utf-8"))
    current_value_law_path = Path(args.current_value_law)
    current_value_law = json.loads(current_value_law_path.read_text(encoding="utf-8")) if current_value_law_path.exists() else None
    forward_path = Path(args.forward)
    forward = json.loads(forward_path.read_text(encoding="utf-8")) if forward_path.exists() else None
    endpoint_ratio_breaker_path = Path(args.endpoint_ratio_breaker)
    endpoint_ratio_breaker = json.loads(endpoint_ratio_breaker_path.read_text(encoding="utf-8")) if endpoint_ratio_breaker_path.exists() else None

    artifact = build_artifact(completion_law, minimal_extension, current_value_law, forward, endpoint_ratio_breaker)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
