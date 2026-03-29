#!/usr/bin/env python3
"""Emit the endpoint-ratio breaker beneath the charged support-extension shell.

Chain role: isolate the unique same-carrier scalar that breaks the fixed
endpoint ratio on the charged support-extension family once eta is known.

Mathematics: endpoint-log subtraction on the existing two-scalar ordered
support-extension carrier.

OPH-derived inputs: the minimal support-extension emitter and the full
two-scalar support-extension completion law.

Output: the smaller endpoint-ratio-breaker primitive beneath the charged
support-extension completion shell.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MINIMAL_EXTENSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"
DEFAULT_COMPLETION_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_endpoint_ratio_breaker.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(minimal_extension: dict, completion_law: dict) -> dict:
    return {
        "artifact": "oph_charged_sector_local_support_extension_endpoint_ratio_breaker",
        "generated_utc": _timestamp(),
        "status": "closed_smaller_primitive",
        "proof_status": "endpoint_ratio_breaker_formula_closed_scalar_open",
        "predictive_promotion_allowed": False,
        "parent_artifacts": {
            "minimal_source_support_extension_emitter": minimal_extension.get("artifact"),
            "support_extension_completion_law": completion_law.get("artifact"),
        },
        "ordered_family_coordinate": minimal_extension.get("ordered_family_coordinate"),
        "linear_basis_vector_centered": minimal_extension.get("linear_basis_vector_centered"),
        "extension_basis_vector_centered": minimal_extension.get("extension_basis_vector_centered"),
        "mu_source_log_per_side": minimal_extension.get("mu_source_log_per_side"),
        "g_active_candidate": minimal_extension.get("g_active_candidate"),
        "precondition_residual_object": "eta_source_support_extension_log_per_side",
        "sigma_source_support_extension_total_log_per_side": None,
        "endpoint_ratio_breaker_invariant_name": "J_endpoint_e_tau_support_extension",
        "endpoint_ratio_breaker_invariant_formula": "log(tau_ext / e_ext)",
        "sigma_equivalence_formula": "sigma_source_support_extension_total_log_per_side = J_endpoint_e_tau_support_extension",
        "endpoint_ratio_formula": "tau_ext / e_ext = exp(sigma_source_support_extension_total_log_per_side)",
        "kappa_ext_formula": completion_law.get("kappa_ext_formula"),
        "centered_extension_formula": completion_law.get("centered_extension_formula"),
        "gamma21_log_per_side_ext_formula": completion_law.get("gamma21_log_per_side_ext_formula"),
        "gamma32_log_per_side_ext_formula": completion_law.get("gamma32_log_per_side_ext_formula"),
        "e_log_centered_ext_formula": completion_law.get("e_log_centered_ext_formula"),
        "mu_log_centered_ext_formula": completion_law.get("mu_log_centered_ext_formula"),
        "tau_log_centered_ext_formula": completion_law.get("tau_log_centered_ext_formula"),
        "e_ext_formula": completion_law.get("e_ext_formula"),
        "mu_ext_formula": completion_law.get("mu_ext_formula"),
        "tau_ext_formula": completion_law.get("tau_ext_formula"),
        "smallest_constructive_missing_object_within_primitive": "sigma_source_support_extension_total_log_per_side",
        "notes": [
            "This is the unique same-carrier endpoint-ratio breaker beneath the eta-only support-extension shell.",
            "It does not enlarge the charged family; it only identifies the remaining span scalar once eta is known.",
            "After eta_source_support_extension_log_per_side is emitted, this scalar fixes tau/e on the support-extension carrier.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged support-extension endpoint-ratio breaker artifact.")
    parser.add_argument("--minimal-extension", default=str(DEFAULT_MINIMAL_EXTENSION))
    parser.add_argument("--completion-law", default=str(DEFAULT_COMPLETION_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    minimal_extension = json.loads(Path(args.minimal_extension).read_text(encoding="utf-8"))
    completion_law = json.loads(Path(args.completion_law).read_text(encoding="utf-8"))
    artifact = build_artifact(minimal_extension, completion_law)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
