#!/usr/bin/env python3
"""Emit the full two-scalar charged support-extension completion law shell.

Chain role: expose the smallest exact charged completion law on the active
ordered support-extension family without inventing the missing source values.

Mathematics: ordered-gap reconstruction on the linear basis plus the canonical
quadratic support-extension direction, with two scalar source slots
`(sigma_ext, eta_ext)`.

OPH-derived inputs: the charged ordered-package value law, the current-support
obstruction certificate, the minimal support-extension emitter shell, and the
shared/local scale binding.

Output: the full two-scalar completion-law shell beneath the charged exactness
audit, together with the live residual scalar order.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VALUE_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"
DEFAULT_OBSTRUCTION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_current_support_obstruction_certificate.json"
DEFAULT_MINIMAL_EXTENSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"
DEFAULT_SCALE_BINDING = ROOT / "particles" / "runs" / "leptons" / "lepton_shared_absolute_scale_binding.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(value_law: dict, obstruction: dict, minimal_extension: dict, scale_binding: dict | None = None) -> dict:
    scale_binding = scale_binding or {}
    x2 = float(value_law["ordered_family_coordinate"][1])
    mu_source = float(value_law["mu_source_log_per_side_readback"])
    linear_basis = [float(value) for value in value_law["linear_basis_vector_centered"]]
    extension_basis = [float(value) for value in value_law["quadratic_basis_vector"]]
    denominator = 2.0 * (1.0 - x2 * x2)
    return {
        "artifact": "oph_charged_sector_local_support_extension_completion_law",
        "generated_utc": _timestamp(),
        "proof_status": "two_scalar_support_extension_completion_law_closed_source_scalars_open",
        "predictive_promotion_allowed": False,
        "source_artifacts": {
            "ordered_package_value_law": value_law.get("artifact"),
            "current_support_obstruction_certificate": obstruction.get("artifact"),
            "minimal_source_support_extension_emitter": minimal_extension.get("artifact"),
            "shared_absolute_scale_binding": scale_binding.get("artifact"),
        },
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "mu_source_log_per_side": mu_source,
        "linear_basis_vector_centered": linear_basis,
        "extension_basis_vector_centered": extension_basis,
        "support_extension_denominator": denominator,
        "g_active_candidate": scale_binding.get("g_e"),
        "sigma_source_support_extension_total_log_per_side": None,
        "eta_source_support_extension_log_per_side": None,
        "kappa_ext": None,
        "centered_extension_formula": "E_e_log_centered_ext = (sigma_source_support_extension_total_log_per_side / 2) * linear_basis_vector_centered + kappa_ext * extension_basis_vector_centered",
        "kappa_ext_formula": "eta_source_support_extension_log_per_side / (2 * (1 - x2^2))",
        "gamma21_log_per_side_ext_formula": "((1 + x2) * sigma_source_support_extension_total_log_per_side - eta_source_support_extension_log_per_side) / 2",
        "gamma32_log_per_side_ext_formula": "((1 - x2) * sigma_source_support_extension_total_log_per_side + eta_source_support_extension_log_per_side) / 2",
        "e_log_centered_ext_formula": "-((3 + x2) * sigma_source_support_extension_total_log_per_side - eta_source_support_extension_log_per_side) / 6",
        "mu_log_centered_ext_formula": "(x2 * sigma_source_support_extension_total_log_per_side - eta_source_support_extension_log_per_side) / 3",
        "tau_log_centered_ext_formula": "((3 - x2) * sigma_source_support_extension_total_log_per_side + eta_source_support_extension_log_per_side) / 6",
        "source_side_ordered_package_log_per_side_ext_formula": {
            "r1": "mu_source_log_per_side - ((3 + x2) * sigma_source_support_extension_total_log_per_side - eta_source_support_extension_log_per_side) / 6",
            "r2": "mu_source_log_per_side + (x2 * sigma_source_support_extension_total_log_per_side - eta_source_support_extension_log_per_side) / 3",
            "r3": "mu_source_log_per_side + ((3 - x2) * sigma_source_support_extension_total_log_per_side + eta_source_support_extension_log_per_side) / 6",
        },
        "shape_singular_values_ext_formula": "exp(E_e_log_centered_ext)",
        "singular_values_abs_ext_formula": "g_active_candidate * exp(E_e_log_centered_ext)",
        "e_ext_formula": "g_active_candidate * exp(e_log_centered_ext)",
        "mu_ext_formula": "g_active_candidate * exp(mu_log_centered_ext)",
        "tau_ext_formula": "g_active_candidate * exp(tau_log_centered_ext)",
        "current_sigma_source_total_log_per_side_readback": float(value_law["sigma_source_total_log_per_side_readback"]),
        "diagnostic_eta_source_support_extension_log_per_side_candidate": minimal_extension.get("eta_source_support_extension_log_per_side_candidate"),
        "diagnostic_kappa_ext_candidate": minimal_extension.get("kappa_ext_candidate"),
        "diagnostic_candidate_status": minimal_extension.get("candidate_support_extension_status"),
        "smallest_constructive_missing_object": "eta_source_support_extension_log_per_side",
        "next_single_residual_object": "eta_source_support_extension_log_per_side",
        "next_residual_after_debug_eta_promotion": "sigma_source_support_extension_total_log_per_side",
        "notes": [
            "The same-support ordered package is exhausted, so exact charged closure requires the full two-scalar support-extension family rather than more same-support tuning.",
            "The live residual order on that family is first eta_source_support_extension_log_per_side and then sigma_source_support_extension_total_log_per_side.",
            "The rigid eta candidate remains diagnostic-only and is not promoted into the emitted scalar slots.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged two-scalar support-extension completion-law shell.")
    parser.add_argument("--value-law", default=str(DEFAULT_VALUE_LAW))
    parser.add_argument("--obstruction", default=str(DEFAULT_OBSTRUCTION))
    parser.add_argument("--minimal-extension", default=str(DEFAULT_MINIMAL_EXTENSION))
    parser.add_argument("--scale-binding", default=str(DEFAULT_SCALE_BINDING))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    value_law = json.loads(Path(args.value_law).read_text(encoding="utf-8"))
    obstruction = json.loads(Path(args.obstruction).read_text(encoding="utf-8"))
    minimal_extension = json.loads(Path(args.minimal_extension).read_text(encoding="utf-8"))
    scale_binding_path = Path(args.scale_binding)
    scale_binding = json.loads(scale_binding_path.read_text(encoding="utf-8")) if scale_binding_path.exists() else None
    artifact = build_artifact(value_law, obstruction, minimal_extension, scale_binding)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
