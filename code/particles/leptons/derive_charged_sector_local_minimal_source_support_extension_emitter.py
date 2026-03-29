#!/usr/bin/env python3
"""Parameterize the first beyond-support charged extension without fitting it.

Chain role: expose the smallest extra charged scalar that can leave the current
linear support and change the `e/mu/tau` hierarchy.

Mathematics: quadratic support-extension basis, ordered-gap reconstruction, and
rigid ratio propagation from the active ordered family.

OPH-derived inputs: the charged ordered-package value law, the support
obstruction certificate, the family coordinate data, and the shared scale lane.

Output: the `eta_source_support_extension_log_per_side` emitter shell together
with candidate formulas and downstream readout consequences.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VALUE_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"
DEFAULT_OBSTRUCTION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_current_support_obstruction_certificate.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "flavor" / "family_excitation_evaluator.json"
DEFAULT_SCALE_BINDING = ROOT / "particles" / "runs" / "leptons" / "lepton_shared_absolute_scale_binding.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    value_law: dict,
    obstruction: dict | None = None,
    family: dict | None = None,
    scale_binding: dict | None = None,
) -> dict:
    obstruction = obstruction or {}
    family = family or {}
    scale_binding = scale_binding or {}
    x2 = float(value_law["ordered_family_coordinate"][1])
    sigma_source = float(value_law["sigma_source_total_log_per_side_readback"])
    linear_basis = [float(value) for value in value_law["linear_basis_vector_centered"]]
    extension_basis = [float(value) for value in value_law["quadratic_basis_vector"]]
    denom = 2.0 * (1.0 - x2 * x2)
    ordered_package = [float(value) for value in value_law["source_side_ordered_package_log_per_side_emitted"]]
    mu_source = float(value_law["mu_source_log_per_side_readback"])
    centered_logs = [value - mu_source for value in ordered_package]
    gamma21_current = centered_logs[1] - centered_logs[0]
    gamma32_current = centered_logs[2] - centered_logs[1]
    rho_ord_current = gamma21_current / gamma32_current
    rho_ord_family = family.get("rho_ord")
    rho_ord_extension = float(rho_ord_family) if rho_ord_family is not None else rho_ord_current
    eta_candidate = (((1.0 + x2) - rho_ord_extension * (1.0 - x2)) / (1.0 + rho_ord_extension)) * sigma_source
    kappa_candidate = eta_candidate / denom
    gamma21_ext = (((1.0 + x2) * sigma_source) - eta_candidate) / 2.0
    gamma32_ext = (((1.0 - x2) * sigma_source) + eta_candidate) / 2.0
    e_log_centered_ext = [
        -((2.0 * gamma21_ext) + gamma32_ext) / 3.0,
        (gamma21_ext - gamma32_ext) / 3.0,
        (gamma21_ext + (2.0 * gamma32_ext)) / 3.0,
    ]
    source_side_ordered_package_ext = [mu_source + value for value in e_log_centered_ext]
    shape_singular_values_ext = [math.exp(value) for value in e_log_centered_ext]
    singular_values_abs_ext = None
    if scale_binding.get("g_e") is not None:
        g_active = float(scale_binding["g_e"])
        singular_values_abs_ext = [g_active * value for value in shape_singular_values_ext]
    endpoint_ratio = math.exp(sigma_source)
    return {
        "artifact": "oph_charged_sector_local_minimal_source_support_extension_emitter",
        "generated_utc": _timestamp(),
        "proof_status": "minimal_support_extension_formula_closed_source_scalar_open",
        "predictive_promotion_allowed": False,
        "source_artifact": obstruction.get("artifact", value_law.get("artifact")),
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "linear_basis_vector_centered": linear_basis,
        "extension_basis_vector_centered": extension_basis,
        "extension_basis_kind": "canonical_quadratic_ordered_three_point_direction",
        "mu_source_log_per_side": value_law.get("mu_source_log_per_side_readback"),
        "sigma_source_total_log_per_side": sigma_source,
        "g_active_candidate": scale_binding.get("g_e"),
        "eta_source_support_extension_log_per_side": None,
        "kappa_ext": None,
        "rho_ord_current": rho_ord_current,
        "rho_ord_support_extension": rho_ord_extension,
        "eta_source_support_extension_log_per_side_candidate": eta_candidate,
        "eta_source_support_extension_log_per_side_candidate_formula": "(((1 + x2) - rho_ord_support_extension * (1 - x2)) / (1 + rho_ord_support_extension)) * sigma_source_total_log_per_side",
        "eta_source_support_extension_candidate_origin": (
            "rigid_ordered_ratio_support_extension_candidate"
            if rho_ord_family is not None
            else "current_package_ratio_support_extension_candidate"
        ),
        "kappa_ext_candidate": kappa_candidate,
        "gamma21_log_per_side_ext_candidate": gamma21_ext,
        "gamma32_log_per_side_ext_candidate": gamma32_ext,
        "E_e_log_centered_ext_candidate": e_log_centered_ext,
        "source_side_ordered_package_log_per_side_ext_candidate": source_side_ordered_package_ext,
        "shape_singular_values_ext_candidate": shape_singular_values_ext,
        "singular_values_abs_ext_candidate": singular_values_abs_ext,
        "candidate_support_extension_status": "eta_candidate_only_total_span_still_open",
        "candidate_next_single_residual_object": "sigma_source_support_extension_total_log_per_side",
        "fixed_current_span_certificate": {
            "sigma_source_total_log_per_side_fixed": sigma_source,
            "sigma_source_total_log_per_side_formula": "current emitted support span",
            "endpoint_ratio_formula": "exp(sigma_source_total_log_per_side)",
            "endpoint_ratio_fixed_value": endpoint_ratio,
            "endpoint_ratio_invariant_under_eta_only_extension": True,
            "why_eta_only_preserves_endpoint_ratio": "the eta-only extension adds +eta/6 to both endpoint centered logs and therefore leaves tau/e unchanged",
        },
        "kappa_ext_formula": "eta_source_support_extension_log_per_side / (2 * (1 - x2^2))",
        "centered_extension_formula": "E_e_log_centered_ext = (sigma_source_total_log_per_side / 2) * linear_basis_vector_centered + kappa_ext * extension_basis_vector_centered",
        "gamma21_log_per_side_ext_formula": "((1 + x2) * sigma_source_total_log_per_side - eta_source_support_extension_log_per_side) / 2",
        "gamma32_log_per_side_ext_formula": "((1 - x2) * sigma_source_total_log_per_side + eta_source_support_extension_log_per_side) / 2",
        "e_log_centered_ext_formula": "-((3 + x2) * sigma_source_total_log_per_side - eta_source_support_extension_log_per_side) / 6",
        "mu_log_centered_ext_formula": "(x2 * sigma_source_total_log_per_side - eta_source_support_extension_log_per_side) / 3",
        "tau_log_centered_ext_formula": "((3 - x2) * sigma_source_total_log_per_side + eta_source_support_extension_log_per_side) / 6",
        "e_ext_formula": "g_active_candidate * exp(e_log_centered_ext)",
        "mu_ext_formula": "g_active_candidate * exp(mu_log_centered_ext)",
        "tau_ext_formula": "g_active_candidate * exp(tau_log_centered_ext)",
        "hierarchy_mover_condition": "eta_source_support_extension_log_per_side != 0",
        "new_beyond_support_scalar_slots_required": 1,
        "support_extension_denominator": denom,
        "smallest_constructive_missing_object": "eta_source_support_extension_log_per_side",
        "next_single_residual_object": "eta_source_support_extension_log_per_side",
        "notes": [
            "The current-support charged package has collapsed to the linear subray.",
            "The first honest hierarchy mover is one beyond-support scalar on the canonical quadratic ordered direction.",
            "At fixed current span, eta-only moves the middle state against a fixed endpoint ratio tau/e = exp(sigma_source_total_log_per_side).",
            "Once eta_source_support_extension_log_per_side is emitted, every downstream charged gap and centered log is algebraically determined.",
            "A rigid ordered-ratio support-extension candidate is recorded for debugging because it uses only current-family data, but it is not promoted as closure because the remaining charged mismatch is then dominated by the still-open total spread scalar.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the minimal support-extension charged emitter artifact.")
    parser.add_argument("--value-law", default=str(DEFAULT_VALUE_LAW))
    parser.add_argument("--obstruction", default=str(DEFAULT_OBSTRUCTION))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--scale-binding", default=str(DEFAULT_SCALE_BINDING))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    value_law = json.loads(Path(args.value_law).read_text(encoding="utf-8"))
    obstruction_path = Path(args.obstruction)
    obstruction = json.loads(obstruction_path.read_text(encoding="utf-8")) if obstruction_path.exists() else None
    family_path = Path(args.family)
    family = json.loads(family_path.read_text(encoding="utf-8")) if family_path.exists() else None
    scale_binding_path = Path(args.scale_binding)
    scale_binding = json.loads(scale_binding_path.read_text(encoding="utf-8")) if scale_binding_path.exists() else None
    artifact = build_artifact(value_law, obstruction, family, scale_binding)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
