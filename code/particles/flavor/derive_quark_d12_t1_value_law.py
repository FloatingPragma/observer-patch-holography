#!/usr/bin/env python3
"""Emit the conditional ray-coordinate D12 quark identity package.

Chain role: record the smallest target-free theorem object left on the emitted
same-family D12 mass ray.

Mathematics: the emitted ray already identifies the unresolved mass-side
coordinate with one scalar, `ray_modulus = t1`, and the full mass-side wrapper
is then induced by
`Delta_ud_overlap = t1 / 5` and
`eta_Q_centered = -((1 - x2^2) / 27) * t1`.

The current source corpus does not emit ``c_d/c_u`` or a numerical ``t1``.
The output therefore records the exact conditional identities while keeping
the scalar-emission route open.  Its retained numerical point is diagnostic
only and cannot promote the package.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
MASS_RAY_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"
OVERLAP_LAW_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_t1_value_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(mass_ray: dict[str, Any], overlap_law: dict[str, Any] | None = None) -> dict[str, Any]:
    sample = dict(mass_ray["sample_same_family_point"])
    formulas = dict(mass_ray["same_family_ray"]["ray_formulas"])
    overlap_law = overlap_law or {}
    exact_transport = dict(overlap_law.get("exact_transport_contract", {}))
    main_builder_branch = dict((overlap_law.get("sigma_branch_contracts") or {}).get("main_builder_sigma_pair") or {})
    main_builder_sigma_status = main_builder_branch.get("provider_status")
    main_builder_sigma_source_emitted = main_builder_sigma_status in {
        "closed_source_emission",
        "source_emitted",
        "source_values_derived_from_source_emission",
    }
    return {
        "artifact": "oph_quark_d12_t1_value_law",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "conditional_mass_ray_identity_value_source_open",
        "public_promotion_allowed": False,
        "exact_missing_object": "source_derived_c_d_over_c_u_or_t1_value",
        "value_emission_status": "open_viable_scalar_emission_route",
        "parent_emitted_object": {
            "artifact": mass_ray["artifact"],
            "id": mass_ray["emitted_object"]["id"],
        },
        "theorem_statement": (
            "Conditional on a source-derived positive coefficient ratio c_d/c_u, the minimal light branch "
            "and emitted same-family D12 mass ray imply t1 = (5/6) * log(c_d / c_u). Equivalently, "
            "t1 = 5 * Delta_ud_overlap and ray_modulus = t1 on D12_ud_mass_ray. "
            "These identities do not emit c_d/c_u or select a numerical point of the ray."
        ),
        "smaller_primitive_equivalent_theorem": {
            "id": "light_quark_overlap_defect_value_law",
            "selector_scalar_name": "Delta_ud_overlap",
            "equivalent_ray_coordinate_identity": "t1 = 5 * Delta_ud_overlap",
        },
        "one_scalar_reduction": {
            "scalar_name": "t1",
            "ray_coordinate_identity": "ray_modulus = t1",
            "induced_formulas": {
                "Delta_ud_overlap": "t1 / 5",
                "eta_Q_centered": "-((1 - x2^2) / 27) * t1",
                "kappa_Q": "-t1 / 54",
            },
            "emitted_ray_formulas": formulas,
        },
        "forced_source_payload_after_t1": {
            "selector_scalar_name": "Delta_ud_overlap",
            "Delta_ud_overlap": "t1 / 5",
            "beta_u_diag_B_source": "t1 / 10",
            "beta_d_diag_B_source": "-t1 / 10",
            "source_readback_u_log_per_side": "(-t1/10, 0, +t1/10)",
            "source_readback_d_log_per_side": "(+t1/10, 0, -t1/10)",
            "J_B_source_u": "t1 / 10",
            "J_B_source_d": "-t1 / 10",
        },
        "transport_side_reduction": {
            "status": "closed_given_fixed_sigma_branch",
            "law_statement": exact_transport.get(
                "law_statement",
                "Once a positive sigma branch is fixed, the D12 odd transport side is an exact affine readback of Delta_ud_overlap and therefore of t1.",
            ),
            "weighted_transport_formulas_after_t1": {
                "tau_u_log_per_side": "sigma_d_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
                "tau_d_log_per_side": "sigma_u_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
                "Lambda_ud_B_transport": "sigma_u_total_log_per_side * sigma_d_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            },
            "closed_identities_after_t1": {
                "tau_u_plus_tau_d": "t1 / 10",
                "tau_d_over_tau_u": "sigma_u_total_log_per_side / sigma_d_total_log_per_side",
            },
            "exact_gap_after_transport_reduction": "supported_sigma_branch_emission",
            "reduced_gap_note": (
                "The transport side no longer contributes a separate exact theorem burden once a sigma branch is fixed. "
                "The remaining nontransport question is the value law itself together with promotion of any candidate sigma branch used downstream."
            ),
        },
        "candidate_public_construction_route": {
            "status": "conditional_identity_value_source_open",
            "minimal_light_branch": {
                "textures": [
                    "y_u = c_u * epsilon^6",
                    "y_d = c_d * epsilon^6",
                ],
                "same_forward_surface_ratio_identification": "c_u / c_d = y_u / y_d",
                "light_isospin_breaking_scalar": "Delta_ud_overlap = (1/6) * log(c_d / c_u)",
            },
            "same_family_d12_hypercharge_law": {
                "C1": "C_1(Y) = (3/5) * Y^2",
                "difference": "C_1(u^c) - C_1(d^c) = 1/5",
                "ray_identity": "Delta_ud_overlap = t1 / 5",
            },
            "candidate_scalar_identities": [
                "log(c_d / c_u) = (6/5) * t1",
                "t1 = 5 * Delta_ud_overlap = (5/6) * log(c_d / c_u)",
            ],
            "closure_note": (
                "The algebraic reduction is internalized, but the source law selecting c_d/c_u or t1 remains open."
            ),
        },
        "current_public_yukawa_frontier": {
            "main_builder_sigma_branch_artifact": main_builder_branch.get(
                "provider_artifact"
            ),
            "main_builder_sigma_branch_status": main_builder_sigma_status,
            "main_builder_sigma_branch_kind": main_builder_branch.get(
                "sigma_source_kind"
            ),
            "main_builder_sigma_source_emitted": main_builder_sigma_source_emitted,
            "transport_side_status_given_that_branch": "closed_conditionally",
            "remaining_exact_blocker_set_if_this_branch_is_used": [
                "source_derived_c_d_over_c_u_or_t1_value",
                *(
                    []
                    if main_builder_sigma_source_emitted
                    else ["supported_source_sigma_branch_emission"]
                ),
            ],
            "target_1_status": "open_value_source_not_emitted",
            "equivalent_ray_coordinate_presentation": "quark_d12_t1_value_law",
            "optional_nonminimal_strengthening": "edge_statistics_mean_surface_compatibility_theorem",
        },
        "sample_same_family_point": {
            "classification": "diagnostic_sample_only_not_value_emission",
            "t1": float(sample["t1_sample"]),
            "ray_modulus": float(sample["ray_modulus"]),
            "x2": float(sample["x2"]),
            "Delta_ud_overlap": float(sample["Delta_ud_overlap"]),
            "eta_Q_centered": float(sample["eta_Q_centered"]),
            "kappa_Q": float(sample["kappa_Q"]),
        },
        "minimal_new_theorem": {
            "id": "quark_d12_t1_value_law",
            "must_emit": "quark_d12_t1_value_law",
            "scalar_name": "t1",
            "unique_intersection_with": "D12_ud_mass_ray",
            "identifies": "ray_modulus = t1",
            "then_emits": ["t1", "ray_modulus", "Delta_ud_overlap", "eta_Q_centered"],
            "must_not_use_target_masses": True,
            "must_not_use_ckm_cp": True,
        },
        "derived_wrapper": {
            "id": "intrinsic_scale_law_D12",
            "derived_from": "quark_d12_t1_value_law",
            "meaning": "the induced mass-side wrapper on D12_ud_mass_ray after t1 is fixed intrinsically",
        },
        "notes": [
            "This artifact proves a conditional coordinate identity; it does not emit the unresolved scalar.",
            "The same-family D12 mass ray remains a viable open derivation route, not a no-go.",
            "The exact transport identity does not upgrade the diagnostic sigma witness consumed by the current main builder.",
            "The larger wrapper intrinsic_scale_law_D12 is only a conditional derived wrapper above quark_d12_t1_value_law.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 quark t1 value-law scaffold.")
    parser.add_argument("--mass-ray", default=str(MASS_RAY_JSON))
    parser.add_argument("--overlap-law", default=str(OVERLAP_LAW_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    overlap_path = Path(args.overlap_law)
    overlap_law = _load_json(overlap_path) if overlap_path.exists() else None
    payload = build_payload(_load_json(Path(args.mass_ray)), overlap_law=overlap_law)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
