#!/usr/bin/env python3
"""Emit the constructive contract for the Ward-projected Thomson endpoint.

This artifact is the implementation target beneath the open low-energy
electromagnetic endpoint.  It does not promote the current continuation model.
It names the exact source-side objects that must be produced next so workers and
local code have constructive deliverables rather than obstruction-only packets.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_contract_current.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_contract() -> dict[str, Any]:
    return {
        "artifact": "oph_ward_projected_thomson_endpoint_contract",
        "generated_utc": _now_utc(),
        "github_issue": 235,
        "closed_blocker_isolation_issue": 223,
        "status": "constructive_contract_emitted_not_endpoint_theorem",
        "promotion_allowed": False,
        "computed_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
        "worker_result_policy": {
            "obstruction_only_result_allowed": False,
            "required_if_primary_route_fails": (
                "emit an implementation-ready replacement object, builder interface, schema, "
                "certificate interface, or runnable local target"
            ),
        },
        "endpoint_formula": {
            "source_anchor": "a0(P) = alpha_em^-1(m_Z^2;P)",
            "target": "alpha_Th^-1(P) = a0(P) + Delta_Th(P)",
            "split": [
                "Delta_lep_src(P)",
                "Delta_had_src(P)",
                "Delta_EW_src(P)",
            ],
        },
        "constructive_objects": [
            {
                "id": "sigma_q_scheme_lock",
                "kind": "certificate_interface",
                "target_status": "source_scheme_shared_by_anchor_and_endpoint",
                "required_fields": [
                    "charge_operator_Q",
                    "ward_projection_condition",
                    "normalization_convention",
                    "scheme_identifier",
                ],
            },
            {
                "id": "delta_lep_source_transport",
                "kind": "builder_interface",
                "target_status": "source_emitted_charged_transport_with_quadrature_bound",
                "required_fields": [
                    "charged_mass_source",
                    "one_loop_kernel",
                    "quadrature_rule",
                    "quadrature_error_bound",
                    "scheme_match_to_a0",
                ],
            },
            {
                "id": "rho_had_spectral_measure",
                "kind": "production_export_contract",
                "target_status": "ward_projected_hadronic_spectral_measure_export",
                "local_contract": "code/particles/hadron/ward_projected_spectral_measure.schema.json",
                "required_fields": [
                    "finite_volume_levels",
                    "ward_projected_residues",
                    "current_normalization",
                    "pushforward_to_rho_had",
                    "continuum_volume_chiral_statistical_matching_budgets",
                ],
            },
            {
                "id": "delta_ew_remainder",
                "kind": "certificate_interface",
                "target_status": "zero_identity_or_source_bound",
                "required_fields": [
                    "scheme_remainder_formula",
                    "zero_identity_status",
                    "source_bound",
                    "bound_interval",
                ],
            },
            {
                "id": "delta_qcd_screening_and_endpoint_remainder",
                "kind": "computed_residual_interface",
                "target_status": "source_emitted_residual_map_for_the_package_scalar",
                "local_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
                "required_fields": [
                    "ward_projected_qcd_screening_map",
                    "electroweak_scheme_remainder_map",
                    "positivity_or_bound_status",
                    "interval_transport_bound",
                ],
            },
            {
                "id": "full_endpoint_interval_certificate",
                "kind": "certificate_interface",
                "target_status": "interval_existence_uniqueness_for_full_map",
                "required_fields": [
                    "alpha_interval",
                    "G_interval_image",
                    "derivative_bound",
                    "transport_error_budget",
                    "unique_fixed_point_statement",
                ],
            },
        ],
        "local_next_steps": [
            "Use thomson_endpoint_package.py to compute the residual endpoint packet from the internal report.",
            "Populate a source-emitted Ward-projected QCD screening and endpoint-remainder map matching that packet.",
            "Upgrade the interval certificate only after the endpoint builder uses source-emitted transport.",
        ],
        "forbidden_solver_inputs": [
            "measured_alpha_0",
            "CODATA_or_NIST_Thomson_endpoint",
            "free_quark_screened_ansatz_as_theorem_input",
            "charged_physical_values_imported_into_the_determinant_line",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the Thomson endpoint constructive contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = build_contract()
    text = json.dumps(contract, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
