#!/usr/bin/env python3
"""Emit the scientific contract for the Ward-projected Thomson endpoint.

The artifact records the established no-go results, the missing source-side
objects, and the exact promotion boundary. It deliberately carries no
project-control or scheduling state.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_contract_current.json"

def build_contract() -> dict[str, Any]:
    return {
        "artifact": "oph_ward_projected_thomson_endpoint_contract",
        "status": "source_spectral_reduction_closed_measure_payload_absent",
        "promotion_allowed": False,
        "promotion_boundary": {
            "promotion_allowed": False,
            "first_missing_scientific_object": (
                "source-emitted same-scheme Ward-projected R_Q(P)"
            ),
            "minimal_required_theorem": "WardProjectedHadronicSpectralEmission_Q",
            "reason": (
                "The source-spectral reduction is established, but the source-emitted "
                "Ward-projected hadronic measure, same-scheme endpoint remainder, and "
                "full interval transport certificate are absent."
            ),
        },
        "computed_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
        "blocking_artifacts": {
            "screening_invariant_no_go": "code/P_derivation/runtime/screening_invariant_no_go_current.json",
            "interval_certificate": "code/P_derivation/runtime/fine_structure_interval_certificate_current.json",
            "residual_map_contract": "code/P_derivation/runtime/r_q_residual_contract_current.json",
            "source_spectral_theorem": "code/P_derivation/runtime/source_spectral_theorem_current.json",
            "spectral_transport_validator": "code/P_derivation/thomson_spectral_transport.py",
        },
        "no_go_results": {
            "detuning_only_bypass": "closed_no_go",
            "screening_scalar_from_current_invariants": "closed_no_go",
            "reason": (
                "The outer detuning coordinate and the current D10 invariant packet do not determine "
                "the zero-momentum Ward-projected current-current spectral transport."
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
                "id": "source_spectral_reduction_theorem",
                "kind": "theorem_artifact",
                "target_status": "reduction_closed_measure_payload_absent",
                "local_artifact": "code/P_derivation/runtime/source_spectral_theorem_current.json",
                "required_fields": [
                    "Ward_projected_transverse_correlator",
                    "subtracted_dispersion_kernel",
                    "source_payload_validation",
                    "nonidentifiability_corollary",
                ],
            },
            {
                "id": "screening_invariant_no_go",
                "kind": "non_identifiability_certificate",
                "target_status": "current_source_invariant_surface_does_not_determine_c_Q",
                "local_artifact": "code/P_derivation/runtime/screening_invariant_no_go_current.json",
                "required_fields": [
                    "source_packet_hash",
                    "compare_only_targets",
                    "lambda_family_witness",
                    "candidate_failures",
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
                "local_artifact": "code/P_derivation/runtime/fine_structure_interval_certificate_current.json",
                "required_fields": [
                    "alpha_interval",
                    "G_interval_image",
                    "derivative_bound",
                    "transport_error_budget",
                    "unique_fixed_point_statement",
                ],
            },
        ],
        "forbidden_solver_inputs": [
            "measured_alpha_0",
            "CODATA_or_NIST_Thomson_endpoint",
            "compare_alpha_inv",
            "P_C_residual",
            "S_required",
            "c_Q_target",
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
