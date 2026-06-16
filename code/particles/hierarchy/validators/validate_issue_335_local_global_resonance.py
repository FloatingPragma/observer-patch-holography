#!/usr/bin/env python3
"""Validate the OPH issue #335 local/global hierarchy-resonance close-out."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str) -> Decimal:
    return Decimal(value)


def main(path: str = "certificates/R_local_global_hierarchy_resonance_closeout_335.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    deps = cert.get("dependencies", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    obstruction = cert.get("obstruction_record", {})
    gates = cert.get("remaining_promotion_gates", [])
    forbidden = cert.get("forbidden_calibrations", [])
    checks = cert.get("checks", {})
    gate_checks = checks.get("remaining_promotion_gates_recorded", {})
    exact_capacity = cert.get("exact_capacity_certificate", {})
    readback = cert.get("finite_readback_resolution_certificate", {})
    round_count = cert.get("round_count_certificate", {})
    screen_sieve = cert.get("screen_sieve_certificate", {})
    chain = cert.get("derivation_chain", [])
    factor_origins = cert.get("factor_origins", {})
    branch_scope = cert.get("branch_scope", {})

    chain_steps = {step.get("step"): step for step in chain}

    validation = {
        "issue_is_335": cert.get("issue") == 335,
        "accepted_closeout": cert.get("accepted") is True,
        "status_is_full_resonance_closeout": (
            cert.get("status") == "closed_full_local_global_hierarchy_resonance"
        ),
        "full_theorem_promoted": cert.get("full_theorem_grade_resonance_promoted") is True,
        "all_dependencies_closed": all(deps.values()),
        "screen_sieve_dependency_present": deps.get(
            "screen_sieve_icosahedral_geometric_strengthening"
        )
        is True,
        "exact_capacity_supplied": (
            acceptance.get("exact_capacity_source_certificate_supplied") is True
            and D(exact_capacity.get("bridge_residual", "1")) == 0
        ),
        "finite_readback_supplied": (
            acceptance.get("finite_readback_resolution_supplied") is True
            and readback.get("status") == "closed_finite_readback_resolution_certificate"
        ),
        "round_count_supplied": (
            acceptance.get("round_count_derivation_supplied") is True
            and round_count.get("status") == "closed_representation_to_spectrum_round_count"
            and round_count.get("m_rep") == 24
        ),
        "screen_sieve_supplied": (
            acceptance.get("screen_sieve_geometric_strengthening_supplied") is True
            and screen_sieve.get("status") == "closed_on_declared_triangulated_screen_branch"
            and screen_sieve.get("orbit_size") == 12
            and screen_sieve.get("gamma_EW") == "(P/12)*log(N/pi)"
        ),
        "full_theorem_status_true": acceptance.get("full_theorem_grade_resonance_proved") is True,
        "exact_bridge_target_recorded": "3.5323546226929906511187512962330547600462" in (
            cert.get("exact_surviving_statement", {}).get("N_EW_public_endpoint", "")
        ),
        "rounded_capacity_rejected": obstruction.get("rounded_N_CRC_status")
        == "diagnostic_only_not_exact_bridge_certificate"
        and abs(D(obstruction.get("rounded_bridge_residual", "0"))) > Decimal("1e-6"),
        "no_promotion_gates_remain": len(gates) == 0
        and gate_checks.get("finite_readback_resolution") is True
        and gate_checks.get("round_count_derivation") is True
        and gate_checks.get("exact_capacity_source_certificate") is True
        and gate_checks.get("screen_sieve_geometric_strengthening") is True,
        "weak_scale_forbidden": any("measured weak scale" in item for item in forbidden),
        "rounded_display_forbidden": any("3.31e122" in item for item in forbidden),
        "derivation_chain_has_nine_steps": len(chain) == 9,
        "step_1_cites_d10_transmutation": "D10 forward transmutation theorem"
        in str(chain_steps.get(1, {}).get("premise", "")),
        "step_5_is_screen_sieve_geometric_strengthening": (
            "icosahedral screen-sieve theorem" in str(chain_steps.get(5, {}).get("premise", ""))
            and "geometric_strengthening_note" in chain_steps.get(5, {})
            and "EW tick-projection certificate"
            in str(chain_steps.get(5, {}).get("geometric_strengthening_note", ""))
            and "(P/12)" in str(chain_steps.get(5, {}).get("conclusion", ""))
        ),
        "step_6_is_ew_projection_bridge": "electroweak tick-projection bridge"
        in str(chain_steps.get(6, {}).get("premise", "")),
        "step_7_is_exact_capacity_contraction": "EW-refined exact-capacity"
        in str(chain_steps.get(7, {}).get("premise", "")),
        "step_8_composes_target_relation": "(P_*/12)*log(N_CRC^EW/pi)"
        in str(chain_steps.get(8, {}).get("conclusion", "")),
        "step_9_is_rg_higgs_compatibility": "RG/Higgs naturality"
        in str(chain_steps.get(9, {}).get("premise", "")),
        "factor_origin_beta_EW_recorded": factor_origins.get("beta_EW", {}).get("value") == "4"
        and "D10" in str(factor_origins.get("beta_EW", {}).get("source_theorem", "")),
        "factor_origin_m_rep_recorded": factor_origins.get("m_rep", {}).get("value") == "24"
        and "representation-to-spectrum"
        in str(factor_origins.get("m_rep", {}).get("source_theorem", "")),
        "factor_origin_icosahedral_orbit_recorded": factor_origins.get(
            "icosahedral_orbit_size_12", {}
        ).get("value")
        == "12"
        and "60 / 5"
        in str(factor_origins.get("icosahedral_orbit_size_12", {}).get("definition", "")),
        "factor_origin_total_charge_recorded": factor_origins.get(
            "total_curvature_charge_12", {}
        ).get("value")
        == "12"
        and "Gauss-Bonnet"
        in str(factor_origins.get("total_curvature_charge_12", {}).get("definition", "")),
        "factor_origin_cell_entropy_scoped": factor_origins.get(
            "cell_entropy_factor_one_over_four", {}
        ).get("value")
        == "1/4"
        and "scope_note"
        in factor_origins.get("cell_entropy_factor_one_over_four", {}),
        "factor_origin_12_in_P_over_12_recorded": factor_origins.get(
            "projection_target_denominator_12_in_P_over_12", {}
        ).get("value")
        == "12"
        and "icosahedral"
        in str(
            factor_origins.get("projection_target_denominator_12_in_P_over_12", {}).get(
                "source_theorem", ""
            )
        ),
        "branch_scope_records_screen_branch": "triangulated S^2"
        in str(branch_scope.get("screen_branch", "")),
        "branch_scope_records_product_gauge_branch": "product adjoint"
        in str(branch_scope.get("oph_product_gauge_branch", "")),
        "branch_scope_includes_scope_note": "cell-entropy"
        in str(branch_scope.get("scope_note", "")),
        "residual_residue_scoped_in_acceptance": acceptance.get(
            "residual_definitional_residue_scoped_as_oph_identification"
        )
        is True
        and "P/beta_EW"
        in str(acceptance.get("residual_definitional_residue_scope_note", "")),
    }
    payload = {"checks": validation, "pass": all(validation.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(
            sys.argv[1]
            if len(sys.argv) > 1
            else "certificates/R_local_global_hierarchy_resonance_closeout_335.json"
        )
    )
