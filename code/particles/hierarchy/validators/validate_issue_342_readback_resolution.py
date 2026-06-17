#!/usr/bin/env python3
"""Validate the OPH issue #342 finite readback-resolution certificate."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str | None, default: str = "0") -> Decimal:
    return Decimal(value if value is not None else default)


def main(path: str = "certificates/R_readback_resolution_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    reg = cert.get("capacity_register", {})
    atoms = [atom for atom in reg.get("atoms", []) if atom.get("selected") is True]
    obs = cert.get("observer_sector", {})
    nf = cert.get("normal_form", {})
    extractor = cert.get("extractor", {})
    readback = cert.get("readback_resolution", {})
    refinement = cert.get("refinement_certificate", {})
    ledger = cert.get("input_ledger", {})
    boundary = cert.get("claim_boundary", {})
    checks = cert.get("checks", {})
    derivation_chain = cert.get("derivation_chain", [])
    factor_origins = cert.get("factor_origins", {})
    branch_scope = cert.get("branch_scope", {})
    obstruction = cert.get("obstruction_record", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    dependencies = cert.get("dependencies", {})
    dependency_artifacts = cert.get("dependency_artifacts", {})
    strict = cert.get("strict_source_capacity_check", {})
    used = set(ledger.get("used_inputs", []))
    forbidden = set(ledger.get("forbidden_inputs", []))

    selected_cap = D(atoms[0].get("capacity")) if len(atoms) == 1 else Decimal(0)
    cap_read_dec = D(readback.get("cap_read"))
    rho_read_dec = D(readback.get("rho_read"))
    rho_star_dec = D(strict.get("rho_star"))
    n_crc_dec = D(strict.get("N_CRC"))
    capacity_residual = D(strict.get("capacity_relative_residual"))
    rho_residual = D(strict.get("rho_relative_residual"))
    tol = D(strict.get("relative_tolerance"))
    rounded_gap = D(obstruction.get("rounded_capacity_relative_gap_vs_N_CRC_EW"))

    chain_step_acceptance_criteria = {
        item.get("acceptance_criterion_closed")
        for item in derivation_chain
        if item.get("acceptance_criterion_closed")
    }

    validation = {
        "issue_is_342": cert.get("issue") == 342,
        "accepted": cert.get("accepted") is True,
        "status_closed_readback_resolution": (
            cert.get("status") == "closed_finite_readback_resolution_certificate"
        ),
        "normal_form_unique": nf.get("unique") is True
        and nf.get("finite_state") is True
        and nf.get("schedule_independent") is True,
        "normal_form_cites_confluence_theorem_artifact": "OBSERVERS_APPENDICES" in str(
            nf.get("source_artifact", "")
        ),
        "observer_sector_central_nonzero": obs.get("nonzero") is True
        and obs.get("central_projector") is True
        and obs.get("stable_self_reading") is True,
        "observer_sector_cites_synthesis_artifact": "OBSERVERS_SYNTHESIS_SECTIONS" in str(
            obs.get("source_artifact", "")
        ),
        "register_is_central_positive": reg.get("central") is True
        and reg.get("positive_spectrum") is True,
        "one_selected_atom": len(atoms) == 1,
        "selected_atom_is_positive": selected_cap > 0,
        "selected_probability_one": len(atoms) == 1 and D(atoms[0].get("probability")) == 1,
        "selected_variance_zero": D(reg.get("selected_variance")) == 0,
        "positive_root_extractor": extractor.get("formula") == "rho_read = sqrt(pi / Cap_read)"
        and extractor.get("positive_root") is True,
        "rho_read_positive": rho_read_dec > 0,
        "display_marked_diagnostic": readback.get("display_only") is True,
        "refinement_bound_recorded": refinement.get("uniform_limit_declared") is True
        and refinement.get("positive_root_closure") is True
        and Decimal("0") <= D(refinement.get("contraction_kappa")) < Decimal("1"),
        "kappa_matches_ew_lambda": refinement.get("kappa_matches_ew_lambda") is True,
        "no_forbidden_inputs_used": not (used & forbidden),
        "rounded_capacity_in_forbidden_inputs": (
            "rounded_3p31e122_capacity_display_as_bridge_witness" in forbidden
        ),
        "no_remaining_boundary": boundary.get("not_closed_here", []) == [],
        "claim_boundary_has_scope": "scope" in boundary
        and "N_CRC^EW" in str(boundary.get("scope", "")),
        "round_count_recorded_elsewhere": any(
            "R_m_rep_24_certificate" in item
            for item in boundary.get("closed_elsewhere", [])
        ),
        "exact_capacity_recorded_elsewhere": any(
            "R_EW_global_capacity_certificate" in item
            for item in boundary.get("closed_elsewhere", [])
        ),
        "verifier_checks_pass": all(checks.values()),
        "derivation_chain_has_eight_steps": len(derivation_chain) == 8,
        "step_2_is_fixed_cutoff_confluence": (
            len(derivation_chain) >= 2
            and "Fixed-cutoff confluence" in str(derivation_chain[1].get("premise", ""))
        ),
        "step_4_closes_single_resolution_criterion": (
            "single effective readback resolution"
            in chain_step_acceptance_criteria
        ),
        "step_7_loads_ew_exact_capacity_certificate": (
            len(derivation_chain) >= 7
            and "R_EW_global_capacity_certificate.json"
            in str(derivation_chain[6].get("source_artifact", ""))
        ),
        "step_8_closes_positive_root_closure_criterion": (
            "positive-root fixed-point closure forces rho_read -> rho_star"
            in chain_step_acceptance_criteria
        ),
        "factor_origin_pi_recorded": "area_law_normalization_pi" in factor_origins
        and factor_origins.get("area_law_normalization_pi", {}).get("value") == "pi",
        "factor_origin_positive_root_one_half": factor_origins.get(
            "positive_root_exponent_one_half", {}
        ).get("value")
        == "1/2",
        "factor_origin_banach_lambda": factor_origins.get(
            "banach_contraction_lambda_one_half", {}
        ).get("value")
        == "1/2"
        and "R_EW_global_capacity_certificate"
        in str(
            factor_origins.get("banach_contraction_lambda_one_half", {}).get(
                "source_artifact", ""
            )
        ),
        "factor_origin_derivative_factor_two": factor_origins.get(
            "derivative_bound_factor_two", {}
        ).get("value")
        == "2",
        "branch_scope_records_d6_branch": "d6_area_law_branch" in branch_scope
        and "area law" in str(branch_scope.get("d6_area_law_branch", "")).lower(),
        "branch_scope_records_ew_branch": "ew_refined_exact_capacity_branch"
        in branch_scope
        and "N_CRC^EW" in str(branch_scope.get("ew_refined_exact_capacity_branch", "")),
        "branch_scope_records_finite_repair_branch": "finite_oph_repair_branch"
        in branch_scope,
        "branch_scope_records_observer_branch": "stable_observer_sector_branch"
        in branch_scope,
        "branch_scope_includes_scope_note": "N_CRC^EW" in str(
            branch_scope.get("scope_note", "")
        ),
        "obstruction_records_rounded_diagnostic": obstruction.get(
            "rounded_N_CRC_status"
        )
        == "diagnostic_only_not_exact_bridge_witness"
        and rounded_gap > tol,
        "acceptance_definitions_emitted": acceptance.get("definitions_emitted") is True,
        "acceptance_single_resolution_proved": acceptance.get(
            "single_effective_readback_resolution_proved_via_central_record_algebra"
        )
        is True,
        "acceptance_finite_to_refinement_proved": acceptance.get(
            "finite_to_refinement_statement_proved_via_banach_and_positive_root_closure"
        )
        is True,
        "acceptance_positive_root_closure_proved": acceptance.get(
            "positive_root_fixed_point_closure_forces_rho_read_to_rho_star"
        )
        is True,
        "acceptance_inputs_and_forbidden_calibrations": acceptance.get(
            "allowed_inputs_and_forbidden_calibrations_stated"
        )
        is True
        and acceptance.get("no_measured_weak_higgs_or_hierarchy_calibration_used")
        is True,
        "acceptance_rounded_display_rejected": acceptance.get(
            "rounded_capacity_display_rejected_as_bridge_witness"
        )
        is True,
        "acceptance_certificate_emitted": acceptance.get(
            "public_certificate_and_verifier_emitted"
        )
        is True,
        "acceptance_ew_dependency_loaded": acceptance.get(
            "ew_exact_capacity_dependency_loaded_and_accepted"
        )
        is True,
        "acceptance_exact_bridge_hypothesis_supplied": acceptance.get(
            "exact_bridge_fixed_point_hypothesis_supplied_for_repair_tick_lemma"
        )
        is True,
        "ew_dependency_recorded": dependencies.get("ew_refined_exact_capacity") is True
        and "R_EW_global_capacity_certificate.json"
        in str(dependency_artifacts.get("ew_refined_exact_capacity", "")),
        "selected_atom_equals_n_crc_ew": selected_cap > 0
        and selected_cap == n_crc_dec,
        "cap_read_equals_n_crc": cap_read_dec == n_crc_dec,
        "rho_read_equals_rho_star": rho_read_dec == rho_star_dec,
        "strict_residuals_at_or_below_tolerance": capacity_residual <= tol
        and rho_residual <= tol,
    }
    payload = {"checks": validation, "pass": all(validation.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(
            sys.argv[1]
            if len(sys.argv) > 1
            else "certificates/R_readback_resolution_certificate.json"
        )
    )
