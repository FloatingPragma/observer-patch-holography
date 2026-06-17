#!/usr/bin/env python3
"""Validate the OPH issue #343 m_rep=24 certificate."""

from __future__ import annotations

import json
import pathlib
import sys


REQUIRED_FACTOR_ORIGIN_KEYS = {
    "dim_su3_adjoint",
    "dim_su2_adjoint",
    "dim_u1",
    "unoriented_total_twelve_curvature_ports",
    "orientation_multiplier",
    "m_rep",
    "exponent_denominator",
}

REQUIRED_BRANCH_SCOPE_KEYS = {
    "oph_realized_compact_gauge_branch",
    "reversible_repair_orientation_branch",
    "cyclic_scheduler_branch",
    "scope_note",
}

REQUIRED_ACCEPTANCE_KEYS = {
    "ac1_define_grammar_sector_spectral_object_observable",
    "ac2_proof_m_rep_24_on_source_side_oph_data",
    "ac3_specialization_to_minus_one_over_48",
    "ac4_negative_controls_or_obstruction_notes",
    "ac5_no_measured_inputs",
    "ac6_emit_public_certificate_and_verifier",
    "ac7_update_only_if_status_changes",
    "all_acceptance_criteria_satisfied",
}

REQUIRED_DEPENDENCY_KEYS = {
    "global_repair_tick_lemma",
    "corpus_realized_product_branch",
    "corpus_orientation_doubling",
    "corpus_twelve_curvature_ports",
}

REQUIRED_CONSUMER_KEYS = {
    "global_repair_tick_lemma",
    "ew_tick_projection_bridge",
    "local_global_resonance_closeout",
}


def _factor_value(cert: dict, key: str) -> object:
    return cert.get("factor_origins", {}).get(key, {}).get("value")


def main(path: str = "certificates/R_m_rep_24_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    rep = cert.get("representation_sector", {})
    result = cert.get("result", {})
    spectral = cert.get("spectral_object", {})
    controls = cert.get("negative_controls", [])
    checks = cert.get("verifier_checks", {})
    forbidden_used = cert.get("forbidden_inputs_used", [])
    dims = [component.get("dimension") for component in rep.get("components", [])]
    derivation_chain = cert.get("derivation_chain", [])
    factor_origins = cert.get("factor_origins", {})
    branch_scope = cert.get("branch_scope", {})
    dependency_artifacts = cert.get("dependency_artifacts", {})
    consumer_artifacts = cert.get("consumer_artifacts", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    claim_boundary = cert.get("claim_boundary", {})

    derivation_premises = [step.get("premise", "") for step in derivation_chain]

    validation = {
        "issue_is_343": cert.get("issue") == 343,
        "accepted": cert.get("accepted") is True,
        "status_closed_round_count": (
            cert.get("status") == "closed_representation_to_spectrum_round_count"
        ),
        "target_relation_recorded": (
            "G_N = g_N^m_rep" in cert.get("target_relation", "")
            and "(N_CRC/pi)^(-1/48)" in cert.get("target_relation", "")
        ),
        "component_dimensions_are_8_3_1": dims == [8, 3, 1],
        "unoriented_dimension_is_12": rep.get("unoriented_adjoint_dimension") == 12,
        "orientation_multiplier_is_2": rep.get("orientation_multiplier") == 2,
        "oriented_support_is_24": rep.get("oriented_support_dimension") == 24,
        "m_rep_is_24": result.get("m_rep") == 24,
        "exponent_is_minus_one_over_48": result.get("specialized_exponent") == "-1/48",
        "full_cycle_decomposition_recorded": (
            result.get("full_cycle_decomposition") == "G_N = g_N^m_rep = g_N^24"
        ),
        "spectral_period_is_24": spectral.get("period") == 24,
        "su5_negative_control_recorded": any(
            item.get("name") == "single-orientation SU(5) adjoint"
            and item.get("status") == "reject_despite_same_integer"
            for item in controls
        ),
        "doubled_su5_rejected": any(
            item.get("name") == "doubled SU(5) adjoint" and item.get("status") == "reject"
            for item in controls
        ),
        "graviton_excluded": any(
            item.get("name") == "include graviton in compact-gauge repair support"
            and item.get("status") == "reject"
            for item in controls
        ),
        "every_negative_control_records_violated_branch": all(
            isinstance(item.get("violated_branch"), str)
            and item.get("violated_branch")
            for item in controls
        ),
        "negative_controls_at_least_seven": len(controls) >= 7,
        "no_forbidden_inputs_used": forbidden_used == [],
        "verifier_checks_pass": all(checks.values()),
        "claim_boundary_closed": claim_boundary.get("not_closed_here") == [],
        "claim_boundary_records_scope": (
            isinstance(claim_boundary.get("scope"), str)
            and "(N_CRC/pi)^(-1/48)" in claim_boundary.get("scope", "")
        ),
        "derivation_chain_has_eight_steps": len(derivation_chain) == 8,
        "derivation_step1_realized_product_branch": (
            len(derivation_premises) >= 1
            and "realized observer-visible product-gauge branch"
            in derivation_premises[0]
        ),
        "derivation_step2_lie_algebra_dimensions": (
            len(derivation_premises) >= 2
            and "Compact Lie algebra adjoint dimensions" in derivation_premises[1]
        ),
        "derivation_step3_twelve_curvature_ports": (
            len(derivation_premises) >= 3
            and "twelve-curvature-port" in derivation_premises[2]
        ),
        "derivation_step4_orientation_doubling": (
            len(derivation_premises) >= 4
            and "orientation-doubling axiom" in derivation_premises[3]
        ),
        "derivation_step5_oriented_support": (
            len(derivation_premises) >= 5
            and "Oriented adjoint support dimension" in derivation_premises[4]
        ),
        "derivation_step6_spectral_object": (
            len(derivation_premises) >= 6
            and "Cyclic repair scheduler" in derivation_premises[5]
        ),
        "derivation_step7_specializes_tick_law": (
            len(derivation_premises) >= 7
            and "Specialization of the parametric global repair-tick law"
            in derivation_premises[6]
        ),
        "derivation_step8_negative_controls": (
            len(derivation_premises) >= 8
            and "Negative-control rejection" in derivation_premises[7]
        ),
        "every_derivation_step_has_source_artifact": all(
            isinstance(step.get("source_artifact"), str) and step.get("source_artifact")
            for step in derivation_chain
        ),
        "every_derivation_step_has_conclusion": all(
            isinstance(step.get("conclusion"), str) and step.get("conclusion")
            for step in derivation_chain
        ),
        "factor_origins_keys_complete": (
            set(factor_origins.keys()) == REQUIRED_FACTOR_ORIGIN_KEYS
        ),
        "factor_origin_dim_su3_is_8": _factor_value(cert, "dim_su3_adjoint") == 8,
        "factor_origin_dim_su2_is_3": _factor_value(cert, "dim_su2_adjoint") == 3,
        "factor_origin_dim_u1_is_1": _factor_value(cert, "dim_u1") == 1,
        "factor_origin_unoriented_is_12": (
            _factor_value(cert, "unoriented_total_twelve_curvature_ports") == 12
        ),
        "factor_origin_orientation_multiplier_is_2": (
            _factor_value(cert, "orientation_multiplier") == 2
        ),
        "factor_origin_m_rep_is_24": _factor_value(cert, "m_rep") == 24,
        "factor_origin_exponent_denominator_is_48": (
            _factor_value(cert, "exponent_denominator") == 48
        ),
        "factor_origin_orientation_cites_corpus": (
            "compact_proof_of_oph.tex"
            in factor_origins.get("orientation_multiplier", {}).get("source_artifact", "")
        ),
        "factor_origin_unoriented_cites_corpus": (
            "compact_proof_of_oph.tex"
            in factor_origins.get("unoriented_total_twelve_curvature_ports", {}).get(
                "source_artifact", ""
            )
        ),
        "factor_origin_exponent_cites_global_repair_tick_cert": (
            "R_N_global_repair_tick_certificate.json"
            in factor_origins.get("exponent_denominator", {}).get("source_artifact", "")
        ),
        "branch_scope_keys_complete": (
            set(branch_scope.keys()) == REQUIRED_BRANCH_SCOPE_KEYS
        ),
        "branch_scope_records_realized_product_branch": (
            "(SU(3) x SU(2) x U(1))/Z6"
            in branch_scope.get("oph_realized_compact_gauge_branch", "")
        ),
        "branch_scope_records_patch_carrier_pipeline": (
            "patch-carrier"
            in branch_scope.get("reversible_repair_orientation_branch", "")
        ),
        "branch_scope_records_cyclic_order_equals_m_rep": (
            "cyclic permutation"
            in branch_scope.get("cyclic_scheduler_branch", "")
            and "m_rep" in branch_scope.get("cyclic_scheduler_branch", "")
        ),
        "branch_scope_includes_scope_note": (
            isinstance(branch_scope.get("scope_note"), str)
            and "m_rep=24" in branch_scope.get("scope_note", "")
        ),
        "dependency_artifacts_keys_complete": (
            set(dependency_artifacts.keys()) == REQUIRED_DEPENDENCY_KEYS
        ),
        "dependency_global_repair_tick_present": (
            "R_N_global_repair_tick_certificate.json"
            in dependency_artifacts.get("global_repair_tick_lemma", "")
        ),
        "consumer_artifacts_keys_complete": (
            set(consumer_artifacts.keys()) == REQUIRED_CONSUMER_KEYS
        ),
        "consumer_local_global_resonance_present": (
            "R_local_global_hierarchy_resonance_closeout_335.json"
            in consumer_artifacts.get("local_global_resonance_closeout", "")
        ),
        "acceptance_criteria_keys_complete": (
            set(acceptance.keys()) == REQUIRED_ACCEPTANCE_KEYS
        ),
        "acceptance_criteria_all_satisfied": (
            acceptance.get("all_acceptance_criteria_satisfied") is True
        ),
        "ac1_string_records_grammar_sector_spectral_observable": (
            isinstance(
                acceptance.get("ac1_define_grammar_sector_spectral_object_observable"),
                str,
            )
            and "repair_grammar"
            in acceptance.get("ac1_define_grammar_sector_spectral_object_observable", "")
        ),
        "ac2_string_records_m_rep_derivation": (
            "2 * (8 + 3 + 1)"
            in acceptance.get("ac2_proof_m_rep_24_on_source_side_oph_data", "")
        ),
        "ac3_string_records_minus_one_over_48": (
            "(N_CRC/pi)^(-1/48)"
            in acceptance.get("ac3_specialization_to_minus_one_over_48", "")
        ),
        "ac5_string_records_no_measured_inputs": (
            "FORBIDDEN_INPUTS"
            in acceptance.get("ac5_no_measured_inputs", "")
        ),
        "used_inputs_cite_compact_proof_corpus": any(
            "compact_proof_of_oph.tex" in entry
            for entry in cert.get("used_inputs", [])
        ),
        "used_inputs_cite_global_repair_tick_cert": any(
            "R_N_global_repair_tick_certificate.json" in entry
            for entry in cert.get("used_inputs", [])
        ),
        "certificate_id_v2": (
            cert.get("certificate_id") == "issue-343-m-rep-24-doubled-sm-adjoint-v2"
        ),
    }
    payload = {"checks": validation, "pass": all(validation.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(
            sys.argv[1]
            if len(sys.argv) > 1
            else "certificates/R_m_rep_24_certificate.json"
        )
    )
