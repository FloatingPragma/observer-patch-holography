#!/usr/bin/env python3
"""Focused guards for the target-free quark-spread obstruction."""

from __future__ import annotations

import json
import math
from pathlib import Path

from derive_quark_sigma_source_nonidentifiability_obstruction import (
    DIRECT_INPUT_PATHS,
    FORBIDDEN_ANCESTORS,
    build_artifact,
)


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def _build() -> dict:
    return build_artifact(
        _load("family_transport_kernel.json"),
        _load("generation_bundle_branch_generator.json"),
        _load("generation_bundle_branch_generator_splitting_obstruction.json"),
        _load("overlap_edge_transport_cocycle.json"),
    )


def test_quark_sigma_obstruction_closes_all_three_backlog_issues_without_promotion() -> None:
    payload = _build()

    assert payload["artifact"] == "oph_quark_sigma_source_nonidentifiability_obstruction"
    assert payload["proof_status"] == "closed_exact_current_corpus_obstruction"
    assert payload["theorem_grade_obstruction"] is True
    assert payload["source_only_sigma_emitted"] is False
    assert payload["public_promotion_allowed"] is False
    assert payload["numeric_quark_rows_allowed"] is False
    assert payload["issue_377_acceptance_met_as_obstruction"] is True
    assert payload["issue_379_acceptance_met_as_obstruction"] is True
    assert payload["issue_380_acceptance_met_as_obstruction"] is True

    composition = payload["issue_composition"]
    assert composition["issue_379_up_type"]["remaining_fiber"] == "R_{>0}"
    assert composition["issue_379_up_type"]["direct_top_codomain_kept_separate"] is True
    assert composition["issue_380_down_type"]["remaining_fiber"] == "R_{>0}"
    assert (
        composition["issue_380_down_type"]["elementary_running_quark_vs_hadron_qcd_distinction_preserved"]
        is True
    )
    assert composition["issue_377_all_quark"]["remaining_fiber"] == "(R_{>0})^2"
    assert composition["issue_377_all_quark"]["composition"] == (
        "issue_379_up_type x issue_380_down_type"
    )


def test_exact_profile_classification_leaves_two_independent_positive_moduli() -> None:
    payload = _build()
    classification = payload["exact_ray_classification"]

    assert classification["fiber"] == "(R_{>0})^2"
    assert classification["fiber_dimension"] == 2
    assert classification["independent_coordinates"] == ["sigma_u", "sigma_d"]
    assert classification["redundant_four_tuple_map"]["rank"] == 2

    rho = float(payload["ordered_source_geometry"]["rho_ord"])
    up = classification["up_type_ray"]
    down = classification["down_type_ray"]
    assert abs(float(up["trace"])) < 1.0e-12
    assert abs(float(down["trace"])) < 1.0e-12
    assert math.isclose(float(up["endpoint_span"]), 1.0, abs_tol=1.0e-12)
    assert math.isclose(float(down["endpoint_span"]), 1.0, abs_tol=1.0e-12)
    assert math.isclose(float(up["gap21_over_gap32"]), rho, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(
        float(down["gap21_over_gap32"]),
        1.0 / rho,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert float(up["max_abs_identity_residual"]) < 1.0e-12
    assert float(down["max_abs_identity_residual"]) < 1.0e-12


def test_free_rescaling_gives_same_source_packet_and_distinct_spreads_and_means() -> None:
    payload = _build()
    action = payload["free_action_certificate"]

    assert action["group"] == "(R_{>0})^2"
    assert action["source_projection_fixed"] is True
    assert action["action_free"] is True
    assert action["action_transitive_on_positive_spread_fiber"] is True
    assert action["requested_tuple_not_invariant"] is True
    assert action["countermodels_share_source_projection"] is True
    assert action["countermodels_have_distinct_requested_tuples"] is True
    assert float(action["max_abs_structural_identity_residual"]) < 1.0e-12

    models = action["formal_countermodels"]
    assert len(models) == 2
    assert models[0]["source_projection_token"] == models[1]["source_projection_token"]
    assert models[0]["sigma_tuple"] != models[1]["sigma_tuple"]
    assert all(
        float(model["sigma_tuple"][key]) > 0.0
        for model in models
        for key in ("sigma_u", "sigma_d")
    )

    affine = payload["affine_downstream_injectivity"]
    assert affine["rank"] == 2
    assert affine["injective"] is True
    assert affine["models_have_distinct_sector_means"] is True
    assert abs(float(affine["determinant"])) > 1.0e-12
    assert math.isclose(
        float(affine["determinant"]),
        -float(affine["A_ud"]) * float(affine["B_ud"]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert abs(float(affine["determinant_residual"])) < 1.0e-12


def test_edge_packet_still_has_two_free_coefficients_without_target_residuals() -> None:
    edge = _build()["edge_statistics_do_not_select_moduli"]

    assert edge["free_coefficients"] == ["c_u", "c_d"]
    assert edge["two_points_share_source_edge_data_and_differ_in_spreads"] is True
    assert edge["unselected_existing_formula_point"]["status"] == (
        "candidate_formula_reconstructed_from_target_free_fields_not_promoted"
    )
    assert edge["zero_coefficient_counterpoint"]["status"] == "formal_countermodel_not_a_prediction"
    assert edge["unselected_existing_formula_point"]["sigma_u"] != edge["zero_coefficient_counterpoint"]["sigma_u"]
    assert edge["unselected_existing_formula_point"]["sigma_d"] != edge["zero_coefficient_counterpoint"]["sigma_d"]


def test_dependency_audit_is_target_free_and_exposes_template_ancestry() -> None:
    payload = _build()
    audit = payload["dependency_audit"]

    expected_inputs = {
        "particles/runs/flavor/family_transport_kernel.json",
        "particles/runs/flavor/generation_bundle_branch_generator.json",
        "particles/runs/flavor/generation_bundle_branch_generator_splitting_obstruction.json",
        "particles/runs/flavor/overlap_edge_transport_cocycle.json",
    }
    assert set(DIRECT_INPUT_PATHS) == expected_inputs
    assert set(audit["direct_input_paths"]) == expected_inputs
    assert set(audit["allowed_ancestors"]).isdisjoint(FORBIDDEN_ANCESTORS)
    assert audit["allowed_forbidden_intersection"] == []
    assert audit["allowed_forbidden_disjoint"] is True
    assert audit["no_target_leak"] is True
    assert audit["loads_contaminated_family_excitation_artifact"] is False
    assert audit["loads_contaminated_edge_comparison_artifact"] is False
    assert audit["loads_spread_or_sector_mean_artifact"] is False
    assert audit["loads_running_quark_reference_rows"] is False
    assert audit["loads_selected_or_current_family_exact_witness"] is False
    assert audit["loads_fitted_sigma_values"] is False
    assert audit["emits_target_sigma_values"] is False

    ancestry = payload["template_ancestry"]
    assert ancestry["family_transport_kernel"]["is_template"] is True
    assert ancestry["family_transport_kernel"]["status"] == "template"
    assert ancestry["generation_bundle_branch_generator"]["is_candidate"] is True
    assert ancestry["generation_bundle_splitting_obstruction"]["verdict"] == "not_forced_by_current_corpus"
    assert ancestry["positive_source_emission_blocked_upstream"] is True
    assert ancestry["conditional_grant_used_for_nonidentifiability_proof"] is True

    serialized = json.dumps(payload, sort_keys=True)
    for forbidden_target_value in (
        "5.579692209267639",
        "3.300314452061615",
        "5.573928426395543",
        "3.296264198808688",
        "172.3523553288311",
        "0.00216",
    ):
        assert forbidden_target_value not in serialized
