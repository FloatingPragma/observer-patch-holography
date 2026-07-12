#!/usr/bin/env python3
"""Focused guards for the fail-closed W/Z/H source-closure receipts."""

from __future__ import annotations

from collections import defaultdict, deque
import hashlib
import json
import pathlib
import sys

import jsonschema
import pytest


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from derive_boson_source_prediction_closure_audit import (  # noqa: E402
    REPO_ROOT,
    build_artifact,
    load_inputs,
)
from derive_d10_ew_quotient_transport_receipt import build_artifact as build_qt_receipt  # noqa: E402
from derive_d10_ew_target_free_repair_value_law import (  # noqa: E402
    evaluate_candidate_from_source_basis,
)


SCHEMA = HERE / "boson_source_prediction_closure_audit.schema.json"
QT_CERTIFICATE_SCHEMA = HERE / "d10_ew_quotient_path_certificate.schema.json"
GENERATED = ROOT / "particles" / "runs" / "calibration" / "boson_source_prediction_closure_audit.json"
QT_GENERATED = ROOT / "particles" / "runs" / "calibration" / "d10_ew_quotient_transport_receipt.json"


def _build() -> dict:
    payloads, manifest, constants, correspondence = load_inputs()
    return build_artifact(
        payloads,
        manifest,
        constants,
        correspondence_manifest=correspondence,
        generated_utc="2026-07-12T00:00:00Z",
    )


def _rows_by_id(rows: list[dict], key: str) -> dict[str, dict]:
    return {row[key]: row for row in rows}


def test_audit_schema_and_fail_closed_claim() -> None:
    report = _build()
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(report)
    jsonschema.Draft202012Validator(schema).validate(
        json.loads(GENERATED.read_text(encoding="utf-8"))
    )

    assert report["claim"]["prediction_promotion_allowed"] is False
    assert report["promotion_decision"]["all_required_gates_passed"] is False
    assert report["promotion_decision"]["wz_status"] == "source_incomplete"
    assert report["promotion_decision"]["higgs_status"] == (
        "conditional_declared_surface_candidate"
    )

    qt_schema = json.loads(QT_CERTIFICATE_SCHEMA.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(qt_schema)
    assert qt_schema["properties"]["evidence_kind"]["const"] == (
        "realized_carrier_enumeration_not_synthetic_fixture"
    )


def test_four_pixel_branches_and_four_mass_surfaces_stay_distinct() -> None:
    report = _build()
    branches = report["source_branches"]
    assert set(branches) == {
        "legacy_d10_P_1p63094",
        "public_endpoint_P_C",
        "source_audit_P_cand",
        "compressed_p_trunk_candidate",
    }
    assert branches["legacy_d10_P_1p63094"]["P"] == "1.63094"
    assert branches["public_endpoint_P_C"]["P"] == (
        "1.630968209403959324879279847782648941"
    )
    assert branches["source_audit_P_cand"]["P"] == (
        "1.63097209569432901817967892561191884270169"
    )
    assert branches["compressed_p_trunk_candidate"]["P"] == (
        "1.63097210492078846050203640439"
    )
    assert all(not row["eligible_for_source_only_mass_prediction"] for row in branches.values())
    assert report["branch_consistency"]["single_end_to_end_branch_selected"] is False

    surfaces = _rows_by_id(report["mass_surfaces"], "surface_id")
    assert set(surfaces) == {
        "d10_selected_current_carrier",
        "d10_running_tree_repair_candidate",
        "d10_freeze_once_reference_adapter",
        "d11_declared_surface_conditional_split",
    }
    assert surfaces["d10_selected_current_carrier"]["values"] == pytest.approx(
        {"W_GeV": 80.38629169244275, "Z_GeV": 91.18290444674243}
    )
    assert surfaces["d10_running_tree_repair_candidate"]["values"] == pytest.approx(
        {"W_GeV": 80.37700001539531, "Z_GeV": 91.18797807794321}
    )
    assert surfaces["d10_freeze_once_reference_adapter"]["values"] == pytest.approx(
        {"W_GeV": 80.3625, "Z_GeV": 91.1879}
    )
    assert surfaces["d11_declared_surface_conditional_split"]["values"] == pytest.approx(
        {"H_GeV": 125.1995304097179, "top_companion_GeV": 172.3523553288312}
    )
    assert all(not row["physical_pole_certified"] for row in surfaces.values())
    assert all(not row["promotion_allowed"] for row in surfaces.values())


def test_raw_neutral_trace_weight_and_normalized_d_are_not_conflated() -> None:
    normal = _build()["repair_coefficient_normalization"]
    n_c = normal["N_c"]
    neutral = normal["neutral_channel"]
    assert neutral["raw_hypercharge_trace_weight"] == n_c
    assert neutral["normalized_chart_coefficient_d"] == n_c / 2.0
    assert neutral["normalization_identity_residual"] == pytest.approx(0.0, abs=5.0e-18)
    assert neutral["uplift_identity_residual"] == pytest.approx(0.0, abs=5.0e-18)
    assert normal["conditionality"]["applies_to_surface"] == (
        "historical_color_balanced_quadratic_descent_only"
    )
    assert normal["conditionality"]["does_not_derive_surface"] == (
        "d10_running_tree_repair_candidate"
    )


def test_qt_receipt_verifies_algebra_but_not_the_path_certificate() -> None:
    report = _build()
    qt = report["d10_quotient_transport_review"]
    assert qt["certificate_status"] == "certificate_assumed_not_enumerated"
    assert qt["certificate_verified"] is False
    assert qt["conditional_theorem_algebra_verified"] is True
    assert qt["unconditional_d10_source_theorem_closed"] is False
    assert qt["promotion_allowed"] is False
    assert qt["domain_checks"]["eta_positive"] is True
    assert {item["id"] for item in qt["premises"]} == {"QT1", "QT2", "QT3", "QT4", "QT5"}
    assert all(not item["source_emitted_by_current_repo"] for item in qt["premises"])
    assert all(
        item["certificate_evidence_status"] == "assumed_not_enumerated"
        for item in qt["premises"]
    )
    by_id = _rows_by_id(qt["premises"], "id")
    assert "coefficients remain free" in by_id["QT3"]["missing"]
    assert "deformations remain unexcluded" in by_id["QT5"]["missing"]
    checks = qt["value_law_consistency_checks"]
    assert checks["canonical_artifact_max_abs_residual"] == 0.0
    assert checks["parallel_fibre_identity_residual"] == pytest.approx(0.0, abs=1.0e-18)
    assert checks["alphaY_recomposition_residual"] == pytest.approx(0.0, abs=1.0e-18)
    assert checks["mass_chart_jacobian_det"] > 0.0
    assert checks["mass_chart_nondegenerate"] is True
    assert qt["relation_to_color_balanced_candidate"]["same_law"] is False
    assert qt["relation_to_color_balanced_candidate"][
        "models_are_alternatives_not_cumulative_gates"
    ] is True
    assert qt["physical_pole_attachment"]["required_convention"] == (
        "s_B=(M_B-i*Gamma_B/2)^2"
    )
    assert qt["physical_pole_attachment"]["analytic_continuation_and_riemann_sheet_fixed"] is False
    contract = qt["certificate_contract"]
    assert contract["schema_present"] is True
    assert contract["certificate_input_present"] is False
    assert contract["certificate_input_status"] == "missing_no_path_enumeration_supplied"
    assert contract["schema_conformance_is_sufficient_for_source_evidence"] is False
    assert contract["synthetic_path_fixture_eligible_as_source_evidence"] is False
    assert contract["simulation_or_exact_enumeration_required"] is True

    source = json.loads(
        (ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json").read_text()
    )
    value_law = json.loads(
        (ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json").read_text()
    )
    assert build_qt_receipt(source, value_law)["value_law_consistency_checks"] == (
        json.loads(QT_GENERATED.read_text())["value_law_consistency_checks"]
    )


def test_pure_candidate_helper_has_no_reference_argument_and_requires_positive_eta() -> None:
    evaluated = evaluate_candidate_from_source_basis(
        alpha_2=0.03377843630219015,
        alpha_y=0.010131601067241624,
        eta_source=0.022147000871961295,
        v_value=246.76711732749683,
    )
    assert evaluated["coherent_emitted_quintet"]["MW_pole"] == pytest.approx(
        80.37700001539531, abs=1.0e-12
    )
    with pytest.raises(ValueError, match="must be positive"):
        evaluate_candidate_from_source_basis(
            alpha_2=0.03377843630219015,
            alpha_y=0.010131601067241624,
            eta_source=0.0,
            v_value=246.76711732749683,
        )


def test_twelve_gate_inventory_treats_d10_models_as_alternatives() -> None:
    report = _build()
    gates = _rows_by_id(report["theorem_gates"], "theorem_id")
    assert len(gates) == 12
    assert gates["T05_d10_source_uniqueness"]["gate_passed_for_full_source_prediction"] is False
    assert gates["T06_conditional_color_weights"]["gate_role"] == "alternative_model"
    assert "T06_conditional_color_weights" not in report["promotion_decision"]["failed_gate_ids"]
    assert "alternatives, not cumulative premises" in report["promotion_decision"][
        "alternative_model_policy"
    ]
    assert gates["T12_full_wzh_source_prediction"]["gate_passed_for_full_source_prediction"] is False


def test_scale_rg_pole_and_full_source_dag_remain_open() -> None:
    report = _build()
    scale = report["absolute_scale_contract"]
    assert scale["source_closed"] is False
    assert scale["legacy_dimensionful_constant"]["name"] == "E_PLANCK_GEV"
    assert scale["legacy_dimensionful_constant"]["value"] == 1.220890e19

    rg = report["rg_matching_contract"]
    assert rg["status"] == "closed_declared_convention_contract_not_rg_matching_theorem"
    assert rg["concrete_transport_receipt_present"] is False
    assert rg["promotion_allowed"] is False

    pole = report["physical_pole_contract"]
    assert pole["required_mass_convention"] == "s_B=(M_B-i*Gamma_B/2)^2"
    assert pole["analytic_continuation_and_riemann_sheet_fixed"] is False
    assert pole["self_energy_functions_present"] is False
    assert pole["BRST_complete_mixing_blocks_present"] is False
    assert pole["Ward_Slavnov_Taylor_and_Nielsen_receipts_present"] is False
    assert pole["physical_residue_certificate_present"] is False
    assert pole["promotion_allowed"] is False

    d11_gate = _rows_by_id(report["theorem_gates"], "theorem_id")[
        "T08_conditional_d11_higgs_closure"
    ]
    assert any("DS1--DS5" in item for item in d11_gate["open_requirements"])

    dag = report["source_dependency_dag"]
    assert dag["single_end_to_end_source_branch_selected"] is False
    assert dag["full_file_level_runtime_ancestry_certified"] is False
    assert dag["human_formula_selection_ancestry_certified"] is False
    graph: dict[str, list[str]] = defaultdict(list)
    indegree: dict[str, int] = {node["id"]: 0 for node in dag["nodes"]}
    for edge in dag["edges"]:
        graph[edge["from"]].append(edge["to"])
        indegree[edge["to"]] += 1
    queue = deque(node for node, degree in indegree.items() if degree == 0)
    visited = 0
    while queue:
        node = queue.popleft()
        visited += 1
        for child in graph[node]:
            indegree[child] -= 1
            if indegree[child] == 0:
                queue.append(child)
    assert visited == len(indegree)


def test_input_hashes_are_deterministic_and_scoped() -> None:
    report = _build()
    repo_entries = report["input_manifest"]
    assert [entry["path"] for entry in repo_entries] == sorted(
        entry["path"] for entry in repo_entries
    )
    for entry in repo_entries:
        assert entry["path"].startswith("code/")
        data = (REPO_ROOT / entry["path"]).read_bytes()
        assert entry["sha256"] == hashlib.sha256(data).hexdigest()
        assert entry["bytes"] == len(data)

    # Reviewed workspace correspondence is integrated semantically, not a
    # runtime dependency of the standalone repository build.
    assert report["reviewed_correspondence_manifest"] == []


def test_regenerated_candidate_and_d11_artifacts_state_the_conditional_boundary() -> None:
    value_law = json.loads(
        (ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json").read_text()
    )
    assert "does not derive or uniquely select" in value_law["theorem"]["statement"]
    assert value_law["promotion_allowed"] is False

    d11 = json.loads(
        (ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json").read_text()
    )
    assert d11["status"] == "candidate_only"
    assert d11["prediction_promotion_allowed"] is False
    assert any("not a full source-only mass prediction" in note for note in d11["notes"])
    assert "full_source_only_W_Z_H_complex_pole_prediction" in d11["strictly_not_claimed"]
