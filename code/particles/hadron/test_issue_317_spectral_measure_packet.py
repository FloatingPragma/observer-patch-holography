#!/usr/bin/env python3
"""Tests for the issue #317 conditional spectral-measure theorem and contract packet.

Coverage, following the re-scoped acceptance list:

(a) the conditional theorem is stated with typed premises and conclusions,
    and the positivity/threshold conclusions require the
    reflection-positivity premise;
(b) the executed witnesses hold, including the gauge witness on the actual
    link-inserted point-split conserved current, and every diagnostic
    witness is labeled non-promoting;
(c) the strict validator fails closed on every adversarial payload class
    the review enumerated (nested measured-HVP provenance, stringly-typed
    target lists, s != E^2, dangling and duplicate identifiers, NaN and
    infinite numerics, unquantified and unordered budgets, covariance
    defects, weight/residue inconsistency, unknown keys, unknown statuses);
(d) the adversarial battery covers every declared semantic requirement;
(e) every gate-approved payload is accepted by the real downstream
    source-transport validator (executed implication);
(f) the three verdicts (contract_certified, physical_source_payload_available,
    physical_promotion_allowed) are computed independently and the physical
    availability predicate fails closed on unknown, renamed, missing, None,
    blocked, and contradictory status combinations;
(g) the stored certificate matches a fresh emission on all
    acceptance-relevant fields.
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import verify_issue_317_spectral_measure_packet as packet_mod
import ward_projected_spectral_measure_validator as strict_validator

STORED = (
    HERE.parent
    / "runs"
    / "hadron"
    / "ward_projected_spectral_measure_proof_packet.json"
)


@pytest.fixture(scope="module")
def packet() -> dict:
    return packet_mod.build_packet()


@pytest.fixture(scope="module")
def stored() -> dict:
    return json.loads(STORED.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def schema() -> dict:
    return strict_validator.load_schema()


def _reject(payload: dict, schema: dict) -> tuple[str, ...]:
    result = strict_validator.validate_production_payload(payload, schema)
    assert result.accepted is False, "payload must be rejected"
    return result.reasons


# ---------------------------------------------------------------------------
# Conditional theorem structure.
# ---------------------------------------------------------------------------


def test_theorem_typed_premises_and_conclusions(packet):
    theorem = packet["theorem"]
    assert theorem["theorem_id"] == "SpectralMeasureQ_ConditionalConstruction"
    premises = theorem["premises"]
    conclusions = theorem["conclusions"]
    assert len(premises) == 10
    assert len(conclusions) == 8
    premise_ids = {p["id"] for p in premises}
    for premise in premises:
        assert premise["premise_type"] == "declared_premise_validated_in_representation_only"
        assert premise["physical_value_produced_here"] is False
        assert premise["contract_fields"], premise["id"]
    for conclusion in conclusions:
        assert conclusion["uses"], conclusion["id"]
        assert set(conclusion["uses"]) <= premise_ids
    assert packet["theorem_structure_check"]["passed"] is True


def test_positivity_conclusion_conditional_on_reflection_positivity(packet):
    """C3/C4 must require premise P4; no sign check may substitute."""
    conclusions = {c["id"]: c for c in packet["theorem"]["conclusions"]}
    assert "P4" in conclusions["C3"]["uses"]
    assert "P4" in conclusions["C4"]["uses"]
    assert conclusions["C3"]["implementation_witness"] is None
    sign_check = packet["machine_witnesses"]["demonstrator_measure"]["rho_had_or_measure"][
        "sign_check_report"
    ]
    assert sign_check["not_a_proof_of_spectral_positivity"] is True
    assert "no statistically significant sign violation" in sign_check["statement"]


def test_packet_accepted_and_all_criteria_machine_checked(packet):
    assert packet["accepted"] is True
    assert "theorem-and-contract layer only" in packet["acceptance_scope"]
    criteria = packet["acceptance_criteria_status"]
    assert len(criteria) == 14
    for name, row in criteria.items():
        assert row["packet_level_passed"] is True, name
        assert row["machine_checks"], name


# ---------------------------------------------------------------------------
# Executed witnesses.
# ---------------------------------------------------------------------------


def test_gauge_witness_exercises_conserved_point_split_current(packet):
    gw = packet["machine_witnesses"]["gauge_and_ward"]
    conserved = gw["conserved_current_gauge_invariance"]
    assert conserved["passed"] is True
    assert conserved["correlator_relative_defect"] < conserved["tolerance"]
    assert "link-inserted point-split conserved current" in conserved["observable"]
    local = gw["local_current_gauge_invariance"]
    assert local["passed"] is True
    assert local["correlator_relative_defect"] < local["tolerance"]
    assert local["plaquette_absolute_defect"] < local["tolerance"]


def test_ward_identity_witness(packet):
    ward = packet["machine_witnesses"]["gauge_and_ward"]["ward_identity"]
    assert ward["passed"] is True
    assert ward["relative_defect"] < ward["tolerance"]


def test_zv_free_field_anchor_scoped_as_convention_check(packet):
    zv = packet["machine_witnesses"]["zv_free_field_anchor"]
    assert zv["passed"] is True
    assert zv["max_abs_deviation_from_one"] < zv["tolerance"]
    assert "not an interacting plateau determination" in zv["statement"]


def test_demonstrator_labeled_diagnostic_with_narrowed_claims(packet):
    demo = packet["machine_witnesses"]["demonstrator_measure"]
    assert demo["checks_passed"] is True
    assert demo["non_promoting_diagnostic"] is True
    assert demo["guards"]["promotion_allowed"] is False
    assert demo["guards"]["production_schema_conformant"] is False
    assert demo["ensemble"]["external_inputs_used"] is False
    # positivity language must be the demoted sign check
    assert demo["rho_had_or_measure"]["positivity_status"].startswith(
        "diagnostic_sign_check_passed"
    )
    # covariance claims are narrowed and the joint covariance is emitted
    assert "correlator covariance and nothing more" in demo["covariance"]["object"]
    joint = demo["joint_extraction_covariance"]
    assert joint["row_basis"] == ["am_v", "z_v", "amplitude_hop", "weight_phys"]
    matrix = np.array(joint["matrix"])
    assert matrix.shape == (joint["dimension"], joint["dimension"])
    assert np.allclose(matrix, matrix.T)
    assert joint["positive_semidefinite"] is True
    assert any("windowed average" in row for row in demo["diagnostic_limitations"])
    corr = np.array(demo["covariance"]["matrix"])
    assert corr.shape == (demo["covariance"]["dimension"], demo["covariance"]["dimension"])
    assert np.allclose(corr, corr.T)
    assert demo["covariance"]["positive_semidefinite"] is True


def test_demonstrator_residue_normalization_chain(packet):
    """w_phys = Z_V^2 (2 kappa)^2 * 24 pi^2 A / m, per conserved_vector.py:
    the hopping-to-physical factor is (2 kappa)^2 and the atom amplitude is
    A = w m / (24 pi^2) (vector_correlator.synthetic_atom_correlator)."""
    demo = packet["machine_witnesses"]["demonstrator_measure"]
    kappa = demo["ensemble"]["kappa"]
    ex = demo["extraction"]
    w_hop = 24.0 * math.pi**2 * ex["amplitude_hop"] / ex["am_vector"]
    expected = ex["z_v"] ** 2 * (2.0 * kappa) ** 2 * w_hop
    assert ex["weight_physical_units"] == pytest.approx(expected, rel=1e-12)


def test_all_diagnostic_witnesses_labeled_non_promoting(packet):
    witnesses = packet["machine_witnesses"]
    assert witnesses["gauge_and_ward"]["non_promoting_diagnostic"] is True
    assert witnesses["zv_free_field_anchor"]["non_promoting_diagnostic"] is True
    assert witnesses["demonstrator_measure"]["non_promoting_diagnostic"] is True
    row = packet["acceptance_criteria_status"][
        "diagnostic_witnesses_explicitly_labeled_non_promoting"
    ]
    assert row["packet_level_passed"] is True


# ---------------------------------------------------------------------------
# Strict validator: adversarial battery (each review-reproduced acceptance
# must now be a rejection).
# ---------------------------------------------------------------------------


def test_conformant_payload_accepted(schema):
    result = strict_validator.validate_production_payload(
        packet_mod.build_conformant_payload(), schema
    )
    assert result.accepted is True
    assert result.reasons == ()


def test_nested_measured_hvp_provenance_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["measured_hvp_input"] = {"source": "EE_TO_HADRONS"}
    reasons = _reject(payload, schema)
    assert any(r.startswith("TARGET_LEAK_DETECTED") for r in reasons)
    assert any(r.startswith("schema_violation") for r in reasons)


def test_external_targets_used_as_string_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["external_targets_used"] = "EE_TO_HADRONS"
    reasons = _reject(payload, schema)
    assert any(r.startswith("TARGET_LEAK_DETECTED") for r in reasons)
    assert "external_targets_used_not_a_list" in reasons


def test_s_not_equal_energy_squared_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["finite_volume_levels"][0]["levels"][0]["s"] = 999.0
    reasons = _reject(payload, schema)
    assert any(r.startswith("s_not_equal_energy_squared") for r in reasons)


def test_dangling_residue_reference_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["ward_projected_residues"][0]["level_id"] = "NO_SUCH_LEVEL"
    reasons = _reject(payload, schema)
    assert "residue_references_unknown_level:NO_SUCH_LEVEL" in reasons
    assert "level_without_residue:L0" in reasons


def test_duplicate_level_ids_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["finite_volume_levels"][0]["levels"].append(
        {"level_id": "L0", "s": 4.0, "energy": 2.0, "weight": 0.25}
    )
    reasons = _reject(payload, schema)
    assert "duplicate_level_id:L0" in reasons


@pytest.mark.parametrize(
    "mutate,expected_fragment",
    [
        (lambda p: p["finite_volume_levels"][0]["levels"][0].__setitem__("energy", float("nan")), "nonfinite_number"),
        (lambda p: p["finite_volume_levels"][0]["levels"][0].__setitem__("weight", float("inf")), "nonfinite_number"),
        (lambda p: p["ward_projected_residues"][0].__setitem__("residue", float("nan")), "nonfinite_number"),
        (lambda p: p["finite_volume_levels"][0]["levels"][0].__setitem__("s", float("inf")), "nonfinite_number"),
        (lambda p: p["covariance"]["matrix"][0].__setitem__(0, float("nan")), "nonfinite_number"),
    ],
)
def test_nonfinite_values_rejected(schema, mutate, expected_fragment):
    payload = packet_mod.build_conformant_payload()
    mutate(payload)
    reasons = _reject(payload, schema)
    assert any(expected_fragment in r for r in reasons)


def test_all_budgets_status_not_quantified_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["systematics"] = {
        name: {"status": "not_quantified"} for name in payload["systematics"]
    }
    reasons = _reject(payload, schema)
    assert any(r.startswith("budget_bound_interval_missing") for r in reasons)


def test_unordered_and_nonfinite_bound_intervals_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["systematics"]["statistical_budget"]["bound_interval"] = {"lo": 0.02, "hi": 0.01}
    reasons = _reject(payload, schema)
    assert any("invalid_interval_order" in r for r in reasons)

    payload = packet_mod.build_conformant_payload()
    payload["systematics"]["continuum_budget"]["bound_interval"] = {"lo": 0.0, "hi": "inf"}
    reasons = _reject(payload, schema)
    assert any("nonfinite_number" in r for r in reasons)


def test_missing_transport_moment_certificate_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    del payload["transport_moment_certificate"]
    reasons = _reject(payload, schema)
    assert "transport_moment_certificate_missing" in reasons


def test_covariance_defects_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["covariance"]["matrix"] = [[1.0e-4, 5.0e-5], [0.0, 1.0e-4]]
    assert "covariance_not_symmetric" in _reject(payload, schema)

    payload = packet_mod.build_conformant_payload()
    payload["covariance"]["matrix"] = [[1.0e-4, 2.0e-4], [2.0e-4, 1.0e-4]]
    assert "covariance_not_positive_semidefinite" in _reject(payload, schema)

    payload = packet_mod.build_conformant_payload()
    payload["covariance"]["dimension"] = 3
    reasons = _reject(payload, schema)
    assert any("covariance" in r and "dimension" in r for r in reasons)


def test_weight_residue_inconsistency_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["ward_projected_residues"][0]["residue"] = 0.4
    assert "weight_residue_inconsistent:L0" in _reject(payload, schema)


def test_negative_and_nonpositive_values_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["ward_projected_residues"][0]["residue"] = -0.5
    payload["finite_volume_levels"][0]["levels"][0]["weight"] = -0.5
    reasons = _reject(payload, schema)
    assert any(r.startswith("negative_residue") for r in reasons)
    assert any(r.startswith("negative_level_weight") for r in reasons)

    payload = packet_mod.build_conformant_payload()
    payload["finite_volume_levels"][0]["levels"][0]["energy"] = 0.0
    payload["finite_volume_levels"][0]["levels"][0]["s"] = 0.0
    reasons = _reject(payload, schema)
    assert any(r.startswith("nonpositive_energy") for r in reasons)


def test_typing_and_guard_violations_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    payload["branch"]["flavors"] = "quenched"
    reasons = _reject(payload, schema)
    assert any("branch_flavors_not_allowlisted" in r for r in reasons)

    payload = packet_mod.build_conformant_payload()
    payload["guards"]["compare_only_external_endpoint"] = "false"
    reasons = _reject(payload, schema)
    assert "guard_not_literal_false:compare_only_external_endpoint" in reasons

    payload = packet_mod.build_conformant_payload()
    payload["projection"]["ward_projected"] = False
    reasons = _reject(payload, schema)
    assert "ward_projection_dropped" in reasons


def test_provenance_violations_rejected(schema):
    payload = packet_mod.build_conformant_payload()
    del payload["provenance"]
    assert "provenance_manifest_missing" in _reject(payload, schema)

    payload = packet_mod.build_conformant_payload()
    payload["provenance"]["source_inputs"].append(
        {"kind": "external_dataset", "identifier": "unknown"}
    )
    reasons = _reject(payload, schema)
    assert any("provenance_source_input_kind_not_allowlisted" in r for r in reasons)

    payload = packet_mod.build_conformant_payload()
    payload["provenance"]["measured_hvp_input_present"] = True
    reasons = _reject(payload, schema)
    assert "provenance_measured_hvp_input_present_not_false" in reasons

    payload = packet_mod.build_conformant_payload()
    payload["provenance"]["source_inputs"] = []
    reasons = _reject(payload, schema)
    assert "provenance_source_inputs_missing_or_empty" in reasons


def test_unknown_keys_rejected_at_closed_boundaries(schema):
    payload = packet_mod.build_conformant_payload()
    payload["undeclared_extra_block"] = {"note": "unknown"}
    reasons = _reject(payload, schema)
    assert any(r.startswith("schema_violation") for r in reasons)

    payload = packet_mod.build_conformant_payload()
    payload["guards"]["undeclared_guard"] = False
    reasons = _reject(payload, schema)
    assert any(r.startswith("schema_violation") for r in reasons)


def test_unknown_positivity_status_rejected(schema):
    for status in ("checked_numerically", "diagnostic_sign_check_passed", "", None):
        payload = packet_mod.build_conformant_payload()
        payload["rho_had_or_measure"]["positivity_status"] = status
        reasons = _reject(payload, schema)
        assert any("positivity_status_not_allowlisted" in r for r in reasons), status


def test_downstream_forbidden_keys_rejected_by_gate(schema):
    payload = packet_mod.build_conformant_payload()
    payload["backend"]["run_id"] = "gate"
    payload["c_Q"] = 0.5
    reasons = _reject(payload, schema)
    assert any("downstream_forbidden_key:c_Q" in r for r in reasons)


def test_adversarial_battery_covers_every_semantic_requirement(packet):
    gate = packet["machine_witnesses"]["acceptance_gate"]
    assert gate["passed"] is True
    assert gate["conformant_payload_accepted"] is True
    assert gate["all_negative_controls_rejected"] is True
    assert gate["adversarial_coverage_complete"] is True
    assert gate["uncovered_semantic_requirements"] == []
    coverage = gate["adversarial_coverage"]
    assert set(coverage) == set(strict_validator.SEMANTIC_REQUIREMENTS)
    for requirement, control_ids in coverage.items():
        assert control_ids, requirement
    for control in gate["negative_controls"]:
        assert control["rejected"] is True, control["control_id"]
        assert control["rejection_reasons"], control["control_id"]


# ---------------------------------------------------------------------------
# Downstream consumer compatibility (executed implication).
# ---------------------------------------------------------------------------


def test_gate_approved_payloads_accepted_by_downstream_validator(packet):
    downstream = packet["machine_witnesses"]["downstream_compatibility"]
    assert downstream["passed"] is True
    variant_ids = {row["variant_id"] for row in downstream["variants"]}
    assert {"conformant_baseline", "multi_level", "two_ensembles", "string_decimal_bounds"} <= variant_ids
    for row in downstream["variants"]:
        assert row["gate_accepted"] is True, row["variant_id"]
        assert row["downstream_source_measure_reasons"] == [], row["variant_id"]
        assert row["downstream_forbidden_keys"] == [], row["variant_id"]
        assert row["implication_holds"] is True, row["variant_id"]


def test_downstream_implication_holds_directly():
    """Re-execute the implication outside the packet against the live modules."""
    sys.path.insert(0, str(packet_mod.P_DERIVATION))
    from thomson_spectral_transport import (
        _validate_source_measure,
        source_payload_forbidden_keys,
    )

    schema = strict_validator.load_schema()
    for variant_id, payload in packet_mod.build_gate_approved_variants():
        result = strict_validator.validate_production_payload(payload, schema)
        assert result.accepted is True, (variant_id, result.reasons)
        reasons: list[str] = []
        _validate_source_measure({"source_measure": payload}, reasons)
        assert reasons == [], (variant_id, reasons)
        assert source_payload_forbidden_keys({"source_measure": payload}) == set()


# ---------------------------------------------------------------------------
# Three verdicts and fail-closed closure predicates.
# ---------------------------------------------------------------------------


def test_certificate_distinguishes_three_verdicts(packet):
    verdicts = packet["verdicts"]
    assert set(verdicts) == {
        "contract_certified",
        "physical_source_payload_available",
        "physical_promotion_allowed",
    }
    assert verdicts["contract_certified"]["value"] is True
    assert verdicts["physical_source_payload_available"]["value"] is False
    assert verdicts["physical_promotion_allowed"]["value"] is False
    availability = verdicts["physical_source_payload_available"]["computed_from"]
    assert availability["reasons"], "unavailability must carry typed reasons"
    assert any("production_payload_absent" in r for r in availability["reasons"])
    promotion = verdicts["physical_promotion_allowed"]["computed_from"]
    assert promotion["reasons"]
    assert "physical_source_payload_not_available" in promotion["reasons"]


def test_availability_predicate_fails_closed_on_unknown_statuses():
    all_success = {
        "export_bundle_status": "complete",
        "closure_grade": "execution_complete",
        "public_unsuppression_ready": True,
        "base_measure_status": "POPULATED",
        "ward_current_status": "SOURCE_CERTIFICATE_VERIFIED",
        "production_payload_path": packet_mod.PRODUCTION_PAYLOAD_PATH,
    }
    adversarial_states = [
        {**all_success, "export_bundle_status": "blocked"},
        {**all_success, "export_bundle_status": "done"},
        {**all_success, "base_measure_status": "MISSING"},
        {**all_success, "ward_current_status": None},
        {**all_success, "closure_grade": "complete_pending_review"},
        {**all_success, "public_unsuppression_ready": False},
        {**all_success, "public_unsuppression_ready": "true"},
        {**all_success, "public_unsuppression_ready": 1},
        {},
    ]
    for state in adversarial_states:
        verdict = packet_mod.physical_source_payload_verdict(state)
        assert verdict["available"] is False, state
        assert verdict["reasons"], state
    # even with every status in its success allowlist, the actual
    # gate-approved payload must exist: absence fails closed
    verdict = packet_mod.physical_source_payload_verdict(dict(all_success))
    assert verdict["available"] is False
    assert any("production_payload_absent" in r for r in verdict["reasons"])


def test_availability_requires_gate_approved_payload(tmp_path):
    """Success statuses plus an existing but non-conforming payload still fail."""
    bad_payload_path = tmp_path / "ward_projected_spectral_measure.production.json"
    bad_payload_path.write_text(json.dumps({"artifact": "wrong"}), encoding="utf-8")
    state = {
        "export_bundle_status": "complete",
        "closure_grade": "execution_complete",
        "public_unsuppression_ready": True,
        "base_measure_status": "POPULATED",
        "ward_current_status": "SOURCE_CERTIFICATE_VERIFIED",
        "production_payload_path": bad_payload_path,
    }
    verdict = packet_mod.physical_source_payload_verdict(state)
    assert verdict["available"] is False
    assert any("production_payload_rejected_by_strict_gate" in r for r in verdict["reasons"])

    good_payload_path = tmp_path / "good.production.json"
    good_payload_path.write_text(
        json.dumps(packet_mod.build_conformant_payload()), encoding="utf-8"
    )
    state["production_payload_path"] = good_payload_path
    verdict = packet_mod.physical_source_payload_verdict(state)
    assert verdict["available"] is True
    assert verdict["reasons"] == []


def test_fail_closed_probe_battery_in_packet(packet):
    probes = packet["machine_witnesses"]["fail_closed_probes"]
    assert probes["passed"] is True
    probe_ids = {row["probe_id"] for row in probes["probes"]}
    assert {
        "status_blocked",
        "status_missing_string",
        "status_none",
        "status_renamed",
        "closure_grade_unknown",
        "public_unsuppression_not_ready",
        "public_unsuppression_stringly_true",
        "all_keys_missing",
        "all_statuses_success_but_payload_absent",
    } <= probe_ids
    for row in probes["probes"]:
        assert row["reported_unavailable"] is True, row["probe_id"]


def test_live_closure_statuses_mirror_artifacts(packet):
    snapshot = packet["verdicts"]["physical_source_payload_available"]["computed_from"][
        "live_status_snapshot"
    ]
    live_base = json.loads(packet_mod.BASE_MEASURE.read_text(encoding="utf-8"))
    assert snapshot["base_measure_status"] == live_base["status"]
    live_ward = json.loads(packet_mod.WARD_CURRENT_DEFINITION.read_text(encoding="utf-8"))
    assert snapshot["ward_current_status"] == live_ward["status"]


# ---------------------------------------------------------------------------
# Typing criteria (live artifact checks).
# ---------------------------------------------------------------------------


def test_typing_criteria_are_live_checks_not_hardcoded(packet):
    empirical = packet["machine_witnesses"]["empirical_typing"]
    assert empirical["passed"] is True
    assert all(empirical["checks"].values()), empirical["checks"]
    assert "comparison_data_manifest.json" in empirical["checked_artifacts"]["comparison_manifest"]
    assert empirical["checks"]["gate_rejects_nested_measured_hvp_provenance"] is True
    higher = packet["machine_witnesses"]["higher_point_typing"]
    assert higher["passed"] is True
    assert all(higher["checks"].values()), higher["checks"]
    assert higher["xi_ledger_source_certificate_status"] == "MISSING_SOURCE_CERTIFICATE"


def test_blinding_and_disclosure_fields_explicit(packet):
    empirical = packet["machine_witnesses"]["empirical_typing"]
    blinding = empirical["blinding_status"]
    assert blinding["comparison_event_performed"] is False
    disclosure = empirical["public_use_disclosure"]
    assert disclosure["usable_for_public_final_values_const"] is True
    assert "comparison-only" in disclosure["confinement"]


# ---------------------------------------------------------------------------
# Claim boundary, dependency edges, stored-vs-fresh agreement.
# ---------------------------------------------------------------------------


def test_claim_boundary_and_dependency_edges(packet):
    boundary = packet["claim_boundary"]
    assert any("production backend execution" in row for row in boundary["not_closed_here"])
    assert any("conditional on" in row for row in boundary["not_closed_here"])
    assert "physical promotion" in " ".join(boundary["not_closed_here"])
    dep = packet["dependency_note"]
    assert "upstream_of_425_not_downstream" not in dep
    assert dep["contract_design_edge"]["direction"].startswith("#317 contract specification")
    assert dep["execution_dependency_edge"]["direction"].startswith("#425 physical payload")
    assert packet["forbidden_inputs"] == packet_mod.FORBIDDEN_TARGETS
    closure = packet["issue_closure_condition"]
    assert closure["theorem_and_contract_acceptance_list_passes"] is True
    assert closure["original_scope_production_condition"]["met_locally"] is False
    assert "resolves this issue under the re-scoped acceptance list" in closure[
        "closing_keyword_policy"
    ]
    assert "production" in closure["closing_keyword_policy"]


def test_stored_packet_matches_fresh_emission(packet, stored):
    assert stored["issue"] == 317
    assert stored["certificate_id"] == packet["certificate_id"]
    assert stored["artifact"] == packet["artifact"]
    assert stored["status"] == packet["status"]
    assert stored["accepted"] is True
    assert stored["theorem"] == packet["theorem"]
    assert stored["theorem_structure_check"] == packet["theorem_structure_check"]
    assert stored["derivation_chain"] == packet["derivation_chain"]
    assert stored["acceptance_criteria_status"] == packet["acceptance_criteria_status"]
    assert stored["claim_boundary"] == packet["claim_boundary"]
    assert stored["dependency_note"] == packet["dependency_note"]
    assert stored["issue_closure_condition"] == packet["issue_closure_condition"]
    for verdict_name, verdict in packet["verdicts"].items():
        assert stored["verdicts"][verdict_name]["value"] == verdict["value"], verdict_name
    stored_avail = stored["verdicts"]["physical_source_payload_available"]["computed_from"]
    fresh_avail = packet["verdicts"]["physical_source_payload_available"]["computed_from"]
    assert stored_avail["reasons"] == fresh_avail["reasons"]
    assert stored_avail["live_status_snapshot"] == fresh_avail["live_status_snapshot"]
    gate_stored = stored["machine_witnesses"]["acceptance_gate"]
    gate_fresh = packet["machine_witnesses"]["acceptance_gate"]
    assert gate_stored == gate_fresh
    ds_stored = stored["machine_witnesses"]["downstream_compatibility"]
    ds_fresh = packet["machine_witnesses"]["downstream_compatibility"]
    assert ds_stored == ds_fresh
    probes_stored = stored["machine_witnesses"]["fail_closed_probes"]
    probes_fresh = packet["machine_witnesses"]["fail_closed_probes"]
    assert probes_stored == probes_fresh
    # numeric witnesses agree to tolerance (identical seeds and inputs)
    demo_stored = stored["machine_witnesses"]["demonstrator_measure"]["extraction"]
    demo_fresh = packet["machine_witnesses"]["demonstrator_measure"]["extraction"]
    assert demo_stored["am_vector"] == pytest.approx(demo_fresh["am_vector"], rel=1e-10)
    assert demo_stored["weight_physical_units"] == pytest.approx(
        demo_fresh["weight_physical_units"], rel=1e-10
    )
    gw_stored = stored["machine_witnesses"]["gauge_and_ward"]
    gw_fresh = packet["machine_witnesses"]["gauge_and_ward"]
    assert gw_stored["conserved_current_gauge_invariance"]["passed"] is True
    assert gw_fresh["conserved_current_gauge_invariance"]["passed"] is True


def test_stored_schema_is_strict(schema):
    """The shared export schema must be closed at every provenance-bearing boundary."""
    assert schema["additionalProperties"] is False
    for key in ("guards", "provenance", "branch", "projection", "covariance"):
        assert schema["properties"][key]["additionalProperties"] is False, key
    assert set(schema["required"]) >= {
        "provenance",
        "external_targets_used",
        "covariance",
        "transport_moment_certificate",
    }
    level_schema = schema["properties"]["finite_volume_levels"]["items"]["properties"][
        "levels"
    ]["items"]
    assert level_schema["additionalProperties"] is False
    assert level_schema["properties"]["energy"]["exclusiveMinimum"] == 0
    assert schema["properties"]["rho_had_or_measure"]["properties"]["positivity_status"][
        "enum"
    ] == ["certified_positive"]
