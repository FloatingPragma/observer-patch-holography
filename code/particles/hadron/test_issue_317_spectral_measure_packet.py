#!/usr/bin/env python3
"""Tests for the issue #317 Ward-projected spectral-measure proof packet.

Coverage: (a) the stored packet certificate matches a fresh emission on all
acceptance-relevant fields, (b) the executed witnesses hold (gauge
invariance, exact Ward identity, free-field Z_V anchor, demonstrator
positivity and covariance), (c) the acceptance gate accepts the conformant
payload and rejects every negative control, (d) the typing criteria are
live artifact checks rather than hard-coded booleans, (e) the author-stated
issue-closure condition is computed from live production-lane artifacts and
reports open, (f) the claim boundary and dependency direction (packet
upstream of issue #425) stay explicit.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import verify_issue_317_spectral_measure_packet as packet_mod

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


def test_packet_accepted_and_all_criteria_machine_checked(packet):
    assert packet["accepted"] is True
    assert "packet construction/contract layer only" in packet["acceptance_scope"]
    criteria = packet["acceptance_criteria_status"]
    assert len(criteria) == 5
    for name, row in criteria.items():
        assert row["packet_level_passed"] is True, name
        assert row["machine_checks"], name


def test_typing_criteria_are_live_checks_not_hardcoded(packet):
    """Criteria 4 and 5 must be derived from live artifact reads."""
    empirical = packet["machine_witnesses"]["empirical_typing"]
    assert empirical["passed"] is True
    assert all(empirical["checks"].values()), empirical["checks"]
    assert "comparison_data_manifest.json" in empirical["checked_artifacts"]["comparison_manifest"]
    higher = packet["machine_witnesses"]["higher_point_typing"]
    assert higher["passed"] is True
    assert all(higher["checks"].values()), higher["checks"]
    assert set(empirical["checks"]) == {
        "comparison_manifest_empty",
        "comparison_manifest_status_no_data",
        "production_schema_rejects_compare_only_const_false",
        "production_schema_rejects_surrogate_const_false",
        "gate_rejects_measured_hvp_comparison_endpoint",
        "gate_rejects_ee_to_hadrons_target_leak",
        "row_class_const_empirical_closure",
        "empirical_surface_not_source_only_const",
        "empirical_use_disclosed_const_true",
        "promotable_as_oph_source_theorem_const_false",
        "satisfies_production_constructive_next_artifact_const_false",
    }
    assert set(higher["checks"]) == {
        "q4_status_two_point_measure_insufficient",
        "q4_promotion_not_allowed",
        "q4_no_external_targets",
        "xi_ledger_distinct_artifact",
        "xi_required_terms_include_qed_and_ew",
        "xi_promotion_not_allowed",
        "xi_no_external_targets",
        "production_branch_excludes_qed",
    }


def test_blinding_and_disclosure_fields_explicit(packet):
    """The blinding state and the empirical public-use disclosure must be
    explicit fields, and the radiative ledger's missing source certificate
    must be disclosed rather than absorbed."""
    empirical = packet["machine_witnesses"]["empirical_typing"]
    blinding = empirical["blinding_status"]
    assert blinding["comparison_event_performed"] is False
    disclosure = empirical["public_use_disclosure"]
    assert disclosure["usable_for_public_final_values_const"] is True
    assert "comparison-only" in disclosure["confinement"]
    higher = packet["machine_witnesses"]["higher_point_typing"]
    assert higher["xi_ledger_source_certificate_status"] == "MISSING_SOURCE_CERTIFICATE"
    criteria = packet["acceptance_criteria_status"]
    radiative_open = criteria["higher_point_and_radiative_corrections_separately_typed"][
        "production_open_items"
    ]
    assert any("Xi_same_scheme" in row for row in radiative_open)


def test_issue_closure_condition_computed_and_open(packet):
    closure = packet["issue_closure_condition"]
    assert closure["met_locally"] is False
    checked = closure["checked_objects"]
    assert checked["production_backend_export_bundle"]["status"] == "open"
    assert checked["production_backend_export_bundle"]["closure_grade"] == "execution_incomplete"
    assert checked["qcd_base_measure"]["status"] == "REQUIRED_NOT_POPULATED"
    assert checked["ward_current_source_certificate"]["status"] == "MISSING_SOURCE_CERTIFICATE"
    # statuses must mirror the live artifacts, not stored strings
    live_base = json.loads(packet_mod.BASE_MEASURE.read_text(encoding="utf-8"))
    assert checked["qcd_base_measure"]["status"] == live_base["status"]


def test_gauge_invariance_witness(packet):
    gi = packet["machine_witnesses"]["gauge_and_ward"]["gauge_invariance"]
    assert gi["passed"] is True
    assert gi["correlator_relative_defect"] < gi["tolerance"]
    assert gi["plaquette_absolute_defect"] < gi["tolerance"]


def test_ward_identity_witness(packet):
    ward = packet["machine_witnesses"]["gauge_and_ward"]["ward_identity"]
    assert ward["passed"] is True
    assert ward["relative_defect"] < ward["tolerance"]


def test_zv_free_field_anchor(packet):
    zv = packet["machine_witnesses"]["zv_free_field_anchor"]
    assert zv["passed"] is True
    assert zv["max_abs_deviation_from_one"] < zv["tolerance"]


def test_demonstrator_measure_positive_with_psd_covariance(packet):
    demo = packet["machine_witnesses"]["demonstrator_measure"]
    assert demo["checks_passed"] is True
    assert demo["rho_had_or_measure"]["positivity_status"].startswith(
        "verified_nonnegative"
    )
    level = demo["finite_volume_levels"][0]["levels"][0]
    assert level["energy"] > 0.0
    assert level["weight"] > 0.0
    assert demo["ward_projected_residues"][0]["residue"] > 0.0
    cov = demo["covariance"]
    assert cov["positive_semidefinite"] is True
    matrix = np.array(cov["matrix"])
    assert matrix.shape == (cov["dimension"], cov["dimension"])
    assert np.allclose(matrix, matrix.T)


def test_demonstrator_uses_no_external_inputs_and_is_non_promoting(packet):
    demo = packet["machine_witnesses"]["demonstrator_measure"]
    assert demo["ensemble"]["external_inputs_used"] is False
    guards = demo["guards"]
    assert guards["promotion_allowed"] is False
    assert guards["production_schema_conformant"] is False
    assert guards["surrogate_hadron_artifact"] is False
    assert guards["target_anchored"] is False


def test_acceptance_gate_accepts_conformant_and_rejects_all_controls(packet):
    gate = packet["machine_witnesses"]["acceptance_gate"]
    assert gate["passed"] is True
    assert gate["conformant_payload_accepted"] is True
    controls = gate["negative_controls"]
    assert len(controls) >= 9
    control_ids = {c["control_id"] for c in controls}
    assert "measured_hvp_comparison_endpoint" in control_ids
    assert "forbidden_target_leak_ee_to_hadrons" in control_ids
    assert "surrogate_hadron_artifact" in control_ids
    assert "quenched_branch" in control_ids
    assert "negative_residue" in control_ids
    for control in controls:
        assert control["rejected"] is True, control


def test_demonstrator_residue_normalization_chain(packet):
    """w_phys = Z_V^2 (2 kappa)^2 * 24 pi^2 A / m, per conserved_vector.py:
    the hopping-to-physical factor is (2 kappa)^2 and the atom amplitude is
    A = w m / (24 pi^2) (vector_correlator.synthetic_atom_correlator)."""
    import math

    demo = packet["machine_witnesses"]["demonstrator_measure"]
    kappa = demo["ensemble"]["kappa"]
    ex = demo["extraction"]
    w_hop = 24.0 * math.pi**2 * ex["amplitude_hop"] / ex["am_vector"]
    expected = ex["z_v"] ** 2 * (2.0 * kappa) ** 2 * w_hop
    assert ex["weight_physical_units"] == pytest.approx(expected, rel=1e-12)


def test_gate_target_leak_fails_closed_directly():
    schema = json.loads(packet_mod.SCHEMA_PATH.read_text(encoding="utf-8"))
    payload = packet_mod.build_conformant_payload()
    payload["external_targets_used"] = ["EE_TO_HADRONS"]
    accepted, reason = packet_mod.accept_payload(payload, schema)
    assert accepted is False
    assert reason.startswith("TARGET_LEAK_DETECTED")


def test_claim_boundary_and_dependency_direction(packet):
    boundary = packet["claim_boundary"]
    assert any("production backend execution" in row for row in boundary["not_closed_here"])
    assert "physical promotion" in " ".join(boundary["not_closed_here"])
    dep = packet["dependency_note"]
    assert dep["upstream_of_425_not_downstream"] is True
    assert "consumed by #317" in dep["direction"]
    assert packet["forbidden_inputs"] == packet_mod.FORBIDDEN_TARGETS


def test_stored_packet_matches_fresh_emission(packet, stored):
    assert stored["issue"] == 317
    assert stored["artifact"] == packet["artifact"]
    assert stored["status"] == packet["status"]
    assert stored["accepted"] is True
    assert stored["issue_closure_condition"] == packet["issue_closure_condition"]
    assert stored["acceptance_criteria_status"] == packet["acceptance_criteria_status"]
    assert stored["theorem"] == packet["theorem"]
    assert stored["derivation_chain"] == packet["derivation_chain"]
    assert stored["claim_boundary"] == packet["claim_boundary"]
    gate_stored = stored["machine_witnesses"]["acceptance_gate"]
    gate_fresh = packet["machine_witnesses"]["acceptance_gate"]
    assert gate_stored["conformant_payload_accepted"] == gate_fresh["conformant_payload_accepted"]
    assert [c["control_id"] for c in gate_stored["negative_controls"]] == [
        c["control_id"] for c in gate_fresh["negative_controls"]
    ]
    # numeric witnesses agree to tolerance (identical seeds and inputs)
    demo_stored = stored["machine_witnesses"]["demonstrator_measure"]["extraction"]
    demo_fresh = packet["machine_witnesses"]["demonstrator_measure"]["extraction"]
    assert demo_stored["am_vector"] == pytest.approx(demo_fresh["am_vector"], rel=1e-10)
    assert demo_stored["weight_physical_units"] == pytest.approx(
        demo_fresh["weight_physical_units"], rel=1e-10
    )
