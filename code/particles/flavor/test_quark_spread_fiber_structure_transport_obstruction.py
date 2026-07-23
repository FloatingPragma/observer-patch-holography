"""Tests for the 2026-07 structure transport into the quark spread fiber."""

from __future__ import annotations

import copy

import pytest

import derive_quark_spread_fiber_structure_transport_obstruction as lane


@pytest.fixture(scope="module")
def documents():
    return {
        "base": lane._load_json(lane.BASE_OBSTRUCTION_JSON),
        "axiom": lane._load_json(lane.AXIOM_LEVEL_JSON),
        "matter_manifest": lane._load_json(lane.MATTER_MANIFEST_JSON),
        "matter_receipt": lane._load_json(lane.MATTER_RECEIPT_JSON),
        "port_manifest": lane._load_json(lane.PORT_MANIFEST_JSON),
        "port_receipt": lane._load_json(lane.PORT_RECEIPT_JSON),
        "candidates": lane._load_json(lane.CANDIDATES_JSON),
    }


def _build(documents):
    return lane.build_artifact(
        documents["base"],
        documents["axiom"],
        documents["matter_manifest"],
        documents["matter_receipt"],
        documents["port_manifest"],
        documents["port_receipt"],
        documents["candidates"],
        {key: "test" for key in documents},
    )


@pytest.fixture(scope="module")
def result(documents):
    return _build(documents)


def test_fork_is_fiber_survives_on_live_corpus(result):
    assert result["fork"] == "ii_fiber_survives"
    assert result["fiber_cut_detected"] is False
    assert result["failed_checks"] == []
    assert result["proof_status"] == "closed_exact_current_corpus_obstruction"
    assert all(result["checks"].values())


def test_blindness_scan_is_not_vacuous_on_matter_documents(result):
    scans = {scan["document"]: scan for scan in result["coefficient_blindness_scans"]}
    assert scans["super_tannakian_matter_manifest"]["matched_leaf_count"] > 0
    assert scans["super_tannakian_matter_receipt"]["matched_leaf_count"] > 0
    assert result["coefficient_blindness_violations"] == []


def test_invariant_line_transport(result):
    lines = result["invariant_line_transport"]
    assert lines["hadronic_lines_one_dimensional"] is True
    assert lines["forbidden_channel_carries_no_line"] is True


def test_no_promotion_under_either_fork(result):
    assert result["public_promotion_allowed"] is False
    assert result["numeric_quark_rows_allowed"] is False
    assert result["source_only_sigma_emitted"] is False
    policy = result["public_prediction_policy"]
    assert policy["running_quark_numeric_rows"] == "withheld"
    assert policy["forward_yukawa_numeric_promotion"] == "withheld"


def test_issue_591_disjunct_recorded(result):
    boundary = result["issue_591_claim_boundary"]
    assert boundary["disjunct_discharged"] == (
        "prove_non_identifiable_at_2026_07_corpus_level"
    )
    assert "#569 family attachment" in boundary["gating_parents"]
    assert boundary["carrier_bullets_remaining_open"]


def test_dependency_audit_has_no_target_leak(result):
    audit = result["dependency_audit"]
    assert audit["no_target_leak"] is True
    assert audit["emits_target_sigma_values"] is False
    assert "PDG_API_quark_rows" in audit["forbidden_ancestors"]


def test_injected_yukawa_coefficient_flips_the_fork(documents):
    doctored = copy.deepcopy(documents)
    doctored["matter_receipt"]["yukawa_sector"]["channels"][0][
        "yukawa_coefficient"
    ] = 1.27
    result = _build(doctored)
    assert result["fork"] == "i_structure_cuts_fiber"
    assert result["fiber_cut_detected"] is True
    assert "matter_documents_coefficient_blind" in result["failed_checks"]
    assert (
        result["proof_status"]
        == "structure_cut_detected_manual_rederivation_required"
    )
    assert result["numeric_quark_rows_allowed"] is False


def test_injected_numeric_string_flips_the_fork(documents):
    doctored = copy.deepcopy(documents)
    doctored["matter_manifest"]["exterior_matter_contract"][
        "yukawa_normalization"
    ] = "1.27e0"
    result = _build(doctored)
    assert result["fork"] == "i_structure_cuts_fiber"
    assert "matter_documents_coefficient_blind" in result["failed_checks"]


def test_widened_invariant_line_flips_the_fork(documents):
    doctored = copy.deepcopy(documents)
    doctored["matter_receipt"]["yukawa_sector"]["channels"][1][
        "invariant_dimension"
    ] = 2
    result = _build(doctored)
    assert result["fork"] == "i_structure_cuts_fiber"
    assert "yukawa_invariant_lines_one_dimensional" in result["failed_checks"]


def test_unblind_selector_menu_flips_the_fork(documents):
    doctored = copy.deepcopy(documents)
    doctored["candidates"]["blindness"]["reads_measured_masses"] = True
    result = _build(doctored)
    assert result["fork"] == "i_structure_cuts_fiber"
    assert "selector_menu_exhausted_without_selection" in result["failed_checks"]
