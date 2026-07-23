"""Tests for the Clebsch register pairing theorem and frozen weight-set scan."""

from __future__ import annotations

import copy

import pytest

import derive_clebsch_register_pairing_selection as lane


@pytest.fixture(scope="module")
def documents():
    return {
        "receipt": lane._load_json(lane.MATTER_RECEIPT_JSON),
        "manifest": lane._load_json(lane.MATTER_MANIFEST_JSON),
    }


def _build(documents):
    return lane.build_artifact(
        documents["receipt"], documents["manifest"], {"matter_receipt": "t", "matter_manifest": "t"}
    )


@pytest.fixture(scope="module")
def result(documents):
    return _build(documents)


def test_selection_lands_on_live_corpus(result):
    assert result["status"] == "pairing_and_weight_set_selected_order_open"
    assert result["proof_status"] == "closed_conditional_selection_theorem"
    assert result["promotion_allowed"] is False


def test_pairing_theorem_checks(result):
    pairing = result["pairing_theorem"]
    assert pairing["passed"] is True
    assert pairing["checks"]["down_couples_through_Sbar"] is True
    assert pairing["checks"]["lepton_couples_through_Sbar"] is True
    assert pairing["checks"]["down_through_up_scalar_has_no_line"] is True
    assert "#314" in pairing["conditional_on"]


def test_weight_generators_read_from_certified_spectrum(result):
    generators = result["weight_generators"]
    assert generators["N_c"] == 3
    assert generators["alphabet"] == ["1/3", "1", "3"]
    assert all(generators["checks"].values())


def test_scan_enumerates_full_space_and_selects_uniquely(result):
    scan = result["weight_set_scan"]
    assert scan["candidate_count"] == 27
    assert scan["survivors_after_F1"]["count"] == 7
    assert scan["survivors_after_F2"]["count"] == 6
    assert scan["unique_unordered_survivor"] is True
    assert scan["surviving_weight_set"] == ["1/3", "1", "3"]
    assert scan["order_assignment"]["status"] == "open"
    assert scan["order_assignment"]["remaining_premise"] == "GENERATION_REGISTER_ORDER"


def test_constraints_frozen_before_evaluation(result):
    scan = result["weight_set_scan"]
    assert scan["constraint_order"] == [
        "F1_measure_balance",
        "F2_register_faithfulness",
    ]
    assert set(scan["frozen_constraints"]) == set(lane.FROZEN_CONSTRAINTS)


def test_solve_path_reads_no_measured_values(result):
    discipline = result["solve_path_discipline"]
    assert discipline["reads_measured_masses"] is False
    assert discipline["reads_measured_ratios"] is False
    assert discipline["reads_down_type_lane_outputs"] is False
    assert "PDG" in discipline["forbidden_inputs"]


def test_premise_bookkeeping(result):
    bookkeeping = result["premise_bookkeeping"]
    assert bookkeeping["replaced_premise"] == "CLEBSCH_REGISTER_SELECTION_THEOREM"
    assert bookkeeping["remaining_open_parts"] == ["GENERATION_REGISTER_ORDER"]
    assert any("MCPR" in item for item in bookkeeping["unchanged_conditions"])


def test_forbidden_channel_regrowth_downgrades_pairing(documents):
    doctored = copy.deepcopy(documents)
    doctored["receipt"]["yukawa_sector"]["forbidden_channel_control"][
        "invariant_dimension"
    ] = 1
    result = _build(doctored)
    assert result["pairing_theorem"]["passed"] is False
    assert result["status"] == "selection_failed_premise_unchanged"
    assert "CLEBSCH_PAIRING" in result["premise_bookkeeping"]["remaining_open_parts"]


def test_doctored_charge_spectrum_fails_closed(documents):
    doctored = copy.deepcopy(documents)
    doctored["receipt"]["realized_package"]["charge_spectrum"]["-2/3"] = 4
    with pytest.raises(SystemExit, match="color multiplicity"):
        _build(doctored)
