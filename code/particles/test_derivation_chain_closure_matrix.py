#!/usr/bin/env python3
"""Smoke tests for the scientific derivation-chain classification matrix."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_derivation_chain_closure_matrix import build_payload  # noqa: E402


def test_derivation_chain_matrix_contains_scientific_boundaries_only() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_particle_derivation_chain_classification_matrix"
    assert payload["classification"] == (
        "scientific_chain_boundaries_and_required_receipts"
    )
    assert payload["classification_summary"]["promotable_chains"] == []
    assert payload["classification_summary"]["source_backend_absent_chains"] == [
        "hadrons"
    ]
    assert payload["provenance"] == {
        "artifact": "oph_blind_prediction_provenance_ledger",
        "promotion_allowed": False,
    }
    assert "pipeline_status" not in payload["source_artifacts"]
    assert "worker_policy" not in payload
    assert "particle_five_gates" not in payload

    rows = {row["chain"]: row for row in payload["rows"]}
    assert set(rows) == {
        "p_closure_root",
        "conditional_classical_carrier_modes",
        "electroweak_massive_bosons",
        "hierarchy_naturality_bridge",
        "higgs_top_declared_surface",
        "charged_leptons",
        "selected_class_quarks",
        "neutrino_absolute_attachment",
        "hadrons",
    }
    for row in rows.values():
        assert "closed_issue_refs" not in row
        assert "open_gates" not in row
        assert "next_artifact" not in row
        assert isinstance(row["required_receipts"], list)
        assert isinstance(row["evidence_artifacts"], list)

    root = rows["p_closure_root"]
    assert root["promotable"] is False
    assert root["status"] == "candidate_nonpromoting_root"
    assert "source_emitted_same_scheme_Ward_projected_R_Q" in root[
        "required_receipts"
    ]

    carrier = rows["conditional_classical_carrier_modes"]
    assert carrier["status"] == (
        "conditional_classical_modes__quantum_particle_receipts_absent"
    )
    assert carrier["outputs"]["photon"]["hard_quadratic_mass_parameter_squared"] == 0
    assert carrier["outputs"]["photon"]["quantum_particle_gate"] == "not_passed"
    assert carrier["required_receipts"]

    hierarchy = rows["hierarchy_naturality_bridge"]
    assert hierarchy["outputs"]["epsilon_H"] == "0"
    assert hierarchy["full_theorem_grade_resonance_promoted"] is False
    assert hierarchy["required_receipts"]

    higgs = rows["higgs_top_declared_surface"]
    assert higgs["status"] == "conditional_declared_surface_higgs_top_candidate"
    assert higgs["outputs"]["higgs"] == 125.1995304097179
    assert higgs["required_receipts"] == ["EWTargetFreeRepairValueLaw_D10"]

    charged = rows["charged_leptons"]
    assert charged["status"] == "corpus_limited_charged_end_to_end_no_go"
    assert charged["promotable"] is False

    quarks = rows["selected_class_quarks"]
    assert quarks["status"] == (
        "quark_source_nonidentifiability__numeric_rows_withheld"
    )
    assert "QUARK_OPERATIONAL_SCHEME_AND_SCALE_SECTION" in quarks[
        "required_receipts"
    ]

    neutrinos = rows["neutrino_absolute_attachment"]
    assert neutrinos["status"] == "rejected_target_informed_weighted_cycle_candidate"
    assert "pre_reference_hash_lock" in neutrinos["required_receipts"]

    hadrons = rows["hadrons"]
    assert hadrons["status"] == "source_backend_absent__empirical_output_class_separate"
    assert hadrons["required_receipts"] == ["source_only_hadron_backend"]
