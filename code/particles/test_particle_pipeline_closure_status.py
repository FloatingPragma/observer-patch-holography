#!/usr/bin/env python3
"""Smoke tests for the simplified particle pipeline closure status."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_particle_pipeline_closure_status import build_status  # noqa: E402


def test_particle_pipeline_closure_status_scope_locks_hadrons_and_workers() -> None:
    status = build_status()

    assert status["scope"]["hadrons_in_current_local_scope"] is False
    assert status["scope"]["chrome_workers_needed_now"] is False
    assert "GLORB/Echosahedron" in status["scope"]["hadron_scope_reason"]
    assert status["finalization_gates"]["obstruction_only_worker_result_allowed"] is True
    gates = {gate["issue"]: gate for gate in status["issue_gates"]}
    assert gates[153]["state"] == "closed_out_of_scope_computationally_blocked"
    assert gates[153]["closable_now"] is True
    assert gates[153]["requires_oph_hardware_backend"] is True
    assert gates[153]["closed_as_out_of_scope"] is True
    assert gates[153]["chrome_workers"] == "do_not_use_for_backend_execution"
    assert gates[157]["state"] == "closed_out_of_scope_computationally_blocked"
    assert gates[157]["closable_now"] is True
    assert gates[157]["requires_oph_hardware_backend"] is True
    assert gates[157]["closed_as_out_of_scope"] is True
    assert gates[157]["chrome_workers"] == "do_not_use_for_backend_execution"
    assert gates[223]["state"] == "closed_blocker_isolated_endpoint_package"
    assert gates[223]["closable_now"] is True
    assert gates[223]["successor_issue"] == 235
    assert gates[223]["chrome_workers"] == "not_needed_for_closed_package"
    assert gates[235]["state"] == "closed_blocker_isolated_source_residual_no_go"
    assert gates[235]["closable_now"] is True
    assert gates[235]["closed_as_first_missing_lemma_isolated"] is True
    assert gates[235]["promotion_allowed"] is False
    assert gates[235]["source_spectral_reduction"] == "source_spectral_reduction_theorem_emitted_measure_payload_absent"
    assert gates[235]["minimal_new_payload"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    assert gates[235]["chrome_workers"] == "not_needed_until_source_spectral_measure_payload_exists"
    companion = {branch["label"]: branch for branch in status["companion_status_branches"]}
    assert companion["Strong CP"]["state"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert gates[224]["state"] == "closed_canonical_guarded_trunk_adoption"
    assert gates[32]["state"] == "closed_declared_convention_contract"
    assert gates[32]["closable_now"] is True
    assert gates[207]["state"] == "closed_current_corpus_codomain_no_go"
    assert gates[207]["closable_now"] is True
    assert gates[201]["state"] == "closed_current_corpus_charged_end_to_end_no_go"
    assert gates[201]["closable_now"] is True
    assert gates[199]["state"] == "closed_current_corpus_global_classification_no_go"
    assert gates[199]["closable_now"] is True
    assert gates[234]["state"] == "closed_provenance_ledger_and_declared_sensitivity_taxonomy"
    assert gates[155]["state"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert gates[155]["closable_now"] is False
    companion = {branch["label"]: branch for branch in status["companion_status_branches"]}
    assert companion["Strong CP"]["state"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert gates[117]["closable_now"] is True
    assert gates[198]["closable_now"] is True
    assert status["latest_nonhadron_predictions"]["higgs"]["value"] == 125.1995304097179
    assert status["latest_nonhadron_predictions"]["higgs"]["unit"] == "GeV"
    assert status["latest_nonhadron_predictions"]["top_quark"]["value"] == 172.35235532883115
    assert status["latest_nonhadron_predictions"]["top_quark"]["unit"] == "GeV"
    assert status["latest_nonhadron_predictions"]["electron_neutrino"]["value"] == 0.017454720257976796
    assert status["latest_nonhadron_predictions"]["electron_neutrino"]["unit"] == "eV"
