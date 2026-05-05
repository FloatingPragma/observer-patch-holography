#!/usr/bin/env python3
"""Smoke tests for the derivation-chain closure matrix."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_derivation_chain_closure_matrix import build_payload  # noqa: E402


def test_derivation_chain_closure_matrix_keeps_stage_gates_explicit() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_particle_derivation_chain_closure_matrix"
    assert payload["status"] == "executable_nonhadron_chain_matrix_emitted"
    assert payload["closure_summary"]["all_derivation_chains_claimed_closed"] is False
    assert payload["closure_summary"]["hardware_gated_chains"] == ["hadrons"]
    assert payload["worker_policy"]["chrome_pro_workers_needed_now"] is False
    assert payload["provenance_status"] == "provenance_ledger_emitted_convention_sensitivity_contract_open"
    rows = {row["chain"]: row for row in payload["rows"]}
    assert set(rows) == {
        "p_closure_root",
        "structural_massless_bosons",
        "electroweak_wz",
        "higgs_top_declared_surface",
        "charged_leptons",
        "selected_class_quarks",
        "neutrino_absolute_attachment",
        "hadrons",
    }
    assert rows["p_closure_root"]["promotable"] is False
    assert rows["p_closure_root"]["open_gates"] == [32, 223, 224]
    assert rows["structural_massless_bosons"]["promotable"] is True
    assert rows["charged_leptons"]["status"] == "current_family_witness_only_end_to_end_nonclosure_theorem"
    assert rows["charged_leptons"]["promotable"] is False
    assert rows["hadrons"]["status"] == "closed_out_of_scope_computationally_blocked"
    assert rows["hadrons"]["open_gates"] == []
    assert rows["hadrons"]["closed_issue_refs"] == [153, 157]
    assert "p_closure_root" in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "hadrons" not in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert payload["closure_summary"]["closed_out_of_scope_chains"] == ["hadrons"]
    assert payload["particle_five_gates"]["153"]["state"] == "closed_out_of_scope_computationally_blocked"
    assert payload["particle_five_gates"]["224"]["state"] == "open_waiting_certified_root"
