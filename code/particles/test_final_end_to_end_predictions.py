#!/usr/bin/env python3
"""Smoke tests for the final current end-to-end prediction bundle."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_final_end_to_end_predictions import build_payload  # noqa: E402


def test_final_end_to_end_predictions_include_particle_five_gates_and_values() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_final_current_end_to_end_particle_predictions"
    assert payload["p_closure"]["may_feed_live_particle_predictions"] is False
    assert payload["hadron_policy"]["predictions_emitted"] is False
    assert payload["hadron_policy"]["github_issues"] == [153, 157]
    gates = {gate["issue"]: gate for gate in payload["particle_five_issue_gates"]}
    assert set(gates) == {32, 153, 207, 223, 224}
    assert gates[153]["state"] == "closed_out_of_scope_computationally_blocked"
    assert gates[224]["state"] == "open_waiting_certified_root"
    predictions = {entry["particle_id"]: entry for entry in payload["predictions"]}
    assert predictions["photon"]["value"] == 0.0
    assert predictions["w_boson"]["value"] == 80.377
    assert predictions["higgs"]["value"] == 125.1995304097179
    assert predictions["electron"]["value"] == 0.0005109989499999994
    assert predictions["top_quark"]["value"] == 172.35235532883115
    assert predictions["electron_neutrino"]["unit"] == "eV"
    assert predictions["electron_neutrino"]["value"] == 0.017454720257976796
    assert payload["direct_top_auxiliary_comparison"]["bridge_status"] == (
        "constructive_conversion_contract_emitted_not_direct_top_theorem"
    )
