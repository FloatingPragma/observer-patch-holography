#!/usr/bin/env python3
"""Tests for the end-to-end exact-PDG quark theorem artifact."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_exact_pdg_end_to_end_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
CHAIN_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"


def test_quark_exact_pdg_end_to_end_theorem_closes() -> None:
    chain = json.loads(CHAIN_JSON.read_text(encoding="utf-8"))
    payload = build_artifact(chain)

    assert payload["proof_status"] == "closed_current_family_exact_pdg_end_to_end_theorem"
    assert payload["target_name"] == (
        "exact_mixed_convention_quark_target_packet_on_declared_current_family_transport_frame"
    )
    assert payload["supporting_chain_artifact"] == "oph_quark_current_family_end_to_end_exact_pdg_derivation_chain"
    assert payload["minimal_exact_blocker_set"] == []
    assert payload["strengthening_above_target"]["status"] == "separate_question"
    masses = payload["exact_running_values_gev"]
    assert abs(float(masses["t"]) - 172.1) < 1.0e-12
    assert abs(float(masses["u"]) - 0.0021600000000000005) < 1.0e-12
