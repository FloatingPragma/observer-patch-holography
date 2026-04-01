#!/usr/bin/env python3
"""Guard the UV fixed-local-collar raw datum artifact."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "uv" / "derive_bw_fixed_local_collar_markov_faithfulness_datum.py"
OUTPUT = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"


def test_bw_fixed_local_collar_markov_faithfulness_datum() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(OUTPUT)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "saved:" in completed.stdout
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_bw_fixed_local_collar_markov_faithfulness_datum"
    assert payload["exact_missing_object"] == "fixed_local_collar_markov_faithfulness_datum"
    assert payload["parent_missing_witness"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    must_emit = payload["contract"]["must_emit"]
    assert must_emit[0].startswith("epsilon_{n,m,delta} = I(")
    assert "lambda_bar_{m,delta}" in must_emit[1]
    assert payload["single_live_missing_clause_artifact"].endswith(
        "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
    )
    closure_lemma = payload["single_live_missing_clause_closure_lemma"]
    assert closure_lemma["id"] == "exact_markov_reference_eventual_common_floor_transfer"
    assert closure_lemma["requires_exact_markov_artifact"].endswith(
        "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
    )
    assert closure_lemma["requires_spectral_floor_artifact"].endswith(
        "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
    )
    assert payload["markov_side_status"] == "latent_from_epsilon_to_zero"
    assert payload["faithfulness_side_status"] == "open"
    assert payload["implies_schedule"]["artifact"].endswith("bw_carried_collar_schedule_scaffold.json")
    assert payload["implies_schedule"]["formula"].startswith("eta_{n,m,delta} = r_FR")
    assert [entry["id"] for entry in payload["schedule_term_frontier"]["missing_emitted_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    decomposition = payload["decomposes_into"]
    assert decomposition["constructive_recovery_witness"]["artifact"].endswith(
        "bw_fixed_local_collar_constructive_recovery_scaffold.json"
    )
    assert decomposition["derived_exact_markov_comparison_witness"]["artifact"].endswith(
        "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
    )
    assert decomposition["faithful_modular_defect_witness"]["artifact"].endswith(
        "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
    )
    assert "projectively_compatible_transported_cap_marginal_family" in payload["already_packaged_below_this_datum"]
    ledger = payload["obligation_ledger"]
    assert [entry["id"] for entry in ledger] == [
        "collarwise_markov_input",
        "collarwise_faithfulness_input",
        "fixed_local_collar_exact_markov_modulus_vanishing",
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
        "vanishing_carried_collar_schedule_on_fixed_local_collars",
    ]
    honesty_gate = payload["honesty_gate"]
    assert honesty_gate["status"] == "open"
    assert honesty_gate["closure_artifact"].endswith("bw_carried_collar_schedule_scaffold.json")
    assert honesty_gate["insufficient_on_their_own"][0]["artifact"].endswith(
        "bw_fixed_local_collar_constructive_recovery_scaffold.json"
    )
    assert honesty_gate["insufficient_on_their_own"][2]["artifact"].endswith(
        "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
    )
    assert honesty_gate["insufficient_on_their_own"][3]["artifact"].endswith(
        "bw_realized_transported_cap_local_system.json"
    )
