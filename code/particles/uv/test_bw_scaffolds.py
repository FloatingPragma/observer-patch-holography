#!/usr/bin/env python3
"""Guard the split UV/BW scaffolds."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONSTRUCTIVE_RECOVERY = ROOT / "particles" / "uv" / "derive_bw_fixed_local_collar_constructive_recovery_scaffold.py"
EXACT_MARKOV = ROOT / "particles" / "uv" / "derive_bw_fixed_local_collar_exact_markov_modulus_scaffold.py"
SPECTRAL_FLOOR = ROOT / "particles" / "uv" / "derive_bw_fixed_local_collar_eventual_spectral_floor_scaffold.py"
MODULAR_DEFECT = ROOT / "particles" / "uv" / "derive_bw_fixed_local_collar_faithful_modular_defect_scaffold.py"
SCHEDULE = ROOT / "particles" / "uv" / "derive_bw_carried_collar_schedule_scaffold.py"
EXTRACTION = ROOT / "particles" / "uv" / "derive_bw_scaling_limit_cap_pair_extraction_scaffold.py"
RIGIDITY = ROOT / "particles" / "uv" / "derive_bw_ordered_cut_pair_rigidity_scaffold.py"
RUNS = ROOT / "particles" / "runs" / "uv"


def _run(script: Path) -> dict:
    output = RUNS / script.name.replace(".py", ".json")
    completed = subprocess.run(
        [sys.executable, str(script), "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "saved:" in completed.stdout
    return json.loads(output.read_text(encoding="utf-8"))


def test_scaling_limit_cap_pair_extraction_scaffold() -> None:
    payload = _run(EXTRACTION)
    assert payload["artifact"] == "oph_bw_scaling_limit_cap_pair_extraction_scaffold"
    assert payload["status"] == "constructive_prelimit_system_two_lower_emitted_witnesses_still_missing"
    assert payload["exact_missing_object"] == "scaling_limit_cap_pair_extraction"
    assert payload["precise_missing_object_name"] == "canonical_scaling_cap_pair_realization_from_transported_cap_marginals"
    assert payload["theorem_contract_name"] == "conditional_scaling_limit_cap_pair_extraction_theorem"
    assert "projectively_compatible_transported_cap_marginal_family" in payload["fills_contract_witnesses"]
    assert payload["remaining_missing_emitted_witness"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    assert payload["remaining_missing_emitted_witness_artifact"].endswith("bw_carried_collar_schedule_scaffold.json")
    assert payload["remaining_missing_emitted_witness_formula"].startswith("eta_{n,m,delta} = r_FR")
    assert payload["smaller_remaining_raw_datum"] == "fixed_local_collar_markov_faithfulness_datum"
    assert payload["smaller_remaining_raw_datum_artifact"].endswith("bw_fixed_local_collar_markov_faithfulness_datum.json")
    assert payload["single_live_missing_clause_artifact"].endswith(
        "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"
    )
    assert payload["markov_side_status"] == "latent_from_epsilon_to_zero"
    assert payload["faithfulness_side_status"] == "open"
    assert [entry["id"] for entry in payload["intermediate_witness_chain"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_exact_markov_modulus_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
        "vanishing_carried_collar_schedule_on_fixed_local_collars",
    ]
    assert [entry["id"] for entry in payload["schedule_term_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    assert [entry["id"] for entry in payload["actual_solver_missing_emitted_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    assert payload["derived_remaining_input_witness"]["id"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    assert payload["missing_input_witnesses"] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    assert payload["follow_on_object"]["id"] == "ordered_null_cut_pair_rigidity"
    gate = payload["remaining_witness_honesty_gate"]
    assert gate["status"] == "open"
    assert gate["insufficient_on_their_own"][0]["artifact"].endswith(
        "bw_fixed_local_collar_constructive_recovery_scaffold.json"
    )
    assert gate["insufficient_on_their_own"][3]["artifact"].endswith("bw_realized_transported_cap_local_system.json")
    ledger = payload["remaining_witness_obligation_ledger"]
    assert ledger[-1]["id"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"


def test_fixed_local_collar_constructive_recovery_scaffold() -> None:
    payload = _run(CONSTRUCTIVE_RECOVERY)
    assert payload["artifact"] == "oph_bw_fixed_local_collar_constructive_recovery_scaffold"
    assert payload["exact_missing_object"] == "constructive_recovery_remainder_vanishing"
    assert payload["contract"]["must_emit"] == "r_FR(epsilon_{n,m,delta}) -> 0"
    assert payload["feeds_parent_schedule"]["artifact"].endswith("bw_carried_collar_schedule_scaffold.json")
    assert payload["feeds_parent_schedule"]["other_term_still_needed_artifact"].endswith(
        "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
    )
    assert [entry["id"] for entry in payload["joint_schedule_term_frontier"]["missing_emitted_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]


def test_fixed_local_collar_exact_markov_modulus_scaffold() -> None:
    payload = _run(EXACT_MARKOV)
    assert payload["artifact"] == "oph_bw_fixed_local_collar_exact_markov_modulus_scaffold"
    assert payload["exact_missing_object"] == "fixed_local_collar_exact_markov_modulus_vanishing"
    assert payload["parent_raw_datum"] == "fixed_local_collar_markov_faithfulness_datum"
    assert payload["contract"]["must_emit"] == "delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0"
    assert payload["feeds_follow_on_modular_defect"]["artifact"].endswith(
        "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
    )


def test_fixed_local_collar_eventual_spectral_floor_scaffold() -> None:
    payload = _run(SPECTRAL_FLOOR)
    assert payload["artifact"] == "oph_bw_fixed_local_collar_eventual_spectral_floor_scaffold"
    assert payload["status"] == "minimal_faithfulness_side_extension"
    assert payload["exact_missing_object"] == "eventual_fixed_local_collar_spectral_floor_for_transported_marginals"
    assert "lambda_bar_{m,delta}" in payload["contract"]["must_emit"]
    assert payload["unlocks"] == [
        "fixed_local_collar_faithful_modular_defect_vanishing",
        "vanishing_carried_collar_schedule_on_fixed_local_collars",
        "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
    ]


def test_fixed_local_collar_faithful_modular_defect_scaffold() -> None:
    payload = _run(MODULAR_DEFECT)
    assert payload["artifact"] == "oph_bw_fixed_local_collar_faithful_modular_defect_scaffold"
    assert payload["exact_missing_object"] == "fixed_local_collar_faithful_modular_defect_vanishing"
    assert payload["smaller_comparison_witness"] == "fixed_local_collar_exact_markov_modulus_vanishing"
    assert payload["smaller_comparison_witness_artifact"].endswith(
        "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
    )
    assert payload["blocking_side_condition_artifact"].endswith(
        "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"
    )
    assert payload["contract"]["must_emit"].startswith("4 * lambda_{*,n,m,delta}^{-1}")
    assert payload["position_inside_carried_schedule"]["other_term_still_needed"] == "r_FR(epsilon_{n,m,delta}) -> 0"
    assert [entry["id"] for entry in payload["joint_schedule_term_frontier"]["missing_emitted_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]


def test_carried_collar_schedule_scaffold() -> None:
    payload = _run(SCHEDULE)
    assert payload["artifact"] == "oph_bw_carried_collar_schedule_scaffold"
    assert payload["exact_missing_object"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    assert payload["theorem_object"] == "canonical_scaling_cap_pair_realization_from_transported_cap_marginals"
    assert payload["derived_missing_witness"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    assert payload["raw_input_frontier"] == "fixed_local_collar_markov_faithfulness_datum"
    assert payload["intermediate_substep"] == "fixed_local_collar_exact_markov_modulus_vanishing"
    assert payload["smaller_raw_datum"] == "fixed_local_collar_markov_faithfulness_datum"
    assert payload["smaller_raw_datum_artifact"].endswith("bw_fixed_local_collar_markov_faithfulness_datum.json")
    assert payload["single_live_missing_clause_artifact"].endswith(
        "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"
    )
    assert payload["schedule_contract"]["formula"].startswith("eta_{n,m,delta} = r_FR")
    assert payload["decomposed_error_terms"]["faithful_modular_defect_remainder"]["artifact"].endswith(
        "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
    )
    recovery = payload["decomposed_error_terms"]["constructive_recovery_remainder"]
    assert recovery["id"] == "constructive_recovery_remainder_vanishing"
    assert recovery["artifact"].endswith("bw_fixed_local_collar_constructive_recovery_scaffold.json")
    assert payload["reduction_from_raw_datum"]["status_on_fill"] == "carried_collar_schedule_closed"
    assert payload["reduction_from_raw_datum"]["constructive_recovery_witness_artifact"].endswith(
        "bw_fixed_local_collar_constructive_recovery_scaffold.json"
    )
    assert [entry["id"] for entry in payload["schedule_term_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    assert payload["termwise_closure_frontier"]["derived_parent_witness"]["id"] == (
        "vanishing_carried_collar_schedule_on_fixed_local_collars"
    )
    assert payload["honesty_gate"]["promotion_rule"].startswith("No theorem promotion is honest")
    assert payload["obligation_ledger"][3]["formula"] == "r_FR(epsilon_{n,m,delta}) -> 0"


def test_ordered_cut_pair_rigidity_scaffold() -> None:
    payload = _run(RIGIDITY)
    assert payload["artifact"] == "oph_bw_ordered_cut_pair_rigidity_scaffold"
    assert payload["exact_missing_object"] == "ordered_null_cut_pair_rigidity"
    assert payload["symbolic_disk_halfline_witness"]["solution_dimension"] == 1
