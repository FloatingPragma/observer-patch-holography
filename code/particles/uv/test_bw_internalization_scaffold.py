#!/usr/bin/env python3
"""Guard the UV/BW internalization scaffold."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "uv" / "derive_bw_internalization_scaffold.py"
OUTPUT = ROOT / "particles" / "runs" / "uv" / "bw_internalization_scaffold.json"


def test_bw_internalization_scaffold_contract() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(OUTPUT)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_uv_bw_internalization_scaffold"
    assert payload["status"] == "minimal_constructive_extension"
    assert payload["public_promotion_allowed"] is False
    assert payload["extension_kind"] == "scaling_limit_cap_pair_plus_ordered_cut_pair_rigidity"
    assert payload["solver_spec"]["output_certificate"]["typeI_required"] is False
    boundary = payload["public_status_boundary"]
    assert boundary["remaining_object"] == "canonical_scaling_cap_pair_realization_from_transported_cap_marginals"
    assert boundary["follow_on_object"] == "independent_bw_rigidity_on_realized_limit"
    assert boundary["dominant_pressure_point"] == "eventual_fixed_local_collar_spectral_floor_for_transported_marginals"
    assert boundary["candidate_extension_status"] == "constructive_prelimit_system_two_lower_emitted_witnesses_still_missing"
    assert (
        boundary["remaining_missing_emitted_witness_artifact"]
        == "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json"
    )
    assert boundary["single_live_missing_clause_artifact"] == (
        "code/particles/runs/uv/bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"
    )
    assert boundary["markov_side_status"] == "latent_from_epsilon_to_zero"
    assert boundary["faithfulness_side_status"] == "open"
    assert [entry["id"] for entry in boundary["actual_solver_missing_emitted_witnesses"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
    ]
    assert boundary["derived_remaining_input_witness"]["id"] == "vanishing_carried_collar_schedule_on_fixed_local_collars"
    assert [entry["id"] for entry in boundary["local_intermediate_witness_chain"]] == [
        "constructive_recovery_remainder_vanishing",
        "fixed_local_collar_exact_markov_modulus_vanishing",
        "fixed_local_collar_faithful_modular_defect_vanishing",
        "vanishing_carried_collar_schedule_on_fixed_local_collars",
    ]
    assert [entry["artifact"] for entry in boundary["local_intermediate_witness_chain"]] == [
        "code/particles/runs/uv/bw_fixed_local_collar_constructive_recovery_scaffold.json",
        "code/particles/runs/uv/bw_fixed_local_collar_exact_markov_modulus_scaffold.json",
        "code/particles/runs/uv/bw_fixed_local_collar_faithful_modular_defect_scaffold.json",
        "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json",
    ]
    assert "code/particles/uv/derive_bw_fixed_local_collar_constructive_recovery_scaffold.py" in boundary[
        "canonical_code_scaffolds"
    ]
    assert "code/particles/uv/derive_bw_fixed_local_collar_exact_markov_modulus_scaffold.py" in boundary[
        "canonical_code_scaffolds"
    ]
    assert "code/particles/uv/derive_bw_fixed_local_collar_eventual_spectral_floor_scaffold.py" in boundary[
        "canonical_code_scaffolds"
    ]
    assert "code/particles/uv/derive_bw_carried_collar_schedule_scaffold.py" in boundary["canonical_code_scaffolds"]
    assert "code/particles/runs/uv/bw_fixed_local_collar_constructive_recovery_scaffold.json" in boundary[
        "canonical_artifacts"
    ]
    assert "code/particles/runs/uv/bw_fixed_local_collar_eventual_spectral_floor_scaffold.json" in boundary[
        "canonical_artifacts"
    ]
    assert "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json" in boundary["canonical_artifacts"]
    rigidity = boundary["symbolic_ordered_cut_pair_rigidity_test"]
    assert rigidity["status"] == "pass"
    assert rigidity["solution_dimension"] == 1
    assert rigidity["surviving_generator_half_line"] == "2*u"
