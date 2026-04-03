#!/usr/bin/env python3
"""Validate the exact three-step quark closure contract artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
T1_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_t1_value_law.py"
PHYSICAL_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
SELECTED_SHEET_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_selected_sheet_closure.py"
EXACT_READOUT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
BACKREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_internal_backread_continuation_closure.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_lane_closure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_lane_closure_contract.json"


def test_quark_lane_closure_contract_records_three_exact_missing_theorems() -> None:
    for script in (
        T1_SCRIPT,
        PHYSICAL_BRANCH_SCRIPT,
        EXACT_READOUT_SCRIPT,
        SELECTED_SHEET_SCRIPT,
        BACKREAD_SCRIPT,
        SCRIPT,
    ):
        subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_lane_closure_contract"
    assert payload["base_theorem_emitted_package_artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert payload["proof_status"] == "maximal_theorem_emitted_package_identified_three_exact_missing_theorems_recorded_exact_sidecar_masses_recorded"
    assert payload["selected_local_sheet_status"]["sigma_id"] == "sigma_ref"
    assert payload["exact_sidecar_mass_surface"]["scope"] == "current_family_only"
    assert payload["exact_sidecar_mass_surface"]["exact_outputs_gev"]["u"] == 0.0021599999999999996
    assert payload["exact_sidecar_mass_surface"]["exact_outputs_gev"]["d"] == 0.004700000000000002
    assert payload["continuation_only_mass_sidecar"]["closed_mass_side_package"]["t1"] == 0.6715870378831591
    theorem_ids = [item["id"] for item in payload["exact_missing_theorems"]]
    assert theorem_ids == [
        "quark_d12_t1_value_law",
        "quark_physical_sigma_ud_lift",
        "quark_absolute_sector_readout_theorem",
    ]
    assert payload["exact_missing_theorems"][0]["must_emit"] == "t1"
    assert payload["exact_missing_theorems"][0]["mathematical_name"] == "Theta_ud^mass"
    assert payload["exact_missing_theorems"][0]["same_family_d12_hypercharge_law"]["equivalent_scalar_identities"] == [
        "log(c_d / c_u) = (6/5) * t1",
        "t1 = 5 * Delta_ud_overlap = (5/6) * log(c_d / c_u)",
    ]
    assert payload["exact_missing_theorems"][1]["must_emit"] == "{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM, ckm_invariants}"
    assert payload["exact_missing_theorems"][1]["mathematical_name"] == "Theta_ud^phys"
    assert payload["exact_missing_theorems"][1]["physical_carrier"]["name"] == "Sigma_ud^phys"
    assert payload["exact_missing_theorems"][2]["must_emit"] == ["g_u", "g_d"]
    assert payload["exact_missing_theorems"][2]["mathematical_name"] == "Theta_ud^abs"
    assert payload["exact_missing_theorems"][2]["affine_mean_split_law"]["A_ud"] == "1 / (2 * (1 + rho_ord - x2^2))"
    assert payload["exact_missing_theorems"][2]["current_exact_obstruction"]["certification_status"] == "placeholder_unpromotable"
