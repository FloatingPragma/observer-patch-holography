#!/usr/bin/env python3
"""Guard the disposable runtime status surface against public-surface drift."""

from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "compute_current_output_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("compute_current_output_table", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_runtime_surface_preserves_repaired_neutrino_rows_and_canonical_refs(tmp_path: pathlib.Path) -> None:
    module = _load_module()
    current_dir = module.build_runtime(tmp_path / "runtime", with_hadrons=False, verbose=False)

    payload = json.loads((current_dir / "results_status.json").read_text(encoding="utf-8"))
    exact_nonhadron = json.loads((current_dir / "exact_nonhadron_masses.json").read_text(encoding="utf-8"))
    exact_fit_surface = json.loads((current_dir / "exact_fits_only.json").read_text(encoding="utf-8"))
    active = payload["surface_state"]["active_local_public_candidates"]
    uv = payload["premise_boundaries"]["uv_bw_internalization"]
    markdown = (current_dir / "RESULTS_STATUS.md").read_text(encoding="utf-8")
    exact_entries = {entry["particle_id"]: entry for entry in exact_nonhadron["entries"]}
    exact_fit_entries = {entry["id"]: entry for entry in exact_fit_surface["entries"]}

    assert active["neutrino_repaired_branch"] is True
    assert payload["comparison_rows"]
    assert payload["inputs"]["hadron_profile"] == "suppressed"
    assert uv["prelimit_system_artifact"] == "code/particles/runs/uv/bw_realized_transported_cap_local_system.json"
    assert uv["remaining_missing_emitted_witness_artifact"] == (
        "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json"
    )
    assert uv["smaller_remaining_raw_datum_artifact"] == (
        "code/particles/runs/uv/bw_fixed_local_collar_markov_faithfulness_datum.json"
    )
    assert uv["neutrino_local_bridge_candidate_context"] == (
        "code/particles/runs/neutrino/neutrino_lambda_nu_bridge_candidate.json"
    )
    assert "## Neutrino Oscillation Comparison" in markdown
    assert exact_nonhadron["status"] == "exact_output_lane_closed_nonhadron_only"
    assert exact_fit_surface["artifact"] == "oph_exact_fits_only_surface"
    assert exact_entries["top_quark"]["mass_gev"] == pytest.approx(172.3523553288312, abs=1.0e-12)
    assert exact_entries["bottom_quark"]["mass_gev"] == pytest.approx(4.182999999999994, abs=1.0e-15)
    assert (
        exact_entries["bottom_quark"]["supporting_scope_closure_artifact"]
        == "code/particles/runs/flavor/quark_current_family_selected_sheet_closure.json"
    )
    assert (
        exact_entries["electron"]["supporting_scope_closure_artifact"]
        == "code/particles/runs/leptons/lepton_current_family_affine_anchor_theorem.json"
    )
    assert exact_fit_entries["charged_current_family_exact_witness"]["supporting_scope_closure_artifact"] == (
        "code/particles/runs/leptons/lepton_current_family_affine_anchor_theorem.json"
    )
    assert exact_fit_entries["quark_current_family_exact_witness"]["supporting_scope_closure_artifact"] == (
        "code/particles/runs/flavor/quark_current_family_selected_sheet_closure.json"
    )
    assert (current_dir / "P_derivation" / "runtime" / "p_closure_trunk_current.json").exists()
    assert (current_dir / "P_derivation" / "runtime" / "thomson_endpoint_contract_current.json").exists()
    assert (current_dir / "P_derivation" / "runtime" / "thomson_endpoint_package_current.json").exists()
    assert (current_dir / "P_derivation" / "runtime" / "rg_matching_threshold_contract_current.json").exists()
    assert (current_dir / "runs" / "status" / "particle_derivation_gap_ledger.json").exists()
    assert (current_dir / "runs" / "status" / "particle_pipeline_closure_status.json").exists()
    provenance = json.loads((current_dir / "runs" / "status" / "blind_prediction_provenance.json").read_text())
    assert provenance["status"] == "provenance_ledger_emitted_convention_sensitivity_contract_open"
    assert (current_dir / "BLIND_PREDICTION_PROVENANCE.md").exists()
    final_predictions = json.loads((current_dir / "runs" / "status" / "final_end_to_end_predictions.json").read_text())
    assert final_predictions["artifact"] == "oph_final_current_end_to_end_particle_predictions"
    assert (current_dir / "FINAL_END_TO_END_PREDICTIONS.md").exists()
    direct_top = json.loads((current_dir / "runs" / "calibration" / "direct_top_bridge_contract.json").read_text())
    assert direct_top["status"] == "hard_no_go_current_corpus_compare_only_direct_top_codomain"
    assert direct_top["promotion_allowed"] is False


def test_runtime_surface_with_hadrons_exercises_intermediate_hadron_lane(tmp_path: pathlib.Path) -> None:
    module = _load_module()
    current_dir = module.build_runtime(tmp_path / "runtime", with_hadrons=True, verbose=False)

    payload = json.loads((current_dir / "results_status.json").read_text(encoding="utf-8"))
    hadron_rows = [row for row in payload["rows"] if row["group"] == "Hadrons"]

    assert payload["inputs"]["with_hadrons"] is True
    assert payload["inputs"]["hadron_profile"] == "serious"
    assert payload["surface_state"]["active_local_public_candidates"]["hadrons_enabled"] is True
    assert hadron_rows
    assert {row["status"] for row in hadron_rows} == {"simulation_dependent"}
    assert (current_dir / "runs" / "status" / "status_table_forward_current.json").exists()
