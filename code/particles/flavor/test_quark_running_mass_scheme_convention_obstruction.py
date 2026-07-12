#!/usr/bin/env python3
"""Guard the structural scheme-convention obstruction for issues #381/#382."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_running_mass_scheme_convention_obstruction.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_running_mass_scheme_convention_obstruction.json"


def _walk_keys(value: Any) -> set[str]:
    if isinstance(value, dict):
        keys = set(value)
        for child in value.values():
            keys.update(_walk_keys(child))
        return keys
    if isinstance(value, list):
        keys: set[str] = set()
        for child in value:
            keys.update(_walk_keys(child))
        return keys
    return set()


def test_quark_running_mass_scheme_convention_obstruction_closes_381_382_only_as_no_go() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_running_mass_scheme_convention_obstruction"
    assert payload["github_issues"] == [381, 382]
    assert payload["closure_kind"] == "theorem_grade_sharper_obstruction"
    assert payload["formal_quotient_obstruction"]["target_action"] == "nontrivial"
    assert "epsilon" in payload["same_source_counterfamily"]["free_parameters"]
    assert "kappa" in payload["same_source_counterfamily"]["free_parameters"]

    dynamics = payload["physical_rg_dynamics_nondefinability"]
    assert dynamics["source_RGI_mass_vector_emitted"] is False
    assert dynamics["source_six_flavor_Yukawa_trajectory_emitted"] is False
    assert dynamics["source_running_v_trajectory_emitted"] is False
    assert dynamics["additional_axiom_used"] is False
    assert "RGI quark mass vector" in dynamics["physical_outputs_changed"]

    chart = payload["external_chart_boundary"]
    assert chart["choosing_MSbar_is_a_physical_axiom"] is False
    assert chart["choosing_a_reporting_scale_is_a_physical_axiom"] is False

    conventions = payload["declared_comparison_conventions"]
    assert conventions["status"] == "external_scheme_metadata_not_OPH_source_output"
    assert conventions["mass_reference_values_consumed"] is False
    assert conventions["light_quarks"] == {
        "particles": ["u", "d", "s"],
        "renormalization_scheme": "MSbar",
        "scale": "mu = 2 GeV",
        "provenance": "declared_external_comparison_convention",
        "source_selected": False,
    }
    assert conventions["heavy_running_quarks"]["particles"] == ["c", "b"]
    assert conventions["heavy_running_quarks"]["renormalization_scheme"] == "MSbar"
    assert conventions["heavy_running_quarks"]["scale"] == "mu = m_q(mu) (self-scale convention)"
    assert conventions["top"]["row_kind"] == "separate_extraction_codomain"
    assert conventions["top"]["included_in_running_mass_sextet"] is False
    assert conventions["top"]["related_issue"] == 383

    assert payload["row_partition"]["single_running_quark_sextet_claim_allowed"] is False
    matrix_audit = payload["stored_matrix_dimensional_audit"]
    assert matrix_audit["stored_entry_dimension"] == "GeV"
    assert matrix_audit["comparison_packet_is_mixed"] is True
    assert matrix_audit["common_renormalization_scale"] is None
    assert matrix_audit["running_higgs_vev_supplied"] is False
    assert matrix_audit["top_pole_to_MSbar_conversion_supplied"] is False
    assert matrix_audit["dimensionless_normalization_supplied"] is False
    assert matrix_audit["current_classification"] == "mixed_scheme_GeV_mass_texture_matrices"
    assert matrix_audit["certified_physical_yukawa_matrices"] is False
    assert matrix_audit["exact_forward_yukawa_claim_allowed"] is False
    assert matrix_audit["physical_dimensionless_relation"] == (
        "y_q(mu_common) = sqrt(2)*m_q^MSbar(mu_common)/v(mu_common)"
    )
    assert payload["reference_data_policy"]["input_artifacts"] == []
    assert payload["reference_data_policy"]["reference_or_target_mass_values_consumed"] is False
    assert payload["reference_data_policy"]["PDG_or_API_quark_rows_consumed"] is False
    assert payload["reference_data_policy"]["current_family_exact_targets_consumed"] is False
    assert payload["reference_data_policy"]["selected_class_exact_target_witnesses_consumed"] is False
    assert payload["reference_data_policy"]["no_target_leak_by_construction"] is True
    assert "mass_gev" not in _walk_keys(payload)

    issue_381 = payload["issue_acceptance"]["381"]
    issue_382 = payload["issue_acceptance"]["382"]
    assert issue_381["acceptance_met_as_sharper_obstruction"] is True
    assert issue_381["positive_source_only_scheme_map_emitted"] is False
    assert issue_382["acceptance_met_as_sharper_obstruction"] is True
    assert issue_382["positive_source_only_scheme_transport_emitted"] is False

    effect = payload["closure_effect"]
    assert effect["issues_381_and_382_closable_as_obstruction"] is True
    assert effect["source_only_running_mass_coordinate_emitted"] is False
    assert effect["numeric_quark_prediction_rows_unblocked"] is False
    assert effect["public_numeric_quark_rows_allowed"] is False
    assert effect["exact_forward_physical_yukawa_claim_unblocked"] is False
    assert "fold_the_top_extraction_coordinate_into_a_running_mass_sextet" in payload["forbidden_promotions"]
    assert "call_a_mixed_scheme_GeV_mass_texture_a_physical_Yukawa_matrix" in payload["forbidden_promotions"]
