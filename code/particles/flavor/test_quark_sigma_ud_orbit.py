#!/usr/bin/env python3
"""Validate the quark sigma_ud orbit frontier artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sigma_ud_orbit.py"


def test_sigma_ud_orbit_closes_to_reference_singleton_by_uniqueness_theorem() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_sigma_ud_orbit_") as tmpdir:
        out = pathlib.Path(tmpdir) / "orbit.json"
        subprocess.run([sys.executable, str(SCRIPT), "--output", str(out)], check=True, cwd=ROOT)
        payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["status"] == "same_label_left_handed_local_orbit_singleton_closed"
    assert payload["exact_missing_object"] is None
    assert payload["elements_origin"] == "reference_sheet_singleton_provider"
    assert payload["selector_status"] == "quark_relative_sheet_selector_closed_to_reference_singleton"
    assert payload["provider_frontier"]["status"] == "same_label_left_handed_local_orbit_singleton_closed"
    assert payload["provider_frontier"]["already_local_diagnostic_orbit_available"] is True
    assert payload["provider_frontier"]["transport_frame_diagnostic_orbit_available"] is True
    assert payload["provider_frontier"]["smallest_missing_runtime_object"] is None
    assert payload["next_exact_object_after_orbit_closure"] == "intrinsic_scale_law_D12"
    assert len(payload["elements"]) == 1
    element = payload["elements"][0]
    assert element["sigma_id"] == "sigma_ref"
    assert element["canonical_token"] == "D12::same_label_left::reference_sheet"
    assert element["coverage_status"] == "reference_sheet_representative_only"
    assert element["selection_proof"]["theorem_grade_select"] is True
    assert payload["selected_sigma"]["sigma_id"] == "sigma_ref"
    assert payload["selected_sigma"]["branch_key"] == ["D12", "sigma_ref"]
    assert payload["selection_gate"]["quark_relative_sheet_selector"]["sigma_id"] == "sigma_ref"
    assert payload["singleton_uniqueness_theorem"]["proof_status"] == "same_label_left_handed_local_orbit_singleton_closed"
    assert payload["provider_frontier"]["emitted_reference_sheet"]["canonical_token"] == element["canonical_token"]
    assert payload["provider_frontier"]["emitted_reference_sheet"]["branch_key"] == ["D12", "sigma_ref"]
    assert payload["provider_frontier"]["emitted_reference_sheet"]["family_branch_key"] == ["D12", None]
    assert element["branch_key"] == ["D12", "sigma_ref"]
    assert element["family_branch_key"] == ["D12", None]
    assert payload["debug_compare_shell_ranking"]["ranked"][0]["canonical_token"] == element["canonical_token"]
    assert payload["already_local_diagnostic_orbit"]["physical_reference_element"]["basis_u"] == "L"
    assert payload["already_local_diagnostic_orbit"]["physical_reference_element"]["basis_d"] == "L"
    assert payload["already_local_diagnostic_orbit"]["best_nonphysical_candidate"]["physical_admissible"] is False
    assert payload["diagnostic_transport_frame_orbit"]["self_overlap_symbol"] == "F0^dagger F1"
    assert payload["diagnostic_transport_frame_orbit"]["artifact"] == "code/particles/runs/flavor/quark_transport_frame_diagnostic_orbit.json"
    assert payload["diagnostic_transport_frame_orbit"]["debug_log_shell_loss"]["transport_frame_self_overlap"] < payload["diagnostic_transport_frame_orbit"]["debug_log_shell_loss"]["current_same_sheet"]
    assert "canonical_token" in payload["diagnostic_transport_frame_orbit"]["missing_sector_attachment"]["cannot_emit"]
