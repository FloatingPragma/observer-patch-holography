#!/usr/bin/env python3
"""Validate the exact split D11 Higgs/top theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
REPAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_free_repair_value_law.py"
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
NO_GO_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_fixed_ray_no_go_theorem.py"
PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_split_pair_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json"


def test_d11_live_exact_split_pair_theorem_closes_exact_pair() -> None:
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REPAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAIR_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_live_exact_split_pair_theorem"
    assert payload["theorem_id"] == "D11DeclaredSurfaceSplitImplication"
    assert payload["legacy_theorem_id"] == "D11SourceSplitForwardExactness"
    assert payload["proof_status"] == "conditional_on_unpromoted_d10_repair_candidate"
    assert payload["status"] == "candidate_only"
    assert payload["public_surface_candidate_allowed"] is False
    assert payload["prediction_promotion_allowed"] is False
    assert payload["upstream_promotion_gate"]["passed"] is False
    assert payload["upstream_promotion_gate"]["actual_status"] == "candidate_only"
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "closed_promotable_EWTargetFreeRepairValueLaw_D10"
    )
    assert payload["closure_logic"]["fixed_ray_blocked"] is True
    assert payload["exact_split_pair"]["mH_gev"] == 125.1995304097179
    assert payload["exact_split_pair"]["mt_pole_gev"] == 172.3523553288312
    assert payload["exact_split_pair"]["w_HT_exact"] == pytest.approx(-0.0003857630977715052, abs=5.0e-18)
    assert payload["shared_split_scalar"]["value"] == pytest.approx(-0.00023118902229730438, abs=1.0e-18)
    assert "promotion_of_the_old_fixed_ray_as_exact_pair" in payload["strictly_not_claimed"]


def test_closed_d10_repair_does_not_promote_declared_d11_surface(tmp_path: pathlib.Path) -> None:
    repair = json.loads(
        (ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json").read_text(
            encoding="utf-8"
        )
    )
    repair["status"] = "closed"
    repair["promotion_allowed"] = True
    repair["proof_status"] = "closed_hypothetical_test_fixture"
    repair_path = tmp_path / "closed_d10_repair.json"
    output_path = tmp_path / "d11_pair.json"
    repair_path.write_text(json.dumps(repair), encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(PAIR_SCRIPT),
            "--d10-repair",
            str(repair_path),
            "--output",
            str(output_path),
        ],
        check=True,
        cwd=ROOT,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["upstream_promotion_gate"]["passed"] is True
    assert payload["status"] == "conditional_theorem_only"
    assert payload["source_surface_promotable"] is False
    assert payload["public_surface_candidate_allowed"] is False
    assert payload["prediction_promotion_allowed"] is False
    assert payload["non_circularity_status"]["promotion_allowed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "source_emitted_D11_core_Jacobian_residual_selectors_scale_and_pole_packet"
    )
