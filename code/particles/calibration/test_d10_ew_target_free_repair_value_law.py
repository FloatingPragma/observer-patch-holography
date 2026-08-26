#!/usr/bin/env python3
"""Validate the current D10 target-free repair candidate artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_free_repair_value_law.py"


def test_d10_target_free_repair_value_law_records_candidate_only_source_only_quintet(
    tmp_path: pathlib.Path,
) -> None:
    output = tmp_path / "target_free.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(output)],
        check=True,
        cwd=ROOT,
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_target_free_repair_value_law"
    assert payload["status"] == "candidate_only"
    assert payload["object_id"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["promotion_allowed"] is False
    assert payload["promotion_blockers"] == ["current_corpus_underdetermination_of_forward_d10_repair_law"]
    chart = payload["repair_chart"]
    assert abs(chart["tau2_tree_exact"] - (-2.311623001746158e-4)) < 1.0e-18
    assert abs(chart["delta_n_tree_exact"] - 2.346358802434819e-4) < 1.0e-18
    quintet = payload["coherent_emitted_quintet"]
    assert abs(quintet["MW_pole"] - 80.37700001539531) < 1.0e-12
    assert abs(quintet["MZ_pole"] - 91.18797807794321) < 1.0e-12
    comparison = payload["compare_only_validation_against_frozen_surface"]
    assert abs(comparison["delta_MW_gev"]) > 1.0e-3
    assert abs(comparison["delta_MZ_gev"]) > 1.0e-5
