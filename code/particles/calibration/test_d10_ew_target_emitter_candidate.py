#!/usr/bin/env python3
"""Validate the D10 target-emitter candidate artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_emitter_candidate.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_emitter_candidate.json"


def test_d10_target_emitter_candidate_records_near_exact_source_only_surface() -> None:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_target_emitter_candidate"
    assert payload["status"] == "new_theorem_candidate_not_yet_in_corpus"
    assert payload["object_id"] == "EWTargetEmitter_D10"
    law = payload["target_emitter_law"]
    assert abs(law["tau2_tree_exact"] - (-2.311623001746158e-4)) < 1.0e-18
    assert abs(law["delta_n_tree_exact"] - 2.346358802434819e-4) < 1.0e-18
    quintet = payload["coherent_emitted_quintet"]
    assert abs(quintet["MW_pole"] - 80.37700001539531) < 1.0e-12
    assert abs(quintet["MZ_pole"] - 91.18797807794321) < 1.0e-12
    comparison = payload["comparison_to_frozen_local_reference_surface"]
    assert abs(comparison["MW_difference_gev"]) < 2.0e-8
    assert abs(comparison["MZ_difference_gev"]) < 2.0e-8
