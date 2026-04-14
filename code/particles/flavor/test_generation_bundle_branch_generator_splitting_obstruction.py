#!/usr/bin/env python3
"""Validate the exact current-corpus obstruction certificate for branch-generator splitting."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_generation_bundle_branch_generator_splitting_obstruction.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator_splitting_obstruction.json"


def test_generation_bundle_branch_generator_splitting_obstruction_certificate() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_generation_bundle_branch_generator_splitting_obstruction_certificate"
    assert payload["proof_status"] == "closed_exact_current_corpus_obstruction_certificate"
    assert payload["target_theorem"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["target_clause"] == (
        "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split"
    )
    obstruction = payload["current_attached_data_obstruction"]
    assert abs(float(obstruction["commutator_operator_norm"]) - 0.04861550547372144) < 1.0e-12
    assert abs(float(obstruction["projector_defect_operator_norm"]) - 0.06363734112184061) < 1.0e-12
    assert abs(float(obstruction["commutator_over_defect_squared"]) - 12.004684594077071) < 1.0e-12
    exact_vanishing = payload["exact_vanishing_not_forced"]
    assert exact_vanishing["min_same_label_overlap_amplitude"] < 1.0
    assert payload["issue_149_resolution_mode"] == "sharp_exact_obstruction"
