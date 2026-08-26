#!/usr/bin/env python3
"""Validate the D10 minimal conditional-promotion artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_minimal_conditional_promotion.py"


def test_d10_minimal_conditional_promotion_records_underdetermination_and_smallest_route(
    tmp_path: pathlib.Path,
) -> None:
    output = tmp_path / "minimal_conditional.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(output)],
        check=True,
        cwd=ROOT,
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_minimal_conditional_promotion"
    assert payload["status"] == "open_split_beneath_target_free_candidate"
    assert payload["candidate_object_id"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["unconditional_theorem"]["name"] == "current_corpus_underdetermination_of_forward_d10_repair_law"
    assert payload["conditional_principle"]["name"] == "ColorBalancedQuadraticRepairDescent_D10"
    assert payload["conditional_theorem"]["name"] == "minimal_conditional_d10_forward_repair_law"
    specialization = payload["n_c_3_specialization"]
    assert abs(specialization["tau2_exact"] - (-4.2477649513449626e-4)) < 1.0e-18
    assert abs(specialization["delta_n_exact"] - 3.3952000966304088e-4) < 1.0e-18
    assert abs(specialization["coherent_quintet"]["MW_pole_gev"] - 80.36921677537332) < 1.0e-12
    assert abs(specialization["coherent_quintet"]["MZ_pole_gev"] - 91.18800674424540) < 1.0e-12
