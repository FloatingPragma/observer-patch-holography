#!/usr/bin/env python3
"""Validate the selected-class public exact quark evaluator."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_selected_class_public_exact_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_selected_class_public_exact_evaluator.json"


def test_quark_selected_class_public_exact_evaluator_is_closed_but_not_off_canonical() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_selected_class_public_exact_evaluator"
    assert payload["proof_status"] == (
        "blocked_source_spread_nonidentifiability_and_physical_yukawa_normalization"
    )
    assert payload["theorem_scope"] == "selected_public_physical_quark_frame_class_only"
    assert payload["claim_tier"] == "selected_class_target_audit_mass_texture_only"
    assert payload["public_promotion_allowed"] is False
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is False
    assert payload["mass_texture_algebra_closed"] is True
    assert payload["physical_yukawa_construction_closed"] is False
    assert payload["arbitrary_P_off_canonical_motion_closed"] is False
    assert payload["selector"]["selected_public_physical_frame_class"]["selected_by"] == "P"
    assert "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM" in payload["closure"]["minimal_exact_blocker_set"]
    assert payload["closure"]["global_frame_classification_claimed"] is False

    assert payload["closure"]["selected_class_exact"] is False
    assert payload["closure"]["selected_class_target_audit"] is True
    masses = payload["target_audit_mass_coordinates_gev"]
    assert math.isclose(masses["u"], 0.00216, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(masses["d"], 0.0047, rel_tol=0.0, abs_tol=1.0e-12)
    assert payload["target_audit_mass_textures"]["forward_certified"] is True
