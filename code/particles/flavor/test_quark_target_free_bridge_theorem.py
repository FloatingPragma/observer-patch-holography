#!/usr/bin/env python3
"""Validate the synthesized quark target-free bridge theorem package."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_target_free_bridge_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_target_free_bridge_theorem.json"


def test_quark_target_free_bridge_theorem_package() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_target_free_bridge_theorem"
    assert payload["proof_status"] == "conditional_bridge_identity_value_source_open"
    assert payload["public_promotion_allowed"] is False
    assert payload["bridge_closure_status"] == "conditional_algebra_closed_value_emission_open"
    assert payload["exact_missing_object"] == "source_derived_c_d_over_c_u_or_t1_value"
    assert payload["theorem_ids"] == [
        "light_quark_overlap_defect_value_law",
        "quark_d12_t1_value_law",
    ]
    assert "Delta_ud_overlap = (1/6) * log(c_d / c_u)" in payload["principal_theorem_statement"]
    assert payload["proof_skeleton"][0]["statement"].endswith("log(y_d / y_u) = log(c_d / c_u).")
    assert "(1/6) * log(y_d / y_u)" in payload["proof_skeleton"][1]["statement"]
    assert payload["single_bridge_gap_to_internalize"] == "source_derived_c_d_over_c_u_or_t1_value"
    assert payload["computed_current_family_target_check"] is None
    assert payload["target_check_attachment"] == {
        "attached": False,
        "classification": "explicit_comparison_only_target_anchored_sidecar",
        "used_as_proof_premise": False,
        "changes_public_promotion": False,
    }


def test_explicit_target_sidecar_cannot_promote_bridge() -> None:
    target = ROOT / "particles" / "runs" / "flavor" / "quark_d12_current_family_target_anchored_value_package.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--current-family-target", str(target)],
        check=True,
    )
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["computed_current_family_target_check"] is not None
    assert payload["target_check_attachment"]["attached"] is True
    assert payload["target_check_attachment"]["used_as_proof_premise"] is False
    assert payload["public_promotion_allowed"] is False

    # Restore the canonical default artifact: no target sidecar is attached.
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
