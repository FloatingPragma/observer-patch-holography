#!/usr/bin/env python3
"""Validate the charged absolute-scale underdetermination theorem artifact."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_scale_underdetermination_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_absolute_scale_underdetermination_theorem.json"


def test_charged_absolute_scale_is_explicitly_underdetermined() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_absolute_scale_underdetermination_theorem"
    assert payload["proof_status"] == "centered_shape_closed_absolute_scale_no_go_common_shift"
    assert payload["public_promotion_allowed"] is False
    assert payload["charged_absolute_equalizer"] == "NO_GO_COMMON_SHIFT"
    assert abs(payload["centered_sum_rule"]["value"]) < 1.0e-12
    assert payload["theorem_forbid_emit"] == ["g_e", "Delta_e_abs"]
    assert payload["next_exact_missing_object"] == "charged_absolute_anchor_A_ch"
    assert payload["minimal_new_theorem"]["required_new_scalar"] == "A_ch"
    assert payload["no_go_theorem"]["id"] == "charged_absolute_shift_invariance_no_go"
    assert payload["no_go_theorem"]["target_transforms"]["g_e"] == "exp(c) * g_e"
    assert payload["future_single_slot_only"]["required_contract"] == "A_ch(logm + c*(1,1,1)) = A_ch(logm) + c"
    assert payload["hard_reject"]["g_e"] == 0.6822819838027987

    compare = payload["compare_only_continuation_target"]
    assert math.isclose(compare["g_e_star"], 0.04577883668665009, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(compare["delta_e_abs_star"], 3.0039867953851016, rel_tol=0.0, abs_tol=1.0e-12)
