#!/usr/bin/env python3
"""Validate the selected-public-class quark D12 pure-B source payload."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_public_source_payload.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_public_source_payload.json"


def test_quark_d12_public_source_payload_emits_selected_class_values() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_d12_public_source_payload"
    assert payload["proof_status"] == (
        "blocked_unpromotable_yukawa_input_and_open_d12_value_source"
    )
    assert payload["public_promotion_allowed"] is False
    assert payload["off_canonical_promotion_allowed"] is False
    assert payload["value_classification"] == (
        "comparison_only_projection_from_unpromotable_inputs"
    )
    assert payload["promotion_blockers"] == [
        "PUBLIC_COMMON_SCALE_DIMENSIONLESS_SOURCE_YUKAWA_THEOREM",
        "SOURCE_DERIVED_C_D_OVER_C_U_OR_T1_VALUE",
    ]
    assert payload["input_status_audit"]["promotion_inherited_fail_closed"] is True
    assert payload["selected_public_physical_frame_class"]["selected_by"] == "P"

    expected_t1 = (5.0 / 6.0) * math.log(
        payload["comparison_only_light_ratio_input"]["m_d_gev"]
        / payload["comparison_only_light_ratio_input"]["m_u_gev"]
    )
    assert math.isclose(payload["d12_scalars"]["t1"], expected_t1, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(payload["beta_u_diag_B_source"], payload["d12_scalars"]["t1"] / 10.0)
    assert math.isclose(payload["beta_d_diag_B_source"], -payload["d12_scalars"]["t1"] / 10.0)
    assert payload["source_readback_u_log_per_side"][1] == 0.0
    assert payload["source_readback_d_log_per_side"][1] == 0.0
    assert math.isclose(payload["pure_B_certificates"]["endpoint_sum_u"], 0.0, abs_tol=1.0e-12)
    assert math.isclose(payload["pure_B_certificates"]["endpoint_sum_d"], 0.0, abs_tol=1.0e-12)
