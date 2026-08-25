#!/usr/bin/env python3
"""Smoke-test the quark diagonal gap-shift scalar-evaluator artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

from derive_quark_diagonal_gap_shift_scalar_evaluator import build_artifact


ROOT = pathlib.Path(__file__).resolve().parents[2]
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
SOURCE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SOURCE_EMISSION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_emission.py"
SOURCE_VALUES_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_values.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_scalar_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_VALUES_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_gap_shift_scalar_evaluator":
        print("wrong quark diagonal scalar-evaluator artifact id", file=sys.stderr)
        return 1
    if payload.get("tau_u_log_per_side") is None or payload.get("tau_d_log_per_side") is None:
        print("quark scalar evaluator should retain the quarantined comparison arithmetic", file=sys.stderr)
        return 1
    if payload.get("source_value_promotion_ready") is not False:
        print("comparison-only scalar values must not be promotion ready", file=sys.stderr)
        return 1
    if payload.get("value_classification") != "comparison_only":
        print("quark scalar values should be classified as comparison-only", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "comparison_only_scalar_projection_value_source_open":
        print("quark scalar evaluator should preserve the open value-source status", file=sys.stderr)
        return 1
    delta = float(payload["beta_u_diag_B_source"]) - float(
        payload["beta_d_diag_B_source"]
    )
    sigma_u = float(payload["sigma_u_total_log_per_side"])
    sigma_d = float(payload["sigma_d_total_log_per_side"])
    expected_tau_u = sigma_d * delta / (2.0 * (sigma_u + sigma_d))
    expected_tau_d = sigma_u * delta / (2.0 * (sigma_u + sigma_d))
    if abs(float(payload["tau_u_log_per_side"]) - expected_tau_u) > 1.0e-12:
        print("tau_u should use the D12 opposite-sector sigma weighting", file=sys.stderr)
        return 1
    if abs(float(payload["tau_d_log_per_side"]) - expected_tau_d) > 1.0e-12:
        print("tau_d should use the D12 opposite-sector sigma weighting", file=sys.stderr)
        return 1
    if payload.get("constructive_missing_objects") != [
        "source_derived_c_d_over_c_u_or_t1_value",
        "supported_source_sigma_branch_emission",
    ]:
        print("scalar evaluator should retain both source-side blockers", file=sys.stderr)
        return 1
    if payload.get("source_values_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_values":
        print("quark scalar evaluator should reference the diagonal common gap-shift source values", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_derived_c_d_over_c_u_or_t1_value":
        print("quark scalar evaluator should point to the source-derived light ratio or t1", file=sys.stderr)
        return 1

    # Closing only the beta-value source is insufficient while the sigma
    # branch remains a hard-coded diagnostic witness.
    diagonal_map = json.loads(
        (ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json").read_text(
            encoding="utf-8"
        )
    )
    spread = json.loads(
        (ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json").read_text(
            encoding="utf-8"
        )
    )
    source_law = json.loads(
        (ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json").read_text(
            encoding="utf-8"
        )
    )
    source_values = json.loads(
        (ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_values.json").read_text(
            encoding="utf-8"
        )
    )
    source_values["proof_status"] = "source_values_derived_from_source_emission"
    source_values["value_classification"] = "source_emitted"
    beta_only = build_artifact(diagonal_map, spread, source_law, source_values)
    if beta_only.get("source_value_promotion_ready") is not False:
        print("source beta values must not upgrade a diagnostic sigma branch", file=sys.stderr)
        return 1
    if beta_only.get("promotion_blockers") != [
        "supported_source_sigma_branch_emission"
    ]:
        print("beta-only closure should retain the sigma-source blocker", file=sys.stderr)
        return 1

    spread["predictive_promotion_allowed"] = True
    spread["spread_emitter_status"] = "source_emitted"
    fully_sourced = build_artifact(diagonal_map, spread, source_law, source_values)
    if fully_sourced.get("source_value_promotion_ready") is not True:
        print("source-emitted beta and sigma inputs should close the local transport map", file=sys.stderr)
        return 1
    if fully_sourced.get("value_classification") != "source_emitted":
        print("fully sourced transport values should be classified as source-emitted", file=sys.stderr)
        return 1
    if fully_sourced.get("constructive_missing_objects") != []:
        print("fully sourced transport map should have no local source blockers", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
