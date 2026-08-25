#!/usr/bin/env python3
"""Validate the light-quark overlap-defect value-law scaffold."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SELECTOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_light_quark_isospin_overlap_defect_selector_law.py"
MASS_RAY_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_ray.py"
OVERLAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_light_quark_overlap_defect_value_law.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "light_quark_overlap_defect_value_law.json"


def main() -> int:
    for script in (SELECTOR_SCRIPT, MASS_RAY_SCRIPT, OVERLAP_SCRIPT, SPREAD_SCRIPT, SCRIPT):
        subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload["artifact"] != "oph_light_quark_overlap_defect_value_law":
        print("wrong selector-value artifact id", file=sys.stderr)
        return 1
    if payload["proof_status"] != "conditional_overlap_identity_selector_value_open":
        print("selector-value identity should remain conditional", file=sys.stderr)
        return 1
    if payload["public_promotion_allowed"] is not False:
        print("selector-value identity must not be publicly promotable", file=sys.stderr)
        return 1
    if payload["exact_missing_object"] != "source_derived_c_d_over_c_u_or_Delta_ud_overlap":
        print("selector-value identity must retain its source-value gap", file=sys.stderr)
        return 1
    if payload["selector_scalar_name"] != "Delta_ud_overlap":
        print("wrong selector scalar name", file=sys.stderr)
        return 1
    if payload["target_free_map"]["formula"] != "Theta_ud^Delta(c_u, c_d) = (1/6) * log(c_d / c_u)":
        print("selector-value identity should expose the conditional bridge formula", file=sys.stderr)
        return 1
    if payload["equivalent_ray_coordinate_presentation"]["theorem_id"] != "quark_d12_t1_value_law":
        print("selector-value scaffold should expose the t1 wrapper as equivalent presentation", file=sys.stderr)
        return 1
    if payload["induced_mass_side_formulas"]["t1"] != "5 * Delta_ud_overlap":
        print("selector-value scaffold should derive t1 algebraically", file=sys.stderr)
        return 1
    if payload["target_1_public_yukawa_consequence"]["does_target_1_close_after_this"] is not False:
        print("selector-value identity must not close target 1 without a value source", file=sys.stderr)
        return 1
    if payload["target_1_public_yukawa_consequence"]["target_1_status"] != "open_value_source_not_emitted":
        print("target 1 should remain open", file=sys.stderr)
        return 1
    if payload["exact_transport_contract_after_selector"]["single_residual_object_on_each_fixed_sigma_branch"] != "Delta_ud_overlap":
        print("transport contract should keep Delta_ud_overlap as the only fixed-branch residual", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
