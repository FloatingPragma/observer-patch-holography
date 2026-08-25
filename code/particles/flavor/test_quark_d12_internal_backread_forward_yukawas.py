#!/usr/bin/env python3
"""Validate the continuation-only D12 internal-backread forward Yukawas."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DESCENT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_internal_backread_descent.py"
FORWARD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_internal_backread_forward_yukawas.py"
DESCENT_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_descent.json"
FORWARD_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_forward_yukawas.json"


def main() -> int:
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)

    descent = json.loads(DESCENT_OUTPUT.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_OUTPUT.read_text(encoding="utf-8"))

    if descent.get("artifact") != "oph_quark_d12_internal_backread_descent":
        print("wrong D12 internal-backread descent artifact id", file=sys.stderr)
        return 1
    if descent.get("exact_missing_object") is not None:
        print("continuation-only backread descent should not leave an exact missing object open on its own sidecar inputs", file=sys.stderr)
        return 1
    if descent.get("even_excitation_proof_status") != "closed":
        print("continuation-only backread descent should close the even excitation slot on its sidecar surface", file=sys.stderr)
        return 1
    if descent.get("J_B_source_u") is None or descent.get("J_B_source_d") is None:
        print("continuation-only backread descent should carry the emitted pure-B source values", file=sys.stderr)
        return 1

    if forward.get("artifact") != "oph_quark_d12_internal_backread_forward_yukawas":
        print("wrong D12 internal-backread forward Yukawas artifact id", file=sys.stderr)
        return 1
    if forward.get("scope") != "D12_continuation_internal_backread_only":
        print("continuation-only backread forward Yukawas should remain on the sidecar scope", file=sys.stderr)
        return 1
    if forward.get("forward_certified") is not True:
        print("continuation-only backread forward Yukawas should be certified on their own sidecar inputs", file=sys.stderr)
        return 1
    if forward.get("certification_status") != "forward_matrix_arithmetic_certified_source_uncertified":
        print("continuation-only backread should certify arithmetic without source provenance", file=sys.stderr)
        return 1
    if forward.get("predictive_promotion_allowed") is not False or forward.get("source_certified") is not False:
        print("continuation-only backread values must not become predictive source outputs", file=sys.stderr)
        return 1
    if forward.get("public_surface_candidate_allowed") is not False:
        print("source-uncertified continuation values must not become a public-surface candidate", file=sys.stderr)
        return 1
    if forward.get("promotion_blockers") != [
        "D12_INTERNAL_BACKREAD_NOT_SOURCE_DERIVED",
        "D12_RAY_VALUE_SOURCE_OPEN",
    ]:
        print("continuation-only backread forward Yukawas should expose their source blockers", file=sys.stderr)
        return 1
    if any(value <= 0.0 for value in forward["singular_values_u"] + forward["singular_values_d"]):
        print("continuation-only backread forward Yukawas should emit positive singular values", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
