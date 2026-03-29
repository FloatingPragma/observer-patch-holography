#!/usr/bin/env python3
"""Smoke-test the charged two-scalar hierarchy-emitter artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
READOUT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_log_spectrum_readout.py"
AUDIT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_exactness_audit.py"
EMITTER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_two_scalar_hierarchy_emitter.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_two_scalar_hierarchy_emitter.json"


def main() -> int:
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EMITTER_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_current_family_two_scalar_hierarchy_emitter":
        print("wrong charged two-scalar emitter artifact id", file=sys.stderr)
        return 1
    if payload.get("hierarchy_emitter_status") != "missing_joint_emission":
        print("charged two-scalar emitter should remain unresolved", file=sys.stderr)
        return 1
    if payload.get("frozen_sigma_branch_impossible") is not True:
        print("charged two-scalar emitter should carry the frozen-sigma impossibility flag", file=sys.stderr)
        return 1
    if payload.get("sigma_e_total_log_per_side_emitted") is not None or payload.get("eta_e_split_log_per_side_emitted") is not None:
        print("predictive sigma/eta slots should remain unset until emitted from OPH inputs", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
