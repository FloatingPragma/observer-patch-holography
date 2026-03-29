#!/usr/bin/env python3
"""Smoke-test the Lambda_MSbar descendant artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "qcd" / "derive_lambda_msbar_descendant.py"
ARTIFACT = ROOT / "particles" / "runs" / "qcd" / "lambda_msbar_descendant.json"


def main() -> int:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_qcd_lambda_msbar3":
        print("wrong Lambda_MSbar artifact id", file=sys.stderr)
        return 1
    outputs = payload.get("outputs", {})
    lambda5 = float(outputs["Lambda_MSbar_5_gev"])
    lambda4 = float(outputs["Lambda_MSbar_4_gev"])
    lambda3 = float(outputs["Lambda_MSbar_3_gev"])
    if not (lambda5 > 0.0 and lambda4 > 0.0 and lambda3 > 0.0):
        print("Lambda_MSbar outputs must be positive", file=sys.stderr)
        return 1
    if not (lambda3 >= lambda4 >= lambda5):
        print("expected Lambda_MSbar threshold chain to grow toward lower n_f", file=sys.stderr)
        return 1
    chain = payload.get("threshold_chain", [])
    if [entry.get("n_f") for entry in chain] != [5, 4, 3]:
        print("threshold chain should run 5 -> 4 -> 3", file=sys.stderr)
        return 1
    if payload.get("law_id") != "oph_qcd_lambda_msbar3_from_d10_alpha3":
        print("unexpected Lambda_MSbar law id", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
