#!/usr/bin/env python3
"""Validate the D10 W-anchor / neutral-shear diagnostics."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_w_anchor_neutral_shear_factorization.py"
FACTORIZATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
BOX = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_box_dominance.json"
SPLIT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_freeze_and_subobject_split.json"


def main() -> int:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    factorization = json.loads(FACTORIZATION.read_text(encoding="utf-8"))
    box = json.loads(BOX.read_text(encoding="utf-8"))
    split = json.loads(SPLIT.read_text(encoding="utf-8"))

    if factorization.get("artifact") != "oph_d10_ew_w_anchor_neutral_shear_factorization":
        print("unexpected factorization artifact", file=sys.stderr)
        return 1
    if factorization.get("status") != "closed_freeze_once_coherent_repair_law":
        print("factorization should close the freeze-once repair law", file=sys.stderr)
        return 1
    exact_missing_law = factorization.get("exact_missing_law") or {}
    if exact_missing_law.get("object_id") != "FreezeOnceCoherentD10ElectroweakRepairLaw_D10":
        print("factorization should expose the freeze-once repair law object id", file=sys.stderr)
        return 1
    central = factorization.get("central_target_point") or {}
    if central.get("delta_MZ_after_exact_W_anchor_mev") is None:
        print("factorization should expose the post-W-anchor Z excess", file=sys.stderr)
        return 1
    if (central.get("neutral_shear_share_of_total") or 0.0) <= 0.5:
        print("central neutral-shear share should dominate", file=sys.stderr)
        return 1
    repaired_quintet = factorization.get("coherent_repaired_quintet") or {}
    if abs((repaired_quintet.get("MW_pole") or 0.0) - 80.377) > 1.0e-12:
        print("repaired quintet should emit exact frozen W", file=sys.stderr)
        return 1
    if abs((repaired_quintet.get("MZ_pole") or 0.0) - 91.18797809193725) > 1.0e-12:
        print("repaired quintet should emit exact frozen Z", file=sys.stderr)
        return 1
    verdict = box.get("verdict") or {}
    if verdict.get("neutral_shear_dominates_total_hypercharge_repair_everywhere") is not True:
        print("neutral shear should dominate over the whole target box", file=sys.stderr)
        return 1
    if split.get("target_freeze_required") is not True:
        print("target-freeze split should require a frozen target surface", file=sys.stderr)
        return 1
    if split.get("status") != "closed_freeze_once_subobject_split":
        print("split should close the freeze-once subobject split", file=sys.stderr)
        return 1
    subobject_split = split.get("subobject_split") or {}
    if subobject_split.get("neutral_shear_object") != "EWNeutralShearLaw_D10":
        print("split should identify the neutral-shear subobject", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
