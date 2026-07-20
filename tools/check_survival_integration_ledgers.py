#!/usr/bin/env python3
"""Verify compact survival-campaign ledgers against the workspace archives."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise AssertionError(f"expected JSON object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_a5(workspace_root: Path) -> None:
    ledger = load_json(
        REPO_ROOT / "code/capacity_readback/runtime/a5_finite_control_status.json"
    )
    source_path = (
        workspace_root
        / "survival-proof-2/artifacts/receipts/issue_closure_packet.json"
    )
    source = load_json(source_path)

    require(
        sha256(source_path) == ledger["source_issue_closure_receipt_sha256"],
        "A5 source receipt hash drift",
    )
    require(source["FINITE_PACKET_EVALUATOR_CONTROL_RECEIPT"] is True, "A5 control receipt failed")
    require(
        source["RAW_CARRIER_MULTIPLICITY_MINIMUM_RECEIPT"] is True,
        "A5 raw-minimum receipt failed",
    )
    require(source["PHYSICAL_N_CLOSURE_RECEIPT"] is False, "A5 packet was retyped as physical closure")
    require(source["status"] == ledger["status"], "A5 status drift")
    require(source["capacity_coordinate"]["public_zero_error_capacity_M0_public"] == 60, "A5 M0 drift")
    require(source["capacity_coordinate"]["intrinsic_capacity_coordinate_D_intrinsic"] is None, "A5 intrinsic D became defined")
    require(ledger["capacity_formula"] == "M_0(k,r)=60", "A5 ledger formula drift")
    require(ledger["raw_dimension_formula"] == "D_raw(k)=60*k", "A5 raw dimension formula drift")
    require(ledger["raw_defect_formula"] == "Delta_raw(k)=60*(k-1)", "A5 raw defect formula drift")


def check_h3_kms(workspace_root: Path) -> None:
    ledger = load_json(REPO_ROOT / "code/geometry/runs/h3_kms_repaired_4k_status.json")
    source_dir = workspace_root / "survival-proof-4/outputs/run_4k_repaired_20260720"
    run_path = source_dir / "campaign_run_receipt.json"
    replay_path = source_dir / "physical_h3_kms_replay_verification.json"
    manifest_path = source_dir / "physical_h3_kms_campaign_manifest.json"
    run = load_json(run_path)
    replay = load_json(replay_path)
    manifest = load_json(manifest_path)

    require(sha256(run_path) == ledger["source_campaign_run_receipt_sha256"], "H3/KMS run receipt hash drift")
    require(sha256(replay_path) == ledger["source_replay_verification_sha256"], "H3/KMS replay receipt hash drift")
    require(run["instrument_status"] == ledger["instrument_status"], "H3/KMS instrument status drift")
    require(run["campaign_complete"] is ledger["campaign_complete"], "H3/KMS campaign completion drift")
    require(run["cell_scientific_status"] == ledger["cell_scientific_status"], "H3/KMS cell status drift")
    require(run["physical_promotion_allowed"] is False, "H3/KMS run permits physical promotion")
    require(run["postrun_scientific_failures"] == [], "H3/KMS run contains a scientific failure")
    require(len(run["postrun_not_evaluated_reasons"]) == 7, "H3/KMS not-evaluated reason count drift")
    require(replay["scientific_status"] == ledger["replay_scientific_status"], "H3/KMS replay status drift")
    require(replay["physical_promotion_allowed"] is False, "H3/KMS replay permits physical promotion")
    require(manifest["rungs"] == ledger["rungs"], "H3/KMS rung ledger drift")
    require(manifest["seeds"] == ledger["seeds"], "H3/KMS seed ledger drift")
    require(manifest["physical_promotion_allowed"] is False, "H3/KMS campaign permits physical promotion")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=REPO_ROOT.parent,
        help="workspace containing survival-proof-2 and survival-proof-4",
    )
    args = parser.parse_args()
    check_a5(args.workspace_root.resolve())
    check_h3_kms(args.workspace_root.resolve())
    print("survival integration ledgers OK: A5 finite control; repaired 4k H3/KMS status")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
