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
    source_dir = workspace_root / "survival-proof-4/outputs/run_4k_acceptance_20260720_v5"
    run_path = source_dir / "campaign_run_receipt.json"
    replay_path = source_dir / "physical_h3_kms_replay_verification.json"
    manifest_path = source_dir / "physical_h3_kms_campaign_manifest.json"
    preflight_path = source_dir / "physical_h3_kms_preflight_report.json"
    geometry_controls_path = source_dir / "physical_h3_kms_geometry_controls_report.json"
    state_path = source_dir / "prime_geometric_cap_state_report.json"
    refinement_path = source_dir / "physical_h3_kms_refinement_report.json"
    replay_manifest_path = source_dir / "replay_bundle/replay_manifest.json"
    numerical_runtime_path = source_dir / "freeze/numerical_runtime.json"
    run = load_json(run_path)
    replay = load_json(replay_path)
    manifest = load_json(manifest_path)
    preflight = load_json(preflight_path)
    geometry_controls = load_json(geometry_controls_path)
    state = load_json(state_path)
    refinement = load_json(refinement_path)

    require(sha256(run_path) == ledger["source_campaign_run_receipt_sha256"], "H3/KMS run receipt hash drift")
    require(sha256(replay_path) == ledger["source_replay_verification_sha256"], "H3/KMS replay receipt hash drift")
    require(sha256(manifest_path) == ledger["source_campaign_manifest_sha256"], "H3/KMS manifest hash drift")
    require(sha256(preflight_path) == ledger["source_physical_preflight_sha256"], "H3/KMS preflight hash drift")
    require(
        sha256(geometry_controls_path) == ledger["source_geometry_controls_sha256"],
        "H3/KMS geometry-controls hash drift",
    )
    require(
        sha256(replay_manifest_path) == ledger["source_replay_manifest_sha256"],
        "H3/KMS replay manifest hash drift",
    )
    require(
        sha256(numerical_runtime_path) == ledger["source_numerical_runtime_sha256"],
        "H3/KMS numerical-runtime receipt hash drift",
    )
    require(
        run["replay_manifest_byte_sha256"].removeprefix("sha256:")
        == ledger["source_replay_manifest_sha256"],
        "H3/KMS replay manifest binding drift",
    )
    require(
        run["source_capture_sha256"].removeprefix("sha256:")
        == ledger["source_capture_semantic_sha256"],
        "H3/KMS source-capture semantic hash drift",
    )
    require(
        run["plan_sha256"].removeprefix("sha256:") == ledger["campaign_plan_sha256"],
        "H3/KMS campaign plan hash drift",
    )
    require(ledger["numerical_runtime_bound"] is True, "H3/KMS runtime binding was dropped")
    require(run["instrument_status"] == ledger["instrument_status"], "H3/KMS instrument status drift")
    require(run["campaign_complete"] is ledger["campaign_complete"], "H3/KMS campaign completion drift")
    require(run["cell_scientific_status"] == ledger["cell_scientific_status"], "H3/KMS cell status drift")
    require(run["physical_promotion_allowed"] is False, "H3/KMS run permits physical promotion")
    require(run["postrun_scientific_failures"] == [], "H3/KMS run contains a scientific failure")
    require(len(run["postrun_not_evaluated_reasons"]) == 7, "H3/KMS not-evaluated reason count drift")
    require(replay["scientific_status"] == ledger["replay_scientific_status"], "H3/KMS replay status drift")
    require(replay["physical_promotion_allowed"] is False, "H3/KMS replay permits physical promotion")
    require(preflight["verdict"] == ledger["physical_preflight_verdict"], "H3/KMS preflight verdict drift")
    require(preflight["PHYSICAL_H3_KMS_PREFLIGHT_RECEIPT"] is False, "H3/KMS preflight was promoted")
    require(preflight["physical_promotion_allowed"] is False, "H3/KMS preflight permits promotion")
    require(preflight["blocker_count"] == ledger["preflight_blocker_count"], "H3/KMS blocker count drift")
    require(preflight["stage_gate_summary"] == ledger["physical_stage_counts"], "H3/KMS stage counts drift")
    replay_receipts = {
        key: value
        for key, value in replay.items()
        if key.endswith("RECEIPT") and isinstance(value, bool)
    }
    require(
        sum(value is True for value in replay_receipts.values())
        == ledger["replay_receipt_pass_count"],
        "H3/KMS replay pass count drift",
    )
    require(
        sum(value is False for value in replay_receipts.values())
        == ledger["replay_receipt_fail_count"],
        "H3/KMS replay failure count drift",
    )
    require(manifest["rungs"] == ledger["rungs"], "H3/KMS rung ledger drift")
    require(manifest["seeds"] == ledger["seeds"], "H3/KMS seed ledger drift")
    require(manifest["campaign_status"] == ledger["campaign_status"], "H3/KMS campaign status drift")
    require(manifest["physical_promotion_allowed"] is False, "H3/KMS campaign permits physical promotion")
    require(manifest["branch_retirement_authorized"] is False, "H3/KMS campaign permits retirement")
    cell_counts = preflight["stages"]["P8_frozen_multiseed_four_rung_campaign"]["evidence"]["cell_status_counts"]
    require(cell_counts == ledger["campaign_cell_counts"], "H3/KMS campaign cell counts drift")
    require(
        refinement["levels"][-1]["patch_count"] == ledger["support_regulator_cell_count"],
        "H3/KMS support-regulator count drift",
    )
    require(
        {
            name: geometry_controls["models"][name]["heldout_score"]
            for name in ("E4", "S2", "H3", "E3")
        }
        == ledger["diagnostic_geometry_scores"],
        "H3/KMS diagnostic geometry scores drift",
    )
    require(
        geometry_controls["physical_gate_eligible"]
        is ledger["diagnostic_geometry_physical_gate_eligible"]
        is False,
        "H3/KMS diagnostic geometry comparison was promoted",
    )
    require(state["state_status"] == ledger["prime_state_status"], "H3/KMS state status drift")
    require(
        state["SOURCE_SELECTION_A5_QUOTIENT_INVARIANCE_RECEIPT"] is False,
        "H3/KMS presentation-bound state selection was promoted",
    )
    require(state["mixed_gns"]["constructed"] is True, "H3/KMS mixed-GNS diagnostic disappeared")
    require(
        refinement["CONSTRUCTED_REPEATED_RHO_IDENTITY_RECEIPT"] is True,
        "H3/KMS repeated-rho diagnostic drift",
    )
    require(
        refinement["PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE"] is False,
        "H3/KMS diagnostic was promoted to MGNS-1",
    )
    require(
        refinement["conditional_expectations_receipt"] is False,
        "H3/KMS conditional expectations were promoted",
    )
    bw_evidence = preflight["stages"]["P4_native_bw01_bw08"]["evidence"]
    require(bw_evidence["native_payload_conformance_receipt"] is True, "native BW payload does not parse")
    require(bw_evidence["native_payload_receipt"] is False, "native BW payload was promoted")
    passed_bw = sorted(
        clause_id
        for clause_id, clause in bw_evidence["clauses"].items()
        if clause["passed"] is True
    )
    require(passed_bw == ledger["native_bw_passing_clauses"], "native BW clause status drift")
    event_evidence = preflight["stages"]["P7_semantic_event_e1_e4_and_frame_fiber_separation"]["evidence"]
    require(event_evidence["event_clause_status"] == ledger["event_clause_status"], "event clause status drift")
    require(event_evidence["heldout_quadratic_cone_passed"] is False, "Lorentzian event cone was promoted")


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
    print(
        "survival integration ledgers OK: A5 finite control; H3/KMS archive status"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
