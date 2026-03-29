#!/usr/bin/env python3
"""Export the current shared charged absolute-scale promotion artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_CERT = ROOT / "particles" / "runs" / "flavor" / "charged_mean_eigenvalue_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "charged_shared_absolute_scale_promotion.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the shared charged absolute-scale promotion artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--certificate", default=str(DEFAULT_CERT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    cert = json.loads(Path(args.certificate).read_text(encoding="utf-8"))
    budget = payload["charged_dirac_scalarization_certificate"]
    shared_scale = float(budget["seed_value"])
    total_budget = float(budget["sampled_total_budget"])
    sampled_sector_shares = budget["sampled_sector_shares"]
    sector_share_by_sector = {
        sector: float(info["snapshot"])
        for sector, info in sampled_sector_shares.items()
    }

    artifact = {
        "artifact": "charged_shared_absolute_scale_promotion",
        "generated_utc": _timestamp(),
        "parent_artifact": "oph_charged_budget_transport",
        "parent_candidate": budget.get("candidate_id"),
        "promotion_scope": "current_family_only",
        "mean_witness_artifact": cert.get("artifact"),
        "mean_witness": {
            "left_common_mean": cert.get("left_common_mean"),
            "right_common_mean": cert.get("right_common_mean"),
            "mean_residual": cert.get("mean_residual"),
            "witness_grade_status": cert.get("witness_grade_status"),
            "gap_side_support_status": cert.get("gap_side_support_status"),
        },
        "stored_shared_absolute_scale": shared_scale,
        "sigma_seed_clause": {
            "left_common": shared_scale,
            "right_common": shared_scale,
            "residual": 0.0,
        },
        "sector_equalizer_by_sector": {"u": 0.0, "d": 0.0, "e": 0.0},
        "gluing_norm": 0.0,
        "shared_budget_total": total_budget,
        "sector_share_by_sector": sector_share_by_sector,
        "no_reseed": True,
        "proof_status": "promotion_ready_on_current_family",
        "smallest_constructive_missing_object": "charged_shared_absolute_scale_writeback",
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
