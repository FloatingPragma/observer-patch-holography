#!/usr/bin/env python3
"""Emit the single end-to-end exact-PDG quark theorem on the explicit chain."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHAIN_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_exact_pdg_end_to_end_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(chain: dict) -> dict:
    return {
        "artifact": "oph_quark_exact_pdg_end_to_end_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_exact_pdg_end_to_end_theorem",
        "target_name": chain["target_name"],
        "theorem_scope": chain["theorem_scope"],
        "public_promotion_allowed": False,
        "target_anchored": chain.get("target_anchored", True),
        "single_running_quark_sextet_claim_allowed": chain.get(
            "single_running_quark_sextet_claim_allowed", False
        ),
        "physical_dimensionless_yukawa_claim_allowed": chain.get(
            "physical_dimensionless_yukawa_claim_allowed", False
        ),
        "comparison_coordinate_partition": chain.get("comparison_coordinate_partition"),
        "supporting_chain_artifact": chain["artifact"],
        "theorem_statement": (
            "On the explicit current-family/common-refinement transport-frame audit chain, the target-attached sigma "
            "readout reconstructs the chosen mixed-convention quark comparison packet and its dimensionful mass "
            "textures. This wrapper does not turn that target packet into a source-only prediction, a common-scale "
            "running sextet, or physical dimensionless Yukawa matrices."
        ),
        "minimal_exact_blocker_set": chain["minimal_exact_blocker_set"],
        "exact_running_values_gev": chain["exact_running_values_gev"],
        "lemma_chain": chain["lemma_chain"],
        "strengthening_above_target": chain["strengthening_above_target"],
        "notes": [
            "This is the single audit wrapper above the explicit target-anchored derivation chain artifact.",
            "The top coordinate remains separate from the five running-mass comparison coordinates.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the end-to-end exact-PDG quark theorem artifact.")
    parser.add_argument("--chain", default=str(CHAIN_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.chain)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
