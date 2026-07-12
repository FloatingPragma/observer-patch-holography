#!/usr/bin/env python3
"""Wrap the legacy-named end-to-end quark mass-texture target audit."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
YUKAWA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_yukawa_theorem.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(yukawa_theorem: dict) -> dict:
    return {
        "artifact": "oph_quark_exact_yukawa_end_to_end_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_target_anchored_mixed_convention_mass_texture_audit",
        "target_name": yukawa_theorem["target_name"],
        "theorem_scope": yukawa_theorem["theorem_scope"],
        "public_promotion_allowed": False,
        "supporting_theorem_artifact": yukawa_theorem["artifact"],
        "theorem_statement": (
            "On the explicit current-family/common-refinement target-audit chain, the target-attached spread lift "
            "and algebraic readout construct GeV-valued mass textures whose singular values reproduce the selected "
            "mixed-convention coordinates. This is not a source-only mass or physical-Yukawa theorem."
        ),
        "forward_yukawa_artifact": yukawa_theorem["forward_yukawa_artifact"],
        "minimal_exact_blocker_set": yukawa_theorem["minimal_exact_blocker_set"],
        "strengthening_above_target": yukawa_theorem["strengthening_above_target"],
        "notes": [
            "This compatibility wrapper preserves the legacy artifact name while classifying the output as a target-audit mass texture.",
            "Physical dimensionless Yukawas require the common-scale normalization blocker recorded here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the end-to-end quark mass-texture target-audit artifact.")
    parser.add_argument("--yukawa-theorem", default=str(YUKAWA_THEOREM_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.yukawa_theorem)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
