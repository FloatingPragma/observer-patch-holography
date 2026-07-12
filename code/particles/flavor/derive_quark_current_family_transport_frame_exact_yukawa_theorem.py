#!/usr/bin/env python3
"""Wrap the legacy-named current-family quark mass-texture audit artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORWARD_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_forward_yukawas.json"
CHAIN_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_yukawa_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(forward: dict, chain: dict) -> dict:
    return {
        "artifact": "oph_quark_current_family_transport_frame_exact_yukawa_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_target_anchored_mixed_convention_mass_texture_audit",
        "target_name": "mixed_convention_quark_mass_textures_on_declared_current_family_transport_frame",
        "theorem_scope": "current_family_common_refinement_transport_frame_only",
        "public_promotion_allowed": False,
        "supporting_artifacts": {
            "exact_forward_yukawas": forward["artifact"],
            "end_to_end_exact_chain": chain["artifact"],
        },
        "theorem_statement": (
            "On the declared current-family/common-refinement transport-frame audit carrier, the target-attached "
            "spread lift and algebraic readout construct GeV-valued matrices whose singular values reproduce the "
            "chosen mixed-convention target coordinates. This certifies a mass-texture calculation, not one "
            "common-scale running sextet or physical dimensionless Yukawa matrices."
        ),
        "forward_yukawa_artifact": {
            "artifact": forward["artifact"],
            "scope": forward["scope"],
            "forward_certified": forward["forward_certified"],
            "certification_status": forward["certification_status"],
            "Y_u": forward["Y_u"],
            "Y_d": forward["Y_d"],
            "V_CKM": forward["V_CKM"],
            "jarlskog": forward["jarlskog"],
            "singular_values_u": forward["singular_values_u"],
            "singular_values_d": forward["singular_values_d"],
            "matrix_classification": forward.get("matrix_classification"),
            "stored_entry_dimension": forward.get("stored_entry_dimension"),
            "physical_yukawa_certified": False,
        },
        "source_chain_artifact": chain["artifact"],
        "minimal_exact_blocker_set": ["QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE"],
        "strengthening_above_target": {
            "status": "separate_question",
            "name": "source_only_common_scale_dimensionless_yukawa_construction",
            "note": (
                "A physical Yukawa theorem would require source-emitted running trajectories, one common scale, "
                "threshold and top conversion, the running Higgs expectation value in the same scheme, and "
                "dimensionless normalization."
            ),
        },
        "notes": [
            "This is a compatibility wrapper for the legacy-named transport-frame matrix artifact.",
            "It certifies only target-audit mass-texture algebra.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current-family transport-frame mass-texture audit artifact.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--chain", default=str(CHAIN_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.forward)),
        _load_json(Path(args.chain)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
