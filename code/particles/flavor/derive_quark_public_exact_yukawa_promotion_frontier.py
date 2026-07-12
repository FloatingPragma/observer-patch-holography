#!/usr/bin/env python3
"""Emit the fail-closed public quark mass/Yukawa promotion frontier."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PUBLIC_SIGMA_FRONTIER_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_strengthened_physical_sigma_lift_frontier.json"
)
PUBLIC_SIGMA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
)
EXACT_YUKAWA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"
)
PUBLIC_EXACT_YUKAWA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_promotion_frontier.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(
    public_sigma_frontier: dict[str, Any],
    public_sigma_theorem: dict[str, Any],
    exact_yukawa_theorem: dict[str, Any],
    public_exact_yukawa_theorem: dict[str, Any],
) -> dict[str, Any]:
    local_forward = dict(exact_yukawa_theorem["forward_yukawa_artifact"])
    promotion_allowed = (
        public_exact_yukawa_theorem.get("public_promotion_allowed") is True
        and public_exact_yukawa_theorem.get("proof_status")
        == "closed_source_only_public_exact_yukawa_end_to_end_theorem"
    )
    return {
        "artifact": "oph_quark_public_exact_yukawa_promotion_frontier",
        "generated_utc": _timestamp(),
        "proof_status": public_exact_yukawa_theorem["proof_status"],
        "target_name": public_exact_yukawa_theorem["target_name"],
        "scope": public_exact_yukawa_theorem["theorem_scope"],
        "public_promotion_allowed": promotion_allowed,
        "non_circularity_status": public_exact_yukawa_theorem.get("non_circularity_status"),
        "resolved_by_theorem_artifact": public_exact_yukawa_theorem["artifact"],
        "final_public_theorem_candidate": {
            "id": public_sigma_theorem["theorem_id"],
            "statement": public_sigma_theorem["theorem_statement"],
            "selected_public_exact_sigma_datum": public_sigma_theorem["descended_physical_sigma_datum"],
            "selected_public_frame_class": public_sigma_theorem["selected_public_physical_frame_class"],
        },
        "alternate_upstream_route": dict(public_sigma_frontier["alternate_upstream_route"]),
        "closed_local_endpoint": {
            "artifact": exact_yukawa_theorem["artifact"],
            "proof_status": exact_yukawa_theorem["proof_status"],
            "theorem_scope": exact_yukawa_theorem["theorem_scope"],
            "forward_yukawa_artifact": local_forward,
            "classification": "target_anchored_dimensionful_mass_texture_audit",
            "certified_physical_yukawa_matrices": False,
        },
        "closed_public_endpoint": {
            "artifact": public_exact_yukawa_theorem["artifact"],
            "proof_status": public_exact_yukawa_theorem["proof_status"],
            "public_exact_outputs": public_exact_yukawa_theorem["public_exact_outputs"],
        },
        "conditional_public_closure_if_final_theorem_proved": {
            "statement": (
                "The final public theorem candidate is now closed, so this conditional closure is realized."
                if promotion_allowed
                else "Public promotion is blocked until a new source map breaks the two-modulus spread action and a common-scale dimensionless Yukawa construction is certified."
            ),
            "induced_exact_outputs": public_exact_yukawa_theorem["public_exact_outputs"]["forward_yukawa_artifact"],
        },
        "minimal_exact_blocker_set": list(
            public_exact_yukawa_theorem.get("minimal_exact_blocker_set", [])
        ),
        "notes": [
            (
                "This frontier is now resolved: a source-only public sigma-datum theorem is closed, and the public exact Yukawa theorem is emitted explicitly."
                if promotion_allowed
                else "This frontier is closed negatively on the current corpus: the source spread has two free positive moduli, and the stored GeV matrices lack common-scale dimensionless normalization."
            ),
            "The alternate upstream route begins from a candidate generator over a hand-written family-kernel template and does not break the two-modulus action.",
            "The local declared-carrier matrix artifact is retained as a target-audit mass texture, not a physical Yukawa endpoint.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the public exact-Yukawa promotion frontier artifact.")
    parser.add_argument("--public-sigma-frontier", default=str(PUBLIC_SIGMA_FRONTIER_JSON))
    parser.add_argument("--public-sigma-theorem", default=str(PUBLIC_SIGMA_THEOREM_JSON))
    parser.add_argument("--exact-yukawa-theorem", default=str(EXACT_YUKAWA_THEOREM_JSON))
    parser.add_argument("--public-exact-yukawa-theorem", default=str(PUBLIC_EXACT_YUKAWA_THEOREM_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.public_sigma_frontier)),
        _load_json(Path(args.public_sigma_theorem)),
        _load_json(Path(args.exact_yukawa_theorem)),
        _load_json(Path(args.public_exact_yukawa_theorem)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
