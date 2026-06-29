#!/usr/bin/env python3
"""Emit the public exact end-to-end Yukawa theorem above the descended sigma datum."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PUBLIC_SIGMA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
EXACT_PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_pdg_end_to_end_theorem.json"
EXACT_YUKAWA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"

MISSING_FOR_PROMOTION = [
    "QUARK_SIGMA_SOURCE_QUOTIENT",
    "QUARK_SIGMA_SOURCE_SELECTOR",
    "QUARK_EDGE_STATISTICS_CORRECTION_THEOREM",
    "QUARK_SIGMA_REFINEMENT_COMPATIBILITY",
    "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(
    public_sigma_theorem: dict[str, Any],
    exact_pdg_theorem: dict[str, Any],
    exact_yukawa_theorem: dict[str, Any],
) -> dict[str, Any]:
    forward = dict(exact_yukawa_theorem["forward_yukawa_artifact"])
    masses = dict(exact_pdg_theorem["supporting_theorem_artifact"] if "supporting_theorem_artifact" in exact_pdg_theorem else {})
    sigma_non_circularity = dict(public_sigma_theorem.get("non_circularity_status") or {})
    promotion_allowed = (
        public_sigma_theorem.get("public_promotion_allowed") is True
        and sigma_non_circularity.get("promotion_allowed", public_sigma_theorem.get("public_promotion_allowed")) is True
        and public_sigma_theorem.get("proof_status") == "closed_source_only_public_physical_sigma_datum_descent"
    )
    proof_status = (
        "closed_source_only_public_exact_yukawa_end_to_end_theorem"
        if promotion_allowed
        else "blocked_target_derived_sigma_source_missing"
    )
    return {
        "artifact": "oph_quark_public_exact_yukawa_end_to_end_theorem",
        "generated_utc": _timestamp(),
        "proof_status": proof_status,
        "target_name": (
            "source_only_public_exact_forward_quark_yukawas"
            if promotion_allowed
            else "selected_class_conditional_forward_quark_yukawas"
        ),
        "theorem_scope": public_sigma_theorem["theorem_scope"],
        "claim_tier": (
            "source_only_public_yukawa_theorem"
            if promotion_allowed
            else "selected_class_conditional_on_source_sigma"
        ),
        "public_promotion_allowed": promotion_allowed,
        "source_only_sigma_emitted": bool(public_sigma_theorem.get("source_only_sigma_emitted")),
        "downstream_algebra_closed": True,
        "display_allowed_as_selected_class_exact_witness": True,
        "non_circularity_status": {
            "promotion_allowed": promotion_allowed,
            "public_sigma_promotion_allowed": public_sigma_theorem.get("public_promotion_allowed"),
            "public_sigma_proof_status": public_sigma_theorem.get("proof_status"),
            "target_derived_sigma_datum_used": sigma_non_circularity.get("target_derived_sigma_datum_used"),
            "source_sigma_selector_closed": sigma_non_circularity.get("source_sigma_selector_closed"),
            "missing_source_object": None
            if promotion_allowed
            else "quark_sigma_source_datum_no_target_leak_required",
            "strict_audit_label": "source_only_public_yukawa_theorem"
            if promotion_allowed
            else "selected_class_conditional_on_source_sigma",
        },
        "supporting_theorem_artifacts": {
            "public_sigma_datum_descent": public_sigma_theorem["artifact"],
            "local_exact_pdg_wrapper": exact_pdg_theorem["artifact"],
            "local_exact_yukawa_wrapper": exact_yukawa_theorem["artifact"],
        },
        "theorem_statement": (
            "Conditional on a no-target source theorem emitting the physical sigma datum on the public quark frame "
            "class selected by P, the selected bridge-fiber descent makes that datum representative-independent, "
            "the affine mean law emits the absolute sector scales algebraically, the ordered three-point readout "
            "yields the running quark sextet, and the exact forward construction emits explicit Yukawa matrices "
            "Y_u and Y_d on that selected class. The current artifact displays the selected-class exact witness but "
            "does not promote it because the sigma datum is inherited from the current-family exact target surface."
        ),
        "selected_public_physical_frame_class": public_sigma_theorem["selected_public_physical_frame_class"],
        "descended_physical_sigma_datum": public_sigma_theorem["descended_physical_sigma_datum"],
        "public_exact_outputs": {
            "exact_running_values_gev": exact_pdg_theorem["exact_running_values_gev"],
            "forward_yukawa_artifact": forward,
        },
        "lemma_chain": [
            public_sigma_theorem["artifact"],
            "oph_quark_absolute_readout_algebraic_collapse",
            exact_pdg_theorem["artifact"],
            exact_yukawa_theorem["artifact"],
        ],
        "minimal_exact_blocker_set": []
        if promotion_allowed
        else MISSING_FOR_PROMOTION,
        "notes": [
            (
                "This is the theorem wrapper that upgrades the closed declared-carrier exact Yukawa chain to the selected public class."
                if promotion_allowed
                else "This artifact displays the selected-class exact witness but is not promotable under the strict non-circularity audit because the public sigma datum descends from an exact target surface."
            ),
            "The exact matrices Y_u and Y_d are the same numerical matrices already emitted on the closed local chain.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the public exact Yukawa end-to-end theorem artifact.")
    parser.add_argument("--public-sigma-theorem", default=str(PUBLIC_SIGMA_JSON))
    parser.add_argument("--exact-pdg-theorem", default=str(EXACT_PDG_JSON))
    parser.add_argument("--exact-yukawa-theorem", default=str(EXACT_YUKAWA_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.public_sigma_theorem)),
        _load_json(Path(args.exact_pdg_theorem)),
        _load_json(Path(args.exact_yukawa_theorem)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
