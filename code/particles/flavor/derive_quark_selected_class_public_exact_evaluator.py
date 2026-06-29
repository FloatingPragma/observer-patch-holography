#!/usr/bin/env python3
"""Emit the selected-class conditional exact quark evaluator.

Chain role: package the evaluator that is algebraically closed on the public
quark frame class selected by P once a physical sigma datum is supplied.

This is intentionally not the arbitrary-P off-canonical moving evaluator. It
records the conditional composition

    P selects f_P + source sigma datum -> absolute readout -> exact masses/Yukawas.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BRIDGE = ROOT / "particles" / "runs" / "flavor" / "quark_target_free_bridge_theorem.json"
DEFAULT_PUBLIC_SIGMA = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
DEFAULT_ABSOLUTE_READOUT = ROOT / "particles" / "runs" / "flavor" / "quark_absolute_readout_algebraic_collapse.json"
DEFAULT_PUBLIC_YUKAWA = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_selected_class_public_exact_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(
    bridge: dict[str, Any],
    public_sigma: dict[str, Any],
    absolute_readout: dict[str, Any],
    public_yukawa: dict[str, Any],
) -> dict[str, Any]:
    yukawas = dict(public_yukawa["public_exact_outputs"]["forward_yukawa_artifact"])
    yukawas.setdefault("promotion_blockers", [])
    return {
        "artifact": "oph_quark_selected_class_public_exact_evaluator",
        "generated_utc": _timestamp(),
        "proof_status": "selected_class_conditional_exact_evaluator_blocked_by_source_sigma_selector",
        "theorem_scope": "selected_public_physical_quark_frame_class_only",
        "claim_tier": "selected_class_conditional_on_source_sigma",
        "public_promotion_allowed": False,
        "source_only_sigma_emitted": bool(public_sigma.get("source_only_sigma_emitted")),
        "downstream_algebra_closed": True,
        "arbitrary_P_off_canonical_motion_closed": False,
        "selector": {
            "input": "P",
            "output": "f_P",
            "selected_public_physical_frame_class": public_sigma["selected_public_physical_frame_class"],
        },
        "provenance": {
            "target_free_bridge_theorem": bridge.get("artifact"),
            "public_physical_sigma_datum_descent": public_sigma.get("artifact"),
            "absolute_readout_algebraic_collapse": absolute_readout.get("artifact"),
            "public_exact_yukawa_end_to_end_theorem": public_yukawa.get("artifact"),
        },
        "bridge_scalars": {
            "source_artifact": bridge.get("artifact"),
            "presentation": bridge.get("equivalent_wrappers"),
            "computed_current_family_target_check": bridge.get("computed_current_family_target_check"),
        },
        "sigma_datum": public_sigma["descended_physical_sigma_datum"],
        "absolute_readout": {
            "source_artifact": absolute_readout.get("artifact"),
            "proof_status": absolute_readout.get("proof_status"),
            "readout_contract": absolute_readout.get("readout_contract"),
            "emitted_values": absolute_readout.get("emitted_values"),
        },
        "masses": public_yukawa["public_exact_outputs"]["exact_running_values_gev"],
        "yukawas": yukawas,
        "closure": {
            "minimal_exact_blocker_set": list(
                public_yukawa.get(
                    "minimal_exact_blocker_set",
                    [
                        "QUARK_SIGMA_SOURCE_SELECTOR",
                        "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
                    ],
                )
            ),
            "selected_class_exact": True,
            "global_frame_classification_claimed": False,
            "off_canonical_P_family_claimed": False,
        },
        "theorem_statement": (
            "Conditional on a no-target source theorem emitting the physical sigma datum, the selected public quark "
            "frame class f_P carries a representative-independent sigma datum. The affine absolute readout then emits "
            "the sector scales, and the exact forward construction emits the running quark sextet and explicit Y_u, "
            "Y_d on that selected class. The current sigma datum is target-derived and is therefore audit/support only."
        ),
        "off_canonical_boundary": {
            "not_closed_here": [
                "edge_statistics_sigma_lift",
                "off_canonical_odd_response_kappa_value_law",
                "arbitrary_P_public_quark_frame_classification",
            ],
            "reason": (
                "The selected-class theorem does not classify all public quark frame classes and does not emit a "
                "continuous arbitrary-P sigma(P) or mass(P) family."
            ),
        },
        "notes": [
            "Use this artifact for selected-class conditional exact masses and Yukawas, not public source-only predictions.",
            "Use quark_p_driven_shared_evaluator_contract only for the candidate off-canonical slider surface.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the selected-class public exact quark evaluator.")
    parser.add_argument("--bridge", default=str(DEFAULT_BRIDGE))
    parser.add_argument("--public-sigma", default=str(DEFAULT_PUBLIC_SIGMA))
    parser.add_argument("--absolute-readout", default=str(DEFAULT_ABSOLUTE_READOUT))
    parser.add_argument("--public-yukawa", default=str(DEFAULT_PUBLIC_YUKAWA))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        _load_json(Path(args.bridge)),
        _load_json(Path(args.public_sigma)),
        _load_json(Path(args.absolute_readout)),
        _load_json(Path(args.public_yukawa)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
