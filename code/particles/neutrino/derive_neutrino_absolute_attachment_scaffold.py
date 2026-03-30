#!/usr/bin/env python3
"""Emit the neutrino absolute-attachment scaffold.

This does not emit ``lambda_nu``. It sharpens the live neutrino frontier:
the remaining theorem object is an internal attachment from the D10 amplitude
sector to the weighted-cycle scale-free normal form.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BRIDGE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_absolute_amplitude_bridge.json"
BRIDGE_CANDIDATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"
THEOREM_OBJECT_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_theorem_object.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(bridge: dict[str, Any], bridge_candidate: dict[str, Any], theorem_object: dict[str, Any]) -> dict[str, Any]:
    diagnostic = dict(bridge["direct_scale_anchor_attachment_diagnostic"])
    return {
        "artifact": "oph_neutrino_absolute_attachment_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_constructive_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": "neutrino_weighted_cycle_absolute_attachment",
        "equivalent_scalar": {
            "name": "lambda_nu",
            "meaning": "positive absolute normalization scalar for the weighted-cycle scale-free normal form",
        },
        "current_no_go": {
            "statement": bridge["no_go_statement"],
            "direct_attachment_diagnostic": diagnostic,
            "current_candidate_interface_artifact": bridge_candidate.get("current_candidate_interface_artifact"),
            "exact_next_theorem_object": bridge_candidate.get("exact_next_theorem_object"),
            "strictly_smaller_missing_clause": bridge_candidate.get("strictly_smaller_missing_clause"),
        },
        "extension_contract": {
            "input_objects": [
                "oph_neutrino_weighted_cycle_theorem_object",
                "oph_neutrino_scale_anchor",
                "oph_neutrino_family_response",
                "oph_forward_majorana_matrix",
            ],
            "forbidden_inputs": [
                "external_oscillation_anchors",
                "PDG_target_backsolve",
                "PMNS_target_seed",
            ],
            "must_emit": "lambda_nu > 0 or an exactly equivalent amplitude attachment A_nu",
            "must_imply": [
                "m_i = lambda_nu * mhat_i",
                "Delta m^2_ij = lambda_nu^2 * Delta_hat_ij",
            ],
            "bridge_statement": (
                "Attach the internal D10 amplitude sector to the weighted-cycle scale-free normal "
                "form without reusing external oscillation anchors."
            ),
            "current_theorem_stack": bridge_candidate.get("bridge_interface_theorem_stack", []),
        },
        "theorem_object_context": {
            "name": theorem_object["theorem_object"]["name"],
            "D_nu_formula": theorem_object["theorem_object"]["D_nu_formula"],
            "p_nu_formula": theorem_object["theorem_object"]["p_nu_formula"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the neutrino absolute-attachment scaffold.")
    parser.add_argument("--bridge", type=Path, default=BRIDGE_JSON)
    parser.add_argument("--bridge-candidate", type=Path, default=BRIDGE_CANDIDATE_JSON)
    parser.add_argument("--theorem-object", type=Path, default=THEOREM_OBJECT_JSON)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    bridge = _load_json(args.bridge)
    bridge_candidate = _load_json(args.bridge_candidate)
    theorem_object = _load_json(args.theorem_object)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(build_artifact(bridge, bridge_candidate, theorem_object), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
