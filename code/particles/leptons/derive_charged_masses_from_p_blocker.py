#!/usr/bin/env python3
"""Emit the current-corpus dependency blocker for charged masses from P."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from charged_absolute_route_common import (
    P_TO_AFFINE_REDUCTION_JSON,
    POST_PROMOTION_ROUTE_JSON,
    artifact_ref,
    load_json,
)


ROOT = Path(__file__).resolve().parents[2]
BRIDGE_NO_GO_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_p_to_affine_anchor_bridge_no_go.json"
BRIDGE_REDUCTION_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_p_to_affine_anchor_reduction.json"
END_TO_END_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_masses_from_p_blocker.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_dict(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key)
    return dict(value) if isinstance(value, dict) else {}


def build_artifact(
    bridge_no_go: dict[str, Any],
    bridge_reduction: dict[str, Any],
    end_to_end: dict[str, Any],
    post_promotion_route: dict[str, Any],
) -> dict[str, Any]:
    exact_chain = list(post_promotion_route.get("exact_irreducible_chain", []))
    operator_gate = exact_chain[0] if len(exact_chain) > 0 else {}
    affine_slot = exact_chain[1] if len(exact_chain) > 1 else {}

    return {
        "artifact": "oph_charged_masses_from_p_blocker",
        "generated_utc": _timestamp(),
        "verdict": "charged_masses_from_P_not_emitted_on_current_corpus",
        "status": "dependency_blocker_certificate",
        "public_promotion_allowed": False,
        "theorem_statement": (
            "No theorem-grade charged masses are emitted from P on the current corpus. The operator-side "
            "promotion theorem is still open, the post-promotion affine slot is still open, and no theorem-grade "
            "bridge from the D10 calibration descendants of P to A_ch has been derived. Therefore A_ch(P), "
            "g_e(P), and the charged mass triple from P remain un-emitted."
        ),
        "failed_dependencies": [
            {
                "id": operator_gate.get("id"),
                "smallest_missing_clause": operator_gate.get("smallest_missing_clause"),
                "effect_on_fill": operator_gate.get("effect_on_fill"),
            },
            {
                "id": affine_slot.get("id"),
                "exact_descended_scalar": affine_slot.get("exact_descended_scalar"),
                "exact_smaller_forcing_object": affine_slot.get("exact_smaller_forcing_object"),
                "effect_on_fill": affine_slot.get("effect_on_fill"),
            },
            {
                "id": bridge_no_go.get("exact_missing_object", {}).get("id"),
                "smallest_honest_contract": bridge_no_go.get("exact_missing_object", {}).get(
                    "smallest_honest_contract"
                ),
                "exact_smallest_bridge_target": bridge_reduction.get("exact_smallest_bridge_target", {}).get(
                    "id"
                ),
                "P_threaded_reformulation": "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu",
                "effect_on_fill": "theorem_grade_charged_determinant_line_section(P)_then_A_ch(P)",
            },
        ],
        "bridge_gap_artifact": {
            "artifact": bridge_no_go.get("artifact"),
            "artifact_ref": artifact_ref(BRIDGE_NO_GO_JSON),
            "verdict": bridge_no_go.get("verdict"),
        },
        "bridge_reduction_artifact": {
            "artifact": bridge_reduction.get("artifact"),
            "artifact_ref": artifact_ref(BRIDGE_REDUCTION_JSON),
            "reduction_theorem_id": bridge_reduction.get("reduction_theorem", {}).get("id"),
            "exact_smallest_bridge_target": bridge_reduction.get("exact_smallest_bridge_target", {}).get("id"),
            "P_threaded_reformulation_id": bridge_reduction.get(
                "P_threaded_source_exactness_reformulation", {}
            ).get("id"),
        },
        "current_forbidden_outputs": [
            "A_ch(P)",
            "g_e(P)",
            "m_e(P)",
            "m_mu(P)",
            "m_tau(P)",
        ],
        "symbolic_readout_if_dependencies_close": {
            "g_e(P)": "exp(A_ch(P))",
            "m_e(P)": "exp(A_ch(P) + e_log_centered(P))",
            "m_mu(P)": "exp(A_ch(P) + mu_log_centered(P))",
            "m_tau(P)": "exp(A_ch(P) + tau_log_centered(P))",
        },
        "alignment_with_global_no_go": {
            "artifact": end_to_end.get("artifact"),
            "verdict": end_to_end.get("verdict"),
            "artifact_ref": artifact_ref(END_TO_END_JSON),
        },
        "notes": [
            "This is a dependency theorem, not a numerical mismatch report.",
            "It isolates exactly what remains before any public charged table can emit masses from P.",
            "The symbolic formulas are already known; what is missing is the theorem-grade production of A_ch(P).",
            "The charged P-bridge can be phrased more sharply as theorem-grade source landing on physical Y_e(P), after which the affine scalar follows from the physical descent package.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current-corpus blocker for charged masses from P.")
    parser.add_argument("--bridge-no-go", default=str(BRIDGE_NO_GO_JSON))
    parser.add_argument("--bridge-reduction", default=str(BRIDGE_REDUCTION_JSON))
    parser.add_argument("--end-to-end", default=str(END_TO_END_JSON))
    parser.add_argument("--post-promotion-route", default=str(POST_PROMOTION_ROUTE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.bridge_no_go)),
        load_json(Path(args.bridge_reduction)),
        load_json(Path(args.end_to_end)),
        load_json(Path(args.post_promotion_route)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
