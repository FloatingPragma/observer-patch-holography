#!/usr/bin/env python3
"""Emit the neutrino closure summary with transitive status gates.

Chain role: collect the weighted-cycle candidate boundary, the exact
irreducibility/no-go results, and the bridge/absolute-attachment status.

Mathematics:
1. The weighted-cycle matrix fixes a scale-free PMNS/hierarchy point once its
   template and declared selectors are supplied.
2. The bridge rigidity artifact contains a compare-only reduced invariant.
3. The absolute attachment remains compare-only.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
WEIGHTED_CYCLE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
BRIDGE_RIGIDITY_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_bridge_rigidity_theorem.json"
ABSOLUTE_ATTACHMENT_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lane_closure_contract.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(
    weighted_cycle: dict[str, Any],
    bridge_rigidity: dict[str, Any],
    absolute_attachment: dict[str, Any],
) -> dict[str, Any]:
    pmns = weighted_cycle["pmns_observables"]
    ratio = float(weighted_cycle["dimensionless_ratio_dm21_over_dm32"])
    scale_free_promotable = (
        weighted_cycle.get("prediction_promotion_allowed") is True
        and weighted_cycle.get("source_only_prediction_eligible") is True
        and weighted_cycle.get("historical_target_exposure") is False
    )
    absolute_promotable = (
        absolute_attachment.get("status") == "theorem_grade_emitted"
        and absolute_attachment.get("public_surface_candidate_allowed") is True
        and (absolute_attachment.get("non_circularity_status") or {}).get("promotion_allowed") is True
    )
    c_nu_value = bridge_rigidity.get("emitted_value")
    if c_nu_value is None:
        c_nu_value = bridge_rigidity.get("display_value")

    return {
        "artifact": "oph_neutrino_lane_closure_contract",
        "generated_utc": _timestamp(),
        "scope": "weighted_cycle_bridge_rigidity_plus_absolute_attachment",
        "proof_status": (
            "source_closed_neutrino_prediction"
            if scale_free_promotable and absolute_promotable
            else "target_informed_weighted_cycle_candidate_with_compare_only_absolute_attachment"
        ),
        "public_promotion_allowed": scale_free_promotable and absolute_promotable,
        "scale_free_prediction_promotion_allowed": scale_free_promotable,
        "source_closure_status": weighted_cycle.get("source_closure_status"),
        "historical_target_exposure": bool(weighted_cycle.get("historical_target_exposure", True)),
        "non_circularity_status": absolute_attachment.get("non_circularity_status"),
        "current_branch_status": {
            "branch": "weighted_cycle_majorana_holonomy",
            "pmns_observables": pmns,
            "dimensionless_ratio_dm21_over_dm32": ratio,
            "no_hidden_discrete_branch": False,
            "status": "target_informed_template_candidate",
            "absolute_family": "m_i = lambda_nu * mhat_i, Delta m^2_ij = lambda_nu^2 * Delta_hat_ij",
        },
        "emitted_bridge_rigidity_theorem": {
            "artifact": bridge_rigidity["artifact"],
            "status": bridge_rigidity["status"],
            "statement": bridge_rigidity["statement"],
            "emitted_formula": bridge_rigidity["emitted_formula"],
            "emitted_value": bridge_rigidity["emitted_value"],
            "display_value": c_nu_value,
            "P_nu": bridge_rigidity["emitted_proxy"]["value"],
        },
        "emitted_absolute_attachment_theorem": {
            "artifact": absolute_attachment["artifact"],
            "status": absolute_attachment["status"],
            "public_surface_candidate_allowed": absolute_attachment["public_surface_candidate_allowed"],
            "B_nu": absolute_attachment["outputs"]["B_nu"],
            "lambda_nu": absolute_attachment["outputs"]["lambda_nu"],
            "masses_eV": absolute_attachment["outputs"]["masses_eV"],
            "delta_m_sq_eV2": absolute_attachment["outputs"]["delta_m_sq_eV2"],
        },
        "closure_chain": [
            "declared weighted-cycle candidate => (mhat_i, Delta_hat_ij, U_PMNS)",
            "compare-only bridge candidate => C_nu = sum_gap^2 * prod_qbar * solar_response_over_mstar^-0.5",
            "B_nu = P_nu * C_nu",
            "lambda_nu = (m_star_eV / q_mean^p_nu) * P_nu * C_nu",
            "m_i = lambda_nu * mhat_i and Delta m^2_ij = lambda_nu^2 * Delta_hat_ij",
        ],
        "notes": [
            (
                "The compare-only continuation adapter is retired from the proof-facing neutrino lane."
                if absolute_promotable
                else "The absolute attachment is displayed as a compare-only candidate; the scale-free weighted-cycle point also fails source and historical-blindness gates."
            ),
            "The bridge corridor and residual correction audits remain diagnostic-only surfaces beneath the weighted-cycle candidate.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact neutrino closure contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        _load_json(WEIGHTED_CYCLE_JSON),
        _load_json(BRIDGE_RIGIDITY_JSON),
        _load_json(ABSOLUTE_ATTACHMENT_JSON),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
