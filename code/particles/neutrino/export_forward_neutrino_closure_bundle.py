#!/usr/bin/env python3
"""Bundle the active neutrino lane with transitive, fail-closed status gates.

Chain role: collect the weighted-cycle PMNS/hierarchy candidate together with
the bridge-rigidity and absolute-attachment artifacts into the forward bundle
used by audits and public-surface gating.

Mathematics: packaging only. Promotion is the conjunction of upstream gates;
missing status fields fail closed.

Inputs: the weighted-cycle candidate, bridge-rigidity artifact,
absolute-attachment artifact, and correlated-profile score.

Output: the forward neutrino closure bundle for downstream reporting.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WEIGHTED_CYCLE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_BRIDGE_RIGIDITY = ROOT / "particles" / "runs" / "neutrino" / "neutrino_bridge_rigidity_theorem.json"
DEFAULT_ABSOLUTE_ATTACHMENT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_theorem.json"
DEFAULT_PROFILE_SCORE = ROOT / "particles" / "runs" / "neutrino" / "nufit61_weighted_cycle_retrospective_score.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the forward neutrino closure-bundle artifact.")
    parser.add_argument("--weighted-cycle", default=str(DEFAULT_WEIGHTED_CYCLE))
    parser.add_argument("--bridge-rigidity", default=str(DEFAULT_BRIDGE_RIGIDITY))
    parser.add_argument("--absolute-attachment", default=str(DEFAULT_ABSOLUTE_ATTACHMENT))
    parser.add_argument("--profile-score", default=str(DEFAULT_PROFILE_SCORE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    weighted_cycle = json.loads(Path(args.weighted_cycle).read_text(encoding="utf-8"))
    bridge_rigidity = json.loads(Path(args.bridge_rigidity).read_text(encoding="utf-8"))
    absolute_attachment = json.loads(Path(args.absolute_attachment).read_text(encoding="utf-8"))
    profile_score_path = Path(args.profile_score)
    profile_score = json.loads(profile_score_path.read_text(encoding="utf-8")) if profile_score_path.exists() else None

    weighted_allowed = bool(weighted_cycle.get("prediction_promotion_allowed", False))
    bridge_allowed = bool(bridge_rigidity.get("prediction_promotion_allowed", False))
    absolute_allowed = bool(absolute_attachment.get("prediction_promotion_allowed", False))
    profile_rejected = bool(
        profile_score
        and profile_score.get("decision", {}).get("current_weighted_cycle_candidate_rejected_by_declared_gate")
    )
    promotion_allowed = weighted_allowed and bridge_allowed and absolute_allowed and not profile_rejected
    absolute_outputs = dict(absolute_attachment.get("outputs") or {})
    absolute_dm2 = dict(absolute_outputs.get("delta_m_sq_eV2") or {})
    absolute_masses = list(absolute_outputs.get("masses_eV") or [])

    payload = {
        "artifact": "oph_forward_neutrino_closure_bundle",
        "generated_utc": _timestamp(),
        "closure_tier": "retrospective_target_informed_weighted_cycle_candidate",
        "prediction_promotion_allowed": promotion_allowed,
        "public_surface_candidate_allowed": promotion_allowed,
        "display_allowed_as_compare_only": True,
        "phase_mode": "weighted_cycle_candidate",
        "selector_law_certified": False,
        "certification_status": "candidate_rejected_and_source_closure_failed" if profile_rejected else "candidate_source_closure_failed",
        "transitive_promotion_gate": {
            "weighted_cycle_prediction_promotion_allowed": weighted_allowed,
            "bridge_prediction_promotion_allowed": bridge_allowed,
            "absolute_attachment_prediction_promotion_allowed": absolute_allowed,
            "correlated_profile_rejected": profile_rejected,
            "all_required": True,
        },
        "weighted_cycle_branch": {
            "D_nu": weighted_cycle["selected_D_nu"],
            "p_nu": weighted_cycle.get("selected_p_nu", weighted_cycle["weight_exponent"]),
            "pmns_observables": weighted_cycle["pmns_observables"],
            "dimensionless_ratio_dm21_over_dm32": weighted_cycle["dimensionless_ratio_dm21_over_dm32"],
            "dimensionless_masses": weighted_cycle["dimensionless_masses"],
            "dimensionless_dm2": weighted_cycle["dimensionless_dm2"],
        },
        "bridge_rigidity": {
            "artifact": bridge_rigidity["artifact"],
            "C_nu": bridge_rigidity["emitted_value"],
            "P_nu": bridge_rigidity["emitted_proxy"]["value"],
            "formula": bridge_rigidity["emitted_formula"],
        },
        "absolute_attachment": absolute_outputs,
        "absolute_attachment_status": "compare_only_blocked",
        "masses_gev_sorted": [value * 1.0e-9 for value in absolute_masses],
        "delta_m21_sq_gev2": float(absolute_dm2.get("21", 0.0)) * 1.0e-18,
        "delta_m31_sq_gev2": float(absolute_dm2.get("31", 0.0)) * 1.0e-18,
        "delta_m32_sq_gev2": float(absolute_dm2.get("32", 0.0)) * 1.0e-18,
        "legacy_absolute_fields_are_compare_only": True,
        "splitting_ratio_r": weighted_cycle["dimensionless_ratio_dm21_over_dm32"],
        "ordering_phase_certified": False,
        "pmns_status": "target_informed_template_candidate_rejected_by_nufit61_profile" if profile_rejected else "target_informed_template_candidate",
        "profile_score": (
            {
                "artifact": profile_score.get("artifact"),
                "current_weighted_cycle_candidate_rejected_by_declared_gate": profile_rejected,
                "TByes_NO_T23_DCP_delta_chi2": profile_score["scores"]["TByes-NO"]["profiles"]["T23/DCP"]["delta_chi2"],
                "TBoff_NO_T23_DCP_delta_chi2": profile_score["scores"]["TBoff-NO"]["profiles"]["T23/DCP"]["delta_chi2"],
            }
            if profile_score
            else None
        ),
        "notes": [
            "The weighted-cycle values are exact linear-algebra outputs conditional on a target-informed template branch.",
            "The NuFIT 6.1 correlated profile rejects the candidate under the declared two-parameter 3-sigma gate.",
            "The absolute attachment and legacy absolute fields are compare-only and cannot feed a prediction table.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
