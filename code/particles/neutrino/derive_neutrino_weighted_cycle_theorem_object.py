#!/usr/bin/env python3
"""Emit the explicit law object behind the weighted-cycle candidate.

Chain role: serialize the retrospective dimensionless weighted-cycle law while
preserving its non-theorem status.

Mathematics: the repaired branch fixes the positive load segment between
`chi = 1 + eps` and `1 + gamma_half`. On that one-dimensional affine segment,
the balanced selector and the least-distortion selector coincide uniquely at
the midpoint, so

  D_nu = (chi + 1 + gamma_half) / 2
  p_nu(gamma, eps) = 1 + gamma + eps / D_nu

The midpoint statement is elementary mathematics conditional on a declared
segment.  It does not derive that segment, the exponent formula, the cycle
operator, or its physical basis placement.  Git history also records target
ranking before this law was promoted.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_SELECTOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json"
DEFAULT_REPAIR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_theorem_object.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Emit the historical neutrino weighted-cycle candidate law audit.")
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE))
    parser.add_argument("--selector", default=str(DEFAULT_SELECTOR))
    parser.add_argument("--repair", default=str(DEFAULT_REPAIR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    cocycle = _load_json(Path(args.cocycle))
    selector = _load_json(Path(args.selector))
    repair = _load_json(Path(args.repair))

    gamma = float(cocycle["theorem_gap_gamma"])
    eps = float(cocycle["defect_gap_ratio"])
    gamma_half = float(cocycle["hermitian_descendant_riesz_margin"]["gamma_half"])
    chi = float(selector["derived_quantities"]["chi"])
    d_nu = float(selector["selected_D_nu"])
    p_nu = float(selector["derived_quantities"]["weight_exponent_value"])

    source_closed = bool((repair.get("source_closure_status") or {}).get("closed", False))
    promotion_allowed = (
        source_closed
        and repair.get("source_only_prediction_eligible") is True
        and repair.get("prediction_promotion_allowed") is True
        and repair.get("historical_target_exposure") is False
    )

    payload = {
        "artifact": "oph_neutrino_weighted_cycle_theorem_object",
        "generated_utc": _timestamp(),
        "status": "retrospective_weighted_cycle_candidate_law",
        "theorem_status": "not_established",
        "source_closure_closed": source_closed,
        "public_surface_candidate_allowed": False,
        "prediction_promotion_allowed": promotion_allowed,
        "source_artifacts": {
            "overlap_edge_transport_cocycle": str(Path(args.cocycle)),
            "transport_load_selector": str(Path(args.selector)),
            "weighted_cycle_repair": str(Path(args.repair)),
        },
        "candidate_law": {
            "name": "p_nu_gamma_eps_midpoint_law",
            "statement": (
                "On the positive affine load segment between chi = 1 + eps and 1 + gamma_half, "
                "the balanced selector and the least-distortion selector coincide uniquely at the midpoint."
            ),
            "selector_family": "D_tau = (1-tau_nu) * chi + tau_nu * (1 + gamma_half)",
            "selected_tau_nu": 0.5,
            "chi_nu_formula": "chi_nu(eps) = 1 + eps",
            "D_nu_formula": "D_nu = (chi_nu + 1 + gamma_half) / 2 = 1 + eps/2 + gamma/4",
            "p_nu_formula": "p_nu(gamma, eps) = 1 + gamma + eps / D_nu",
            "w_e_formula": "w_e = q_e ** p_nu",
        },
        "live_inputs": {
            "gamma": gamma,
            "eps": eps,
            "gamma_half": gamma_half,
            "chi_nu": chi,
            "D_nu": d_nu,
            "p_nu": p_nu,
        },
        "live_outputs": {
            "dimensionless_ratio_dm21_over_dm32": float(repair["dimensionless_ratio_dm21_over_dm32"]),
            "pmns_observables": dict(repair["pmns_observables"]),
        },
        "remaining_open_object": {
            "name": "source_derived_weighted_cycle_operator_and_basis_map",
            "status": "open",
            "statement": (
                "A source-emitted operator, exponent law, basis placement, and pre-reference lock are required "
                "before this construction can become a neutrino prediction."
            ),
        },
        "audit": {
            "midpoint_fact_valid_conditional_on_declared_segment": True,
            "segment_endpoints_derived_for_neutrino_transport": False,
            "exponent_formula_derived": False,
            "cycle_operator_derived": False,
            "physical_basis_placement_derived": False,
            "historical_target_exposure": bool(repair.get("historical_target_exposure", True)),
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
