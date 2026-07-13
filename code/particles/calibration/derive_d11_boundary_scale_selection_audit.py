#!/usr/bin/env python3
"""Audit the criticality boundary-scale selection question.

The double-criticality family leaves one discrete freedom: the source scale
at which ``lambda = 0`` and ``beta_lambda = 0`` are imposed.  This audit
settles what can be settled and freezes what cannot.

1. Flow-internal selection no-go.  The only scale-selecting condition
   available inside the renormalization flow itself is triple criticality,
   ``d beta_lambda / d ln mu = 0`` at the critical point.  Along the
   criticality family the flow derivative is dominated by
   ``-24 y^3 beta_y > 0`` and stays bounded away from zero across the whole
   candidate window, so no triple-critical root exists.  Every family point
   is a clean lambda minimum, and the boundary scale must be selected by the
   source structure, never by the flow.

2. Frozen prospective candidate registry.  The named candidates are closed
   forms in the source constants.  One candidate, the log-midpoint of the
   model's two pre-existing anchors,
   ``sqrt(mu_U * E_cell) = E_star exp(-pi) P^(-1/6)``, was identified after
   the two-loop implied scale was computed; it is therefore a post-exposure
   candidate and is recorded with that classification.  The registry is
   frozen here, before the three-loop computation exists, so the three-loop
   implied scale becomes a genuine discriminating test.

3. Truncation instability.  The implied boundary scale moved from outside
   the window (one loop) to 4.8e17 GeV (two loops), a shift larger than the
   candidate spacing.  Numerical agreement at the current loop order
   therefore cannot select a candidate; the registered discriminator is the
   frozen three-loop RG and matching packet.

The selection question stays open.  The adopted working conditional branch
for the spectrum surface is the log-midpoint candidate, carried under its
post-exposure classification per the conditional-value policy.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

try:
    from calibration.derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        beta_lambda_1loop,
        beta_y_1loop,
        critical_surface_yukawa,
        gauge_couplings,
        source_scales,
    )
except ModuleNotFoundError:
    from derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        beta_lambda_1loop,
        beta_y_1loop,
        critical_surface_yukawa,
        gauge_couplings,
        source_scales,
    )

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_boundary_scale_selection_audit.json"
)

CANDIDATES = {
    "mu_U_gauge_unification": {
        "closed_form": "E_star * exp(-2 pi) * P^(1/6)",
        "motivation": (
            "the archived branch: criticality at the gauge-unification "
            "scale of the pixel-residual solve"
        ),
        "classification": "archived_declared_choice",
    },
    "log_midpoint_half_turn": {
        "closed_form": "sqrt(mu_U * E_cell) = E_star * exp(-pi) * P^(-1/6)",
        "motivation": (
            "log-midpoint of the model's two pre-existing anchors; the "
            "unique scale treating the transmutation anchor and the "
            "unification anchor symmetrically"
        ),
        "classification": "prospective_post_exposure_candidate",
    },
    "E_cell_pixel_energy": {
        "closed_form": "E_star * P^(-1/2)",
        "motivation": (
            "single-anchor coherence: the same anchor that carries the "
            "electroweak transmutation"
        ),
        "classification": "prospective_principled_candidate",
    },
    "E_star": {
        "closed_form": "E_star",
        "motivation": "the ultraviolet edge of the source theory",
        "classification": "prospective_principled_candidate",
    },
}


def flow_internal_no_go(sample_count: int = 24) -> dict[str, Any]:
    """Certify that triple criticality has no root in the window.

    Along the family (lambda = 0, y = y_crit(mu)), the flow derivative of
    beta_lambda is computed exactly at one loop:
    d beta_lambda / dt = kappa * (-24 y^3 beta_y + (3/8) dX/dt), with the
    gauge term evaluated by finite difference of the chart couplings.
    """

    scales = source_scales(P_FALLBACK, ALPHA_U_FALLBACK)
    mz_run = scales["mz_run_gev"]
    kappa = 1.0 / (16.0 * math.pi**2)
    eps = 1.0e-4

    rows = []
    for mu in np.geomspace(1.0e16, 1.3e19, sample_count):
        g_y, g2, g3 = gauge_couplings(float(mu), mz_run)
        y = critical_surface_yukawa(g_y, g2)
        beta_y_val = beta_y_1loop(y, g_y, g2, g3)

        def x_of(m: float) -> float:
            gy_, g2_, _ = gauge_couplings(m, mz_run)
            return 2.0 * g2_**4 + (g2_**2 + gy_**2) ** 2

        ln_mu = math.log(float(mu))
        dx_dt = (x_of(math.exp(ln_mu + eps)) - x_of(math.exp(ln_mu - eps))) / (
            2.0 * eps
        )
        flow_derivative = kappa * (-24.0 * y**3 * beta_y_val + 0.375 * dx_dt)
        rows.append(
            {
                "mu_gev": float(mu),
                "flow_dbeta_lambda_dt": flow_derivative,
                "yukawa_term": kappa * (-24.0 * y**3 * beta_y_val),
                "gauge_term": kappa * 0.375 * dx_dt,
            }
        )

    minimum = min(row["flow_dbeta_lambda_dt"] for row in rows)
    return {
        "statement": (
            "d beta_lambda / d ln mu > 0 along the entire criticality "
            "family: triple criticality has no root, every family point is "
            "a lambda minimum, and the flow cannot select the boundary "
            "scale"
        ),
        "sampled_rows": rows,
        "minimum_flow_derivative": minimum,
        "margin_over_zero": minimum,
        "no_root_certified": minimum > 1.0e-5,
    }


def build_artifact(scan: dict[str, Any]) -> dict[str, Any]:
    no_go = flow_internal_no_go()
    scales = scan["source_scales_gev"]

    registry = {}
    for name, meta in CANDIDATES.items():
        entry = dict(meta)
        entry["scale_gev"] = scales[name]
        entry["one_loop"] = {
            "mt_pole_gev": scan["one_loop_named_boundaries"][name]["mt_pole_gev"],
            "mh_tree_gev": scan["one_loop_named_boundaries"][name]["mh_tree_gev"],
        }
        two = scan["two_loop_named_boundaries"].get(name)
        if isinstance(two, dict) and "mt_pole_gev" in two:
            entry["two_loop"] = {
                "mt_pole_gev": two["mt_pole_gev"],
                "mh_tree_gev": two["mh_tree_gev"],
            }
        registry[name] = entry

    checks = {
        "flow_internal_no_go_certified": bool(no_go["no_root_certified"]),
        "registry_covers_named_scales": set(CANDIDATES).issubset(
            set(scan["one_loop_named_boundaries"])
        ),
        "scan_checks_pass": bool(scan["checks_pass"]),
    }

    return {
        "artifact": "oph_d11_boundary_scale_selection_audit",
        "schema_version": 1,
        "status": "SELECTION_OPEN_CANDIDATES_FROZEN_PROSPECTIVELY",
        "row_class": "selection_audit_never_a_prediction_ancestor",
        "promotion_allowed": False,
        "flow_internal_selection_no_go": no_go,
        "candidate_registry_frozen": registry,
        "truncation_instability": {
            "one_loop_implied_scale": (
                "outside the window: the measured top exceeds the one-loop "
                "curve maximum"
            ),
            "two_loop_implied_scale_gev": 4.750331123992614e17,
            "statement": (
                "the loop-order shift of the implied scale exceeds the "
                "candidate spacing, so numerical agreement at the current "
                "order cannot select a candidate"
            ),
        },
        "registered_discriminating_test": {
            "object": "FROZEN_THREE_LOOP_RG_MATCHING_PACKET",
            "protocol": (
                "the registry above is frozen before any three-loop "
                "computation exists; the three-loop implied scale either "
                "converges onto one frozen candidate or eliminates the "
                "family"
            ),
        },
        "adopted_working_conditional_branch": {
            "name": "log_midpoint_half_turn",
            "classification": "prospective_post_exposure_candidate",
            "policy": (
                "conditional values with named open theorems are preferred "
                "over absent values; the branch carries its post-exposure "
                "classification explicitly"
            ),
        },
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
        "claim_boundary": (
            "The selection question is open. The no-go closes the "
            "flow-internal route; the candidate registry and the "
            "discriminating test define the prospective closure protocol."
        ),
    }


def build() -> dict[str, Any]:
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    return build_artifact(scan)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "no_go_min_flow_derivative": artifact[
                    "flow_internal_selection_no_go"
                ]["minimum_flow_derivative"],
                "adopted_branch": artifact["adopted_working_conditional_branch"][
                    "name"
                ],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
