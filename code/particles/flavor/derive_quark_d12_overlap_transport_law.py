#!/usr/bin/env python3
"""Expose the one-scalar D12 overlap transport law for the quark pure-B branch.

Chain role: collapse the D12 continuation `u/d` pure-B payload pair to one scalar
`Delta_ud_overlap` once the emitted spread totals `(sigma_u, sigma_d)` are fixed.

Mathematics: on the one-scalar D12 continuation family, the odd transport pair
is affine in `Delta_ud_overlap` with coefficients set exactly by the already
emitted spread totals:
`tau_u = sigma_d * Delta / (2 (sigma_u + sigma_d))`,
`tau_d = sigma_u * Delta / (2 (sigma_u + sigma_d))`,
`Lambda = sigma_u sigma_d * Delta / (2 (sigma_u + sigma_d))`.

OPH-derived inputs: the closed quark spread map and the current D12 continuation
mass-branch artifact.

Output: a smaller D12 continuation primitive showing that the remaining D12
mass-side scalar is `Delta_ud_overlap`, not an unconstrained tau-pair.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
D12_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _evaluate_branch(delta_value: float, sigma_u: float, sigma_d: float) -> dict[str, float]:
    denom = 2.0 * (sigma_u + sigma_d)
    tau_u = sigma_d * delta_value / denom
    tau_d = sigma_u * delta_value / denom
    lam = sigma_u * sigma_d * delta_value / denom
    return {
        "Delta_ud_overlap": float(delta_value),
        "tau_u_log_per_side": float(tau_u),
        "tau_d_log_per_side": float(tau_d),
        "Lambda_ud_B_transport": float(lam),
        "tau_sum_half_delta_identity": float((tau_u + tau_d) - (0.5 * delta_value)),
        "tau_ratio_minus_sigma_ratio": float((tau_d / tau_u) - (sigma_u / sigma_d)),
        "lambda_minus_sigma_u_tau_u": float(lam - sigma_u * tau_u),
        "lambda_minus_sigma_d_tau_d": float(lam - sigma_d * tau_d),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 quark overlap transport law artifact.")
    parser.add_argument("--spread-map", default=str(SPREAD_MAP_JSON))
    parser.add_argument("--d12-branch", default=str(D12_BRANCH_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    spread_map = _load_json(Path(args.spread_map))
    d12_branch = _load_json(Path(args.d12_branch))

    sigma_u = float(spread_map["sigma_u_total_log_per_side"])
    sigma_d = float(spread_map["sigma_d_total_log_per_side"])
    candidate_delta = float(d12_branch["candidate_mass_branch_from_t1_over_5"]["Delta_ud_overlap"])
    best_delta = float(d12_branch["best_honest_one_scalar_mass_point_on_same_family"]["Delta_ud_overlap"])

    artifact = {
        "artifact": "oph_quark_d12_overlap_transport_law",
        "generated_utc": _timestamp(),
        "status": "closed_smaller_primitive",
        "theorem_tier": "D12_continuation_only",
        "strictly_smaller_than": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "next_single_residual_object": "Delta_ud_overlap",
        "input_kind": "closed_spread_totals_plus_one_scalar_overlap_selector",
        "spread_totals": {
            "sigma_u_total_log_per_side": sigma_u,
            "sigma_d_total_log_per_side": sigma_d,
        },
        "transport_formulas": {
            "tau_u_log_per_side": "sigma_d_total_log_per_side * Delta_ud_overlap / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            "tau_d_log_per_side": "sigma_u_total_log_per_side * Delta_ud_overlap / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            "Lambda_ud_B_transport": "sigma_u_total_log_per_side * sigma_d_total_log_per_side * Delta_ud_overlap / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            "tau_u_plus_tau_d": "Delta_ud_overlap / 2",
            "tau_d_over_tau_u": "sigma_u_total_log_per_side / sigma_d_total_log_per_side",
            "Lambda_from_tau_u": "sigma_u_total_log_per_side * tau_u_log_per_side",
            "Lambda_from_tau_d": "sigma_d_total_log_per_side * tau_d_log_per_side",
        },
        "candidate_branch_from_t1_over_5": _evaluate_branch(candidate_delta, sigma_u, sigma_d),
        "best_same_family_mass_point": _evaluate_branch(best_delta, sigma_u, sigma_d),
        "notes": [
            "On the D12 one-scalar overlap family, the odd quark payload pair is not free once the spread totals are fixed.",
            "The remaining D12 mass-side scalar is therefore Delta_ud_overlap, with tau_u, tau_d, and Lambda all determined affinely from it.",
            "This does not close the D12 selector value itself, and it does not solve the CKM/CP residual generator.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
