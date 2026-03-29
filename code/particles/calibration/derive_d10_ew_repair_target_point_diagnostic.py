#!/usr/bin/env python3
"""Emit frozen-target diagnostics on the D10 repair chart.

Chain role: sharpen the D10 electroweak repair burden after the current-carrier
obstruction is known.

Mathematics: inverse chart evaluation only. This script freezes a target
`(M_W, M_Z)` pair and computes the unique repair-chart point
`(tau2_tree_exact, delta_n_tree_exact)` that would emit it on the already-fixed
carrier basis.

OPH-derived inputs: the local D10 source-pair carrier basis and the present
public/readout surface.

Output: diagnostic-only target-point artifacts showing that once the target spec
is frozen, the remaining D10 burden is the value law that emits one unique
repair point rather than an unconstrained chart search.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_target_point_diagnostic.json"

DEFAULT_REPO_PINNED = {
    "spec_id": "repo_pinned_particles_public_surface_2026_03_28",
    "source_kind": "repo_pinned_reference_values",
    "description": "Pinned pair currently shown on the local /particles public surface.",
    "MW_pole_gev": 80.377,
    "MZ_pole_gev": 91.18797809193725,
}

DEFAULT_OFFICIAL_CURRENT = {
    "spec_id": "official_pdg_current_surface_2026_03_28",
    "source_kind": "official_pdg_listing_surface",
    "description": "Current official PDG-like target pair used for D10 exactness audits.",
    "MW_pole_gev": 80.3692,
    "MZ_pole_gev": 91.188,
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _enrich_carrier(carrier: dict[str, float]) -> dict[str, float]:
    payload = dict(carrier)
    payload["beta_EW"] = (payload["alpha2_mz"] - payload["alphaY_mz"]) / (
        payload["alpha2_mz"] + payload["alphaY_mz"]
    )
    payload["alpha2_star"] = payload["alpha2_mz"]
    payload["alphaY_star"] = payload["alphaY_mz"] * (1.0 - 2.0 * payload["eta_source"])
    return payload


def n_ew_fiber(tau2: float, carrier: dict[str, float]) -> float:
    beta = carrier["beta_EW"]
    eta = carrier["eta_source"]
    return 1.0 + (
        beta * tau2 + 2.0 * (1.0 + beta) * tau2**3 - (1.0 - beta) * eta
    ) / (1.0 + 4.0 * tau2**2)


def _carrier_from_source_pair(source_pair: dict[str, Any]) -> dict[str, float]:
    slots = dict(source_pair.get("source_pair") or source_pair.get("source_slots") or {})
    return _enrich_carrier(
        {
            "alpha2_mz": float(slots["alpha2_mz"]),
            "alphaY_mz": float(slots["alphaY_mz"]),
            "eta_source": float(source_pair["eta_source"]),
            "v_inherited": float(slots["v_inherited"]),
        }
    )


def compute_target_point(spec: dict[str, Any], carrier: dict[str, float]) -> dict[str, Any]:
    mw = float(spec["MW_pole_gev"])
    mz = float(spec["MZ_pole_gev"])
    alpha2 = carrier["alpha2_mz"]
    alpha_y = carrier["alphaY_mz"]
    v_value = carrier["v_inherited"]

    tau2 = (mw * mw) / (math.pi * v_value * v_value * alpha2) - 1.0
    delta_n = (mz * mz) / (math.pi * v_value * v_value * (alpha_y + alpha2)) - n_ew_fiber(tau2, carrier)

    delta_alpha2 = alpha2 * tau2
    delta_alphaY = (
        alpha_y * (8.0 * carrier["eta_source"] * tau2 * tau2 - tau2) / (1.0 + 4.0 * tau2 * tau2)
        + (alpha_y + alpha2) * delta_n
    )
    alpha2_prime = carrier["alpha2_star"] + delta_alpha2
    alpha_y_prime = carrier["alphaY_star"] + delta_alphaY

    return {
        "artifact": "oph_d10_ew_repair_target_point_diagnostic",
        "generated_utc": _timestamp(),
        "object_id": "EWRepairTargetPointDiagnostic_D10",
        "status": "diagnostic_only_inverse_target_point",
        "spec_id": spec["spec_id"],
        "source_kind": spec["source_kind"],
        "description": spec["description"],
        "MW_pole_target_gev": mw,
        "MZ_pole_target_gev": mz,
        "tau2_tree_exact_target": tau2,
        "delta_n_tree_exact_target": delta_n,
        "delta_alpha2_tree_target": delta_alpha2,
        "delta_alphaY_tree_target": delta_alphaY,
        "alpha2_prime_target": alpha2_prime,
        "alphaY_prime_target": alpha_y_prime,
        "alpha_em_eff_inv_target": (alpha_y_prime + alpha2_prime) / (alpha_y_prime * alpha2_prime),
        "sin2w_eff_target": alpha_y_prime / (alpha_y_prime + alpha2_prime),
        "v_report_target_gev": v_value,
        "carrier_basis_scalar": carrier,
        "formulas": {
            "tau2_tree_exact_target": "MW_target^2 / (pi * v_inherited^2 * alpha2_mz) - 1",
            "n_EW_fiber": "1 + (beta_EW*tau2 + 2*(1+beta_EW)*tau2^3 - (1-beta_EW)*eta_source) / (1 + 4*tau2^2)",
            "delta_n_tree_exact_target": "MZ_target^2 / (pi * v_inherited^2 * (alphaY_mz + alpha2_mz)) - n_EW_fiber(tau2_tree_exact_target)",
            "delta_alpha2_tree_target": "alpha2_mz * tau2_tree_exact_target",
            "delta_alphaY_tree_target": "alphaY_mz * (8*eta_source*tau2^2 - tau2)/(1+4*tau2^2) + (alphaY_mz + alpha2_mz) * delta_n_tree_exact_target",
        },
        "notes": [
            "This artifact is intentionally diagnostic only. It does not derive the target point from OPH.",
            "Once the D10 target specification is frozen, the repair chart determines one unique target point.",
            "The missing predictive object is still the OPH repair value law that emits this point without inverse fitting.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute the frozen-target diagnostic point on the D10 repair chart.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--target-spec", choices=["repo_pinned", "official_current"], default="official_current")
    parser.add_argument("--target-json", default="")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = _load_json(Path(args.source_pair))
    carrier = _carrier_from_source_pair(source_pair)
    spec = (
        _load_json(Path(args.target_json))
        if args.target_json
        else DEFAULT_REPO_PINNED
        if args.target_spec == "repo_pinned"
        else DEFAULT_OFFICIAL_CURRENT
    )

    payload = compute_target_point(spec, carrier)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
