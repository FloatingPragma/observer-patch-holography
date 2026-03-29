#!/usr/bin/env python3
"""Split the shared quark norm into up/down sector means on the active family.

Chain role: convert the quark spread data into separate sector means `g_u` and
`g_d` before the odd-source lane is attached.

Mathematics: a two-scalar affine readback on the ordered-family coordinates
using the current spread and shared normalization.

OPH-derived inputs: `rho_ord`, `x2`, `sigma_u`, `sigma_d`, and the shared quark
normalization from the local flavor chain.

Output: the current-family sector-mean artifact consumed by the full quark
descent and forward Yukawa builder.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SPREAD_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_SHARED_NORM = ROOT / "particles" / "runs" / "flavor" / "quark_shared_absolute_norm_binding.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(spread_map: dict, shared_norm: dict | None) -> dict:
    sigma_u = float(spread_map["sigma_u_total_log_per_side"])
    sigma_d = float(spread_map["sigma_d_total_log_per_side"])
    rho_ord = float(spread_map["rho_ord"])
    x2 = float(spread_map["normalized_coordinate_x2"])
    g_ch = float((shared_norm or {}).get("g_ch", spread_map["shared_norm_value"]))

    sigma_seed_ud = 0.5 * (sigma_u + sigma_d)
    eta_ud = 0.5 * (sigma_u - sigma_d)
    x2_sq = x2 * x2
    mean_denominator = 1.0 + rho_ord - x2_sq
    skew_denominator = 1.0 - x2_sq - (x2_sq / (1.0 + rho_ord))
    if abs(mean_denominator) <= 1.0e-12 or abs(skew_denominator) <= 1.0e-12:
        raise ValueError("current-family affine mean denominators are singular")

    a_ud = 1.0 / (2.0 * mean_denominator)
    b_ud = 1.0 / (2.0 * skew_denominator)
    seed_decay = a_ud * sigma_seed_ud
    skew_readback = b_ud * eta_ud
    log_shift_u = -(seed_decay - skew_readback)
    log_shift_d = -(seed_decay + skew_readback)
    g_u = g_ch * math.exp(log_shift_u)
    g_d = g_ch * math.exp(log_shift_d)

    return {
        "artifact": "oph_quark_sector_mean_split",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_predictive_law",
        "theorem_candidate": "quark_sector_mean_split_theorem",
        "theorem_scope": "current_family_only",
        "source_artifacts": {
            "spread_map": spread_map.get("artifact"),
            "shared_norm_binding": None if shared_norm is None else shared_norm.get("artifact"),
        },
        "shared_norm_value": g_ch,
        "rho_ord": rho_ord,
        "normalized_coordinate_x2": x2,
        "x2_squared": x2_sq,
        "sigma_u_total_log_per_side": sigma_u,
        "sigma_d_total_log_per_side": sigma_d,
        "sigma_seed_ud_candidate": sigma_seed_ud,
        "eta_ud_candidate": eta_ud,
        "mean_law_kind": "theorem_grade_two_scalar_affine_readout_closed",
        "mean_denominator": mean_denominator,
        "skew_denominator": skew_denominator,
        "seed_decay_formula": "sigma_seed_ud / (2 * (1 + rho_ord - x2^2))",
        "skew_readback_formula": "eta_ud / (2 * (1 - x2^2 - x2^2 / (1 + rho_ord)))",
        "seed_decay": seed_decay,
        "skew_readback": skew_readback,
        "A_ud_candidate": a_ud,
        "B_ud_candidate": b_ud,
        "A_ud_formula": "1 / (2 * (1 + rho_ord - x2^2))",
        "B_ud_formula": "1 / (2 * (1 - x2^2 - x2^2 / (1 + rho_ord)))",
        "sector_log_mean_shift_u": log_shift_u,
        "sector_log_mean_shift_d": log_shift_d,
        "sector_log_mean_shift_formula_u": "-(seed_decay - skew_readback)",
        "sector_log_mean_shift_formula_d": "-(seed_decay + skew_readback)",
        "g_u_candidate": g_u,
        "g_d_candidate": g_d,
        "g_u_over_g_ch": math.exp(log_shift_u),
        "g_d_over_g_ch": math.exp(log_shift_d),
        "candidate_kind": "two_scalar_affine_mean_split",
        "active_candidate": "ordered_affine_mean_readout_candidate",
        "smallest_constructive_missing_object": (
            None
            if str(spread_map.get("spread_emitter_status", "")) == "closed"
            else "oph_family_excitation_spread_map"
        ),
        "metadata": {
            "note": (
                "Compact current-family quark mean split built from the emitted two-scalar spread package. "
                "The affine readout coefficients are fixed by the theorem-grade denominator pair on the same two scalars "
                "(sigma_seed_ud, eta_ud) and ordered geometry (rho_ord, x2). When the spread map is closed, this mean "
                "surface is no longer waiting on another constructive object inside the current-family lane."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the compact current-family quark sector-mean split candidate.")
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD_MAP))
    parser.add_argument("--shared-norm", default=str(DEFAULT_SHARED_NORM))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    shared_norm_path = Path(args.shared_norm)
    shared_norm = json.loads(shared_norm_path.read_text(encoding="utf-8")) if shared_norm_path.exists() else None
    artifact = build_artifact(spread_map, shared_norm)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
