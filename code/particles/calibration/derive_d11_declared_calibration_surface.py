#!/usr/bin/env python3
"""Emit the declared D10/D11 calibration surface for the Higgs/top lane.

Chain role: make the D11 core and Jacobian payload explicit on the declared
running, matching, and threshold surface used by the forward-seed theorem.

Mathematics: fix one declared D11 calibration surface carrying the core
coordinates and linear Jacobian for the compact Higgs/top readout.

OPH-derived inputs: the closed D10 gauge source family plus the declared
D10/D11 running, matching, and threshold conventions.

Output: a theorem-facing calibration-surface artifact consumed by the live
forward seed and the compare-only inverse adapter.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_D10_SOURCE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d10_source: Path) -> dict[str, object]:
    return {
        "artifact": "oph_d11_declared_calibration_surface",
        "generated_utc": _timestamp(),
        "proof_status": "declared_d10_d11_calibration_surface_fixed",
        "surface_kind": "declared_running_matching_threshold_surface",
        "predictive_promotion_allowed": False,
        "d10_source_artifact": str(d10_source),
        "scope": [
            "compact_d11_forward_seed",
            "d11_jacobian_mass_readout",
            "compare_only_inverse_slice",
        ],
        "declared_surface_inputs": {
            "running_package": "printed_d10_d11_running_surface",
            "matching_package": "printed_d10_d11_matching_surface",
            "threshold_package": "printed_d10_d11_threshold_surface",
        },
        "transport_family": "sm_below_sync__mssm_like_gauge_above_sync",
        "mu_sync_gev": 6.8e11,
        "mu_eval_gev": 160.61247,
        "core": {
            "mt_ms_gev": 160.61247,
            "mt_pole_core_gev": 170.26125,
            "mH_core_gev": 126.62263,
            "y_t_core_mt": 0.92046435,
            "lambda_core_mt": 0.13164915,
            "alpha_s_mt": 0.11018777,
            "pole_ratio_core": 1.06007492,
        },
        "jacobian": {
            "d_mt_pole_d_y_t": 184.97,
            "d_mH_d_lambda": 480.0,
        },
        "fixed_ray_candidate": {
            "kappa_HT": 16.0 / 9.0,
            "seed_symbol": "sigma_D11_HT",
        },
        "notes": [
            "This artifact makes the declared D10/D11 running, matching, and threshold surface explicit for the compact Higgs/top lane.",
            "The live forward seed and the compare-only inverse slice both consume this surface rather than a diagnostic-sidecar payload.",
            "The artifact fixes the D11 core coordinates and linear Jacobian on the declared surface; it does not by itself claim reference fitting or a larger recovered-core theorem.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the declared D11 calibration surface.")
    parser.add_argument("--d10-source", default=str(DEFAULT_D10_SOURCE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(Path(args.d10_source))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
