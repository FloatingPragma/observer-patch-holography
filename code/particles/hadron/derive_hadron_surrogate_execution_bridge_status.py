#!/usr/bin/env python3
"""Record the surrogate hadron execution bridge as a diagnostic artifact.

Chain role: preserve the strongest currently available end-to-end hadron
execution bridge without confusing surrogate execution with production lattice
QCD closure.

Mathematics: no new hadron theorem is derived here; this script records the
surrogate HMC/RHMC execution bridge completed in the March 28, 2026
intermediate follow-up bundle.

OPH-derived inputs: the emitted seeded 2+1-family schedule/writeback contract
and the surrogate execution summary values extracted from the bundle review.

Output: a diagnostic-only surrogate execution status artifact consumed by the
hadron audit and paper/status summaries.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "hadron_surrogate_execution_bridge_status.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict:
    return {
        "artifact": "oph_hadron_surrogate_execution_bridge_status",
        "generated_utc": _timestamp(),
        "status": "surrogate_hmc_execution_bridge_complete",
        "theorem_tier": "diagnostic_execution_bridge_only",
        "public_promotion_allowed": False,
        "surrogate_execution": True,
        "runtime_receipt_frozen": {
            "N_therm": 2048,
            "N_sep": 512,
        },
        "surrogate_finest_ensemble": {
            "ensemble_id": "qcd_2p1_seed_n0",
            "aLambda_msbar3": 0.0550900879921322,
            "pi_iso_mass_gev_candidate": 0.13503938383599978,
            "N_iso_mass_gev_candidate": 0.9389602105777344,
            "pi_iso_forward_window_limit_exists": True,
            "N_iso_forward_window_limit_exists": True,
        },
        "surrogate_error_summary": {
            "worst_abs_error_gev_vs_reference": 6.255606831506721e-05,
            "worst_rel_error_vs_reference": 4.634578345753951e-04,
            "all_six_channel_limits_exist": True,
        },
        "surrogate_completion_scope": [
            "non-null runtime receipt",
            "executed cfg/source writeback on the seeded 2+1 family",
            "populated published_statistical_error",
            "populated published_systematics",
            "forward-window convergence from realized arrays",
        ],
        "honest_remaining_physical_blockers": [
            "production unquenched RHMC/HMC on SU(3) gauge fields",
            "real Dirac solves and baryon contractions",
            "production autocorrelation/statistics studies",
            "production continuum / finite-volume / chiral systematics",
        ],
        "notes": [
            "This artifact records the strongest current hadron execution bridge, but it is surrogate-only and not promotable.",
            "The surrogate bridge is useful because it closes the runtime-receipt -> writeback -> evaluation -> convergence -> systematics path on the emitted schema.",
            "Public hadron rows remain blocked on replacing the surrogate kernel with production unquenched execution and production systematics.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the surrogate hadron execution bridge status artifact.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
