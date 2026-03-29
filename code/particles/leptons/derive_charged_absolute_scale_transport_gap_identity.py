#!/usr/bin/env python3
"""Expose the current-family charged absolute-scale transport-gap identity.

Chain role: make the typed charged absolute-scale restore candidate less ad hoc
by tying its common gap subtraction to the emitted overlap-edge theorem gap.

Mathematics: current-family identity check plus exact one-step coordinate
conversion between the linear shared seed `g_e_raw` and the log-coordinate
candidate `mu_e_abs`.

OPH-derived inputs: the overlap-edge transport cocycle gap, the charged shared
absolute-scale writeback, and the lepton absolute-scale binding.

Output: a current-family identity artifact showing that the typed absolute-scale
candidate is `mu_e_abs = log(g_e_raw) - gamma_gap`.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_WRITEBACK = ROOT / "particles" / "runs" / "flavor" / "charged_shared_absolute_scale_writeback.json"
DEFAULT_BINDING = ROOT / "particles" / "runs" / "leptons" / "lepton_shared_absolute_scale_binding.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_absolute_scale_transport_gap_identity.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged absolute-scale transport-gap identity artifact.")
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE))
    parser.add_argument("--writeback", default=str(DEFAULT_WRITEBACK))
    parser.add_argument("--binding", default=str(DEFAULT_BINDING))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    cocycle = _load_json(Path(args.cocycle))
    writeback = _load_json(Path(args.writeback))
    binding = _load_json(Path(args.binding))

    gamma_gap = float(cocycle["theorem_gap_gamma"])
    gamma_min = float(binding["gamma_min_common"])
    g_raw = float(writeback["stored_shared_absolute_scale_raw"])
    mu_seed = float(writeback["stored_shared_absolute_scale_log_seed"])
    mu_candidate = float(mu_seed - gamma_gap)
    g_effective = float(math.exp(mu_candidate))

    artifact = {
        "artifact": "oph_charged_absolute_scale_transport_gap_identity",
        "generated_utc": _timestamp(),
        "scope": "current_family_only",
        "proof_status": "current_family_gap_identity_closed_hierarchy_still_open",
        "public_promotion_allowed": False,
        "theorem_gap_gamma": gamma_gap,
        "gamma_min_common": gamma_min,
        "identity_formula": "gamma_min_common = theorem_gap_gamma",
        "identity_residual": gamma_min - gamma_gap,
        "typed_restore_formulas": {
            "mu_e_absolute_log_seed": "log(stored_shared_absolute_scale_raw)",
            "mu_e_absolute_log_candidate": "log(stored_shared_absolute_scale_raw) - theorem_gap_gamma",
            "g_e_effective_candidate": "exp(mu_e_absolute_log_candidate)",
            "g_e_effective_candidate_alt": "stored_shared_absolute_scale_raw * exp(-theorem_gap_gamma)",
        },
        "typed_restore_values": {
            "stored_shared_absolute_scale_raw": g_raw,
            "mu_e_absolute_log_seed": mu_seed,
            "mu_e_absolute_log_candidate": mu_candidate,
            "g_e_effective_candidate": g_effective,
        },
        "lane_context": {
            "remaining_lane_frontier": "eta_source_support_extension_log_per_side",
            "next_single_residual_object_after_eta": "sigma_source_support_extension_total_log_per_side",
            "why_not_closure": "The current-family absolute-scale coordinate is cleaner, but charged hierarchy closure still depends first on the same-carrier eta-then-sigma support pair.",
        },
        "parent_artifacts": {
            "overlap_edge_transport_cocycle": cocycle.get("artifact"),
            "charged_shared_absolute_scale_writeback": writeback.get("artifact"),
            "lepton_shared_absolute_scale_binding": binding.get("artifact"),
        },
        "notes": [
            "This is a current-family identity, not a recovered-core charged-lepton mass theorem.",
            "It removes one source of arbitrariness from the typed absolute-scale restore candidate by identifying the subtracted common gap with the emitted transport-gap scalar.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
