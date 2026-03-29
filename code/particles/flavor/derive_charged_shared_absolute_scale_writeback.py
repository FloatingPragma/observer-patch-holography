#!/usr/bin/env python3
"""Build the current-family charged shared absolute-scale writeback bundle.

Chain role: record the shared charged absolute-scale writeback in one place and
keep the linear-vs-log coordinate split explicit for downstream consumers.

Mathematics: representation bookkeeping only; the emitted writeback remains a
current-family shared-budget artifact and does not close the charged hierarchy.

OPH-derived inputs: the shared charged-budget promotion surface together with
the active charged ordered-gap readout.

Output: a writeback artifact plus lepton/quark bindings that expose both the
linear seed `g_e` and the log-coordinate audit shell consistently.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PROMOTION = ROOT / "particles" / "runs" / "flavor" / "charged_shared_absolute_scale_promotion.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "charged_shared_absolute_scale_writeback.json"
DEFAULT_LEPTON_BINDING = ROOT / "particles" / "runs" / "leptons" / "lepton_shared_absolute_scale_binding.json"
DEFAULT_QUARK_BINDING = ROOT / "particles" / "runs" / "flavor" / "quark_shared_absolute_norm_binding.json"
DEFAULT_BUNDLE = ROOT / "particles" / "runs" / "flavor" / "charged_shared_absolute_scale_writeback_bundle.json"
DEFAULT_LEPTON_LOG = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current-family charged shared absolute-scale writeback.")
    parser.add_argument("--input", default=str(DEFAULT_PROMOTION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--lepton-binding", default=str(DEFAULT_LEPTON_BINDING))
    parser.add_argument("--quark-binding", default=str(DEFAULT_QUARK_BINDING))
    parser.add_argument("--bundle-output", default=str(DEFAULT_BUNDLE))
    parser.add_argument("--lepton-log-readout", default=str(DEFAULT_LEPTON_LOG))
    args = parser.parse_args()

    promotion = json.loads(Path(args.input).read_text(encoding="utf-8"))
    lepton_log_path = Path(args.lepton_log_readout)
    lepton_log = json.loads(lepton_log_path.read_text(encoding="utf-8")) if lepton_log_path.exists() else None
    shared_scale = float(promotion["stored_shared_absolute_scale"])
    sector_share_by_sector = dict(promotion.get("sector_share_by_sector", {}))
    common_log_restore_candidate = None
    if lepton_log and isinstance(lepton_log.get("gap_pair_e"), dict):
        gap_pair = dict(lepton_log["gap_pair_e"])
        gamma21 = float(gap_pair["gamma21_log_per_side"])
        gamma32 = float(gap_pair["gamma32_log_per_side"])
        gamma_min = min(gamma21, gamma32)
        mu_abs_log_seed = math.log(shared_scale)
        mu_common_log = mu_abs_log_seed - gamma_min
        common_log_restore_candidate = {
            "candidate_status": "partial_common_log_restore_only",
            "gamma_min_common": gamma_min,
            "g_e_linear_seed": shared_scale,
            "mu_e_absolute_log_seed": mu_abs_log_seed,
            "mu_common_log": mu_common_log,
            "mu_e_absolute_log_candidate": mu_common_log,
            "shared_absolute_scale_effective": math.exp(mu_common_log),
            "effective_scale_formula": "exp(log(stored_shared_absolute_scale_raw) - gamma_min_common)",
            "why_not_closure": (
                "This fixes the linear-vs-log type mismatch in the restore candidate, "
                "but it does not fix the remaining charged-lepton hierarchy/open-value burden."
            ),
        }

    writeback = {
        "artifact": "charged_shared_absolute_scale_writeback",
        "generated_utc": _timestamp(),
        "parent_artifact": promotion.get("artifact"),
        "parent_candidate": promotion.get("parent_candidate"),
        "writeback_scope": "current_family_only",
        "theorem_scope": "shared_budget_only",
        "proof_status": "current_family_writeback_complete",
        "stored_shared_absolute_scale_raw": shared_scale,
        "stored_shared_absolute_scale": shared_scale,
        "stored_shared_absolute_scale_log_seed": math.log(shared_scale),
        "sigma_seed_clause": {
            "left_common": shared_scale,
            "right_common": shared_scale,
            "residual": 0.0,
        },
        "sector_equalizer_by_sector": {"u": 0.0, "d": 0.0, "e": 0.0},
        "gluing_norm": 0.0,
        "shared_budget_total": promotion.get("shared_budget_total"),
        "sector_share_by_sector": sector_share_by_sector,
        "mean_witness": promotion.get("mean_witness"),
        "common_log_restore_candidate": common_log_restore_candidate,
        "remaining_mass_mover": "oph_family_excitation_spread_map",
        "immediate_consumer_bindings": [
            "lepton_shared_absolute_scale_binding",
            "quark_shared_absolute_norm_binding",
        ],
        "metadata": {
            "no_reseed": bool(promotion.get("no_reseed", True)),
            "promotion_scope": promotion.get("promotion_scope"),
            "note": (
                "Current-family writeback of the stored shared charged absolute scale. "
                "This unblocks the lepton absolute-scale input and pins the shared quark "
                "absolute norm without promoting the outer shared-budget route to a global theorem. "
                "A later partial memo also suggested a common-log restore candidate after the centered-log "
                "consumer; that candidate is recorded here for auditability but is not promoted as charged-lepton closure."
            ),
        },
    }

    lepton_binding = {
        "artifact": "lepton_shared_absolute_scale_binding",
        "generated_utc": _timestamp(),
        "parent_artifact": writeback["artifact"],
        "writeback_scope": "current_family_only",
        "theorem_scope": "shared_budget_only",
        "proof_status": "current_family_writeback_complete",
        "closure_route": "shared_absolute_scale_writeback",
        "g_e": shared_scale,
        "g_e_raw_seed": shared_scale,
        "mu_e_absolute_log_seed": math.log(shared_scale),
        "shared_budget_share_e": sector_share_by_sector.get("e"),
    }
    if common_log_restore_candidate is not None:
        lepton_binding.update(
            {
                "g_e_effective_candidate": common_log_restore_candidate["shared_absolute_scale_effective"],
                "g_e_effective_formula": common_log_restore_candidate["effective_scale_formula"],
                "g_e_linear_seed": common_log_restore_candidate["g_e_linear_seed"],
                "mu_e_absolute_log_candidate": common_log_restore_candidate["mu_e_absolute_log_candidate"],
                "mu_e_absolute_log_seed": common_log_restore_candidate["mu_e_absolute_log_seed"],
                "mu_common_log": common_log_restore_candidate["mu_common_log"],
                "gamma_min_common": common_log_restore_candidate["gamma_min_common"],
                "candidate_status": common_log_restore_candidate["candidate_status"],
            }
        )

    quark_binding = {
        "artifact": "quark_shared_absolute_norm_binding",
        "generated_utc": _timestamp(),
        "parent_artifact": writeback["artifact"],
        "writeback_scope": "current_family_only",
        "theorem_scope": "shared_budget_only",
        "proof_status": "current_family_writeback_complete",
        "closure_route": "shared_absolute_scale_writeback",
        "g_ch": shared_scale,
        "g_u": shared_scale,
        "g_d": shared_scale,
        "shared_budget_share_u": sector_share_by_sector.get("u"),
        "shared_budget_share_d": sector_share_by_sector.get("d"),
        "remaining_mass_mover": "oph_family_excitation_spread_map",
    }

    bundle = {
        "artifact": "charged_shared_absolute_scale_writeback_bundle",
        "generated_utc": _timestamp(),
        "writeback": writeback,
        "lepton_binding_artifact": lepton_binding["artifact"],
        "quark_binding_artifact": quark_binding["artifact"],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(writeback, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lepton_path = Path(args.lepton_binding)
    lepton_path.parent.mkdir(parents=True, exist_ok=True)
    lepton_path.write_text(json.dumps(lepton_binding, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    quark_path = Path(args.quark_binding)
    quark_path.parent.mkdir(parents=True, exist_ok=True)
    quark_path.write_text(json.dumps(quark_binding, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    bundle_path = Path(args.bundle_output)
    bundle_path.parent.mkdir(parents=True, exist_ok=True)
    bundle_path.write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"saved: {out_path}")
    print(f"saved: {lepton_path}")
    print(f"saved: {quark_path}")
    print(f"saved: {bundle_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
