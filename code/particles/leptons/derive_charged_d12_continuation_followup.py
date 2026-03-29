#!/usr/bin/env python3
"""Record the strongest current D12 charged-lepton continuation bridge.

Chain role: keep the honest charged theorem lane separate from a strong but
assumption-bearing D12 continuation branch.

Mathematics: evaluate the balanced circulant continuation branch, compute the
same-carrier source pair `(eta_ext, sigma_ext)`, and measure the compare-only
common-scale gap against charged-lepton references.

OPH-derived inputs: the current charged support-extension readback shell,
especially `x2` and the active charged scale candidate `g_active_candidate`.

Output: a diagnostic artifact that preserves the best available D12 continuation
pair and its assumptions without promoting it onto the theorem-grade public
surface.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
PAIR_READBACK_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_source_scalar_pair_readback.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_d12_continuation_followup.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged-lepton D12 continuation follow-up artifact.")
    parser.add_argument("--pair-readback", default=str(PAIR_READBACK_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    pair_readback = _load_json(Path(args.pair_readback))
    references = _load_json(REFERENCE_JSON)["entries"]

    x2 = float(pair_readback["ordered_family_coordinate"][1])
    g_active_candidate = float(pair_readback["g_active_candidate"])

    delta = 2.0 / 9.0
    a_value = 1.0
    abs_b = a_value / math.sqrt(2.0)
    roots = sorted(
        a_value + 2.0 * abs_b * math.cos(delta + 2.0 * math.pi * idx / 3.0)
        for idx in range(3)
    )
    r_e, r_mu, r_tau = roots

    sigma = 2.0 * math.log(r_tau / r_e)
    eta = 2.0 * ((1.0 - x2) * math.log(r_e) + (1.0 + x2) * math.log(r_tau) - 2.0 * math.log(r_mu))
    gamma21 = (((1.0 + x2) * sigma) - eta) / 2.0
    gamma32 = (((1.0 - x2) * sigma) + eta) / 2.0
    kappa = eta / (2.0 * (1.0 - x2 * x2))
    centered_logs = [
        -((2.0 * gamma21) + gamma32) / 3.0,
        (gamma21 - gamma32) / 3.0,
        (gamma21 + 2.0 * gamma32) / 3.0,
    ]
    shape = [math.exp(value) for value in centered_logs]

    target = [
        float(references["electron"]["value_gev"]),
        float(references["muon"]["value_gev"]),
        float(references["tau"]["value_gev"]),
    ]
    reference_logs = [math.log(value) for value in target]
    reference_mean = sum(reference_logs) / 3.0
    reference_centered = [value - reference_mean for value in reference_logs]
    centered_residual = [centered_logs[idx] - reference_centered[idx] for idx in range(3)]
    centered_residual_norm = math.sqrt(sum(value * value for value in centered_residual))
    g_e_needed = math.exp(reference_mean)

    payload = {
        "artifact": "oph_charged_d12_continuation_followup",
        "generated_utc": _timestamp(),
        "status": "diagnostic_only_d12_continuation",
        "public_promotion_allowed": False,
        "theorem_tier": "D12_continuation_only",
        "d12_continuation_assumptions": {
            "A1_uniform_Z6_center_label_ensemble": "epsilon = 1/6",
            "A2_balanced_circulant_branch": "|b|/a = 1/sqrt(2)",
            "A3_oph_phase_choice": "delta = 2/9",
            "scale_free_normalization_for_ratio_invariants": "a = 1",
        },
        "current_theorem_lane_status": {
            "forward_closure_state": "support_extension_missing",
            "next_single_residual_object": "eta_source_support_extension_log_per_side",
            "next_single_residual_object_after_eta": "sigma_source_support_extension_total_log_per_side",
            "channel_norm_proof_status": "shared_budget_only",
            "current_g_active_candidate": g_active_candidate,
        },
        "d12_continuation_pair": {
            "x2": x2,
            "delta": delta,
            "a": a_value,
            "abs_b": abs_b,
            "roots_sorted": {"r_e": r_e, "r_mu": r_mu, "r_tau": r_tau},
            "eta_source_support_extension_log_per_side": eta,
            "sigma_source_support_extension_total_log_per_side": sigma,
            "gamma21_log_per_side": gamma21,
            "gamma32_log_per_side": gamma32,
            "kappa_ext": kappa,
            "E_e_log_centered_emitted": centered_logs,
            "shape_singular_values_emitted": shape,
        },
        "compare_only_shape_check_against_reference_masses": {
            "reference_masses": {"e": target[0], "mu": target[1], "tau": target[2]},
            "reference_centered_logs": reference_centered,
            "continuation_centered_logs": centered_logs,
            "centered_log_residual": centered_residual,
            "centered_log_residual_norm": centered_residual_norm,
            "g_e_compare_only_needed_for_exact_absolute_masses": g_e_needed,
            "mean_scale_ratio": g_e_needed / g_active_candidate,
            "log_shift_needed_relative_to_current_g_active_candidate": math.log(g_e_needed / g_active_candidate),
            "shape_to_reference_scale_ratios": {
                "e": target[0] / shape[0],
                "mu": target[1] / shape[1],
                "tau": target[2] / shape[2],
            },
        },
        "immediate_followup_required": {
            "live_scalar_gap": "eta_source_support_extension_log_per_side then sigma_source_support_extension_total_log_per_side",
            "absolute_scale_gap": "upgrade channel norm from shared_budget_only to shared_budget_closed or sector_local_closed so g_e is emitted on the charged lane",
            "core_theorem_gap": "derive or replace A1-A3 if charged leptons are to become recovered-core outputs",
            "paper_status_decision": "either keep the live theorem lane open or add an explicit D12 continuation theorem/corollary with A1-A3 marked as extra assumptions",
        },
        "notes": [
            "This artifact records the strongest current charged-lepton D12 continuation bridge.",
            "It is not theorem-grade OPH closure because the continuation assumptions A1-A3 are extra assumptions and the charged absolute scale g_e is still compare-only here.",
            "The charged public rows should remain suppressed until the same-carrier eta/sigma pair and the charged absolute scale are emitted on the live theorem lane.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
