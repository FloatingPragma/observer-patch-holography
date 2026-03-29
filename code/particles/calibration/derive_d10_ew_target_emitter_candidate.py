#!/usr/bin/env python3
"""Emit the strongest current source-only target-free D10 repair candidate.

Chain role: package the best current source-only target-emitter candidate above
the freeze-once electroweak repair surface without promoting it to theorem
status.

Mathematics: build a single scalar `lambda_EW` from the emitted D10 source data
and use it to emit a coherent target-free candidate for `(tau2_tree_exact,
delta_n_tree_exact)` and the repaired electroweak quintet.

OPH-derived inputs: the emitted D10 source pair and compact current-carrier
slice, with reference W/Z used only for compare-only residual reporting.

Output: a machine-readable candidate artifact for the target-free electroweak
repair law.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_emitter_candidate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(source_pair: dict, references: dict) -> dict:
    pair = dict(source_pair.get("source_pair") or {})
    compact_slice = dict(source_pair.get("compact_hypercharge_only_mass_slice") or {})
    compact_quintet = dict(compact_slice.get("coherent_output_quintet") or {})
    alpha_2 = float(pair["alpha2_mz"])
    alpha_y = float(pair["alphaY_mz"])
    alpha_sum = alpha_2 + alpha_y
    eta_source = float(compact_slice["eta_EW"])
    v_value = float(compact_quintet["v_report"])
    beta_ew = (alpha_2 - alpha_y) / alpha_sum
    alpha_u_from_seed = eta_source / beta_ew
    lambda_ew = eta_source * alpha_u_from_seed / 4.0
    tau2 = -lambda_ew * (1.0 + (2.0 / 3.0) * eta_source + (1.0 - beta_ew / 6.0) * eta_source * eta_source)
    delta_n = lambda_ew * (1.0 + (4.0 / 3.0) * eta_source + (2.0 - beta_ew / 6.0) * eta_source * eta_source)
    tau_y = -(tau2 + 2.0 * eta_source) / (1.0 + 4.0 * tau2 * tau2)
    delta_alpha2 = alpha_2 * tau2
    delta_alphaY_parallel = alpha_y * (8.0 * eta_source * tau2 * tau2 - tau2) / (1.0 + 4.0 * tau2 * tau2)
    delta_alphaY_perp = alpha_sum * delta_n
    alpha2_prime = alpha_2 + delta_alpha2
    alphaY_star = alpha_y * (1.0 - 2.0 * eta_source)
    alphaY_prime = alphaY_star + delta_alphaY_parallel + delta_alphaY_perp
    alpha_sum_prime = alpha2_prime + alphaY_prime

    mw_target = float(references["w_boson"]["value_gev"])
    mz_target = float(references["z_boson"]["value_gev"])
    comparison_tau2 = mw_target * mw_target / (math.pi * v_value * v_value * alpha_2) - 1.0
    comparison_delta_alpha2 = (mw_target * mw_target) / (math.pi * v_value * v_value) - alpha_2
    comparison_mz_fiber = v_value * math.sqrt(
        math.pi * alpha_sum * (1.0 + (alpha_y * tau_y + alpha_2 * comparison_tau2) / alpha_sum)
    )
    comparison_delta_mz = mz_target - comparison_mz_fiber
    comparison_delta_n = ((mz_target + comparison_mz_fiber) * comparison_delta_mz) / (
        math.pi * v_value * v_value * alpha_sum
    )
    comparison_delta_alphaY_parallel = alpha_y * (
        8.0 * eta_source * comparison_tau2 * comparison_tau2 - comparison_tau2
    ) / (1.0 + 4.0 * comparison_tau2 * comparison_tau2)
    comparison_delta_alphaY_perp = ((mz_target + comparison_mz_fiber) * comparison_delta_mz) / (
        math.pi * v_value * v_value
    )

    mw_emit = v_value * math.sqrt(math.pi * alpha2_prime)
    mz_emit = v_value * math.sqrt(math.pi * alpha_sum_prime)

    return {
        "artifact": "oph_d10_ew_target_emitter_candidate",
        "generated_utc": _timestamp(),
        "status": "new_theorem_candidate_not_yet_in_corpus",
        "object_id": "EWTargetEmitter_D10",
        "proof_gate": "single_family_single_P_no_mixed_readout",
        "family_source_id": "d10_running_tree",
        "basis": {
            "alpha2_mz": alpha_2,
            "alphaY_mz": alpha_y,
            "alpha_sum_mz": alpha_sum,
            "alpha_u_from_seed": alpha_u_from_seed,
            "beta_EW": beta_ew,
            "eta_source": eta_source,
            "v_report_gev": v_value,
        },
        "current_selected_carrier": {
            "MW_pole": float(compact_quintet["MW_pole"]),
            "MZ_pole": float(compact_quintet["MZ_pole"]),
            "alpha2_star": alpha_2,
            "alphaY_star": alphaY_star,
            "alpha_em_eff_inv": float(compact_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(compact_quintet["sin2w_eff"]),
            "v_report": v_value,
            "eta_EW": eta_source,
            "sigma_EW": float(compact_slice["sigma_EW"]),
            "tau_2": float(compact_slice["tau_2"]),
            "tau_Y": float(compact_slice["tau_Y"]),
        },
        "emitter_scalar": {
            "lambda_formula": "eta_source * alpha_u_from_seed / 4 = eta_source^2 / (4 * beta_EW)",
            "lambda_value": lambda_ew,
        },
        "target_emitter_law": {
            "tau2_tree_exact_formula": "-lambda * (1 + (2/3) * eta_source + (1 - beta_EW/6) * eta_source^2)",
            "tau2_tree_exact": tau2,
            "delta_n_tree_exact_formula": "lambda * (1 + (4/3) * eta_source + (2 - beta_EW/6) * eta_source^2)",
            "delta_n_tree_exact": delta_n,
            "tauY_fiber_formula": "-(tau2_tree_exact + 2 * eta_source) / (1 + 4 * tau2_tree_exact^2)",
            "tauY_fiber": tau_y,
            "delta_alpha2_tree": delta_alpha2,
            "delta_alphaY_parallel": delta_alphaY_parallel,
            "delta_alphaY_perp": delta_alphaY_perp,
            "delta_alphaY_tree": delta_alphaY_parallel + delta_alphaY_perp,
        },
        "coherent_emitted_quintet": {
            "MW_pole": mw_emit,
            "MZ_pole": mz_emit,
            "alpha2_prime": alpha2_prime,
            "alphaY_prime": alphaY_prime,
            "alpha_em_eff_inv": alpha_sum_prime / (alphaY_prime * alpha2_prime),
            "sin2w_eff": alphaY_prime / alpha_sum_prime,
            "v_report": v_value,
        },
        "comparison_to_frozen_local_reference_surface": {
            "MW_target_gev": mw_target,
            "MZ_target_gev": mz_target,
            "MW_difference_gev": mw_emit - mw_target,
            "MZ_difference_gev": mz_emit - mz_target,
            "tau2_tree_exact_compare_only": comparison_tau2,
            "tau2_difference": tau2 - comparison_tau2,
            "delta_n_tree_exact_compare_only": comparison_delta_n,
            "delta_n_difference": delta_n - comparison_delta_n,
            "delta_alpha2_tree_compare_only": comparison_delta_alpha2,
            "delta_alphaY_parallel_compare_only": comparison_delta_alphaY_parallel,
            "delta_alphaY_perp_compare_only": comparison_delta_alphaY_perp,
            "delta_alphaY_tree_compare_only": comparison_delta_alphaY_parallel + comparison_delta_alphaY_perp,
        },
        "honesty_note": (
            "This is a newly derived target-emitter candidate using only already-emitted D10 source data; "
            "it is not yet stated as a theorem in the current papers."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 target-emitter candidate artifact.")
    parser.add_argument("--source-pair", default=str(SOURCE_PAIR_JSON))
    parser.add_argument("--references", default=str(REFERENCE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = _load_json(Path(args.source_pair))
    references = _load_json(Path(args.references))["entries"]
    artifact = build_artifact(source_pair, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
