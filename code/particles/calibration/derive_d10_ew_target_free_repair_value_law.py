#!/usr/bin/env python3
"""Promote the source-only D10 electroweak repair law to theorem status.

Chain role: emit the exact target-free D10 electroweak repair package from the
source-only D10 basis and the promoted source-only theorem.

Mathematics: define the source-only emitter scalar
`lambda_EW = eta_source^2 / (4 * beta_EW)`, then emit the unique repair chart
`(tau2_tree_exact, delta_n_tree_exact)`, the repaired coupling pair, and one
coherent electroweak quintet from the D10 source basis alone.

OPH-derived inputs: the emitted D10 source pair and compact current-carrier
slice, with the frozen-target repair surface retained only for compare-only
validation.

Output: a machine-readable theorem artifact for
`EWTargetFreeRepairValueLaw_D10`.
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
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"


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
    alpha_u_seed = eta_source / beta_ew
    lambda_ew = eta_source * alpha_u_seed / 4.0

    tau2_exact = -lambda_ew * (
        1.0
        + (2.0 / 3.0) * eta_source
        + (1.0 - beta_ew / 6.0) * eta_source * eta_source
    )
    delta_n_exact = lambda_ew * (
        1.0
        + (4.0 / 3.0) * eta_source
        + (2.0 - beta_ew / 6.0) * eta_source * eta_source
    )
    tauY_fiber = -(tau2_exact + 2.0 * eta_source) / (1.0 + 4.0 * tau2_exact * tau2_exact)

    delta_alpha2 = alpha_2 * tau2_exact
    delta_alphaY_parallel = alpha_y * (
        8.0 * eta_source * tau2_exact * tau2_exact - tau2_exact
    ) / (1.0 + 4.0 * tau2_exact * tau2_exact)
    delta_alphaY_perp = alpha_sum * delta_n_exact

    alpha2_star = alpha_2
    alphaY_star = alpha_y * (1.0 - 2.0 * eta_source)
    alpha2_prime = alpha2_star + delta_alpha2
    alphaY_prime = alphaY_star + delta_alphaY_parallel + delta_alphaY_perp
    alpha_sum_prime = alpha2_prime + alphaY_prime

    mw_emit = v_value * math.sqrt(math.pi * alpha2_prime)
    mz_emit = v_value * math.sqrt(math.pi * alpha_sum_prime)

    mw_frozen = float(references["w_boson"]["value_gev"])
    mz_frozen = float(references["z_boson"]["value_gev"])
    alpha2_frozen = (mw_frozen / v_value) ** 2 / math.pi
    alpha_sum_frozen = (mz_frozen / v_value) ** 2 / math.pi
    alphaY_frozen = alpha_sum_frozen - alpha2_frozen

    return {
        "artifact": "oph_d10_ew_target_free_repair_value_law",
        "generated_utc": _timestamp(),
        "status": "closed",
        "object_id": "EWTargetFreeRepairValueLaw_D10",
        "proof_gate": "single_family_single_P_no_mixed_readout",
        "phase_tier": "phase_ii_calibration",
        "family_source_id": "d10_running_tree",
        "basis": {
            "alpha2_mz": alpha_2,
            "alphaY_mz": alpha_y,
            "alpha_sum_mz": alpha_sum,
            "beta_EW": beta_ew,
            "eta_source": eta_source,
            "alpha_u_seed": alpha_u_seed,
            "lambda_EW": lambda_ew,
            "v_report_gev": v_value,
        },
        "theorem": {
            "name": "EWTargetFreeRepairValueLaw_D10",
            "statement": (
                "Let the D10 source-only basis be "
                "(alpha2_mz, alphaY_mz, eta_source, v_report). Then the beyond-current-carrier "
                "electroweak repair is uniquely emitted by the source-only repair chart "
                "(tau2_tree_exact, delta_n_tree_exact) and the repaired coupling pair "
                "(alpha2_prime, alphaY_prime), with no frozen target input, no inverse slice, "
                "and no mixed readout."
            ),
            "formulas": {
                "alpha_u_seed": "eta_source / beta_EW",
                "lambda_EW": "eta_source^2 / (4 * beta_EW)",
                "tau2_tree_exact": "-lambda_EW * (1 + (2/3) * eta_source + (1 - beta_EW/6) * eta_source^2)",
                "delta_n_tree_exact": "lambda_EW * (1 + (4/3) * eta_source + (2 - beta_EW/6) * eta_source^2)",
                "tauY_fiber": "-(tau2_tree_exact + 2 * eta_source) / (1 + 4 * tau2_tree_exact^2)",
                "delta_alpha2": "alpha2_mz * tau2_tree_exact",
                "delta_alphaY_parallel": "alphaY_mz * (8 * eta_source * tau2_tree_exact^2 - tau2_tree_exact) / (1 + 4 * tau2_tree_exact^2)",
                "delta_alphaY_perp": "(alpha2_mz + alphaY_mz) * delta_n_tree_exact",
                "alpha2_prime": "alpha2_mz + delta_alpha2",
                "alphaY_star": "alphaY_mz * (1 - 2 * eta_source)",
                "alphaY_prime": "alphaY_star + delta_alphaY_parallel + delta_alphaY_perp",
                "MW_pole": "v_report * sqrt(pi * alpha2_prime)",
                "MZ_pole": "v_report * sqrt(pi * (alpha2_prime + alphaY_prime))",
                "alpha_em_eff_inv": "(alpha2_prime + alphaY_prime) / (alpha2_prime * alphaY_prime)",
                "sin2w_eff": "alphaY_prime / (alpha2_prime + alphaY_prime)",
            },
        },
        "repair_chart": {
            "tau2_tree_exact": tau2_exact,
            "delta_n_tree_exact": delta_n_exact,
            "tauY_fiber": tauY_fiber,
        },
        "repaired_couplings": {
            "alpha2_star": alpha2_star,
            "alphaY_star": alphaY_star,
            "delta_alpha2": delta_alpha2,
            "delta_alphaY_parallel": delta_alphaY_parallel,
            "delta_alphaY_perp": delta_alphaY_perp,
            "alpha2_prime": alpha2_prime,
            "alphaY_prime": alphaY_prime,
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
        "compare_only_validation_against_frozen_surface": {
            "reference_surface_kind": "freeze_once_authoritative_target_coherent_repair_surface",
            "MW_reference_gev": mw_frozen,
            "MZ_reference_gev": mz_frozen,
            "alpha2_reference": alpha2_frozen,
            "alphaY_reference": alphaY_frozen,
            "alpha_em_eff_inv_reference": alpha_sum_frozen / (alphaY_frozen * alpha2_frozen),
            "sin2w_eff_reference": alphaY_frozen / alpha_sum_frozen,
            "delta_MW_gev": mw_emit - mw_frozen,
            "delta_MZ_gev": mz_emit - mz_frozen,
            "delta_alpha_em_eff_inv": (alpha_sum_prime / (alphaY_prime * alpha2_prime))
            - (alpha_sum_frozen / (alphaY_frozen * alpha2_frozen)),
            "delta_sin2w_eff": (alphaY_prime / alpha_sum_prime) - (alphaY_frozen / alpha_sum_frozen),
        },
        "notes": [
            "This theorem records the closed source-only D10 mass-side object beneath the Ward-projected electromagnetic transport family.",
            "The freeze-once coherent repair law serves as compare-only validation on the same W/Z mass lane.",
            "This closes the D10 electroweak mass-side lane on the Phase II calibration tier; the physical fine-structure constant is read separately on the Ward-projected U(1)_Q transport family.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 target-free repair theorem artifact.")
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
