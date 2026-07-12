#!/usr/bin/env python3
"""Record the current runtime target-free-input D10 repair candidate.

Chain role: evaluate the current D10 electroweak repair formula from the
source-side D10 basis while keeping its source-selection theorem open.

Mathematics: define the candidate emitter scalar
`lambda_EW = eta_source^2 / (4 * beta_EW)`, then deterministically evaluate
the selected repair chart `(tau2_tree_exact, delta_n_tree_exact)`, the repaired
coupling pair, and one coherent electroweak quintet.  The current corpus does
not yet derive or uniquely select this chart from the finite quotient carrier.

Inputs: the emitted D10 source pair and compact current-carrier slice, with the
frozen-target repair surface retained only for compare-only validation.  A
runtime formula with no target argument is not by itself a prospective
source-entailment proof.

Output: a machine-readable candidate artifact for
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


def evaluate_candidate_from_source_basis(
    *,
    alpha_2: float,
    alpha_y: float,
    eta_source: float,
    v_value: float,
) -> dict:
    """Evaluate the selected candidate law without loading reference data.

    This is the canonical numerical evaluator shared by the candidate artifact
    and the conditional quotient-transport receipt.  Purity of this function
    establishes only that its runtime evaluation has no target argument; it
    does not prove that the formula itself is selected by the OPH source laws.
    """

    if not (alpha_2 > 0.0 and alpha_y > 0.0 and eta_source > 0.0 and v_value > 0.0):
        raise ValueError("alpha_2, alpha_y, eta_source, and v_value must be positive")
    alpha_sum = alpha_2 + alpha_y
    beta_ew = (alpha_2 - alpha_y) / alpha_sum
    if not 0.0 < beta_ew < 1.0:
        raise ValueError("candidate domain requires 0 < beta_EW < 1")
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
    tau_y_fiber = -(tau2_exact + 2.0 * eta_source) / (
        1.0 + 4.0 * tau2_exact * tau2_exact
    )

    delta_alpha2 = alpha_2 * tau2_exact
    delta_alpha_y_parallel = alpha_y * (
        8.0 * eta_source * tau2_exact * tau2_exact - tau2_exact
    ) / (1.0 + 4.0 * tau2_exact * tau2_exact)
    delta_alpha_y_perp = alpha_sum * delta_n_exact

    alpha2_star = alpha_2
    alpha_y_star = alpha_y * (1.0 - 2.0 * eta_source)
    alpha2_prime = alpha2_star + delta_alpha2
    alpha_y_prime = alpha_y_star + delta_alpha_y_parallel + delta_alpha_y_perp
    if not (alpha2_prime > 0.0 and alpha_y_prime > 0.0):
        raise ValueError("candidate repaired couplings must remain positive")
    alpha_sum_prime = alpha2_prime + alpha_y_prime

    return {
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
        "repair_chart": {
            "tau2_tree_exact": tau2_exact,
            "delta_n_tree_exact": delta_n_exact,
            "tauY_fiber": tau_y_fiber,
        },
        "repaired_couplings": {
            "alpha2_star": alpha2_star,
            "alphaY_star": alpha_y_star,
            "delta_alpha2": delta_alpha2,
            "delta_alphaY_parallel": delta_alpha_y_parallel,
            "delta_alphaY_perp": delta_alpha_y_perp,
            "alpha2_prime": alpha2_prime,
            "alphaY_prime": alpha_y_prime,
        },
        "coherent_emitted_quintet": {
            "MW_pole": v_value * math.sqrt(math.pi * alpha2_prime),
            "MZ_pole": v_value * math.sqrt(math.pi * alpha_sum_prime),
            "alpha2_prime": alpha2_prime,
            "alphaY_prime": alpha_y_prime,
            "alpha_em_eff_inv": alpha_sum_prime / (alpha_y_prime * alpha2_prime),
            "sin2w_eff": alpha_y_prime / alpha_sum_prime,
            "v_report": v_value,
        },
    }


def build_artifact(source_pair: dict, references: dict) -> dict:
    pair = dict(source_pair.get("source_pair") or {})
    compact_slice = dict(source_pair.get("compact_hypercharge_only_mass_slice") or {})
    compact_quintet = dict(compact_slice.get("coherent_output_quintet") or {})
    alpha_2 = float(pair["alpha2_mz"])
    alpha_y = float(pair["alphaY_mz"])
    eta_source = float(compact_slice["eta_EW"])
    v_value = float(compact_quintet["v_report"])
    evaluated = evaluate_candidate_from_source_basis(
        alpha_2=alpha_2,
        alpha_y=alpha_y,
        eta_source=eta_source,
        v_value=v_value,
    )
    basis = evaluated["basis"]
    repair_chart = evaluated["repair_chart"]
    repaired_couplings = evaluated["repaired_couplings"]
    emitted_quintet = evaluated["coherent_emitted_quintet"]
    alpha_sum = float(basis["alpha_sum_mz"])
    alpha2_prime = float(repaired_couplings["alpha2_prime"])
    alphaY_prime = float(repaired_couplings["alphaY_prime"])
    alpha_sum_prime = alpha2_prime + alphaY_prime
    mw_emit = float(emitted_quintet["MW_pole"])
    mz_emit = float(emitted_quintet["MZ_pole"])

    mw_frozen = float(references["w_boson"]["value_gev"])
    mz_frozen = float(references["z_boson"]["value_gev"])
    alpha2_frozen = (mw_frozen / v_value) ** 2 / math.pi
    alpha_sum_frozen = (mz_frozen / v_value) ** 2 / math.pi
    alphaY_frozen = alpha_sum_frozen - alpha2_frozen

    return {
        "artifact": "oph_d10_ew_target_free_repair_value_law",
        "generated_utc": _timestamp(),
        "status": "candidate_only",
        "object_id": "EWTargetFreeRepairValueLaw_D10",
        "proof_gate": "blocked_by_current_corpus_underdetermination",
        "phase_tier": "phase_ii_calibration",
        "promotion_allowed": False,
        "promotion_blockers": [
            "current_corpus_underdetermination_of_forward_d10_repair_law",
        ],
        "family_source_id": "d10_running_tree",
        "basis": basis,
        "theorem": {
            "name": "EWTargetFreeRepairValueLaw_D10",
            "statement": (
                "Given the current selected D10 candidate rule and the basis "
                "(alpha2_mz, alphaY_mz, eta_source, v_report), the formulas below "
                "deterministically evaluate (tau2_tree_exact, delta_n_tree_exact), the "
                "repaired coupling pair (alpha2_prime, alphaY_prime), and one coherent "
                "quintet without a frozen target argument in the forward calculation. "
                "The present corpus does not derive or uniquely select the candidate rule "
                "from the finite quotient transport, so this statement is not a promoted "
                "source theorem."
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
        "repair_chart": repair_chart,
        "repaired_couplings": repaired_couplings,
        "coherent_emitted_quintet": emitted_quintet,
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
            "This artifact records the current runtime target-free-input D10 mass-side candidate beneath the Ward-projected electromagnetic transport family.",
            "The unconditional source-only underdetermination theorem remains the support boundary on the present corpus, so this candidate is not promoted to theorem status.",
            "The exact conditional QT1-QT5 quotient-path theorem reproduces this formula, but the current corpus does not yet emit QT1-QT5 or bind them to one strict source pixel branch.",
            "The freeze-once coherent repair law serves as compare-only validation on the same W/Z mass lane.",
            "At the repo default P this candidate nearly coincides numerically with the frozen benchmark surface, but that numerical coincidence is not by itself a proof of target-free closure.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 target-free-input repair candidate artifact.")
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
