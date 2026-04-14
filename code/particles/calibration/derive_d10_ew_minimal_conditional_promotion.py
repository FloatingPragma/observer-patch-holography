#!/usr/bin/env python3
"""Record the historical sharpest target-free D10 repair split.

Chain role: preserve the sharpest pre-theorem split beneath the promoted
target-free D10 repair law.

Mathematics: two statements are recorded. First, the current emitted source data
only determine a two-parameter quadratic repair family, so the target-free
repair law is underdetermined on the present corpus. Second, once one extra
color-balanced quadratic descent principle is assumed, the repair law closes
uniquely.

OPH-derived inputs: the emitted D10 source pair and compact current-carrier
slice.

Output: a machine-readable artifact containing the unconditional
underdetermination theorem and the smallest historical conditional route
beneath the promoted target-free repair theorem.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_minimal_conditional_theorem.json"
DEFAULT_COLOR_COUNT = 3


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(source_pair: dict, *, color_count: int) -> dict:
    pair = dict(source_pair.get("source_pair") or {})
    compact_slice = dict(source_pair.get("compact_hypercharge_only_mass_slice") or {})
    compact_quintet = dict(compact_slice.get("coherent_output_quintet") or {})
    alpha_y = float(pair["alphaY_mz"])
    alpha_2 = float(pair["alpha2_mz"])
    eta_source = float(compact_slice["eta_EW"])
    v_value = float(compact_quintet["v_report"])
    beta_ew = (alpha_2 - alpha_y) / (alpha_2 + alpha_y)

    tau2_exact = -(math.sqrt(color_count) / 2.0) * eta_source**2
    delta_n_exact = (color_count / 2.0) * (1.0 - beta_ew) * eta_source**2
    delta_alpha2_tree = alpha_2 * tau2_exact
    delta_alphaY_parallel = alpha_y * (8.0 * eta_source * tau2_exact * tau2_exact - tau2_exact) / (
        1.0 + 4.0 * tau2_exact * tau2_exact
    )
    delta_alphaY_perp = color_count * alpha_y * eta_source**2
    alpha_y_star = alpha_y * (1.0 - 2.0 * eta_source)
    alpha_2_prime = alpha_2 + delta_alpha2_tree
    alpha_y_prime = alpha_y_star + delta_alphaY_parallel + delta_alphaY_perp
    alpha_sum_prime = alpha_y_prime + alpha_2_prime

    return {
        "artifact": "oph_d10_ew_minimal_conditional_promotion",
        "generated_utc": _timestamp(),
        "status": "historical_split_superseded_by_target_free_repair_theorem",
        "superseded_by": "EWTargetFreeRepairValueLaw_D10",
        "source_artifact": source_pair.get("artifact"),
        "unconditional_theorem": {
            "name": "current_corpus_underdetermination_of_forward_d10_repair_law",
            "statement": (
                "For arbitrary real (c,d), the source-only family tau2_exact = -c*eta_source^2 and "
                "delta_n_exact = d*(1-beta_EW)*eta_source^2 defines a coherent runtime-forward D10 repair "
                "law. Therefore the current emitted corpus does not determine a unique forward repair coefficient pair."
            ),
            "free_family": {
                "tau2_exact": "-c * eta_source^2",
                "delta_n_exact": "d * (1 - beta_EW) * eta_source^2",
            },
        },
        "conditional_principle": {
            "name": "ColorBalancedQuadraticRepairDescent_D10",
            "statement": (
                "The first nonzero beyond-current-carrier D10 repair is quadratic in eta_source, with "
                "charged contraction delta_alpha2_tree = -(sqrt(N_c)/2)*alpha2_mz*eta_source^2 and "
                "neutral color-balanced uplift delta_alphaY_perp = N_c*alphaY_mz*eta_source^2."
            ),
        },
        "conditional_theorem": {
            "name": "minimal_conditional_d10_forward_repair_law",
            "statement": (
                "Assuming ColorBalancedQuadraticRepairDescent_D10, the repair law is uniquely "
                "tau2_exact = -(sqrt(N_c)/2)*eta_source^2 and delta_n_exact = (N_c/2)*(1-beta_EW)*eta_source^2."
            ),
            "formulas": {
                "tau2_exact": "-(sqrt(N_c)/2) * eta_source^2",
                "delta_n_exact": "(N_c/2) * (1 - beta_EW) * eta_source^2",
                "delta_alpha2_tree": "-(sqrt(N_c)/2) * alpha2_mz * eta_source^2",
                "delta_alphaY_perp": "N_c * alphaY_mz * eta_source^2",
            },
        },
        "realized_color_count": int(color_count),
        "basis": {
            "alphaY_mz": alpha_y,
            "alpha2_mz": alpha_2,
            "beta_EW": beta_ew,
            "eta_source": eta_source,
            "v_inherited_gev": v_value,
        },
        "n_c_3_specialization": {
            "tau2_exact": tau2_exact,
            "delta_n_exact": delta_n_exact,
            "delta_alpha2_tree": delta_alpha2_tree,
            "delta_alphaY_parallel": delta_alphaY_parallel,
            "delta_alphaY_perp": delta_alphaY_perp,
            "delta_alphaY_total": delta_alphaY_parallel + delta_alphaY_perp,
            "coherent_quintet": {
                "MW_pole_gev": v_value * math.sqrt(math.pi * alpha_2_prime),
                "MZ_pole_gev": v_value * math.sqrt(math.pi * alpha_sum_prime),
                "alpha_em_eff_inv": alpha_sum_prime / (alpha_y_prime * alpha_2_prime),
                "sin2w_eff": alpha_y_prime / alpha_sum_prime,
                "v_report_gev": v_value,
            },
        },
        "notes": [
            "This artifact is retained as historical scaffolding beneath the promoted target-free D10 repair theorem.",
            "It records the sharpest earlier split between source-only underdetermination and the smallest honest conditional closure route.",
            "The W/Z mass lane and the Ward-projected electromagnetic theorem sit above this split on separate D10 surfaces.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 minimal conditional-promotion artifact.")
    parser.add_argument("--source-pair", default=str(SOURCE_PAIR_JSON))
    parser.add_argument("--color-count", type=int, default=DEFAULT_COLOR_COUNT)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = _load_json(Path(args.source_pair))
    artifact = build_artifact(source_pair, color_count=args.color_count)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
