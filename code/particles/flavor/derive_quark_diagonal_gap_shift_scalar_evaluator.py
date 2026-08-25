#!/usr/bin/env python3
"""Emit the quark diagonal gap-shift scalar-evaluator artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"
DEFAULT_SPREAD = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_SOURCE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"
DEFAULT_SOURCE_VALUES = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_values.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    diagonal_gap_shift_map: dict,
    spread_map: dict,
    source_law: dict | None = None,
    source_values: dict | None = None,
) -> dict:
    b_ord = [float(value) for value in diagonal_gap_shift_map["B_ord"]]
    x2 = float(spread_map["normalized_coordinate_x2"])
    source_law = source_law or {}
    source_values = source_values or {}
    beta_u = source_values.get("beta_u_diag_B_source")
    beta_d = source_values.get("beta_d_diag_B_source")
    values_present = beta_u is not None and beta_d is not None
    source_values_closed = (
        values_present
        and source_values.get("proof_status")
        == "source_values_derived_from_source_emission"
        and source_values.get("value_classification") == "source_emitted"
    )
    sigma_u = float(spread_map["sigma_u_total_log_per_side"])
    sigma_d = float(spread_map["sigma_d_total_log_per_side"])
    sigma_source_closed = (
        spread_map.get("predictive_promotion_allowed") is True
        and spread_map.get("spread_emitter_status")
        in {
            "closed_source_emission",
            "source_emitted",
            "source_values_derived_from_source_emission",
        }
    )
    source_value_promotion_ready = source_values_closed and sigma_source_closed
    comparison_only = values_present and not source_value_promotion_ready
    delta_ud_overlap = (
        None if not values_present else float(beta_u) - float(beta_d)
    )
    denominator = 2.0 * (sigma_u + sigma_d)
    tau_u = (
        None
        if delta_ud_overlap is None
        else sigma_d * delta_ud_overlap / denominator
    )
    tau_d = (
        None
        if delta_ud_overlap is None
        else sigma_u * delta_ud_overlap / denominator
    )
    value_classification = (
        "source_emitted"
        if source_value_promotion_ready
        else "comparison_only"
        if comparison_only
        else "absent"
    )
    missing_objects: list[str] = []
    if not source_values_closed:
        missing_objects.append(
            source_values.get(
                "smallest_constructive_missing_object",
                "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
            )
        )
    if not sigma_source_closed:
        missing_objects.append("supported_source_sigma_branch_emission")
    missing_object = missing_objects[0] if missing_objects else None
    return {
        "artifact": "oph_family_excitation_diagonal_gap_shift_scalar_evaluator",
        "generated_utc": _timestamp(),
        "proof_status": (
            "current_family_scalar_law_closed_with_source_derived_values"
            if source_value_promotion_ready
            else "comparison_only_scalar_projection_value_source_open"
            if comparison_only
            else "current_family_scalar_law_closed_waiting_source_readback_pair"
        ),
        "predictive_promotion_allowed": False,
        "source_value_promotion_ready": source_value_promotion_ready,
        "value_classification": value_classification,
        "promotion_blockers": missing_objects,
        "input_kind": "ordered_branch_generator_spectral_package_plus_diagonal_gap_shift_carrier",
        "output_kind": "sector_diagonal_gap_shift_scalar",
        "source_artifact": diagonal_gap_shift_map.get("artifact"),
        "source_law_artifact": source_law.get("artifact"),
        "source_law_status": source_law.get("proof_status"),
        "source_values_artifact": source_values.get("artifact"),
        "source_values_status": source_values.get("proof_status"),
        "source_emission_artifact": source_values.get("source_emission_artifact"),
        "source_emission_status": source_values.get("source_emission_status"),
        "spread_artifact": spread_map.get("artifact"),
        "spread_emitter_status": spread_map.get("spread_emitter_status"),
        "spread_predictive_promotion_allowed": spread_map.get(
            "predictive_promotion_allowed"
        ),
        "sigma_source_closed": sigma_source_closed,
        "sigma_u_total_log_per_side": sigma_u,
        "sigma_d_total_log_per_side": sigma_d,
        "B_ord": b_ord,
        "B_ord_norm_sq": sum(value * value for value in b_ord),
        "normalized_coordinate_x2": x2,
        "beta_u_diag_B_source": beta_u,
        "beta_d_diag_B_source": beta_d,
        "Delta_ud_overlap": delta_ud_overlap,
        "Delta_ud_overlap_formula": "beta_u_diag_B_source - beta_d_diag_B_source",
        "delta_gamma21_u_log_per_side": tau_u,
        "delta_gamma32_u_log_per_side": tau_u,
        "delta_gamma21_d_log_per_side": tau_d,
        "delta_gamma32_d_log_per_side": tau_d,
        "diagonal_common_gap_shift_u_log_per_side": tau_u,
        "diagonal_common_gap_shift_d_log_per_side": tau_d,
        "tau_u_log_per_side": tau_u,
        "tau_d_log_per_side": tau_d,
        "tau_u_formula": "sigma_d_total_log_per_side * Delta_ud_overlap / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
        "tau_d_formula": "sigma_u_total_log_per_side * Delta_ud_overlap / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
        "tau_projection_u_formula": "dot(delta_E_u_diag_log_per_side, B_ord) / B_ord_norm_sq",
        "tau_projection_d_formula": "dot(delta_E_d_diag_log_per_side, B_ord) / B_ord_norm_sq",
        "delta_sigma_u_diag_total_log_per_side_formula": "delta_gamma21_u_log_per_side + delta_gamma32_u_log_per_side",
        "delta_sigma_d_diag_total_log_per_side_formula": "delta_gamma21_d_log_per_side + delta_gamma32_d_log_per_side",
        "delta_a_u_formula": "tau_u_log_per_side",
        "delta_a_d_formula": "tau_d_log_per_side",
        "delta_b_u_formula": "normalized_coordinate_x2 * tau_u_log_per_side / (1 - normalized_coordinate_x2^2)",
        "delta_b_d_formula": "normalized_coordinate_x2 * tau_d_log_per_side / (1 - normalized_coordinate_x2^2)",
        "equal_gap_guard_u_formula": "delta_gamma21_u_log_per_side - delta_gamma32_u_log_per_side",
        "equal_gap_guard_d_formula": "delta_gamma21_d_log_per_side - delta_gamma32_d_log_per_side",
        "coefficient_compatibility_u_formula": "(1 - normalized_coordinate_x2^2) * delta_b_u - normalized_coordinate_x2 * delta_a_u",
        "coefficient_compatibility_d_formula": "(1 - normalized_coordinate_x2^2) * delta_b_d - normalized_coordinate_x2 * delta_a_d",
        "smallest_constructive_missing_object": missing_object,
        "constructive_missing_objects": missing_objects,
        "notes": [
            "The diagonal lift beneath the current closed surface is one-dimensional on B_ord = [-1, 0, 1].",
            "The live builder consumes the diagonal lift first as the tau-pair tau_u and tau_d on the B-ordered branch.",
            "The source-side beta pair fixes Delta_ud_overlap = beta_u - beta_d; it is not itself the transport tau pair.",
            "The exact D12 overlap transport law weights tau_u and tau_d by the opposite-sector sigma totals on a fixed positive sigma branch.",
            "Numerical values remain comparison-only unless both the beta pair and the sigma branch are source-emitted; their presence alone cannot close this evaluator.",
            "No PDG quark values are consumed here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal gap-shift scalar-evaluator artifact.")
    parser.add_argument("--map", default=str(DEFAULT_MAP))
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD))
    parser.add_argument("--source-law", default=str(DEFAULT_SOURCE_LAW))
    parser.add_argument("--source-values", default=str(DEFAULT_SOURCE_VALUES))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    diagonal_gap_shift_map = json.loads(Path(args.map).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    source_law_path = Path(args.source_law)
    source_law = json.loads(source_law_path.read_text(encoding="utf-8")) if source_law_path.exists() else None
    source_values_path = Path(args.source_values)
    source_values = json.loads(source_values_path.read_text(encoding="utf-8")) if source_values_path.exists() else None
    artifact = build_artifact(diagonal_gap_shift_map, spread_map, source_law, source_values)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
