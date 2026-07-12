#!/usr/bin/env python3
"""Compare and provenance-audit the retrospective RSCC quark candidate."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp

from quark_rscc_completion_candidate import CLAIM_CLASS, evaluate


CODE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_rscc_completion_candidate_audit.json"
)
PREDICTOR = Path(__file__).resolve().with_name(
    "quark_rscc_completion_candidate.py"
)
OLD_TEMPLATE = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_s3_d12_template_postdiction.json"
)
MODULE_ARITHMETIC = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_rscc_module_arithmetic.json"
)
EXPECTED_BUNDLE_SHA256 = (
    "9e87957c7cbd8031fef0819f2bbb42dac82862776b33b747ca7dac3c61530792"
)

REFERENCE = {
    "u": {"value": 0.002160, "sigma": 0.000065, "chart": "MSbar at 2 GeV"},
    "d": {"value": 0.004700, "sigma": 0.000070, "chart": "MSbar at 2 GeV"},
    "s": {"value": 0.092900, "sigma": 0.000700, "chart": "MSbar at 2 GeV"},
    "c": {"value": 1.272900, "sigma": 0.004500, "chart": "MSbar at self-scale"},
    "b": {"value": 4.186000, "sigma": 0.006000, "chart": "MSbar at self-scale"},
    "t": {"value": 172.100000, "sigma": 0.600000, "chart": "separate top extraction"},
}


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _runtime_literal_firewall() -> dict[str, Any]:
    tree = ast.parse(PREDICTOR.read_text(encoding="utf-8"))
    constants = {
        str(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, (str, int, float))
    }
    target_tokens = {
        str(row[field]) for row in REFERENCE.values() for field in ("value", "sigma")
    }
    hits = sorted(constants & target_tokens)
    return {
        "predictor_sha256": _sha256(PREDICTOR),
        "exact_target_literal_hits": hits,
        "runtime_comparison_separation_passes": not hits,
        "scope": (
            "Direct runtime separation only; it does not clean Stage-5 ancestry, "
            "retrospective formula selection, or mixed-branch D10 scale use."
        ),
    }


def _implied_denominator_audit(
    candidate: dict[str, Any], old: dict[str, Any]
) -> dict[str, Any]:
    mp.mp.dps = 60
    P = mp.mpf(str(candidate["repository_inputs"]["P"]))
    width = mp.pi * mp.mpf(str(candidate["repository_inputs"]["alpha_U"]))
    old_coordinates = old["derived_scalars"]
    old_mean_u = mp.mpf(old_coordinates["mean_u"])
    old_mean_d = mp.mpf(old_coordinates["mean_d"])
    old_a_u = mp.mpf(old_coordinates["a_u"])
    old_a_d = mp.mpf(old_coordinates["a_d"])

    base_mean_u = 3 * P + (5 + mp.mpf(4) / 15) * width
    base_mean_d = 2 * P + mp.mpf(4) / 15 * width
    base_a_u = 3 * P + (5 + mp.mpf(4) / 5) * width
    base_a_d = 2 * P + width
    implied = {
        "mu_up": width**2 / (old_mean_u - base_mean_u),
        "mu_down": width**2 / (base_mean_d - old_mean_d),
        "a_up": width**2 / (old_a_u - base_a_u),
        "a_down_quadratic_given_linear_1_over_32": (
            2 * width**2 / (old_a_d - base_a_d - width / 32)
        ),
    }
    return {
        "target_informed_effective_denominators": {
            key: float(value) for key, value in implied.items()
        },
        "selected_ledger_denominators": {
            "mu_up": 29,
            "mu_down": 432,
            "a_up": 22,
            "a_down_quadratic": 840,
        },
        "postdeclared_module_expression_grammar_enumerated": False,
        "competitor_module_ledgers_excluded": False,
        "interpretation": (
            "The declared composite dimensions closely track the already visible "
            "effective coordinates. Direct sums, tensor products, ranks, signs, and "
            "sector assignments were not frozen in a prospective grammar."
        ),
    }


def build_audit(bundle: Path | None = None) -> dict[str, Any]:
    candidate = evaluate(include_d10_gev_display=True)
    display = candidate["conditional_D10_gev_display"]
    if display is None:
        raise AssertionError("conditional D10 display was not generated")
    coordinates = {
        key: float(value) for key, value in display["coordinates_gev"].items()
    }

    rows: dict[str, Any] = {}
    residual_sum = 0.0
    maximum_relative = 0.0
    for key in ("u", "d", "s", "c", "b", "t"):
        predicted = coordinates[key]
        reference = float(REFERENCE[key]["value"])
        sigma = float(REFERENCE[key]["sigma"])
        relative = predicted / reference - 1.0
        nominal_pull = (predicted - reference) / sigma
        residual_sum += nominal_pull**2
        maximum_relative = max(maximum_relative, abs(relative))
        rows[key] = {
            "conditional_D10_coordinate_gev": predicted,
            "reference_gev": reference,
            "reference_sigma_gev": sigma,
            "reference_chart": REFERENCE[key]["chart"],
            "relative_error_percent": 100 * relative,
            "nominal_pull": nominal_pull,
        }

    bundle_receipt: dict[str, Any] = {
        "expected_sha256": EXPECTED_BUNDLE_SHA256,
        "path_checked": None,
        "matches": None,
    }
    if bundle is not None:
        observed = _sha256(bundle)
        bundle_receipt.update(
            {
                "path_checked": str(bundle),
                "observed_sha256": observed,
                "matches": observed == EXPECTED_BUNDLE_SHA256,
            }
        )

    old = _read_json(OLD_TEMPLATE)
    module_arithmetic = _read_json(MODULE_ARITHMETIC)
    provenance = candidate["provenance"]
    ablation = evaluate(
        include_d10_gev_display=True,
        correction_mode="lower_order_ablation",
    )
    ablation_display = ablation["conditional_D10_gev_display"]
    if ablation_display is None:
        raise AssertionError("lower-order D10 ablation was not generated")
    ablation_coordinates = {
        key: float(value)
        for key, value in ablation_display["coordinates_gev"].items()
    }
    ablation_residual_sum = 0.0
    ablation_maximum_relative = 0.0
    for key, row in REFERENCE.items():
        reference = float(row["value"])
        sigma = float(row["sigma"])
        predicted = ablation_coordinates[key]
        ablation_maximum_relative = max(
            ablation_maximum_relative,
            abs(predicted / reference - 1.0),
        )
        ablation_residual_sum += ((predicted - reference) / sigma) ** 2
    return {
        "artifact": "oph_quark_rscc_completion_candidate_audit_v1",
        "claim_class": CLAIM_CLASS,
        "promotion_allowed": False,
        "submitted_bundle_receipt": bundle_receipt,
        "runtime_separation_audit": _runtime_literal_firewall(),
        "descriptive_mixed_chart_comparison": {
            "rows": rows,
            "max_abs_relative_error_percent": 100 * maximum_relative,
            "raw_diagonal_residual_sum": residual_sum,
            "rms_nominal_pull": math.sqrt(residual_sum / 6),
            "statistical_interpretation_allowed": False,
            "comparison_packet_note": (
                "The submitted RSCC packet changes the u-row nominal sigma from "
                "0.000070 GeV in the preceding bundle to 0.000065 GeV."
            ),
        },
        "target_informed_ledger_audit": _implied_denominator_audit(candidate, old),
        "negative_control": {
            "description": (
                "Remove every RSCC w^2 correction and the whole delta_g correction "
                "while retaining the lower-order P,w scaffold and inherited readout."
            ),
            "max_abs_relative_error_percent": 100 * ablation_maximum_relative,
            "raw_diagonal_residual_sum": ablation_residual_sum,
            "beats_full_rscc_maximum_error": (
                ablation_maximum_relative < maximum_relative
            ),
            "beats_full_rscc_raw_residual_sum": (
                ablation_residual_sum < residual_sum
            ),
            "interpretation": (
                "The submitted quadratic module ledger is not selected even by its "
                "own retrospective comparison metric."
            ),
        },
        "exact_content_retained": {
            "artifact": module_arithmetic["artifact"],
            "proof_status": module_arithmetic["proof_status"],
            "dimension_arithmetic_passes": module_arithmetic["exact_checks"][
                "dimension_arithmetic_passes"
            ],
            "gaussian_identity_is_conditional_on_assumed_variance": True,
        },
        "ancestry_and_branch_audit": {
            "pixel_branch_status": provenance["pixel_branch_status"],
            "pixel_source_uses_internal_stage5_quark_model": provenance[
                "pixel_source_uses_internal_stage5_quark_model"
            ],
            "d10_active_readout_selector_status": provenance[
                "d10_active_readout_selector_status"
            ],
            "d10_mixed_sources_detected": provenance["d10_mixed_sources_detected"],
            "rscc_P_matches_d10_P": provenance["rscc_P_matches_d10_P"],
            "rscc_alpha_U_matches_d10_alpha_U": provenance[
                "rscc_alpha_U_matches_d10_alpha_U"
            ],
            "source_only_ancestry_passes": False,
        },
        "precision_audit": {
            "alpha_U_input_decimal_digits": len(
                str(candidate["repository_inputs"]["alpha_U"]).replace(".", "")
            ),
            "D10_v_input_decimal_digits": len(
                str(candidate["repository_inputs"]["conditional_D10_v_gev"]).replace(
                    ".", ""
                )
            ),
            "printed_coordinate_digits_are_uncertainty_certified": False,
            "status": (
                "High-precision arithmetic makes the declared formulas deterministic; "
                "it does not create input accuracy or a theory uncertainty interval."
            ),
        },
        "physical_nonclosure": {
            "F1_clean_source_root": False,
            "F2_refinement_natural_family_carrier": False,
            "F3_physical_channel_functor": False,
            "F4_affine_mean_law": False,
            "F5_common_scale_physical_readout": False,
            "F6_RG_threshold_scheme_packet": False,
            "module_minimality_or_uniqueness_proved": False,
            "effect_projectors_ranks_and_signs_constructed": False,
            "full_unitary_isotropy_derived_from_S3": False,
            "MAR_two_cumulant_truncation_proved": False,
            "count_only_24_slot_register_emits_family_non_singlet": False,
            "rescaling_residual_is_a_derived_OPH_minimization_law": False,
        },
        "rescaling_orbit_correction": (
            "For an orbit through the declared E_q, ||lambda E_q-E_q||^2 has "
            "its minimum at lambda=1. For a generic base spectrum v, "
            "||lambda v-E_q||^2 is minimized at <v,E_q>/||v||^2, not generally 1. "
            "No current OPH theorem requires minimizing the postulated residual."
        ),
        "verdict": (
            "RSCC is a useful explicit, falsifiable parameterization of part of the "
            "missing flavor-source packet. It removes five legacy decimals from direct "
            "runtime inputs but replaces them with a target-informed module/effect/sign "
            "ledger and inherits the unresolved carrier, affine, scale, and RG laws. It "
            "does not discharge any physical Flavor Source Closure receipt."
        ),
        "candidate": candidate,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    artifact = build_audit(args.bundle)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
