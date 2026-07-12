#!/usr/bin/env python3
"""Compare and provenance-audit the quarantined S3/D12 quark ansatz.

All target values are confined to this compare-only module.  The reported
residual sum is descriptive: the rows mix schemes/scales, the formula was
discovered from the same targets, and neither source covariance nor model
selection is included.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from quark_s3_d12_template_postdiction import (
    CHARGED_BUDGET,
    CLAIM_CLASS,
    FAMILY_KERNEL,
    evaluate,
)


CODE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_s3_d12_template_postdiction_audit.json"
)
PREDICTOR = Path(__file__).resolve().with_name("quark_s3_d12_template_postdiction.py")
EXPECTED_BUNDLE_SHA256 = "bd283ee0d7cccd9a9592e69f74598a3ea1590274d7c5e6144f25f7920880e6e4"

# Frozen comparison packet supplied with the retrospective bundle.  These are
# not one common-scale running-mass sextet.
REFERENCE = {
    "u": {"value": 0.002160, "sigma": 0.000070, "chart": "MSbar at 2 GeV"},
    "d": {"value": 0.004700, "sigma": 0.000070, "chart": "MSbar at 2 GeV"},
    "s": {"value": 0.092900, "sigma": 0.000700, "chart": "MSbar at 2 GeV"},
    "c": {"value": 1.272900, "sigma": 0.004500, "chart": "MSbar at self-scale"},
    "b": {"value": 4.186000, "sigma": 0.006000, "chart": "MSbar at self-scale"},
    "t": {"value": 172.100000, "sigma": 0.600000, "chart": "separate top extraction"},
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _literal_firewall() -> dict[str, Any]:
    source = PREDICTOR.read_text(encoding="utf-8")
    tree = ast.parse(source)
    constants = {
        str(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, (str, int, float))
    }
    target_tokens = {
        str(row[field]) for row in REFERENCE.values() for field in ("value", "sigma")
    }
    target_hits = sorted(target_tokens & constants)
    path_hits = sorted(token for token in constants if "particle_reference_values" in token)
    return {
        "predictor_sha256": _sha256(PREDICTOR),
        "exact_target_literal_hits": target_hits,
        "reference_path_literal_hits": path_hits,
        "runtime_comparison_separation_passes": not target_hits and not path_hits,
        "scope": (
            "This checks direct evaluator separation only. It does not establish clean ancestry, "
            "blind discovery, or source-only status."
        ),
    }


def _ancestry_checks(prediction: dict[str, Any]) -> dict[str, Any]:
    kernel = _read_json(FAMILY_KERNEL)
    charged = _read_json(CHARGED_BUDGET)
    last_refinement = kernel["refinements"][-1]
    eigenvalues = [float(value) for value in last_refinement["eigenvalues"]]
    gaps = [float(value) for value in last_refinement["spectral_gaps"]]
    reconstructed_gch = sum(eigenvalues) / len(eigenvalues) + min(gaps)
    stored_gch = float(prediction["repository_inputs"]["g_ch_template"])
    stored_status = charged["charged_dirac_scalarization_certificate"][
        "stored_shared_absolute_scale_status"
    ]

    return {
        "family_kernel_is_handwritten_template": kernel.get("status") == "template",
        "family_kernel_metadata_requests_replacement": "Replace this template" in str(
            kernel.get("metadata", {}).get("note", "")
        ),
        "g_ch_reconstructed_as_template_mean_plus_min_gap": reconstructed_gch,
        "g_ch_stored": stored_gch,
        "g_ch_identity_residual": reconstructed_gch - stored_gch,
        "g_ch_is_dimensionless": True,
        "stored_scale_status": stored_status,
        "five_template_descendants_consumed": [
            "S_13",
            "S_23",
            "delta_21",
            "DeltaS_13",
            "g_ch_template",
        ],
        "p_source_uses_stage5_quark_model": prediction["provenance"][
            "p_source_uses_stage5_quark_model"
        ],
        "source_only_ancestry_passes": False,
    }


def build_audit(bundle: Path | None = None) -> dict[str, Any]:
    prediction = evaluate()
    coordinates = {
        name: float(value)
        for name, value in prediction["dimensionless_output_coordinates"].items()
    }

    rows: dict[str, Any] = {}
    residual_sum = 0.0
    max_relative = 0.0
    for name in ("u", "d", "s", "c", "b", "t"):
        coordinate = coordinates[name]
        reference = float(REFERENCE[name]["value"])
        sigma = float(REFERENCE[name]["sigma"])
        relative = coordinate / reference - 1.0
        nominal_z = (coordinate - reference) / sigma
        residual_sum += nominal_z * nominal_z
        max_relative = max(max_relative, abs(relative))
        rows[name] = {
            "dimensionless_coordinate": coordinate,
            "reference_gev": reference,
            "reference_sigma_gev": sigma,
            "reference_chart": REFERENCE[name]["chart"],
            "relative_error_percent_after_unproved_1GeV_unit_insertion": 100 * relative,
            "nominal_z_after_unproved_1GeV_unit_insertion": nominal_z,
        }

    up_gap_ratio = math.log(coordinates["c"] / coordinates["u"]) / math.log(
        coordinates["t"] / coordinates["c"]
    )
    down_gap_ratio = math.log(coordinates["s"] / coordinates["d"]) / math.log(
        coordinates["b"] / coordinates["s"]
    )

    bundle_receipt: dict[str, Any] = {
        "expected_sha256": EXPECTED_BUNDLE_SHA256,
        "path_checked": None,
        "matches": None,
    }
    if bundle is not None:
        bundle_receipt.update(
            {
                "path_checked": str(bundle),
                "observed_sha256": _sha256(bundle),
                "matches": _sha256(bundle) == EXPECTED_BUNDLE_SHA256,
            }
        )

    return {
        "artifact": "oph_quark_s3_d12_template_postdiction_audit_v2",
        "claim_class": CLAIM_CLASS,
        "promotion_allowed": False,
        "submitted_bundle_receipt": bundle_receipt,
        "runtime_separation_audit": _literal_firewall(),
        "ancestry_audit": _ancestry_checks(prediction),
        "descriptive_mixed_chart_comparison": {
            "rows": rows,
            "max_abs_relative_error_percent": 100 * max_relative,
            "raw_diagonal_residual_sum": residual_sum,
            "rms_nominal_z": math.sqrt(residual_sum / 6),
            "statistical_interpretation_allowed": False,
            "reasons": [
                "formula and selectors were discovered from the comparison spectrum",
                "rows do not share one renormalization scale or mass convention",
                "the dimensionless template scale has no derived GeV bridge",
                "no source-theory uncertainty or covariance matrix is supplied",
                "the postdeclared denominator grammar omits continuous and structural choices",
            ],
        },
        "shape_diagnostic": {
            "mixed_chart_up_adjacent_log_gap_ratio": up_gap_ratio,
            "mixed_chart_down_adjacent_log_gap_ratio": down_gap_ratio,
            "common_scale_validation_performed": False,
            "note": (
                "Adding the universal Q mode abandons the old reciprocal-ray constraint; "
                "without RG transport this produces no common-scale shape prediction."
            ),
        },
        "local_grammar_audit_correction": {
            "postdeclared_models": 219615,
            "chosen_is_unique_best_raw_residual_sum": True,
            "models_tied_for_best_maximum_error": 8,
            "models_below_0p3_percent_maximum_error": 9,
            "models_below_1_percent_maximum_error": 252,
            "global_look_elsewhere_correction": False,
        },
        "typing_and_invariance_failures": [
            "S_ij log-overlap suppressions are added to a linearly scaling TT-dagger eigenvalue gap",
            "DeltaS_13 is a selected basis-dependent matrix entry, not a conjugation-invariant scalar",
            "the dimensionless g_ch template statistic is used as though it carried GeV units",
            "no theorem maps S3 heat isotypic values of multiplicities 1,4,1 to three generations",
        ],
        "verdict": (
            "The arithmetic is reproducible as a retrospectively designed diagnostic, but the "
            "bundle is not an OPH quark-mass closure, postdiction, or physical Yukawa evaluation."
        ),
        "prediction": prediction,
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
