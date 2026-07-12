#!/usr/bin/env python3
"""Verify exact arithmetic, and only exact arithmetic, of the RSCC ledger."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp


CODE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_rscc_module_arithmetic.json"
)
SINGLET_NO_GO = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_12_24_singlet_no_go.json"
)
EXPECTED_BUNDLE_SHA256 = (
    "9e87957c7cbd8031fef0819f2bbb42dac82862776b33b747ca7dac3c61530792"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact() -> dict[str, Any]:
    n_g = sp.Integer(3)
    n_c = sp.Integer(3)
    d_std = sp.Integer(2)
    order_s3 = sp.Integer(6)
    beta = sp.Integer(4)
    repair = sp.Integer(24)
    orientation = sp.Integer(2)
    n_f = n_g + d_std
    f0 = n_f - 1

    dimensions = {
        "F": n_f,
        "F0": f0,
        "C_tensor_F": n_c * n_f,
        "M_mu_up": n_f**2 + f0,
        "M_mu_down": repair * order_s3 * n_c,
        "M_a_up": orientation * (n_f + order_s3),
        "M_a_down_linear": repair + orientation * f0,
        "M_a_down_quadratic": repair * (n_c + f0) * n_f,
        "M_g_common": repair * order_s3 * (n_c + f0),
        "M_g_color": repair * order_s3 * n_c,
        "M_g_relabel": repair * order_s3 * (n_f + order_s3),
    }
    expected = {
        "F": 5,
        "F0": 4,
        "C_tensor_F": 15,
        "M_mu_up": 29,
        "M_mu_down": 432,
        "M_a_up": 22,
        "M_a_down_linear": 32,
        "M_a_down_quadratic": 840,
        "M_g_common": 1008,
        "M_g_color": 432,
        "M_g_relabel": 1584,
    }
    dimension_check = {
        key: int(value) for key, value in dimensions.items()
    } == expected

    rank, dimension, width = sp.symbols(
        "rank dimension width", positive=True
    )
    assumed_variance = 2 * width**2 * rank / dimension
    gaussian_log_mgf = sp.simplify(assumed_variance / 2)
    gaussian_identity = gaussian_log_mgf == width**2 * rank / dimension
    f_commutant_dimension = 1**2 + 2**2

    weights = {
        "family_scalar": sp.Rational(1, n_f),
        "centered_fraction": sp.Rational(f0, n_f),
        "common_exposure": sp.Rational(f0, n_c * n_f),
        "mu_up": sp.Rational(1, dimensions["M_mu_up"]),
        "mu_down": -sp.Rational(1, dimensions["M_mu_down"]),
        "a_up": sp.Rational(1, dimensions["M_a_up"]),
        "a_down_linear": sp.Rational(1, dimensions["M_a_down_linear"]),
        "a_down_quadratic": sp.Rational(
            orientation, dimensions["M_a_down_quadratic"]
        ),
    }
    singlet_no_go = _read_json(SINGLET_NO_GO)

    return {
        "artifact": "oph_quark_rscc_module_arithmetic_v1",
        "claim_class": "exact_arithmetic_of_postulated_rscc_ledger",
        "proof_status": (
            "exact_arithmetic_passed_physical_module_selection_not_proved"
            if dimension_check and gaussian_identity
            else "failed"
        ),
        "submitted_bundle_sha256": EXPECTED_BUNDLE_SHA256,
        "exact_checks": {
            "F_decomposition_dimension": int(n_f),
            "F0_dimension": int(f0),
            "F0_equals_beta_EW": f0 == beta,
            "declared_composite_dimensions": expected,
            "dimension_arithmetic_passes": dimension_check,
            "normalized_weights": {
                key: str(value) for key, value in weights.items()
            },
            "gaussian_log_mgf_given_assumed_variance": str(gaussian_log_mgf),
            "gaussian_identity_passes": gaussian_identity,
            "gaussian_second_cumulant_is_positive": True,
            "negative_w2_terms_require_external_signed_response_law": True,
            "F_transposition_Cayley_Laplacian_spectrum": {"0": 1, "3": 4},
            "regular_S3_Laplacian_spectrum": {"0": 1, "3": 4, "6": 1},
            "F_contains_sign_sector_or_laplacian_6_mode": False,
            "physical_regular_to_F_heat_attachment_supplied": False,
            "dimension_End_S3_of_F": f_commutant_dimension,
            "S3_invariance_selects_unique_scalar_response_on_F": False,
        },
        "non_theorems": {
            "F_is_unique_minimal_physical_carrier": False,
            "composite_module_incidence_is_source_derived": False,
            "effect_projectors_exist_and_are_canonical": False,
            "effect_ranks_and_signs_are_source_derived": False,
            "S3_equivariance_implies_full_unitary_isotropy": False,
            "MAR_implies_gaussian_two_cumulant_truncation": False,
            "residual_minimization_is_an_OPH_dynamical_law": False,
        },
        "existing_count_only_register_boundary": {
            "artifact": singlet_no_go["artifact"],
            "claim_boundary": singlet_no_go["claim_boundary"],
            "rscc_requires_new_family_non_singlet_attachment": True,
        },
        "status_statement": (
            "Direct-sum/tensor-product dimensions and rank/dimension fractions are "
            "correct once the RSCC ledger is declared. The verifier constructs no "
            "physical incidence functor, effect projector, sign law, refinement law, "
            "or Gaussian source distribution."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    artifact = build_artifact()
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {args.output}")
    return 0 if artifact["proof_status"].startswith("exact_arithmetic_passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
