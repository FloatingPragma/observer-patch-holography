#!/usr/bin/env python3
"""Emit the screening-invariant non-identifiability certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
from typing import Any

from paper_math import _dec, to_serializable


DEFAULT_PACKAGE = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_package_current.json"
DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "screening_invariant_no_go_current.json"


def _canonical_hash(payload: dict[str, Any]) -> str:
    data = json.dumps(to_serializable(payload), sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def _exact_packet(package: dict[str, Any]) -> dict[str, Any]:
    return package["codata_mapped_endpoint_packet"]["exact_one_loop_package"]


def _d10_packet(package: dict[str, Any]) -> dict[str, Any]:
    return package["codata_mapped_endpoint_packet"]["d10_source_point"]


def build_no_go(package: dict[str, Any]) -> dict[str, Any]:
    """Build the source-invariant no-go certificate.

    The proof is by non-identifiability.  At a fixed D10 source point the
    current source packet fixes the anchor, lepton term, naive quark term, and
    first-order screen.  A free lambda in S_lambda=1-x+lambda x^2 leaves that
    source packet unchanged while changing the endpoint.
    """
    exact = _exact_packet(package)
    d10 = _d10_packet(package)
    screen = exact["screening_scalar"]

    with localcontext() as ctx:
        ctx.prec = 96
        lepton_delta = _dec(exact["lepton_delta_alpha_inv"])
        quark_naive = _dec(exact["quark_delta_alpha_inv_naive"])
        implemented_delta = _dec(exact["implemented_transport_delta_alpha_inv"])
        required_delta = _dec(exact["required_transport_delta_alpha_inv"])
        missing_delta = _dec(exact["missing_source_transport_delta_alpha_inv"])
        implemented_screen = _dec(screen["implemented_screening_factor"])
        required_screen = _dec(screen["required_screening_factor"])
        qcd_x = _dec(screen["qcd_x"])
        c_q = _dec(screen["residual_second_order_coefficient"])
        alpha_u = _dec(d10["alpha_U"])

        beta_candidate = Decimal(2) / Decimal(3)
        beta_candidate_gap = quark_naive * qcd_x * qcd_x * (beta_candidate - c_q)
        alpha_u_gap = missing_delta - alpha_u

        lambda_a = Decimal(0)
        lambda_b = beta_candidate
        endpoint_delta_a = lepton_delta + quark_naive * (Decimal(1) - qcd_x + lambda_a * qcd_x * qcd_x)
        endpoint_delta_b = lepton_delta + quark_naive * (Decimal(1) - qcd_x + lambda_b * qcd_x * qcd_x)
        target_equivalence_tolerance = Decimal(10) ** Decimal(-68)
        screen_equivalent = quark_naive * (required_screen - implemented_screen)
        cq_equivalent = quark_naive * qcd_x * qcd_x * c_q

    source_packet = {
        "P": d10["P"],
        "source_anchor_alpha_inv_mz": exact["source_anchor_alpha_inv_mz"],
        "alpha_U": d10["alpha_U"],
        "alpha1_mz": d10["alpha1_mz"],
        "alpha2_mz": d10["alpha2_mz"],
        "alpha3_mz": d10["alpha3_mz"],
        "alphaY_mz": d10["alphaY_mz"],
        "mZ_run": d10["mZ_run"],
        "v": d10["v"],
        "lepton_delta_alpha_inv": exact["lepton_delta_alpha_inv"],
        "quark_delta_alpha_inv_naive": exact["quark_delta_alpha_inv_naive"],
        "implemented_screening_factor": screen["implemented_screening_factor"],
        "qcd_x": screen["qcd_x"],
        "transport_kernel": exact["transport_kernel"],
        "kernel_evaluation": exact.get("kernel_evaluation", "closed_form_one_loop"),
    }

    return to_serializable(
        {
            "artifact": "oph_screening_invariant_no_go_certificate",
            "claim_status": "source_invariant_surface_does_not_determine_screening_scalar",
            "promotion_allowed": False,
            "source_packet_hash": _canonical_hash(source_packet),
            "source_packet": source_packet,
            "compare_only_targets": {
                "compare_alpha_inv": package["codata_mapped_endpoint_packet"]["compare_alpha_inv"],
                "implemented_transport_delta_alpha_inv": +implemented_delta,
                "required_transport_delta_alpha_inv": +required_delta,
                "missing_source_transport_delta_alpha_inv": +missing_delta,
                "required_screening_factor": +required_screen,
                "residual_second_order_coefficient": +c_q,
            },
            "target_equivalences": {
                "missing_delta": +missing_delta,
                "quark_naive_times_screen_gap": +screen_equivalent,
                "quark_naive_times_x2_cq": +cq_equivalent,
                "tolerance": +target_equivalence_tolerance,
                "all_equal_to_precision": (
                    abs(missing_delta - screen_equivalent) <= target_equivalence_tolerance
                    and abs(missing_delta - cq_equivalent) <= target_equivalence_tolerance
                ),
            },
            "candidate_failures": {
                "c_Q_beta_EW_over_2Nc": {
                    "candidate": +beta_candidate,
                    "candidate_minus_required": +(beta_candidate - c_q),
                    "inverse_alpha_overshoot": +beta_candidate_gap,
                    "passes": False,
                },
                "missing_delta_equals_alpha_U": {
                    "candidate": +alpha_u,
                    "missing_minus_candidate": +alpha_u_gap,
                    "passes": False,
                },
            },
            "nonidentifiability_witness": {
                "family": "S_lambda(P)=1-x(P)+lambda*x(P)^2",
                "source_hash_independent_of_lambda": True,
                "lambda_a": +lambda_a,
                "lambda_b": +lambda_b,
                "endpoint_delta_lambda_a": +endpoint_delta_a,
                "endpoint_delta_lambda_b": +endpoint_delta_b,
                "endpoint_delta_difference": +(endpoint_delta_b - endpoint_delta_a),
                "formula_for_difference": "quark_delta_alpha_inv_naive*x^2*(lambda_b-lambda_a)",
            },
            "required_source_object": {
                "id": "ward_projected_qcd_screening_and_endpoint_remainder",
                "required_primitives": [
                    "rho_had(s;P)",
                    "Delta_EW(P)",
                    "same_scheme_as_a0(P)",
                    "certified_quadrature_or_tail_bound",
                    "interval_fixed_point_certificate",
                ],
            },
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the screening-invariant no-go certificate.")
    parser.add_argument("--package", default=str(DEFAULT_PACKAGE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    package = json.loads(Path(args.package).read_text(encoding="utf-8"))
    payload = build_no_go(package)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
