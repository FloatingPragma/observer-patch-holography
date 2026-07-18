#!/usr/bin/env python3
"""Emit the blocked interval-certificate package for issue #235."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import json
from pathlib import Path
from typing import Any

from interval_backend import DecimalIntervalBackend, Interval
from paper_math import PaperMathContext, _dec, to_serializable
from thomson_spectral_transport import validate_source_transport_payload


DEFAULT_PACKAGE = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_package_current.json"
DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "fine_structure_interval_certificate_current.json"
DEFAULT_R_Q_OUT = Path(__file__).resolve().parent / "runtime" / "r_q_residual_contract_current.json"
DEFAULT_SOURCE_SPECTRAL_THEOREM = (
    Path(__file__).resolve().parent / "runtime" / "source_spectral_theorem_current.json"
)

DEFAULT_I_A_LO = Decimal("137.035999156")
DEFAULT_I_A_HI = Decimal("137.035999198")
DEFAULT_L_B_TARGET = Decimal("100")
DEFAULT_L_R_TARGET = Decimal("1000")


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _exact_packet(package: dict[str, Any]) -> dict[str, Any]:
    return package["codata_mapped_endpoint_packet"]["exact_one_loop_package"]


def _make_pixel_interval(ctx: PaperMathContext, backend: DecimalIntervalBackend, i_a: Interval) -> Interval:
    phi = backend.interval(ctx.phi)
    sqrt_pi = backend.interval(ctx.sqrt_pi)
    # A larger inverse alpha gives a smaller P.
    lower = backend.add(phi, backend.div(sqrt_pi, backend.interval(i_a.hi)))
    upper = backend.add(phi, backend.div(sqrt_pi, backend.interval(i_a.lo)))
    return backend.interval(lower.lo, upper.hi)


def build_interval_certificate(
    package: dict[str, Any],
    *,
    i_a_lo: Decimal = DEFAULT_I_A_LO,
    i_a_hi: Decimal = DEFAULT_I_A_HI,
    l_b_target: Decimal = DEFAULT_L_B_TARGET,
    l_r_target: Decimal = DEFAULT_L_R_TARGET,
    precision: int = 96,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the interval theorem scaffold and the residual-map contract.

    The output is intentionally blocked.  It computes the interval target that
    a source-only residual map must satisfy, plus the Banach derivative budget.
    It does not certify D10 roots or supply the residual map.
    """
    backend = DecimalIntervalBackend(precision=precision)
    ctx = PaperMathContext(precision=80, su2_cutoff=80, su3_cutoff=60)
    exact = _exact_packet(package)
    compare = package["codata_mapped_endpoint_packet"]

    with localcontext() as dec_ctx:
        dec_ctx.prec = precision
        i_a = backend.interval(i_a_lo, i_a_hi)
        i_p = _make_pixel_interval(ctx, backend, i_a)
        p_center = _dec(compare["compare_p_from_outer_equation"])
        b_center = _dec(exact["implemented_endpoint_alpha_inv"])
        residual_center = _dec(exact["missing_source_transport_delta_alpha_inv"])
        p_radius = max(abs(p_center - i_p.lo), abs(i_p.hi - p_center))
        b_radius_target = l_b_target * p_radius
        b_image_target = backend.interval(b_center - b_radius_target, b_center + b_radius_target)
        r_required = backend.interval(i_a.lo - b_image_target.hi, i_a.hi - b_image_target.lo)
        a_t_image_if_r_supplied = backend.add(b_image_target, r_required)
        kappa_target = ctx.sqrt_pi * (l_b_target + l_r_target) / (i_a.lo * i_a.lo)
        r_prime_ceiling_if_lb = (i_a.lo * i_a.lo / ctx.sqrt_pi) - l_b_target

    r_q_contract = to_serializable(
        {
            "artifact": "oph_r_q_residual_map_contract",
            "generated_utc": _now_utc(),
            "github_issue": 235,
            "status": "closed_blocker_isolated_missing_source_R_Q",
            "promotion_allowed": False,
            "issue_235_resolution": {
                "close_recommendation": "close_as_first_missing_lemma_isolated",
                "exact_alpha_promoted": False,
                "first_missing_lemma": "source-emitted same-scheme Ward-projected R_Q(P)",
                "minimal_new_theorem": "WardProjectedHadronicSpectralEmission_Q",
            },
            "source_only_guard": {
                "codata_used_to_construct_R_Q": False,
                "measured_endpoint_used_to_construct_R_Q": False,
                "comparison_values_allowed_only_for_target_interval": True,
            },
            "required_map": "R_Q(P)=Delta_had_src(P)+Delta_EW_src(P)-implemented_free_quark_screen_component(P)",
            "source_spectral_theorem_artifact": str(
                DEFAULT_SOURCE_SPECTRAL_THEOREM.relative_to(DEFAULT_OUT.parent.parent)
            ),
            "source_spectral_theorem_status": "source_spectral_reduction_theorem_emitted_measure_payload_absent",
            "required_image_on_I_P": r_required.to_json(),
            "required_center_at_compare_pixel": +residual_center,
            "required_prime_abs_target": +l_r_target,
            "same_branch_requirements": [
                "same_D10_source_family_as_a0(P)",
                "same_Ward_projected_U1_Q_current",
                "same_matching_scheme_as_B(P)",
                "certified_quadrature_or_tail_bound",
            ],
            "forbidden_inputs": [
                "CODATA/NIST alpha",
                "compare_alpha_inv",
                "P_C residual",
                "S_required",
                "c_Q target",
            ],
        }
    )

    certificate = to_serializable(
        {
            "artifact": "oph_fine_structure_interval_certificate",
            "generated_utc": _now_utc(),
            "github_issue": 235,
            "status": "closed_blocker_isolated_missing_R_Q",
            "promotion_allowed": False,
            "issue_235_resolution": {
                "close_recommendation": "close_as_first_missing_lemma_isolated",
                "exact_alpha_promoted": False,
                "first_missing_lemma": "source-emitted same-scheme Ward-projected R_Q(P)",
                "minimal_new_theorem": "WardProjectedHadronicSpectralEmission_Q",
            },
            "interval_backend": {
                "library": "DecimalIntervalBackend",
                "precision_bits": None,
                "precision_decimal_digits": precision,
                "rounding": "high_precision_decimal_with_explicit_padding_not_theorem_grade",
                "promotion_backend_required": "arb_or_mpfi_directed_outward",
            },
            "constants": {
                "phi": backend.interval(ctx.phi).to_json(),
                "sqrt_pi": backend.interval(ctx.sqrt_pi).to_json(),
                "I_A": i_a.to_json(),
                "I_P": i_p.to_json(),
            },
            "B_certificate": {
                "status": "target_bound_not_interval_AD_certificate",
                "definition": "a0 + Delta_lep_closed_form + (1 - Nc alpha3/pi) Delta_q_naive_closed_form",
                "kernel": "closed_form_one_loop",
                "quadrature_error": backend.interval(0).to_json(),
                "B_center_at_compare_pixel": +b_center,
                "B_image_target_using_L_B": b_image_target.to_json(),
                "B_prime_abs_target": +l_b_target,
                "requires_interval_D10_root_boxes": True,
            },
            "R_Q_certificate": {
                "status": "missing_source_artifact",
                "contract_artifact": DEFAULT_R_Q_OUT.relative_to(DEFAULT_OUT.parent.parent).as_posix(),
                "R_Q_image_required": r_required.to_json(),
                "R_Q_prime_abs_target": +l_r_target,
                "source_only": True,
                "codata_used": False,
                "same_scheme_as_B_required": True,
            },
            "composition": {
                "A_T_image_if_R_Q_contract_is_met": a_t_image_if_r_supplied.to_json(),
                "self_map_conditional_on_R_Q": a_t_image_if_r_supplied.subset_of(i_a),
                "contraction_kappa_with_target_bounds": +kappa_target,
                "banach_conditional_pass": kappa_target < Decimal(1),
                "R_Q_prime_ceiling_if_L_B_target_holds": +r_prime_ceiling_if_lb,
            },
            "blocking_items": [
                "populated source spectral measure payload",
                "source-emitted R_Q(P) map",
                "directed-rounding interval backend",
                "interval D10 root boxes and derivative bounds",
                "same-scheme Ward-projected spectral quadrature bound",
            ],
            "conclusion": {
                "unique_fixed_point_in_I_P": False,
                "reason": "The Banach theorem is conditional until R_Q(P) and theorem-grade interval arithmetic are supplied. The current corpus isolates that source object but does not emit it.",
            },
        }
    )
    return certificate, r_q_contract


def build_source_payload_interval_certificate(source_transport_payload: dict[str, Any]) -> dict[str, Any]:
    """Build the source-only interval theorem certificate for a supplied payload.

    Unlike ``build_interval_certificate``, this path does not take a CODATA or
    compare-alpha interval.  The source payload must carry the endpoint image,
    derivative bound, directed-rounding backend certificate, and fixed-point
    enclosure.
    """
    validation = validate_source_transport_payload(source_transport_payload).to_json()
    interval_certificate = validation.get("interval_certificate", {})
    promotion_allowed = bool(validation["promotion_allowed"])
    return to_serializable(
        {
            "artifact": "oph_source_fine_structure_interval_certificate",
            "generated_utc": _now_utc(),
            "github_issue": 235,
            "status": (
                "source_interval_theorem_closed_for_supplied_payload"
                if promotion_allowed
                else "blocked_source_interval_payload_incomplete"
            ),
            "promotion_allowed": promotion_allowed,
            "external_inputs_used": False,
            "codata_or_nist_input_used": False,
            "source_payload_validation": validation,
            "exact_interval_theorem": interval_certificate.get(
                "theorem",
                {
                    "map": "G(P)=phi+sqrt(pi)/A_T(P)",
                    "self_map_condition": "G(I_P) subset I_P",
                    "uniqueness_condition": "sqrt(pi) * sup|A_T'(P)| / inf(A_T(I_P))^2 < 1",
                },
            ),
            "certified_intervals": interval_certificate.get("certified_intervals"),
            "bounds": interval_certificate.get("bounds"),
            "conclusion": {
                "exact_alpha_emitted_from_supplied_payload": promotion_allowed,
                "reason": (
                    "The supplied source payload carries a theorem-grade interval fixed-point "
                    "certificate."
                    if promotion_allowed
                    else "Exact alpha still requires the populated source spectral measure payload, same-scheme remainder, and interval certificate fields."
                ),
            },
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the blocked fine-structure interval certificate.")
    parser.add_argument("--package", default=str(DEFAULT_PACKAGE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--r-q-output", default=str(DEFAULT_R_Q_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    package = json.loads(Path(args.package).read_text(encoding="utf-8"))
    certificate, r_q_contract = build_interval_certificate(package)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    r_q_path = Path(args.r_q_output)
    r_q_path.parent.mkdir(parents=True, exist_ok=True)
    r_q_path.write_text(json.dumps(r_q_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.print_json:
        print(json.dumps(certificate, indent=2, sort_keys=True))
    else:
        print(f"saved: {out_path}")
        print(f"saved: {r_q_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
