#!/usr/bin/env python3
"""Compute the Ward-projected Thomson endpoint package.

The package is deliberately conditional.  It computes the endpoint residual
required by the CODATA/NIST comparison and expresses that residual in the
same variables used by the OPH P-closure code.  It does not let the comparison
constant enter the solver path.  The package isolates the missing
source-emitted residual map and interval-certificate theorem without treating
project state as scientific evidence.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
from typing import Any

from alpha_gap_audit import (
    CODATA_2022_ALPHA_INV,
    CODATA_2022_ALPHA_INV_UNCERTAINTY,
    DEFAULT_INTERVAL_CERTIFICATE,
    build_alpha_gap_audit,
    certified_source_point,
)
from paper_math import PaperMathContext, _dec, to_serializable


DEFAULT_REPORT = Path(__file__).resolve().parent / "runtime" / "full_p_alpha_report_current.json"
DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_package_current.json"
DEFAULT_ENDPOINT_PRECISION = 80


def _screening_payload(
    *,
    ctx: PaperMathContext,
    d10_alpha3_mz: Decimal,
    running: dict[str, Any],
    required_delta: Decimal,
) -> dict[str, Any]:
    lepton_delta = _dec(running["lepton_delta_alpha_inv"])
    quark_naive = _dec(running["quark_delta_alpha_inv_naive"])
    implemented_screen = _dec(running["quark_screening_factor"])
    qcd_x = Decimal(ctx.n_c) * d10_alpha3_mz / ctx.pi
    required_screen = (required_delta - lepton_delta) / quark_naive
    residual_screen = required_screen - implemented_screen
    residual_c2 = residual_screen / (qcd_x * qcd_x)
    return {
        "implemented_screening_factor": +implemented_screen,
        "required_screening_factor": +required_screen,
        "residual_screening_factor": +residual_screen,
        "qcd_x": +qcd_x,
        "residual_second_order_coefficient": +residual_c2,
        "definition": (
            "If x=N_c alpha_3(m_Z;P)/pi and the implemented screen is 1-x, "
            "then c_Q=(S_req-(1-x))/x^2 at the CODATA-mapped pixel point."
        ),
    }


def _endpoint_at_compare_point(
    *,
    ctx: PaperMathContext,
    compare_alpha_inv: Decimal,
    compare_alpha_inv_uncertainty: Decimal,
) -> dict[str, Any]:
    with localcontext() as dec_ctx:
        dec_ctx.prec = ctx.work_precision
        compare_alpha = ctx.one / compare_alpha_inv
        compare_p = ctx.outer_p_from_alpha(compare_alpha)
        compare_p_uncertainty = ctx.sqrt_pi * compare_alpha_inv_uncertainty / (
            compare_alpha_inv * compare_alpha_inv
        )
        d10 = ctx.build_d10_from_p(compare_p)
        exact_running = ctx.structured_thomson_running(d10)
        asymptotic_running = ctx.structured_thomson_running_asymptotic(d10)

        def package_for_running(running: dict[str, Any]) -> dict[str, Any]:
            implemented_delta = _dec(running["total_delta_alpha_inv"])
            implemented_endpoint = d10.alpha_em_inv_mz + implemented_delta
            required_delta = compare_alpha_inv - d10.alpha_em_inv_mz
            residual = compare_alpha_inv - implemented_endpoint
            return {
                "transport_kernel": running["transport_kernel"],
                "kernel_evaluation": running.get("kernel_evaluation"),
                "quadrature_error_bound": running.get("quadrature_error_bound"),
                "source_anchor_alpha_inv_mz": +d10.alpha_em_inv_mz,
                "implemented_transport_delta_alpha_inv": +implemented_delta,
                "implemented_endpoint_alpha_inv": +implemented_endpoint,
                "required_transport_delta_alpha_inv": +required_delta,
                "missing_source_transport_delta_alpha_inv": +residual,
                "lepton_delta_alpha_inv": +_dec(running["lepton_delta_alpha_inv"]),
                "quark_delta_alpha_inv_naive": +_dec(running["quark_delta_alpha_inv_naive"]),
                "quark_delta_alpha_inv_screened": +_dec(running["quark_delta_alpha_inv_screened"]),
                "screening_scalar": _screening_payload(
                    ctx=ctx,
                    d10_alpha3_mz=d10.alpha3_mz,
                    running=running,
                    required_delta=required_delta,
                ),
            }

        return {
            "compare_alpha_inv": +compare_alpha_inv,
            "compare_alpha_inv_standard_uncertainty": +compare_alpha_inv_uncertainty,
            "compare_alpha": +compare_alpha,
            "compare_p_from_outer_equation": +compare_p,
            "compare_p_standard_uncertainty": +compare_p_uncertainty,
            "d10_source_point": {
                "P": +d10.p,
                "M_U": +d10.mu_u,
                "mZ_run": +d10.mz_run,
                "v": +d10.v,
                "alpha_U": +d10.alpha_u,
                "alpha1_mz": +d10.alpha1_mz,
                "alpha2_mz": +d10.alpha2_mz,
                "alpha3_mz": +d10.alpha3_mz,
                "alphaY_mz": +d10.alpha_y_mz,
                "alpha_em_inv_mz": +d10.alpha_em_inv_mz,
                "sin2w_mz": +d10.sin2w_mz,
            },
            "exact_one_loop_package": package_for_running(exact_running),
            "asymptotic_package": package_for_running(asymptotic_running),
        }


def build_endpoint_package(
    report: dict[str, Any],
    *,
    source_point: dict[str, Any] | None = None,
    compare_alpha_inv: Decimal = CODATA_2022_ALPHA_INV,
    compare_alpha_inv_uncertainty: Decimal = CODATA_2022_ALPHA_INV_UNCERTAINTY,
    precision: int = DEFAULT_ENDPOINT_PRECISION,
    su2_cutoff: int = 80,
    su3_cutoff: int = 60,
) -> dict[str, Any]:
    """Build the conditional endpoint package for the P/alpha lane."""
    ctx = PaperMathContext(precision=precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff)
    implemented_root = build_alpha_gap_audit(
        report,
        source_point=source_point,
        compare_alpha_inv=compare_alpha_inv,
        compare_alpha_inv_uncertainty=compare_alpha_inv_uncertainty,
    )
    compare_point = _endpoint_at_compare_point(
        ctx=ctx,
        compare_alpha_inv=compare_alpha_inv,
        compare_alpha_inv_uncertainty=compare_alpha_inv_uncertainty,
    )
    exact_packet = compare_point["exact_one_loop_package"]
    return to_serializable(
        {
            "artifact": "oph_ward_projected_thomson_endpoint_package",
            "classification": "conditional_endpoint_package__source_transport_map_absent",
            "promotion_allowed": False,
            "compare_input_role": (
                "The CODATA/NIST value is used after the internal report is built, only to compute "
                "the residual endpoint packet that a source-only theorem would have to emit."
            ),
            "source_only_guard": {
                "codata_enters_solver": False,
                "hidden_external_alpha_allowed": False,
                "measured_endpoint_allowed_as_transport_input": False,
                "free_quark_screening_promotable": False,
            },
            "conditional_theorem": {
                "statement": (
                    "OPH axioms plus the realized Standard Model branch, the D10 source map, "
                    "a source-only Ward-projected Thomson transport map A_T(P), and an interval "
                    "contraction certificate for G(P)=phi+sqrt(pi)/A_T(P) imply a unique pixel "
                    "fixed point and alpha(P*)=1/A_T(P*)."
                ),
                "endpoint_map": "A_T(P)=a0(P)+Delta_lep(P)+Delta_had(P)+Delta_EW(P)",
                "fixed_point_map": "G(P)=phi+sqrt(pi)/A_T(P)",
                "codata_consequence": (
                    "If the source-only endpoint map evaluates to 137.035999177... at its fixed "
                    "point, the outer equation gives P=1.63096820940395932487927984778..."
                ),
            },
            "implemented_fixed_point_gap": implemented_root,
            "codata_mapped_endpoint_packet": compare_point,
            "first_non_internalized_object": {
                "id": "ward_projected_qcd_screening_and_endpoint_remainder",
                "missing_source_transport_delta_alpha_inv_at_codata_p": exact_packet[
                    "missing_source_transport_delta_alpha_inv"
                ],
                "required_screening_factor_at_codata_p": exact_packet["screening_scalar"][
                    "required_screening_factor"
                ],
                "residual_second_order_screening_coefficient_at_codata_p": exact_packet[
                    "screening_scalar"
                ]["residual_second_order_coefficient"],
                "blocker_statement": (
                    "The charged-lepton one-loop kernel is populated. The remaining scalar is in "
                    "the Ward-projected hadronic spectral transport plus electroweak scheme "
                    "remainder. A theorem must emit that map from OPH data; a fitted residual or "
                    "external endpoint comparison does not supply the source-only theorem."
                ),
            },
            "scientific_result": {
                "theorem_grade_object_defined": True,
                "first_non_internalized_scalar_isolated": True,
                "charged_spectrum_is_limiting_missing_object": False,
                "screening_or_endpoint_remainder_is_limiting_missing_object": True,
                "imported_thomson_endpoint_removed_from_solver_path": True,
                "measured_alpha_derivation_supplied": False,
                "source_residual_map_supplied": False,
                "result_scope": "conditional_missing_object_isolation",
            },
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the conditional Thomson endpoint package.")
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--interval-certificate", default=str(DEFAULT_INTERVAL_CERTIFICATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--compare-alpha-inv", default=str(CODATA_2022_ALPHA_INV))
    parser.add_argument("--compare-alpha-inv-uncertainty", default=str(CODATA_2022_ALPHA_INV_UNCERTAINTY))
    parser.add_argument("--precision", type=int, default=DEFAULT_ENDPOINT_PRECISION)
    parser.add_argument("--su2-cutoff", type=int, default=80)
    parser.add_argument("--su3-cutoff", type=int, default=60)
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    interval_certificate = json.loads(
        Path(args.interval_certificate).read_text(encoding="utf-8")
    )
    payload = build_endpoint_package(
        report,
        source_point=certified_source_point(interval_certificate),
        compare_alpha_inv=_dec(args.compare_alpha_inv),
        compare_alpha_inv_uncertainty=_dec(args.compare_alpha_inv_uncertainty),
        precision=args.precision,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
    )
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
