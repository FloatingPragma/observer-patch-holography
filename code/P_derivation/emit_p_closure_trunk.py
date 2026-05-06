#!/usr/bin/env python3
"""Emit the compressed OPH P-closure trunk artifact.

This is the code-side version of the five-equation simplification:

    P -> M_U -> alpha_U -> running couplings -> a0(P) -> alpha_in(P) -> P

The artifact is intentionally claim-safe.  The structured Thomson layer is a
continuation model until the Ward-projected zero-momentum transport theorem,
the hadronic spectral object, and the interval certificate are closed.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

from paper_math import build_report


DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "p_closure_trunk_current.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _structured_status(mode: str, structured_running: dict[str, Any] | None) -> str:
    if mode == "mz_anchor":
        return "source_anchor_debug_path"
    if structured_running is None:
        return "no_low_energy_transport_layer"
    if structured_running.get("transport_kernel") == "exact_1loop":
        return "structured_thomson_continuation_not_endpoint_theorem"
    return "paper_compression_asymptotic_structured_running_audit"


def build_p_closure_trunk(report: dict[str, Any]) -> dict[str, Any]:
    """Build a compact, machine-readable trunk from a paper-math report."""
    d10 = dict(report["d10"])
    structured_running = report.get("structured_running")
    mode = str(report["mode"])
    trunk_status = _structured_status(mode, structured_running)

    layers: list[dict[str, Any]] = [
        {
            "id": "layer_1_unification_scale",
            "formula": "M_U(P) = E_P * exp(-2*pi) * P^(1/6)",
            "operation_type": "exp",
            "output_keys": ["M_U"],
            "outputs": {"M_U": d10["mu_u"]},
            "claim_status": "implemented_paper_equation",
        },
        {
            "id": "layer_2_pixel_closure",
            "formula": "ellbar_SU2(4*pi^2*alpha2) + ellbar_SU3(4*pi^2*alpha3) = P/4",
            "operation_type": "heat_kernel_exp_sums",
            "output_keys": ["alpha_U"],
            "outputs": {"alpha_U": d10["alpha_u"]},
            "claim_status": "implemented_paper_equation",
        },
        {
            "id": "layer_3_running_couplings",
            "formula": "alpha_i^-1(m_Z) = alpha_U^-1 + b_i/(2*pi) * log(M_U/m_Z)",
            "operation_type": "log_linear",
            "output_keys": ["alpha1_mz", "alpha2_mz", "alpha3_mz"],
            "outputs": {
                "alpha1_mz": d10["alpha1_mz"],
                "alpha2_mz": d10["alpha2_mz"],
                "alpha3_mz": d10["alpha3_mz"],
            },
            "claim_status": "implemented_with_declared_running_matching_conventions",
        },
        {
            "id": "layer_4_electroweak_mixing",
            "formula": "alpha_em^-1(m_Z^2;P) = alpha2^-1 + (5/3)*alpha1^-1",
            "operation_type": "linear",
            "output_keys": ["a0_alpha_em_inv_mz", "sin2w_mz"],
            "outputs": {
                "a0_alpha_em_inv_mz": report["source_anchor_alpha_inv"],
                "sin2w_mz": d10["sin2w_mz"],
            },
            "claim_status": "source_locked_d10_anchor",
        },
        {
            "id": "layer_5_structured_thomson_running",
            "formula": "alpha_in^-1(P) = a0(P) + Delta_Th(P)",
            "operation_type": "structured_log_linear_continuation",
            "output_keys": ["alpha_inv", "Delta_Th_impl"],
            "outputs": {
                "alpha_inv": report["alpha_inv"],
                "Delta_Th_impl": (
                    structured_running.get("total_delta_alpha_inv")
                    if isinstance(structured_running, dict)
                    else None
                ),
            },
            "claim_status": trunk_status,
        },
    ]

    return {
        "artifact": "oph_p_closure_trunk_current",
        "generated_utc": _now_utc(),
        "source_report_mode": mode,
        "source_report_precision": report["precision"],
        "claim_status": "compressed_candidate_trunk_not_final_particle_root",
        "claim_boundary": (
            "This artifact is the canonical compressed P-trunk candidate for the code path. "
            "It is not the certified particle-pipeline root: promotion requires the Ward-projected "
            "Thomson endpoint theorem, OPH-internal RG/matching/scheme closure, and an interval-level "
            "fixed-point certificate."
        ),
        "closed_form_candidate": {
            "name": "golden_ratio_electromagnetic_width",
            "formula": "P = phi + alpha_in(P) * sqrt(pi)",
            "phi": report["phi"],
            "sqrt_pi": report["sqrt_pi"],
        },
        "fixed_point_candidate": {
            "P": report["p"],
            "alpha": report["alpha"],
            "alpha_inv": report["alpha_inv"],
            "outer_equation_residual": report["god_equation_residual"],
            "alpha_fixed_point_residual": report["alpha_fixed_point_residual"],
            "source_anchor_alpha_inv_mz": report["source_anchor_alpha_inv"],
        },
        "d10_source_point": {
            "M_U": d10["mu_u"],
            "mZ_run": d10["mz_run"],
            "v": d10["v"],
            "alpha_U": d10["alpha_u"],
            "alpha1_mz": d10["alpha1_mz"],
            "alpha2_mz": d10["alpha2_mz"],
            "alpha3_mz": d10["alpha3_mz"],
            "alphaY_mz": d10["alpha_y_mz"],
            "alpha_em_inv_mz": d10["alpha_em_inv_mz"],
            "sin2w_mz": d10["sin2w_mz"],
        },
        "five_layer_chain": layers,
        "structured_running": structured_running,
        "promotion_gates": [
            {
                "id": "ward_projected_thomson_endpoint",
                "github_issue": 235,
                "required_status": "source_emitted_zero_momentum_transport_theorem",
            },
            {
                "id": "rg_matching_threshold_scheme",
                "github_issue": 32,
                "required_status": "OPH_internal_running_matching_and_scheme_theorem",
            },
            {
                "id": "interval_fixed_point_certificate",
                "github_issue": 235,
                "required_status": "interval_existence_uniqueness_certificate_for_the_full_map",
            },
            {
                "id": "particle_codepath_adoption",
                "github_issue": 224,
                "required_status": "all_particle_consumers_read_the_certified_trunk_artifact",
            },
        ],
        "consumer_policy": {
            "may_feed_live_particle_predictions": False,
            "may_feed_compare_or_audit_surfaces": True,
            "hidden_external_alpha_allowed": False,
            "default_thomson_endpoint_allowed": False,
            "torus_mode_interpretation_status": "speculative_not_pipeline_input",
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the compressed OPH P-closure trunk artifact.")
    parser.add_argument(
        "--mode",
        choices=("thomson_structured_running", "thomson_structured_running_asymptotic", "mz_anchor"),
        default="thomson_structured_running_asymptotic",
        help="Alpha readout used by the closure report.",
    )
    parser.add_argument("--precision", type=int, default=18, help="Decimal precision for this artifact build.")
    parser.add_argument("--su2-cutoff", type=int, default=24, help="Representation cutoff for the SU(2) edge sum.")
    parser.add_argument("--su3-cutoff", type=int, default=16, help="Representation cutoff for the SU(3) edge sum.")
    parser.add_argument("--scan-points", type=int, default=24, help="Alpha-space scan points used to bracket closure.")
    parser.add_argument("--max-iterations", type=int, default=12, help="Maximum outer fixed-point bisection iterations.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON artifact path.")
    parser.add_argument("--print-json", action="store_true", help="Print the artifact JSON after writing it.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report(
        precision=args.precision,
        mode=args.mode,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        scan_points=args.scan_points,
        max_iterations=args.max_iterations,
    )
    artifact = build_p_closure_trunk(report)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
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
