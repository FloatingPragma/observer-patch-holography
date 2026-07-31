#!/usr/bin/env python3
"""Emit the empirical fine-structure endpoint display surface."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
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
from thomson_endpoint_package import build_endpoint_package


DEFAULT_REPORT = Path(__file__).resolve().parent / "runtime" / "full_p_alpha_report_current.json"
DEFAULT_SOURCE_THEOREM = Path(__file__).resolve().parent / "runtime" / "source_spectral_theorem_current.json"
DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "measured_endpoint_calibration_current.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_measured_endpoint_calibration(
    report: dict[str, Any],
    *,
    source_point: dict[str, Any] | None = None,
    source_spectral_theorem: dict[str, Any] | None = None,
    compare_alpha_inv: Decimal = CODATA_2022_ALPHA_INV,
    compare_alpha_inv_uncertainty: Decimal = CODATA_2022_ALPHA_INV_UNCERTAINTY,
    precision: int = 80,
) -> dict[str, Any]:
    """Build a non-promoting empirical endpoint artifact from a measured endpoint.

    This object is explicit about its external input. It may feed the OPH plus
    empirical hadron closure display surface, but it cannot satisfy the
    source-only theorem gate.
    """
    ctx = PaperMathContext(precision=precision, su2_cutoff=80, su3_cutoff=60)
    gap = build_alpha_gap_audit(
        report,
        source_point=source_point,
        compare_alpha_inv=compare_alpha_inv,
        compare_alpha_inv_uncertainty=compare_alpha_inv_uncertainty,
    )
    endpoint = build_endpoint_package(
        report,
        source_point=source_point,
        compare_alpha_inv=compare_alpha_inv,
        compare_alpha_inv_uncertainty=compare_alpha_inv_uncertainty,
        precision=precision,
        su2_cutoff=80,
        su3_cutoff=60,
    )
    exact_packet = endpoint["codata_mapped_endpoint_packet"]["exact_one_loop_package"]
    source_obstruction = (source_spectral_theorem or {}).get("current_corpus_obstruction", {})

    with localcontext() as dec_ctx:
        dec_ctx.prec = precision
        alpha = +(Decimal(1) / compare_alpha_inv)
        p_value = +ctx.p_from_inverse_alpha(compare_alpha_inv)
        p_sigma = +(ctx.sqrt_pi * compare_alpha_inv_uncertainty / (compare_alpha_inv * compare_alpha_inv))

    return to_serializable(
        {
            "artifact": "oph_measured_fine_structure_endpoint_calibration",
            "generated_utc": _now_utc(),
            "status": "oph_plus_empirical_hadron_closure_endpoint",
            "row_class": "oph_plus_empirical_hadron_closure",
            "promotion_allowed": False,
            "exact_alpha_promoted": False,
            "external_input_used": True,
            "external_input_role": (
                "The measured Thomson-limit inverse fine-structure constant supplies the empirical "
                "hadron-completed endpoint for numeric tables and displays. It is excluded from "
                "the source-only solver and from theorem promotion."
            ),
            "empirical_hadron_closure": {
                "measured_thomson_endpoint_used": True,
                "raw_cross_section_data_used_in_this_calibration": False,
                "external_input_scope": "measured_thomson_endpoint_only",
                "source_registry": "code/particles/hadron/empirical_ee_hadrons_sources.yaml",
                "empirical_payload_schema": (
                    "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json"
                ),
                "dispersion_payload_status": (
                    "separate_empirical_spectral_endpoint_not_consumed_by_this_"
                    "measured_endpoint_calibration"
                ),
                "row_class": "oph_plus_empirical_hadron_closure",
                "source_only_theorem_status": "not_promoted",
            },
            "source_only_guard": {
                "codata_enters_solver": False,
                "measured_endpoint_allowed_as_transport_input": False,
                "may_satisfy_source_spectral_payload_gate": False,
                "must_be_labeled_external_calibration": True,
            },
            "calibrated_values": {
                "alpha_inv_0": compare_alpha_inv,
                "alpha_inv_0_standard_uncertainty": compare_alpha_inv_uncertainty,
                "alpha_0": alpha,
                "P_from_outer_equation": p_value,
                "P_standard_uncertainty": p_sigma,
                "definition": "P = phi + sqrt(pi)/alpha_inv_0",
            },
            "current_source_candidate": {
                "alpha_inv": (source_point or report)["alpha_inv"],
                "alpha": (source_point or report)["alpha"],
                "P": (source_point or report).get("P", report["p"]),
                "missing_inverse_alpha_units": gap["missing_transport_delta_alpha_inv"],
                "p_gap_implemented_minus_calibrated": gap["p_gap_implemented_minus_compare"],
            },
            "codata_mapped_endpoint_requirement": {
                "a0_alpha_inv_mz": exact_packet["source_anchor_alpha_inv_mz"],
                "required_transport_delta_alpha_inv": exact_packet["required_transport_delta_alpha_inv"],
                "implemented_transport_delta_alpha_inv": exact_packet[
                    "implemented_transport_delta_alpha_inv"
                ],
                "missing_source_transport_delta_alpha_inv": exact_packet[
                    "missing_source_transport_delta_alpha_inv"
                ],
                "required_screening_factor": exact_packet["screening_scalar"][
                    "required_screening_factor"
                ],
                "residual_second_order_coefficient": exact_packet["screening_scalar"][
                    "residual_second_order_coefficient"
                ],
            },
            "why_source_computation_is_blocked": {
                "summary": (
                    "The missing object is a nonperturbative Ward-projected hadronic spectral "
                    "measure. It is outside the scalar algebra emitted by the current corpus. "
                    "The current corpus has a schema and contract, with no production "
                    "finite-volume vector levels, current residues, normalization, continuum "
                    "pushforward, systematics budgets, or endpoint quadrature certificate."
                ),
                "current_corpus_obstruction": source_obstruction,
                "local_hardware_boundary": (
                    "A source computation requires real finite-volume gauge-field production "
                    "and spectral reconstruction on the theorem-emitted family, followed by "
                    "continuum, chiral, finite-volume, current-matching, tail, and directed "
                    "interval error certificates. That production payload is absent from this "
                    "workspace and is outside the current local symbolic/code path."
                ),
                "nonidentifiability": (
                    "Distinct positive spectral-measure completions can agree on every invariant "
                    "emitted by the current corpus while giving different Thomson moments. The "
                    "endpoint therefore cannot be recovered by higher precision or more scalar "
                    "algebra alone."
                ),
            },
            "consumer_policy": {
                "may_feed_numeric_prediction_tables": True,
                "may_feed_public_plots_and_svg_values": True,
                "may_feed_compare_or_audit_surfaces": True,
                "may_feed_source_theorem_claim": False,
                "must_show_caveat": (
                    "fine-structure endpoint uses OPH plus empirical hadron closure; source-only "
                    "spectral payload is absent"
                ),
            },
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the empirical alpha endpoint artifact.")
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--interval-certificate", default=str(DEFAULT_INTERVAL_CERTIFICATE))
    parser.add_argument("--source-spectral-theorem", default=str(DEFAULT_SOURCE_THEOREM))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--compare-alpha-inv", default=str(CODATA_2022_ALPHA_INV))
    parser.add_argument("--compare-alpha-inv-uncertainty", default=str(CODATA_2022_ALPHA_INV_UNCERTAINTY))
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_theorem_path = Path(args.source_spectral_theorem)
    source_theorem = _load_json(source_theorem_path) if source_theorem_path.exists() else None
    interval_certificate = _load_json(Path(args.interval_certificate))
    payload = build_measured_endpoint_calibration(
        _load_json(Path(args.report)),
        source_point=certified_source_point(interval_certificate),
        source_spectral_theorem=source_theorem,
        compare_alpha_inv=_dec(args.compare_alpha_inv),
        compare_alpha_inv_uncertainty=_dec(args.compare_alpha_inv_uncertainty),
        precision=args.precision,
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
