#!/usr/bin/env python3
"""Audit the P/alpha closure result against a compare-only alpha target."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
from typing import Any


DEFAULT_REPORT = Path(__file__).resolve().parent / "runtime" / "full_p_alpha_report_current.json"
DEFAULT_INTERVAL_CERTIFICATE = (
    Path(__file__).resolve().parent
    / "runtime"
    / "p_interval_contraction_certificate_2026-07-14.json"
)
CODATA_2022_FIXTURE = Path(__file__).resolve().parent / "codata_2022_alpha_fixture.json"


def _load_codata_2022_fixture() -> tuple[Decimal, Decimal]:
    """Load the pinned compare-only CODATA values from their data artifact."""

    payload = json.loads(CODATA_2022_FIXTURE.read_text(encoding="utf-8"))
    if set(payload) != {
        "artifact",
        "claim_status",
        "inverse_fine_structure_constant",
        "license",
        "source",
    }:
        raise ValueError("CODATA alpha fixture has an unexpected top-level shape")
    if payload["artifact"] != "codata_2022_fine_structure_compare_fixture":
        raise ValueError("CODATA alpha fixture identity drifted")
    if payload["claim_status"] != "compare_only_empirical_input":
        raise ValueError("CODATA alpha fixture must remain compare-only")
    if payload["source"].get("edition") != "CODATA 2022":
        raise ValueError("CODATA alpha fixture edition drifted")
    values = payload["inverse_fine_structure_constant"]
    return Decimal(values["value"]), Decimal(values["standard_uncertainty"])


CODATA_2022_ALPHA_INV, CODATA_2022_ALPHA_INV_UNCERTAINTY = (
    _load_codata_2022_fixture()
)


def _dec(value: Any) -> Decimal:
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def certified_source_point(interval_certificate: dict[str, Any]) -> dict[str, Any]:
    """Extract the certified structured-running point estimate."""
    return dict(
        interval_certificate["modes"]["thomson_structured_running"][
            "fixed_point_point_estimate_display_only"
        ]
    )


def build_alpha_gap_audit(
    report: dict[str, Any],
    *,
    source_point: dict[str, Any] | None = None,
    compare_alpha_inv: Decimal = CODATA_2022_ALPHA_INV,
    compare_alpha_inv_uncertainty: Decimal = CODATA_2022_ALPHA_INV_UNCERTAINTY,
) -> dict[str, Any]:
    """Return the compare-only gap between the implemented closure and a target.

    The target is not a solver input. It is used only after the report exists,
    so the output can identify which term a full transport derivation would need
    to supply.
    """
    with localcontext() as ctx:
        ctx.prec = 80
        point = source_point or report
        alpha_inv = _dec(point["alpha_inv"])
        source_anchor = _dec(report["source_anchor_alpha_inv"])
        p = _dec(point.get("P", point.get("p")))
        phi = _dec(report["phi"])
        sqrt_pi = _dec(report["sqrt_pi"])
        structured_running = report.get("structured_running") or {}
        implemented_delta = _dec(
            structured_running.get("total_delta_alpha_inv", alpha_inv - source_anchor)
        )
        required_delta = compare_alpha_inv - source_anchor
        missing_delta = compare_alpha_inv - alpha_inv
        compare_alpha = Decimal(1) / compare_alpha_inv
        compare_p = phi + compare_alpha * sqrt_pi

        return {
            "claim_status": "open_transport_gap_not_full_derivation",
            "claim_boundary": (
                "The P closure and current structured-running continuation emit a reproducible "
                "numerical witness. A measured-alpha derivation is not closed until the missing "
                "same-family Thomson transport contribution is derived without using the compare target."
            ),
            "source_report_mode": report.get("mode"),
            "source_report_precision": report.get("precision"),
            "source_point_certificate": (
                "p_interval_contraction_certificate_2026-07-14.json::"
                "modes.thomson_structured_running.fixed_point_point_estimate_display_only"
                if source_point is not None
                else None
            ),
            "source_anchor_alpha_inv_mz": str(+source_anchor),
            "implemented_alpha_inv": str(+alpha_inv),
            "implemented_alpha": str(+(Decimal(1) / alpha_inv)),
            "implemented_p": str(+p),
            "implemented_transport_delta_alpha_inv": str(+implemented_delta),
            "implemented_transport_fraction": str(+(alpha_inv / source_anchor - Decimal(1))),
            "compare_alpha_inv": str(+compare_alpha_inv),
            "compare_alpha_inv_standard_uncertainty": str(+compare_alpha_inv_uncertainty),
            "compare_alpha": str(+compare_alpha),
            "compare_p_from_outer_equation": str(+compare_p),
            "required_transport_delta_alpha_inv": str(+required_delta),
            "missing_transport_delta_alpha_inv": str(+missing_delta),
            "missing_transport_fraction": str(+(missing_delta / source_anchor)),
            "missing_fraction_of_implemented_transport": str(+(missing_delta / implemented_delta)),
            "alpha_inv_gap_ppm": str(+(missing_delta / compare_alpha_inv * Decimal(1_000_000))),
            "alpha_inv_gap_sigma": str(+(missing_delta / compare_alpha_inv_uncertainty)),
            "p_gap_implemented_minus_compare": str(+(p - compare_p)),
            "diagnosis": (
                "The gap is in the low-energy Thomson transport/readout layer, not in the outer "
                "fixed-point algebra. The implemented structured-running term is too small by "
                f"{+missing_delta} in inverse-alpha units."
            ),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit a P/alpha closure report against a compare-only inverse-alpha target."
    )
    parser.add_argument("--report", default=str(DEFAULT_REPORT), help="Path to a JSON report from derive_p.py.")
    parser.add_argument(
        "--compare-alpha-inv",
        default=str(CODATA_2022_ALPHA_INV),
        help="Compare-only inverse-alpha target. Default is CODATA 2022.",
    )
    parser.add_argument(
        "--compare-alpha-inv-uncertainty",
        default=str(CODATA_2022_ALPHA_INV_UNCERTAINTY),
        help="One-sigma uncertainty for the compare-only inverse-alpha target.",
    )
    parser.add_argument("--json", action="store_true", help="Print the full audit as JSON.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    audit = build_alpha_gap_audit(
        report,
        compare_alpha_inv=_dec(args.compare_alpha_inv),
        compare_alpha_inv_uncertainty=_dec(args.compare_alpha_inv_uncertainty),
    )
    if args.json:
        print(json.dumps(audit, indent=2, sort_keys=True))
        return 0

    print(f"claim status                         = {audit['claim_status']}")
    print(f"implemented alpha^-1                 = {audit['implemented_alpha_inv']}")
    print(f"compare alpha^-1                     = {audit['compare_alpha_inv']}")
    print(f"source anchor alpha^-1(mZ^2)         = {audit['source_anchor_alpha_inv_mz']}")
    print(f"implemented Delta alpha^-1 transport = {audit['implemented_transport_delta_alpha_inv']}")
    print(f"required Delta alpha^-1 transport    = {audit['required_transport_delta_alpha_inv']}")
    print(f"missing Delta alpha^-1 transport     = {audit['missing_transport_delta_alpha_inv']}")
    print(f"alpha^-1 gap                         = {audit['alpha_inv_gap_ppm']} ppm")
    print(f"alpha^-1 gap                         = {audit['alpha_inv_gap_sigma']} sigma")
    print(f"implemented P                        = {audit['implemented_p']}")
    print(f"compare P from outer equation        = {audit['compare_p_from_outer_equation']}")
    print(f"P gap                                = {audit['p_gap_implemented_minus_compare']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
