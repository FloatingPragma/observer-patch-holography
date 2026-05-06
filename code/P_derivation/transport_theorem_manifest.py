#!/usr/bin/env python3
"""Emit the theorem-status manifest for the missing Thomson transport layer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from alpha_gap_audit import DEFAULT_REPORT, build_alpha_gap_audit


ENDPOINT_PACKAGE = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_package_current.json"


THEOREM_STATUSES = [
    {
        "id": "ward_projected_source_lock",
        "title": "Ward-projected source lock",
        "status": "closed_criterion_not_populated_transport",
        "promotable_to_measured_alpha": False,
        "depends_on": [],
        "blocking_gap": "Computes the correct lane but does not compute t_Q(0;P).",
    },
    {
        "id": "leptonic_one_loop_source_transport",
        "title": "Leptonic one-loop source transport",
        "status": "implemented_numerical_conditional",
        "promotable_to_measured_alpha": False,
        "depends_on": ["stage5_source_only_lepton_mass_theorem", "quadrature_error_bound"],
        "blocking_gap": "Stage-5 lepton masses and quadrature are not yet theorem-certified on this branch.",
    },
    {
        "id": "hadronic_spectral_transport",
        "title": "Hadronic spectral transport",
        "status": "residual_packet_computed_source_map_open",
        "promotable_to_measured_alpha": False,
        "depends_on": ["source_emitted_rho_had_spectral_density", "threshold_support", "ope_tail_match"],
        "blocking_gap": (
            "The endpoint package isolates the residual screening/remainder scalar, but no OPH-emitted "
            "hadronic electromagnetic spectral map rho_had(s;P) or equivalent source remainder map exists."
        ),
    },
    {
        "id": "electroweak_matching_remainder",
        "title": "Electroweak matching remainder",
        "status": "open",
        "promotable_to_measured_alpha": False,
        "depends_on": ["scheme_identity_or_source_formula", "matching_error_bound"],
        "blocking_gap": "The relation between the D10 anchor scheme and Thomson scheme is not bounded.",
    },
    {
        "id": "full_thomson_endpoint",
        "title": "Full Thomson endpoint",
        "status": "endpoint_package_closed_source_residual_open",
        "promotable_to_measured_alpha": False,
        "depends_on": [
            "leptonic_one_loop_source_transport",
            "hadronic_spectral_transport",
            "electroweak_matching_remainder",
        ],
        "blocking_gap": (
            "The missing endpoint packet is numerically isolated. Delta_Th(P) is not source-emitted "
            "until the residual screening/remainder map is emitted on an interval."
        ),
    },
    {
        "id": "fixed_point_closure_with_certified_transport",
        "title": "Fixed-point closure with certified transport",
        "status": "open",
        "promotable_to_measured_alpha": False,
        "depends_on": ["full_thomson_endpoint", "interval_contraction_or_uniqueness_certificate"],
        "blocking_gap": "The current certificate samples the implemented map only, not a closed Delta_Th(P) map.",
    },
]


def build_manifest(report: dict[str, Any]) -> dict[str, Any]:
    audit = build_alpha_gap_audit(report)
    return {
        "artifact": "oph_p_alpha_thomson_transport_theorem_manifest",
        "claim_status": "missing_theorems_scaffold_not_measured_alpha_derivation",
        "source_report_mode": audit["source_report_mode"],
        "source_report_precision": audit["source_report_precision"],
        "implemented_alpha_inv": audit["implemented_alpha_inv"],
        "compare_alpha_inv": audit["compare_alpha_inv"],
        "implemented_transport_delta_alpha_inv": audit["implemented_transport_delta_alpha_inv"],
        "required_transport_delta_alpha_inv": audit["required_transport_delta_alpha_inv"],
        "missing_transport_delta_alpha_inv": audit["missing_transport_delta_alpha_inv"],
        "endpoint_package": {
            "path": str(ENDPOINT_PACKAGE.relative_to(Path(__file__).resolve().parents[1])),
            "exists": ENDPOINT_PACKAGE.exists(),
        },
        "theorems": THEOREM_STATUSES,
        "promotion_rule": {
            "codata_may_enter_solver": False,
            "requires_source_emitted_hadronic_spectral_density": True,
            "requires_single_scheme_lock": True,
            "requires_transport_error_bound": True,
            "requires_interval_fixed_point_certificate": True,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the P/alpha transport theorem manifest.")
    parser.add_argument("--report", default=str(DEFAULT_REPORT), help="Path to a JSON report from derive_p.py.")
    parser.add_argument("--output", help="Optional path for the manifest JSON.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    manifest = build_manifest(report)
    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
