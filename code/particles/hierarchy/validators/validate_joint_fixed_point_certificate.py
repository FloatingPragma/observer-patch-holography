#!/usr/bin/env python3
"""Validate the OPH issue #338 joint fixed-point boundary certificate."""

from __future__ import annotations

import json
import pathlib
import sys


def main(path: str = "certificates/R_PN_joint_fixed_point_certificate_report.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    product = cert.get("product_contraction_certificate", {})
    coupled = cert.get("coupled_contraction_certificate", {})
    boundary = cert.get("claim_boundary", {})
    coupled_required = coupled.get("required", "").replace(" ", "")
    checks = {
        "issue_is_338": cert.get("issue") == 338,
        "product_theorem_status": cert.get("status") == "closed_product_branch_theorem_with_explicit_coupled_branch_boundary",
        "N_display_is_warning_only": "rounded" in cert.get("N_display_warning", ""),
        "backsolve_is_diagnostic_only": "CIRCULAR_DIAGNOSTIC_ONLY" in cert.get("N_backsolved_warning", ""),
        "product_component_condition_recorded": product.get("status") == "conditional_on_component_contractions",
        "coupled_residual_boundary_recorded": coupled.get("status") == "residual_coupled_branch_boundary",
        "weak_inputs_forbidden": "weak scale" in boundary.get("forbidden_input", ""),
        "coupled_metric_condition_recorded": "max(a+b/r,d+r*c)<1" in coupled_required,
    }
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "certificates/R_PN_joint_fixed_point_certificate_report.json"))
