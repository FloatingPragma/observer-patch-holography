#!/usr/bin/env python3
"""Validate the OPH issue #344 exact-capacity certificate."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str) -> Decimal:
    return Decimal(value)


def main(path: str = "certificates/R_EW_global_capacity_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    exact = cert.get("exact_capacity_fixed_point", {})
    contraction = cert.get("contraction_certificate", {})
    rounded = cert.get("rounded_capacity_diagnostic", {})
    forbidden = cert.get("forbidden_calibrations", [])
    boundary = cert.get("claim_boundary", {})

    checks = {
        "issue_is_344": cert.get("issue") == 344,
        "accepted": cert.get("accepted") is True,
        "status_closed_exact_capacity": (
            cert.get("status") == "closed_bridge_refined_global_capacity_fixed_point_certificate"
        ),
        "bridge_residual_zero": D(exact.get("bridge_residual", "1")) == 0,
        "fixed_point_residual_zero": D(exact.get("fixed_point_residual_x", "1")) == 0,
        "projection_exponent_matches_4P": (
            abs(D(exact.get("projection_exponent", "0")) - D(exact.get("target_exponent_4P", "1")))
            <= Decimal("1e-40")
        ),
        "v_identity_matches": abs(D(exact.get("v_identity_error", "1"))) <= Decimal("1e-40"),
        "contraction_factor_half": D(contraction.get("lipschitz_constant", "0")) == Decimal("0.5"),
        "residual_contracts": (
            abs(
                D(contraction.get("sample_residual_ratio", "0"))
                - D(contraction.get("residual_contracts_by", "1"))
            )
            <= Decimal("1e-40")
        ),
        "rounded_capacity_is_diagnostic": (
            rounded.get("status") == "diagnostic_only_not_exact_bridge_certificate"
        ),
        "rounded_capacity_fails_bridge": abs(D(rounded.get("bridge_residual", "0"))) > Decimal("1e-6"),
        "weak_scale_forbidden": any("measured weak scale" in item for item in forbidden),
        "lambda_forbidden": any("measured Lambda" in item for item in forbidden),
        "rounded_display_forbidden": any("3.31e122" in item for item in forbidden),
        "remaining_gates_are_two": len(boundary.get("not_closed_here", [])) == 2,
    }
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(sys.argv[1] if len(sys.argv) > 1 else "certificates/R_EW_global_capacity_certificate.json")
    )
