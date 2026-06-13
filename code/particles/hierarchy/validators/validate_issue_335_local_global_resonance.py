#!/usr/bin/env python3
"""Validate the OPH issue #335 local/global hierarchy-resonance close-out."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str) -> Decimal:
    return Decimal(value)


def main(path: str = "certificates/R_local_global_hierarchy_resonance_closeout_335.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    deps = cert.get("dependencies", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    obstruction = cert.get("obstruction_record", {})
    gates = cert.get("remaining_promotion_gates", [])
    forbidden = cert.get("forbidden_calibrations", [])
    checks = cert.get("checks", {})
    gate_checks = checks.get("remaining_promotion_gates_recorded", {})
    exact_capacity = cert.get("exact_capacity_certificate", {})

    validation = {
        "issue_is_335": cert.get("issue") == 335,
        "accepted_conditional_closeout": cert.get("accepted") is True,
        "status_is_conditional_closeout": cert.get("status")
        == "closed_as_exact_surviving_conditional_statement",
        "full_theorem_not_promoted": cert.get("full_theorem_grade_resonance_promoted") is False,
        "all_dependencies_closed": all(deps.values()),
        "exact_capacity_supplied": (
            acceptance.get("exact_capacity_source_certificate_supplied") is True
            and D(exact_capacity.get("bridge_residual", "1")) == 0
        ),
        "finite_readback_and_round_count_recorded": (
            acceptance.get("finite_readback_and_round_count_remain_for_full_resonance") is True
        ),
        "full_theorem_status_false": acceptance.get("full_theorem_grade_resonance_proved") is False,
        "exact_bridge_target_recorded": "3.5323546226929906511187512962330547600462" in (
            cert.get("exact_surviving_statement", {}).get("N_EW_public_endpoint", "")
        ),
        "rounded_capacity_rejected": obstruction.get("rounded_N_CRC_status")
        == "diagnostic_only_not_exact_bridge_certificate"
        and abs(D(obstruction.get("rounded_bridge_residual", "0"))) > Decimal("1e-6"),
        "two_promotion_gates_recorded": len(gates) == 2
        and gate_checks.get("finite_readback_resolution") is True
        and gate_checks.get("round_count_derivation") is True
        and gate_checks.get("exact_capacity_source_certificate") is True,
        "weak_scale_forbidden": any("measured weak scale" in item for item in forbidden),
        "rounded_display_forbidden": any("3.31e122" in item for item in forbidden),
    }
    payload = {"checks": validation, "pass": all(validation.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(
            sys.argv[1]
            if len(sys.argv) > 1
            else "certificates/R_local_global_hierarchy_resonance_closeout_335.json"
        )
    )
