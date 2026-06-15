#!/usr/bin/env python3
"""Validate the OPH issue #337 electroweak tick-projection certificate."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str) -> Decimal:
    return Decimal(value)


def main(path: str = "certificates/R_EW_tick_projection_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    exact = cert.get("exact_bridge", {})
    rounded = cert.get("rounded_capacity_diagnostic", {})
    forbidden = cert.get("forbidden_calibrations", [])
    boundary = cert.get("claim_boundary", {})
    derivation_chain = cert.get("derivation_chain", [])
    factor_origins = cert.get("factor_origins", {})
    acceptance = cert.get("acceptance_criteria_status", {})

    chain_steps = {step.get("step"): step for step in derivation_chain if isinstance(step, dict)}
    closed_elsewhere = boundary.get("closed_elsewhere", [])

    checks = {
        "issue_is_337": cert.get("issue") == 337,
        "accepted": cert.get("accepted") is True,
        "status_closed_projection": cert.get("status") == "closed_projection_map_with_exact_bridge_condition",
        "bridge_residual_zero": D(exact.get("bridge_residual", "1")) == 0,
        "projection_exponent_matches_4P": abs(D(exact.get("projection_exponent_error", "1"))) <= Decimal("1e-40"),
        "v_bridge_matches_source": abs(D(exact.get("v_bridge_error", "1"))) <= Decimal("1e-40"),
        "rounded_N_is_diagnostic": rounded.get("status") == "diagnostic_only_not_exact_bridge_certificate",
        "rounded_N_fails_exact_bridge": abs(D(rounded.get("bridge_residual", "0"))) > Decimal("1e-6"),
        "weak_scale_forbidden": any("measured weak scale" in item for item in forbidden),
        "rounded_display_forbidden": any("rounded N_CRC" in item for item in forbidden),
        "exact_source_boundary_recorded": "B_EW" in boundary.get("source_certificate_required", ""),
        "derivation_chain_has_seven_steps": len(chain_steps) == 7,
        "step_3_derives_projection_formula": (
            "Pi_EW" in str(chain_steps.get(3, {}).get("conclusion", ""))
            and "24*pi" in str(chain_steps.get(3, {}).get("conclusion", ""))
        ),
        "step_4_records_resonance_target_scope": "scope_note" in chain_steps.get(4, {}),
        "step_6_cites_capacity_certificate": "R_EW_global_capacity" in str(chain_steps.get(6, {}).get("source", "")),
        "step_7_derives_P_over_12_form": "(P_star/12)" in str(chain_steps.get(7, {}).get("conclusion", "")) and "48/4" in str(chain_steps.get(7, {}).get("conclusion", "")),
        "factor_origin_beta_EW_recorded": factor_origins.get("beta_EW", {}).get("value") == "4"
            and "D10" in str(factor_origins.get("beta_EW", {}).get("source_theorem", "")),
        "factor_origin_m_rep_recorded": factor_origins.get("m_rep", {}).get("value") == "24"
            and "representation-to-spectrum" in str(factor_origins.get("m_rep", {}).get("source_theorem", "")),
        "factor_origin_48_recorded": factor_origins.get("tick_exponent_denominator_48", {}).get("value") == "48"
            and "global repair-tick" in str(factor_origins.get("tick_exponent_denominator_48", {}).get("source_theorem", "")),
        "factor_origin_12_recorded": factor_origins.get("projection_target_denominator_12_in_P_over_12", {}).get("value") == "12"
            and "2 * m_rep / beta_EW" in str(factor_origins.get("projection_target_denominator_12_in_P_over_12", {}).get("definition", "")),
        "acceptance_projection_map_defined": acceptance.get("projection_map_defined") is True,
        "acceptance_4P_proved_under_resonance_target": acceptance.get("sampling_exponent_4P_proved_under_resonance_target") is True,
        "acceptance_factor_4_origin_documented": acceptance.get("factor_4_origin_documented") is True,
        "acceptance_factor_12_origin_documented": acceptance.get("factor_12_origin_documented") is True,
        "acceptance_compatible_with_local_D10": acceptance.get("compatibility_with_local_D10_transmutation_certificate") is True,
        "acceptance_no_measured_weak_inputs": acceptance.get("no_measured_weak_scale_inputs") is True,
        "acceptance_resonance_target_scoped_as_oph_condition": acceptance.get("resonance_target_scoped_as_oph_condition") is True,
        "boundary_records_closed_elsewhere": (
            any("D10" in item for item in closed_elsewhere)
            and any("representation-to-spectrum" in item or "R_m_rep_24" in item for item in closed_elsewhere)
            and any("global repair-tick" in item or "R_N_global_repair_tick" in item for item in closed_elsewhere)
            and any("R_EW_global_capacity" in item or "EW-refined" in item for item in closed_elsewhere)
        ),
        "boundary_includes_scope_note": "resonance target" in str(boundary.get("scope", "")).lower()
            or "ew-refined capacity" in str(boundary.get("scope", "")).lower(),
    }
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "certificates/R_EW_tick_projection_certificate.json"))
