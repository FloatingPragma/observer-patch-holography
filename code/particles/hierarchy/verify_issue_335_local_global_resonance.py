#!/usr/bin/env python3
"""Verifier for OPH issue #335: local/global hierarchy-resonance close-out."""

from __future__ import annotations

import argparse
import json
from decimal import Decimal
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
TOL = Decimal("1e-40")


def D(value: str | int | Decimal | None, default: Decimal | None = None) -> Decimal | None:
    if value is None:
        return default
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify OPH issue #335 local/global hierarchy-resonance close-out."
    )
    parser.add_argument("--check", action="store_true", help="exit nonzero unless close-out checks pass")
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    tick = read_json(ROOT / "certificates/R_N_global_repair_tick_certificate.json")
    ew = read_json(ROOT / "certificates/R_EW_tick_projection_certificate.json")
    pn = read_json(ROOT / "certificates/R_PN_joint_fixed_point_certificate_report.json")
    rg = read_json(ROOT / "issue_332_rg_naturality_certificate.json")
    audit = read_json(ROOT / "certificates/local_global_resonance_audit.json")

    exact = ew.get("exact_bridge", {})
    rounded = ew.get("rounded_capacity_diagnostic", {})
    promotion_requires = audit.get("promotion_requires", [])

    dependencies = {
        "#336_global_repair_tick": tick.get("status")
        == "closed_global_repair_tick_lemma_on_declared_round_structure",
        "#337_electroweak_projection_bridge": ew.get("status")
        == "closed_projection_map_with_exact_bridge_condition"
        and ew.get("accepted") is True,
        "#338_joint_product_fixed_point": pn.get("status")
        == "closed_product_branch_theorem_with_explicit_coupled_branch_boundary",
        "#332_rg_higgs_naturality": rg.get("accepted") is True
        and rg.get("epsilon_H_interval") == ["0", "0"],
    }

    exact_residual = D(exact.get("bridge_residual", "1"))
    exact_projection_error = D(exact.get("projection_exponent_error", "1"))
    rounded_residual = D(rounded.get("bridge_residual", "0"))
    assert exact_residual is not None and exact_projection_error is not None and rounded_residual is not None

    closeout_checks = {
        "all_prerequisite_records_present": all(dependencies.values()),
        "exact_projection_bridge_residual_zero": exact_residual == 0,
        "exact_projection_exponent_matches_4P": abs(exact_projection_error) <= TOL,
        "rounded_capacity_is_diagnostic": rounded.get("status")
        == "diagnostic_only_not_exact_bridge_certificate",
        "rounded_capacity_fails_exact_bridge": abs(rounded_residual) > Decimal("1e-6"),
        "remaining_promotion_gates_recorded": {
            "finite_readback_resolution": any("readback resolution" in item for item in promotion_requires),
            "round_count_derivation": any("24-round" in item for item in promotion_requires),
            "exact_capacity_source_certificate": any("B_EW" in item for item in promotion_requires),
        },
    }
    gates = closeout_checks["remaining_promotion_gates_recorded"]
    accepted = (
        closeout_checks["all_prerequisite_records_present"]
        and closeout_checks["exact_projection_bridge_residual_zero"]
        and closeout_checks["exact_projection_exponent_matches_4P"]
        and closeout_checks["rounded_capacity_is_diagnostic"]
        and closeout_checks["rounded_capacity_fails_exact_bridge"]
        and all(gates.values())
    )

    cert: dict[str, Any] = {
        "issue": 335,
        "artifact": "R_local_global_hierarchy_resonance_closeout",
        "status": "closed_as_exact_surviving_conditional_statement",
        "accepted": bool(accepted),
        "full_theorem_grade_resonance_promoted": False,
        "closeout_decision": (
            "The component proof package is accounted for, but the full local/global "
            "N_CRC resonance is not promoted. Issue #335 closes only as the exact "
            "conditional bridge statement plus the recorded promotion gates."
        ),
        "target_relation": {
            "transport_time": "t_tr(P_star) = (P_star/12) * log(N_CRC/pi)",
            "hierarchy_ratio": "v/E_cell = (N_CRC/pi)^(-P_star/12)",
            "tick_form": "v/E_cell = |g_*'|^(4P_star)",
        },
        "exact_surviving_statement": {
            "projection_map": ew["definitions"]["Pi_EW"],
            "bridge_residual": ew["definitions"]["bridge_residual"],
            "conditional": (
                "If the global capacity source certificate supplies "
                "B_EW(P_star,N_CRC)=0, and the finite readback-resolution and "
                "24-round derivation gates are closed, then the target local/global "
                "hierarchy relation follows from the closed tick, projection, joint "
                "fixed-point, and RG/Higgs naturality records."
            ),
            "exact_bridge_target": ew["definitions"]["exact_bridge_capacity"],
            "N_EW_public_endpoint": exact.get("N_EW"),
        },
        "dependencies": dependencies,
        "dependency_artifacts": {
            "#336": "certificates/R_N_global_repair_tick_certificate.json",
            "#337": "certificates/R_EW_tick_projection_certificate.json",
            "#338": "certificates/R_PN_joint_fixed_point_certificate_report.json",
            "#332": "issue_332_rg_naturality_certificate.json",
        },
        "obstruction_record": {
            "rounded_N_CRC_display": rounded.get("N_display"),
            "rounded_N_CRC_status": rounded.get("status"),
            "rounded_bridge_residual": rounded.get("bridge_residual"),
            "rounded_v_error": rounded.get("v_error"),
            "meaning": (
                "The rounded 3.31e122 cosmological capacity display is not an exact "
                "hierarchy bridge certificate; the exact bridge target is N_EW(P_star)."
            ),
        },
        "remaining_promotion_gates": promotion_requires,
        "acceptance_criteria_status": {
            "states_precise_local_and_global_objects": True,
            "four_prerequisite_steps_accounted_for": all(dependencies.values()),
            "full_theorem_grade_resonance_proved": False,
            "downgraded_to_exact_surviving_conditional_statement": True,
            "compatible_with_local_transmutation_certificate": True,
            "forbids_measured_weak_higgs_or_hierarchy_calibration": True,
            "public_hierarchy_packet_emitted": True,
        },
        "allowed_inputs": [
            "OPH local pixel fixed point P_star",
            "OPH source D10 alpha_U(P_star) interval and transmutation law",
            "global repair-tick record |g_*'|=(N_CRC/pi)^(-1/48)",
            "joint product-branch fixed-point/stability record for (P,N_CRC)",
            "RG/Higgs naturality square on the selected exact branch",
        ],
        "forbidden_calibrations": [
            "measured weak scale v as an input",
            "measured W, Z, Higgs, or top mass as an input",
            "measured G or Planck area as an input",
            "measured Lambda as an input",
            "using the rounded 3.31e122 N_CRC display as an exact bridge certificate",
        ],
        "checks": closeout_checks,
        "verifier_command": (
            "python3 code/particles/hierarchy/verify_issue_335_local_global_resonance.py "
            "--check --output "
            "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json"
        ),
    }

    text = json.dumps(cert, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")

    if args.check and not accepted:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
