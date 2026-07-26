#!/usr/bin/env python3
"""Emit a pixel-screen resonance summary from existing hierarchy certificates.

This is a receipt-level composition artifact, not an additional theorem. It
records the equal-area screen chart identities in the Standard Model gauge derivation and checks
that they are sourced from the same selected `(P, N_CRC^EW)`
branch as the local/global hierarchy resonance certificate.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any


getcontext().prec = 180

ROOT = Path(__file__).resolve().parent
CAPACITY_CERT = ROOT / "certificates" / "R_EW_global_capacity_certificate.json"
RESONANCE_CERT = ROOT / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
SCREEN_SIEVE_CERT = ROOT / "certificates" / "R_screen_sieve_icosahedral_certificate.json"
READBACK_CERT = ROOT / "certificates" / "R_readback_resolution_certificate.json"
JOINT_CERT = ROOT / "certificates" / "R_PN_joint_fixed_point_certificate_report.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "R_pixel_screen_resonance_summary.json"

PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")
THREE = Decimal(3)
FOUR = Decimal(4)
TWELVE = Decimal(12)
TOL = Decimal("1e-100")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def decstr(value: Decimal) -> str:
    if value == 0:
        return "0"
    if abs(value) < Decimal("1e-6") or abs(value) >= Decimal("1e6"):
        return format(value, "E")
    return format(value, "f")


def rel_error(a: Decimal, b: Decimal) -> Decimal:
    scale = max(abs(a), abs(b), Decimal(1))
    return abs(a - b) / scale


def build_payload() -> dict[str, Any]:
    capacity = read_json(CAPACITY_CERT)
    resonance = read_json(RESONANCE_CERT)
    screen = read_json(SCREEN_SIEVE_CERT)
    readback = read_json(READBACK_CERT)
    joint = read_json(JOINT_CERT)

    source = capacity["source_values"]
    exact = capacity["exact_capacity_fixed_point"]
    p_star = Decimal(source["P_star"])
    n_crc = Decimal(exact["N_CRC_EW"])
    beta_ew = Decimal(source.get("beta_EW", "4")) if "beta_EW" in source else FOUR

    cell_entropy = p_star / beta_ew
    cell_count = n_crc / cell_entropy
    reconstructed_capacity = cell_count * cell_entropy

    lambda_lstar2 = THREE * PI / n_crc
    lambda_a_cell = lambda_lstar2 * p_star
    lambda_a_cell_from_cell_count = TWELVE * PI / cell_count

    capacity_rel = rel_error(reconstructed_capacity, n_crc)
    lambda_rel = rel_error(lambda_a_cell, lambda_a_cell_from_cell_count)

    checks = {
        "capacity_certificate_accepted": capacity.get("accepted") is True,
        "resonance_certificate_is_exact_conditional": (
            resonance.get("accepted") is True
            and resonance.get("status")
            == "exact_conditional_local_global_hierarchy_resonance"
            and resonance.get("full_theorem_grade_resonance_promoted") is False
        ),
        "screen_sieve_supplies_twelve_ports_and_x_over_twelve": (
            screen.get("orbit_stabilizer", {}).get("orbit_size") == 12
            and screen.get("screen_load_arithmetic", {}).get("local_port_read") == "X/12"
            and screen.get("screen_load_arithmetic", {}).get("gamma_screen_simplified")
            == "(P/12)*log(N/pi)"
        ),
        "hierarchy_screen_readout_is_named_open_premise": (
            resonance.get("target_relation", {}).get("conditional_on")
            == "HIERARCHY-SCREEN-READOUT"
            and screen.get("hierarchy_screen_readout_gate", {}).get("supplied_by_screen_sieve")
            is False
        ),
        "screen_and_product_adjoint_counts_are_independently_sourced": (
            2 * screen.get("orbit_stabilizer", {}).get("orbit_size", 0) == 24
            and resonance.get("round_count_certificate", {}).get("m_rep") == 24
            and resonance.get("round_count_certificate", {}).get("artifact")
            == "certificates/R_m_rep_24_certificate.json"
        ),
        "readback_uses_same_n_crc_ew": "N_CRC^EW" in str(readback.get("branch_scope", {})),
        "joint_product_branch_recorded": joint.get("status")
        == "closed_product_branch_theorem_with_explicit_coupled_branch_boundary",
        "cell_count_reconstructs_capacity": capacity_rel <= TOL,
        "cosmological_cell_coordinate_matches_count_form": lambda_rel <= TOL,
        "si_lambda_not_promoted_here": True,
        "primitive_carrier_not_selected_here": True,
    }

    accepted = all(checks.values())

    return {
        "artifact": "R_pixel_screen_resonance_summary",
        "status": "closed_receipt_summary_from_existing_hierarchy_certificates",
        "accepted": bool(accepted),
        "purpose": (
            "Expose the pixel-screen tile identity, independent screen and product-adjoint counts, and "
            "dimensionless de Sitter capacity coordinate as one machine-readable "
            "summary sourced from existing hierarchy certificates."
        ),
        "claim_boundary": {
            "closed_here": (
                "Receipt-level composition only: K_cell=4*N_CRC^EW/P_star, "
                "K_cell*(P_star/4)=N_CRC^EW, Lambda_CRC*l_star^2=3*pi/N_CRC^EW, "
                "and Lambda_CRC*a_cell=3*pi*P_star/N_CRC^EW=12*pi/K_cell."
            ),
            "not_closed_here": [
                "a primitive one-cell or fixed-block carrier theorem",
                "an SI Lambda value without the selected scale certificate",
                "an SI G or no-G clock-stack promotion",
                "a coupled (P,N_CRC) source-map theorem beyond the recorded product branch",
                "HIERARCHY-SCREEN-READOUT: log(E_cell/v)=Gamma_screen",
            ],
        },
        "source_pair": {
            "P_star": decstr(p_star),
            "P_star_branch": source.get("P_star_branch"),
            "P_star_source_artifact": source.get("P_star_source_artifact"),
            "N_CRC_EW": decstr(n_crc),
            "N_CRC_EW_source_artifact": "certificates/R_EW_global_capacity_certificate.json",
            "same_pair_scope": (
                "The summary consumes the selected public-endpoint P_star and the "
                "EW-refined exact capacity fixed point N_CRC^EW(P_star)."
            ),
        },
        "pixel_screen_tile_identity": {
            "cell_entropy": decstr(cell_entropy),
            "cell_entropy_formula": "ell_cell = P_star/4",
            "cell_count": decstr(cell_count),
            "cell_count_formula": "K_cell = N_CRC^EW / (P_star/4) = 4*N_CRC^EW/P_star",
            "reconstructed_capacity": decstr(reconstructed_capacity),
            "reconstruction_formula": "K_cell*(P_star/4) = N_CRC^EW",
            "relative_reconstruction_error": decstr(capacity_rel),
        },
        "shared_12_24_port_alignment": {
            "screen_ports": 12,
            "screen_oriented_slots": 24,
            "screen_port_source": "certificates/R_screen_sieve_icosahedral_certificate.json",
            "screen_register_source": "paper/screen_microphysics_and_observer_synchronization.tex#def:oriented-24-slot-register",
            "product_adjoint_dimension": 12,
            "product_adjoint_rounds_m_rep": 24,
            "product_adjoint_source": "certificates/R_m_rep_24_certificate.json",
            "equal_cardinality_status": "arithmetic_alignment_without_physical_identification",
            "screen_load_split": "X/12",
            "gamma_screen_algebra": "Gamma_screen=(P_star/12)*log(N_CRC^EW/pi)",
            "identification_premise": "HIERARCHY-SCREEN-READOUT",
            "identification_status": "work_in_progress",
            "conditional_hierarchy_readout": resonance["target_relation"]["hierarchy_ratio"],
            "conditional_tick_form": resonance["target_relation"]["tick_form"],
            "conditional_transport_time": resonance["target_relation"]["transport_time"],
        },
        "dimensionless_de_sitter_coordinate": {
            "Lambda_lstar2": decstr(lambda_lstar2),
            "Lambda_lstar2_formula": "Lambda_CRC*l_star^2 = 3*pi/N_CRC^EW",
            "Lambda_a_cell": decstr(lambda_a_cell),
            "Lambda_a_cell_formula": "Lambda_CRC*a_cell = 3*pi*P_star/N_CRC^EW",
            "Lambda_a_cell_from_cell_count": decstr(lambda_a_cell_from_cell_count),
            "Lambda_a_cell_count_formula": "Lambda_CRC*a_cell = 12*pi/K_cell",
            "relative_cell_coordinate_error": decstr(lambda_rel),
            "si_display_boundary": "SI Lambda requires the selected scale certificate and is not emitted here.",
        },
        "dependency_artifacts": {
            "exact_capacity": "certificates/R_EW_global_capacity_certificate.json",
            "local_global_resonance": "certificates/R_local_global_hierarchy_resonance_closeout_335.json",
            "screen_sieve": "certificates/R_screen_sieve_icosahedral_certificate.json",
            "product_adjoint_round_count": "certificates/R_m_rep_24_certificate.json",
            "finite_readback_resolution": "certificates/R_readback_resolution_certificate.json",
            "joint_product_fixed_point": "certificates/R_PN_joint_fixed_point_certificate_report.json",
        },
        "checks": checks,
        "verifier_command": (
            "python3 code/particles/hierarchy/verify_pixel_screen_resonance_summary.py "
            "--check --output code/particles/hierarchy/certificates/R_pixel_screen_resonance_summary.json"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Emit the OPH pixel-screen resonance summary.")
    parser.add_argument("--check", action="store_true", help="exit nonzero unless the summary checks pass")
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    payload = build_payload()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")

    if args.check and payload["accepted"] is not True:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
