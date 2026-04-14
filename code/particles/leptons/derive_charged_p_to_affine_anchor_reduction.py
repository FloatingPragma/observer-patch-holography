#!/usr/bin/env python3
"""Emit the exact reduction from a D10 landing theorem to A_ch(P).

This does not close the bridge. It fixes the smallest honest theorem target for
the charged P-driven absolute lane: a theorem-grade landing on the charged
determinant line or theorem-grade physical Y_e(P), from which A_ch(P) is
algebraic.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    ANCHOR_SECTION_JSON,
    DETERMINANT_LINE_JSON,
    P_TO_AFFINE_REDUCTION_JSON,
    artifact_ref,
    load_json,
)


ROOT = Path(__file__).resolve().parents[2]
D10_FORWARD_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json"
D10_REPAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
DEFAULT_OUT = P_TO_AFFINE_REDUCTION_JSON


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d10_forward: dict, d10_repair: dict, determinant: dict, anchor: dict) -> dict:
    return {
        "artifact": "oph_charged_p_to_affine_anchor_reduction",
        "generated_utc": _timestamp(),
        "proof_status": "closed_bridge_reduction_to_determinant_line_or_physical_Y_e",
        "public_promotion_allowed": False,
        "theorem_statement": (
            "Let B_D10(P) denote the theorem-grade D10 calibration descendant package emitted from P. "
            "If B_D10(P) admits a theorem-grade landing either on theorem-grade physical Y_e(P) or on a "
            "charged determinant-line section s_det(P), then the charged affine anchor is forced algebraically by "
            "A_ch(P) = (1/3) * log(det(Y_e(P))) = (1/3) * s_det(P). Consequently g_e(P), Delta_e_abs(P), and the "
            "charged mass triple are algebraic descendants. Therefore the exact open bridge target above the D10 "
            "family is a theorem-grade landing on physical Y_e(P) or the charged determinant line, not an additional "
            "independent scalar theorem for A_ch(P)."
        ),
        "closed_d10_inputs": {
            "forward_transmutation_certificate_artifact": d10_forward.get("artifact"),
            "forward_transmutation_certificate_artifact_ref": artifact_ref(D10_FORWARD_JSON),
            "target_free_repair_artifact": d10_repair.get("artifact"),
            "target_free_repair_artifact_ref": artifact_ref(D10_REPAIR_JSON),
            "emitted_from_P": [
                "alpha_U(P)",
                "t_unified(P)",
                "t_transmutation(P)",
                "v_report(P)",
                "alpha2_prime(P)",
                "alphaY_prime(P)",
                "MW_pole(P)",
                "MZ_pole(P)",
                "alpha_em_eff_inv(P)",
                "sin2w_eff(P)",
            ],
        },
        "reduction_theorem": {
            "id": "charged_P_to_A_ch_reduces_to_determinant_line_bridge",
            "statement": (
                "A theorem-grade bridge B_D10(P) -> theorem-grade physical Y_e(P) induces det(Y_e(P)); "
                "a theorem-grade bridge B_D10(P) -> charged determinant-line section s_det(P) induces "
                "A_ch(P) directly. In either case no extra independent affine theorem remains: "
                "A_ch(P) = (1/3) * log(det(Y_e(P))) = (1/3) * s_det(P)."
            ),
        },
        "exact_smallest_bridge_target": {
            "id": "d10_to_charged_determinant_line_bridge",
            "equivalent_landings": [
                "theorem_grade_physical_Y_e(P)",
                "charged_determinant_line_section(P)",
            ],
            "why_smaller_than_A_ch": [
                "The determinant-line section already carries the correct common-shift covariance.",
                "The affine anchor is a readout of that line rather than an extra free scalar once the line is emitted.",
                "The charged determinant-line extension theorem records A_ch = (1/3) * s_det = (1/3) * log(det(Y_e)).",
            ],
        },
        "P_threaded_source_exactness_reformulation": {
            "id": "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu",
            "statement": (
                "Assume the D10 calibration descendants of P land theorem-grade on physical charged data "
                "Y_e(P). Assume theorem-grade C_hat_e is promoted on that same physical Y_e(P) surface. "
                "Then the canonical determinant-line lift theorem fixes "
                "mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P))), defines "
                "C_tilde_e(Y_e(P)) = C_hat_e(Y_e(P)) + mu_phys(Y_e(P)) I, and the "
                "determinant-line section and affine anchor follow canonically by "
                "s_det(Y_e(P)) = 3 mu_phys(Y_e(P)) and A_ch(P) = mu_phys(Y_e(P))."
            ),
            "role": (
                "This reformulates the charged P-bridge burden as theorem-grade source landing on physical Y_e(P), "
                "after which the affine scalar and determinant-line data follow canonically rather than through an extra independent lift theorem."
            ),
        },
        "induced_objects_once_bridge_closes": {
            "determinant_line_section": {
                "artifact": determinant.get("artifact"),
                "artifact_ref": artifact_ref(DETERMINANT_LINE_JSON),
                "formula": "s_det(P) = log(det(Y_e(P)))",
            },
            "affine_anchor": {
                "artifact": anchor.get("artifact"),
                "artifact_ref": artifact_ref(ANCHOR_SECTION_JSON),
                "formula": "A_ch(P) = (1/3) * s_det(P) = (1/3) * log(det(Y_e(P)))",
            },
            "absolute_outputs": {
                "g_e(P)": "exp(A_ch(P))",
                "Delta_e_abs(P)": "log(g_ch_shared) - A_ch(P)",
                "m_e(P)": "exp(A_ch(P) + e_log_centered(P))",
                "m_mu(P)": "exp(A_ch(P) + mu_log_centered(P))",
                "m_tau(P)": "exp(A_ch(P) + tau_log_centered(P))",
            },
        },
        "do_not_claim_now": [
            "theorem_grade_physical_Y_e(P)",
            "charged_determinant_line_section(P)",
            "A_ch(P)",
            "g_e(P)",
            "m_e(P)",
            "m_mu(P)",
            "m_tau(P)",
        ],
        "notes": [
            "This is a reduction theorem, not a bridge-closure theorem.",
            "It fixes the smallest honest bridge target for the D10-to-charged absolute lane.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged P-to-affine-anchor reduction theorem.")
    parser.add_argument("--d10-forward", default=str(D10_FORWARD_JSON))
    parser.add_argument("--d10-repair", default=str(D10_REPAIR_JSON))
    parser.add_argument("--determinant-line", default=str(DETERMINANT_LINE_JSON))
    parser.add_argument("--anchor-section", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.d10_forward)),
        load_json(Path(args.d10_repair)),
        load_json(Path(args.determinant_line)),
        load_json(Path(args.anchor_section)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
