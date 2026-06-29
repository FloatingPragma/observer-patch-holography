#!/usr/bin/env python3
"""Emit the required charged determinant trace-lift theorem target.

This is not a theorem closure artifact. It is the machine-readable promotion
gate for the charged-lepton absolute-mass lane: the downstream readout from
the affine anchor is algebraic, but the live corpus does not emit the
source-only determinant attachment that would produce that anchor from P.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from charged_absolute_route_common import (
    DETERMINANT_CHARACTER_FRONTIER_JSON,
    DETERMINANT_TRACE_LIFT_ATTACHMENT_REQUIRED_JSON,
    DETERMINANT_TRACE_NORMALIZATION_NO_GO_JSON,
    artifact_ref,
    load_json,
)


ROOT = Path(__file__).resolve().parents[2]
CURRENT_FAMILY_EXACT_READOUT_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"
)
CURRENT_FAMILY_AFFINE_ANCHOR_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"
)

FORBIDDEN_ANCESTORS = [
    "particle_reference_values",
    "lepton_current_family_exact_readout",
    "charged_lepton_reference_targets",
    "target_centered_log_shape_exact",
    "m_e",
    "m_mu",
    "m_tau",
    "log_m_e_m_mu_m_tau",
    "g_e_star",
    "Delta_e_abs_star",
    "PDG",
    "CODATA",
    "compare_only",
]

MISSING_FOR_PROMOTION = [
    "charged_branch_generator_splitting_promotion",
    "sector_isolated_charged_determinant_exponent_vector_M_ch",
    "source_side_same_label_q_psi_readout_certificate",
    "charged_determinant_trace_lift_attachment",
    "NO_TARGET_LEAK_DAG_CHARGED_A_CH",
]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _formula(symbol: str, centered_log: float) -> str:
    sign = "+" if centered_log >= 0.0 else "-"
    return f"{symbol}(P)=exp(A_ch(P){sign}{abs(centered_log):.15f})"


def build_artifact(
    exact_readout: dict[str, Any],
    affine_anchor: dict[str, Any],
    frontier: dict[str, Any],
    normalization_no_go: dict[str, Any],
) -> dict[str, Any]:
    centered_logs = [float(value) for value in exact_readout["centered_log_shape_exact"]]
    same_family_anchor = dict(affine_anchor["current_family_affine_anchor"])
    return {
        "artifact": "oph_charged_determinant_trace_lift_attachment_required",
        "generated_utc": _timestamp(),
        "status": "missing_theorem",
        "proof_status": "missing_theorem",
        "promotion_allowed": False,
        "public_promotion_allowed": False,
        "claim_tier": "exact_target_anchored_current_family_witness / continuation_gap",
        "source_only": False,
        "promotable": False,
        "public_theorem_value": None,
        "current_closed_chain": {
            "A_ch_to_charged_masses": True,
            "P_to_A_ch": False,
            "downstream_readout": {
                "m_e(P)": "exp(A_ch(P)+ell_e(P))",
                "m_mu(P)": "exp(A_ch(P)+ell_mu(P))",
                "m_tau(P)": "exp(A_ch(P)+ell_tau(P))",
            },
        },
        "required_identity": "3*A_ch(P)=sum_psi M_ch[psi]*log(q_psi(P))",
        "equivalent_defect": "N_det(P)=s_det(P)-sum_psi M_ch[psi]*log(q_psi(P))=0",
        "required_objects": {
            "charged_branch_generator": "C_hat_e, not only C_hat_e^{cand}",
            "sector_isolated_exponent_vector": "M_ch emitted before charged-mass comparison",
            "same_label_source_readout": "q_psi(P) emitted source-only",
            "determinant_trace_lift": "s_det^phys(P)=sum_psi M_ch[psi]*log(q_psi(P))",
            "no_target_leak_receipt": "NO_TARGET_LEAK_DAG_CHARGED_A_CH",
        },
        "missing_for_promotion": MISSING_FOR_PROMOTION,
        "forbidden_ancestors": FORBIDDEN_ANCESTORS,
        "conditional_readout_if_closed": {
            "centered_log_source_boundary": (
                "The displayed centered logs are the exact same-family witness checksum, "
                "not a source-only charged-mass prediction."
            ),
            "centered_logs": {
                "ell_e": centered_logs[0],
                "ell_mu": centered_logs[1],
                "ell_tau": centered_logs[2],
            },
            "formulas": {
                "electron": _formula("m_e", centered_logs[0]),
                "muon": _formula("m_mu", centered_logs[1]),
                "tau": _formula("m_tau", centered_logs[2]),
            },
        },
        "target_anchored_checksum_only": {
            "artifact_ref": artifact_ref(CURRENT_FAMILY_AFFINE_ANCHOR_JSON),
            "A_ch_current_family": same_family_anchor.get("value"),
            "g_e_current_family": affine_anchor.get("current_family_geometric_mean", {}).get("value"),
            "status": "checksum_not_promotion",
            "why": (
                "A_ch_current_family is computed from the same-family charged determinant. "
                "It is useful as an audit checksum but is a forbidden ancestor for a "
                "source-only A_ch(P) theorem."
            ),
        },
        "unit_convention": {
            "mass_space_anchor": "A_ch^GeV=(1/3)*log((m_e*m_mu*m_tau)/GeV^3)",
            "dimensionless_yukawa_anchor": "a_ch=(1/3)*log(det(Y_e))",
            "conversion": "A_ch^GeV=log(v(P)/(sqrt(2)*GeV))+a_ch",
            "preference": "derive dimensionless a_ch source-only, then display GeV masses after electroweak-scale convention is fixed",
        },
        "supporting_artifacts": {
            "determinant_character_frontier": {
                "artifact": frontier.get("artifact"),
                "artifact_ref": artifact_ref(DETERMINANT_CHARACTER_FRONTIER_JSON),
                "proof_status": frontier.get("proof_status"),
            },
            "trace_normalization_no_go": {
                "artifact": normalization_no_go.get("artifact"),
                "artifact_ref": artifact_ref(DETERMINANT_TRACE_NORMALIZATION_NO_GO_JSON),
                "proof_status": normalization_no_go.get("proof_status"),
            },
            "same_family_exact_readout_checksum": {
                "artifact": exact_readout.get("artifact"),
                "artifact_ref": artifact_ref(CURRENT_FAMILY_EXACT_READOUT_JSON),
                "proof_status": exact_readout.get("proof_status"),
            },
        },
        "non_promotion_rule": [
            "centered charged logs",
            "charged log ratios",
            "Koide-like shape identities",
            "target-centered exact same-family witnesses",
            "same-family determinant anchors",
            "compare-only g_e_star",
            "compare-only Delta_e_abs_star",
            "debug singular values",
            "fitted M_ch vectors",
            "PDG/CODATA charged-lepton rows",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the required charged determinant trace-lift attachment gate."
    )
    parser.add_argument("--exact-readout", default=str(CURRENT_FAMILY_EXACT_READOUT_JSON))
    parser.add_argument("--affine-anchor", default=str(CURRENT_FAMILY_AFFINE_ANCHOR_JSON))
    parser.add_argument("--frontier", default=str(DETERMINANT_CHARACTER_FRONTIER_JSON))
    parser.add_argument("--normalization-no-go", default=str(DETERMINANT_TRACE_NORMALIZATION_NO_GO_JSON))
    parser.add_argument("--output", default=str(DETERMINANT_TRACE_LIFT_ATTACHMENT_REQUIRED_JSON))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.exact_readout)),
        load_json(Path(args.affine_anchor)),
        load_json(Path(args.frontier)),
        load_json(Path(args.normalization_no_go)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
