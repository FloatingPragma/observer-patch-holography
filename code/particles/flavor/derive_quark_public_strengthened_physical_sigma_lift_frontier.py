#!/usr/bin/env python3
"""Emit the public sigma-lift frontier and its source non-identifiability boundary."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
STRENGTHENED_LOCAL_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
ABSOLUTE_COLLAPSE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_absolute_readout_algebraic_collapse.json"
EXACT_PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_pdg_end_to_end_theorem.json"
EXACT_YUKAWA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"
LINE_LIFT_JSON = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
GENERATOR_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
PUBLIC_SIGMA_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_public_strengthened_physical_sigma_lift_frontier.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _upstream_attached_data_obstruction(line_lift: dict[str, Any], generator: dict[str, Any]) -> dict[str, float]:
    lamb = np.asarray(generator["simple_spectrum_certificate"]["eigenvalues"], dtype=float)
    diagnostics = list(line_lift["same_refinement_edge_diagnostic_by_refinement"])
    if len(diagnostics) < 2:
        return {
            "commutator_operator_norm": 0.0,
            "projector_defect_operator_norm": 0.0,
            "commutator_over_defect": 0.0,
            "commutator_over_defect_squared": 0.0,
        }

    projectors_by_level = []
    for entry in diagnostics[:2]:
        projectors_by_level.append([_decode_complex_matrix(projector) for projector in entry["projectors"]])

    centered = []
    for projectors in projectors_by_level:
        current = sum(lamb[idx] * projectors[idx] for idx in range(3))
        current = current - np.trace(current) / 3.0 * np.eye(3, dtype=complex)
        centered.append(current)

    commutator = centered[1] @ centered[0] - centered[0] @ centered[1]
    defect = max(np.linalg.norm(projectors_by_level[1][idx] - projectors_by_level[0][idx], ord=2) for idx in range(3))
    commutator_norm = float(np.linalg.norm(commutator, ord=2))
    return {
        "commutator_operator_norm": commutator_norm,
        "projector_defect_operator_norm": float(defect),
        "commutator_over_defect": float(commutator_norm / defect) if defect > 0.0 else 0.0,
        "commutator_over_defect_squared": float(commutator_norm / (defect * defect)) if defect > 0.0 else 0.0,
    }


def build_artifact(
    strengthened_local: dict[str, Any],
    absolute_collapse: dict[str, Any],
    exact_pdg: dict[str, Any],
    exact_yukawa: dict[str, Any],
    line_lift: dict[str, Any],
    generator: dict[str, Any],
    public_sigma_theorem: dict[str, Any],
) -> dict[str, Any]:
    sigma = dict(strengthened_local["theorem_grade_physical_sigma_datum"])
    masses = dict(exact_pdg["exact_running_values_gev"])
    exact_forward = dict(exact_yukawa["forward_yukawa_artifact"])
    upstream_obstruction = _upstream_attached_data_obstruction(line_lift, generator)
    public_promotion_allowed = public_sigma_theorem.get("public_promotion_allowed") is True
    return {
        "artifact": "oph_quark_public_strengthened_physical_sigma_lift_frontier",
        "generated_utc": _timestamp(),
        "proof_status": public_sigma_theorem["proof_status"],
        "frontier_role": (
            "resolved_direct_public_descent_theorem_plus_historical_upstream_alternative"
            if public_promotion_allowed
            else "blocked_direct_public_descent_candidate_plus_historical_upstream_alternative"
        ),
        "public_promotion_allowed": public_promotion_allowed,
        "non_circularity_status": public_sigma_theorem.get("non_circularity_status"),
        "resolved_by_theorem_artifact": public_sigma_theorem["artifact"],
        "final_public_theorem_candidate": {
            "id": public_sigma_theorem["theorem_id"],
            "kind": "public_bridge_descent_of_exact_sigma_readout",
            "statement": public_sigma_theorem["theorem_statement"],
            "induces_global_contract": public_sigma_theorem["induces_global_contract"],
            "must_emit": {
                "sigma_ud_phys_element": "{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM, ckm_invariants}",
                "physical_sigma_datum": public_sigma_theorem["descended_physical_sigma_datum"],
            },
            "selected_public_exact_sigma_datum": public_sigma_theorem["descended_physical_sigma_datum"],
            "selected_public_frame_class": public_sigma_theorem["selected_public_physical_frame_class"],
            "why_sufficient": public_sigma_theorem["why_sufficient"],
        },
        "algebraic_consequence_after_closure": {
            "supporting_algebraic_collapse_artifact": absolute_collapse["artifact"],
            "statement": (
                "If a new no-target source theorem emits both positive spread moduli, selected-fiber descent makes "
                "them representative-independent and the affine mean law emits the conditional mass coordinates. "
                "The present target packet remains mixed-convention, and physical Yukawa matrices additionally "
                "require common-scale RG transport and dimensionless Higgs normalization."
            ),
            "target_audit_mass_coordinates_gev": masses,
            "target_audit_mass_textures": {
                "artifact": exact_forward["artifact"],
                "forward_certified": exact_forward["forward_certified"],
                "certification_status": exact_forward["certification_status"],
                "physical_yukawa_certified": False,
                "singular_values_u": exact_forward["singular_values_u"],
                "singular_values_d": exact_forward["singular_values_d"],
                "V_CKM": exact_forward["V_CKM"],
                "jarlskog": exact_forward["jarlskog"],
            },
        },
        "algebraic_consequence_if_closed": {
            "supporting_algebraic_collapse_artifact": absolute_collapse["artifact"],
            "statement": (
                "This consequence is unrealized on the current corpus. Its source equations leave an exact "
                "(R_{>0})^2 spread fiber, and selected-class descent is representative independence rather than source selection."
            ),
            "target_audit_mass_coordinates_gev": masses,
            "target_audit_mass_textures": {
                "artifact": exact_forward["artifact"],
                "forward_certified": exact_forward["forward_certified"],
                "certification_status": exact_forward["certification_status"],
                "physical_yukawa_certified": False,
                "singular_values_u": exact_forward["singular_values_u"],
                "singular_values_d": exact_forward["singular_values_d"],
                "V_CKM": exact_forward["V_CKM"],
                "jarlskog": exact_forward["jarlskog"],
            },
        },
        "alternate_upstream_route": {
            "id": "oph_generation_bundle_branch_generator_splitting",
            "status": "upstream_alternative_route_currently_deprioritized",
            "why_relevant": (
                "If the persistent simple-spectrum splitting theorem closes, the overlap-edge line lift upgrades from "
                "candidate readout to public same-label transport class, providing an upstream route into the public "
                "same-label physical sigma carrier."
            ),
            "statement": (
                "On the realized generation bundle, the centered compressed generation-bundle branch generator has a "
                "persistent simple spectrum with label-stable one-dimensional Riesz projectors functorial under every "
                "refinement intertwiner, and the same-label transport classes are presentation-independent modulo "
                "objectwise U(1)."
            ),
            "smaller_exact_missing_clause": generator["actual_generator_transfer_candidate"]["smaller_exact_missing_clause"],
            "upstream_object": line_lift["upstream_missing_object"],
            "current_attached_data_obstruction": {
                "reason": (
                    "On the displayed level-0/level-1 projector systems, the natural centered operators built from the "
                    "proxy ordered spectrum have a commutator that is linear in the projector defect rather than "
                    "quadratic-small, so the displayed current artifacts do not themselves force the transfer lemma."
                ),
                **upstream_obstruction,
            },
        },
        "supporting_local_artifacts": {
            "restricted_strengthened_physical_sigma_lift": strengthened_local["artifact"],
            "absolute_readout_algebraic_collapse": absolute_collapse["artifact"],
            "exact_pdg_end_to_end_theorem": exact_pdg["artifact"],
            "exact_yukawa_end_to_end_theorem": exact_yukawa["artifact"],
            "overlap_edge_line_lift": line_lift["artifact"],
            "generation_bundle_branch_generator": generator["artifact"],
        },
        "notes": [
            (
                "This artifact records a source-only public sigma theorem as closed on the selected physical quark frame class."
                if public_promotion_allowed
                else "This artifact records selected-class target-audit data and a theorem-grade two-modulus source-spread obstruction; it does not promote numeric quark rows."
            ),
            "The upstream generator route remains candidate-only and begins from a hand-written family-kernel template.",
            "Mixed-convention target mass coordinates and dimensionful mass textures are retained only for contract-level audit.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the final public-frontier theorem package above the local exact-PDG quark chain.")
    parser.add_argument("--strengthened-local", default=str(STRENGTHENED_LOCAL_JSON))
    parser.add_argument("--absolute-collapse", default=str(ABSOLUTE_COLLAPSE_JSON))
    parser.add_argument("--exact-pdg", default=str(EXACT_PDG_JSON))
    parser.add_argument("--exact-yukawa", default=str(EXACT_YUKAWA_JSON))
    parser.add_argument("--line-lift", default=str(LINE_LIFT_JSON))
    parser.add_argument("--generator", default=str(GENERATOR_JSON))
    parser.add_argument("--public-sigma-theorem", default=str(PUBLIC_SIGMA_THEOREM_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.strengthened_local)),
        _load_json(Path(args.absolute_collapse)),
        _load_json(Path(args.exact_pdg)),
        _load_json(Path(args.exact_yukawa)),
        _load_json(Path(args.line_lift)),
        _load_json(Path(args.generator)),
        _load_json(Path(args.public_sigma_theorem)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
