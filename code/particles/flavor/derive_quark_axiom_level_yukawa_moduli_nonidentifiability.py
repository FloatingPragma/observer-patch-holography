#!/usr/bin/env python3
"""Emit the Axioms-1--5 quark-Yukawa non-definability theorem.

This artifact answers the strongest possible source-only question.  It does
not merely audit the present builders: it records an explicit counterfamily
of physically inequivalent one-Higgs, three-generation Yukawa packages that
obey the same OPH axioms and have the same MAR complexity vector.

No quark reference value, fitted spread, selected-family target, or numerical
flavor template is loaded.  The result uses no additional axiom.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, Any]:
    return {
        "artifact": "oph_quark_axiom_level_yukawa_moduli_nonidentifiability",
        "generated_utc": _timestamp(),
        "proof_status": "closed_axiom_level_nondefinability_theorem",
        "claim_tier": "axioms_1_5_source_only_obstruction",
        "scope": "OPH_Axioms_1_5_plus_structural_SM_branch_and_fixed_P",
        "additional_axioms_used": False,
        "theorem_grade_obstruction": True,
        "source_only_numeric_quark_spectrum_emitted": False,
        "public_numeric_quark_rows_allowed": False,
        "theorem_statement": (
            "OPH Axioms 1--5, a fixed pixel closure P, and the MAR-selected structural Standard-Model "
            "package do not define the quark Yukawa eigenvalues. For every generic admissible one-Higgs "
            "three-generation package there is a continuous positive-rescaling family with the same screen, "
            "overlap, recovery, gauge, anomaly, hypercharge, Higgs, generation, CKM, CP-capability, weak-UV, "
            "and MAR-complexity data but different Yukawa singular values and quark masses. Hence a unique "
            "quark mass spectrum does not factor through the stated axioms and source signature."
        ),
        "counterfamily": {
            "baseline": (
                "Y_q = U_qL diag(exp(mu_q + sigma_q*v_q1), exp(mu_q + sigma_q*v_q2), "
                "exp(mu_q + sigma_q*v_q3)) U_qR^dagger, q in {u,d}, sum_i v_qi = 0"
            ),
            "family": (
                "Y_q(lambda_q) = U_qL diag(exp(mu_q + lambda_q*sigma_q*v_q1), "
                "exp(mu_q + lambda_q*sigma_q*v_q2), exp(mu_q + lambda_q*sigma_q*v_q3)) "
                "U_qR^dagger"
            ),
            "parameter_space": "(lambda_u,lambda_d) in (R_{>0})^2",
            "free_action": (
                "(lambda_u,lambda_d).(sigma_u,sigma_d) = "
                "(lambda_u*sigma_u,lambda_d*sigma_d)"
            ),
            "simple_spectrum_preserved": True,
            "left_and_right_frames_preserved": True,
            "CKM_matrix_preserved": True,
            "CP_capability_preserved": True,
            "normalized_determinants_preserved": (
                "det(exp(lambda_q*sigma_q*v_q)) = "
                "exp(lambda_q*sigma_q*sum_i(v_qi)) = 1"
            ),
            "quark_spreads_changed": True,
            "quark_mass_singular_values_changed": True,
        },
        "axiom_invariance_audit": {
            "Axiom_1_screen_net": "unchanged by Yukawa-eigenvalue rescaling",
            "Axiom_2_overlap_consistency": "unchanged by Yukawa-eigenvalue rescaling",
            "Axiom_3_MaxEnt_and_refinement": {
                "gauge_invariant_local_Yukawa_densities_allowed": True,
                "constraint_values_or_multipliers_numerically_fixed_by_axiom": False,
                "map_from_P_to_Yukawa_multipliers_supplied": False,
                "homogeneous_refinement_selects_initial_modulus": False,
            },
            "Axiom_4_recoverable_generalized_entropy": (
                "contains no numerical flavor-eigenvalue selector"
            ),
            "Axiom_5_MAR": {
                "complexity_vector": "(chi_cpl,N_nonab,N_c,N_g)",
                "Yukawa_completability_preserved": True,
                "intrinsic_CKM_CP_capability_preserved": True,
                "weak_sector_UV_counting_clause_preserved": True,
                "Yukawa_eigenvalues_are_components_of_complexity_vector": False,
                "different_Yukawa_invariants_are_physically_equivalent": False,
                "counterfamily_members_have_equal_MAR_score": True,
                "counterfamily_members_remain_distinct_MAR_minima": True,
            },
            "fixed_P": "unchanged across the counterfamily",
        },
        "why_minimal_Yukawas_does_not_follow_from_MAR": {
            "current_order_contains_no_continuous_Yukawa_coordinate": True,
            "positive_rescaling_has_no_smallest_positive_element": True,
            "infimum": "lambda_u=lambda_d=0",
            "zero_limit_effect": "massless or degenerate quarks rather than the observed spectrum",
            "adding_a_norm_entropy_cost_description_length_or_RG_functional": (
                "would change the MAR selection rule and therefore add physical selection content"
            ),
        },
        "proof_steps": [
            "Fix any generic admissible one-Higgs three-generation Yukawa package with simple spectra.",
            "Apply independent positive rescalings to the two centered log-spectrum profiles while fixing their frames.",
            "Gauge representations, anomaly cancellation, hypercharges, Higgs content, CKM data, and CP capability are unchanged.",
            "The MAR complexity vector is unchanged because it contains only structural discrete entries.",
            "The packages are physically inequivalent because their Yukawa singular values differ and physical equivalence preserves Yukawa invariants.",
            "Therefore the same OPH axioms and fixed P admit multiple quark spectra, contradicting uniqueness of any proposed source-only mass map.",
        ],
        "corollaries": {
            "unique_source_map_P_to_sigma_u_sigma_d_exists": False,
            "unique_source_map_P_to_six_quark_masses_exists": False,
            "selected_frame_representative_independence_breaks_counterfamily": False,
            "refinement_naturality_breaks_counterfamily": False,
            "MAR_breaks_counterfamily_under_current_definition": False,
            "artifact_level_two_modulus_obstruction_is_axiomatically_expected": True,
        },
        "reference_data_policy": {
            "direct_input_artifacts": [],
            "quark_reference_values_consumed": False,
            "PDG_or_API_rows_consumed": False,
            "current_family_targets_consumed": False,
            "selected_class_target_witnesses_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "no_target_leak_by_construction": True,
        },
        "public_policy": {
            "numeric_prediction_status": "not_defined_by_current_axioms",
            "allowed_public_result": "axiom_level_nondefinability_theorem",
            "target_reconstruction_tables": "audit_only_not_predictions",
        },
        "notes": [
            "This is a non-definability theorem under the current axioms, not a claim that every conceivable future theory must leave Yukawa values free.",
            "No additional axiom or numerical normalization principle is introduced.",
            "A positive numerical theorem would have to refute this counterfamily using a consequence already present in Axioms 1--5; no such consequence is currently stated or emitted.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the OPH Axioms-1--5 quark-Yukawa non-definability theorem."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
