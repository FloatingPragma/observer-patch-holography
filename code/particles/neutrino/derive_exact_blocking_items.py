#!/usr/bin/env python3
"""Summarize the exact remaining neutrino blockers on the live canonical tree.

Chain role: consolidate the builder-facing isotropic neutrino status, the
proof-facing scalar-certificate sufficiency theorem, and the remaining PMNS
blockers into one machine-readable audit.

Mathematics: no new neutrino formulas are introduced here. This is a status
consolidation layer over the already-emitted isotropic forward bundle, the
same-label scalar certificate shell, and the intrinsic eta-chain validation.

OPH-derived inputs: the live forward neutrino closure bundle, the same-label
scalar certificate shell, the proof-facing intrinsic eta demo, and the PMNS
shared-basis status file.

Output: an exact blocker audit plus a smaller current-snapshot summary.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORWARD_JSON = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
CERTIFICATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
PMNS_JSON = ROOT / "particles" / "runs" / "neutrino" / "pmns_from_shared_basis.json"
CHARGED_LEFT_JSON = ROOT / "particles" / "runs" / "neutrino" / "shared_charged_lepton_left_basis.json"
ETA_DEMO_JSON = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_eta_demo_payload.json"
INTRINSIC_VALIDATION_JSON = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_mixing_law_validation.json"
REPAIR_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
AMPLITUDE_BRIDGE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_absolute_amplitude_bridge.json"
BRIDGE_CANDIDATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"
IRREDUCIBILITY_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_attachment_irreducibility_theorem.json"
BRIDGE_SCALAR_CORRIDOR_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_attachment_bridge_scalar_corridor.json"
DEFAULT_EXACT_OUT = ROOT / "particles" / "runs" / "neutrino" / "exact_blocking_items.json"
DEFAULT_SUMMARY_OUT = ROOT / "particles" / "runs" / "neutrino" / "current_snapshot_blocker_summary.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _intrinsic_builder_mass_splittings(forward: dict, pmns: dict) -> dict[str, float]:
    return {
        "delta_m21_sq_gev2": pmns.get("intrinsic_delta_m21_sq_gev2", forward.get("delta_m21_sq_gev2")),
        "delta_m31_sq_gev2": pmns.get("intrinsic_delta_m31_sq_gev2", forward.get("delta_m31_sq_gev2")),
    }


def build_exact_blockers(
    forward: dict,
    certificate: dict,
    pmns: dict,
    charged_left: dict,
    eta_demo: dict,
    intrinsic_validation: dict,
    repair: dict,
    amplitude_bridge: dict | None,
    bridge_candidate: dict | None,
    irreducibility: dict | None,
    bridge_scalar_corridor: dict | None,
) -> tuple[dict, dict]:
    same_label_present = bool(certificate.get("sufficient_for_intrinsic_mass_eigenstates"))
    charged_basis_present = charged_left.get("status") == "closed"
    pmns_present = pmns.get("status") == "closed"
    repair_present = repair.get("artifact") == "oph_neutrino_weighted_cycle_repair"
    repair_shape_closed = repair.get("physical_window_status") == "pmns_and_hierarchy_repaired"
    absolute_normalization_open = repair.get("absolute_normalization_status") == "open_one_positive_scale"
    repair_scale_free_mass = dict(repair.get("scale_free_mass_normal_form") or {})
    repair_scale_free_dm2 = dict((repair.get("scale_free_dm2_normal_form") or {}).get("dm2") or {})
    repair_symbolic_family = dict(repair.get("symbolic_absolute_family") or {})
    intrinsic_builder_mass_splittings = _intrinsic_builder_mass_splittings(forward, pmns)
    reduced_bridge_object = dict(
        (irreducibility or {}).get("reduced_remaining_object")
        or (bridge_candidate or {}).get("smallest_exact_missing_object")
        or {}
    )
    physical_branch_closed = bool(pmns.get("physical_branch_closed", False)) or (
        repair_present and repair_shape_closed and not absolute_normalization_open
    )
    eta_payload = dict(eta_demo.get("eta_e") or {})
    exact_blockers = []
    if not same_label_present:
        exact_blockers.append(
            {
                "name": "live_same_label_scalar_certificate",
                "kind": "proof_facing_source_object",
                "current_snapshot_status": "present" if same_label_present else "absent",
                "required_contract": "oph_same_label_scalar_certificate_required_contract",
            }
        )
    if not charged_basis_present:
        exact_blockers.append(
            {
                "name": "shared_charged_lepton_left_basis",
                "kind": "pmns_and_public_flavor_row_object",
                "current_snapshot_status": "present" if charged_basis_present else "absent",
                "required_contract": "oph_shared_charged_left_basis_required_contract",
            }
        )

    branch_repair_required = (
        same_label_present
        and charged_basis_present
        and not repair_shape_closed
        and pmns_present
        and not physical_branch_closed
    )
    if branch_repair_required:
        exact_blockers.append(
            {
                "name": "physical_neutrino_branch_repair",
                "kind": "branch_selection_or_minimal_repair_theorem",
                "current_snapshot_status": "open",
                "required_contract": "emit_a_physically_correct_flavor_branch_or_prove_a_no_go_for_the_current_continuation_branch",
            }
        )
    if same_label_present and charged_basis_present and repair_shape_closed and absolute_normalization_open:
        exact_blockers.append(
            {
                "name": reduced_bridge_object.get("name", "one_positive_neutrino_bridge_correction_invariant"),
                "kind": "reduced_bridge_correction_invariant",
                "symbol": reduced_bridge_object.get("symbol", "C_nu"),
                "definition": reduced_bridge_object.get(
                    "definition",
                    "C_nu = B_nu / (I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1)",
                ),
                "bridge_reconstruction": reduced_bridge_object.get(
                    "bridge_reconstruction",
                    "B_nu = (I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1) * C_nu",
                ),
                "equivalence_theorem": reduced_bridge_object.get("equivalence_theorem"),
                "current_snapshot_status": reduced_bridge_object.get("status", "irreducible_on_current_corpus"),
                "required_contract": "emit_one_positive_neutrino_bridge_correction_invariant_above_the_emitted_proxy",
                "insufficiency_theorem": "neutrino_weighted_cycle_absolute_scale_no_go",
            }
        )

    fully_completed = same_label_present and charged_basis_present and pmns_present and physical_branch_closed
    reason_not_fully_completed = ""
    if not fully_completed:
        if same_label_present and charged_basis_present and repair_shape_closed and absolute_normalization_open:
            reason_not_fully_completed = (
                "The old isotropic branch has been repaired at the physical-pattern level: the weighted-cycle branch lands in the observed PMNS window and gives the right splitting hierarchy. "
                "The remaining exact attachment object is the reduced bridge-correction invariant C_nu above the already-emitted positive proxy P_nu. "
                "Until that reduced correction scalar is emitted, the branch still carries one positive absolute orbit, so absolute masses and absolute delta m^2 values remain compare-only unless an external atmospheric anchor is supplied."
            )
        elif branch_repair_required:
            reason_not_fully_completed = (
                "The continuation branch closes numerically from the same-label scalar certificate and shared charged basis, "
                "but the emitted PMNS angles and mass splittings are quantitatively wrong. The honest remaining object is a "
                "branch-repair or minimal-repair theorem at the flavor-to-neutrino interface or in the intrinsic neutrino ansatz."
            )
        else:
            reason_not_fully_completed = (
                "The intrinsic chain is exact once the same-label scalar certificate exists, but the current snapshot still "
                + (
                    "lacks a live physical certificate"
                    if not same_label_present
                    else "waits on PMNS writeback from the closed intrinsic bundle"
                )
                + (
                    " and still lacks the shared charged-lepton left basis required for PMNS/public flavor rows."
                    if not charged_basis_present
                    else "."
                )
            )

    exact_payload = {
        "artifact": "oph_exact_neutrino_blocker_audit_v8",
        "generated_utc": _timestamp(),
        "fully_completed": fully_completed,
        "reason_not_fully_completed": reason_not_fully_completed,
        "closed_theorem_chain": [
            "oph_native_scale_anchor_m_star_equals_v2_over_mu_u",
            "oph_fixed_cutoff_trace_pullback_metric",
            "neutrino_only_isotropy_obstruction",
            "same_label_scalar_certificate_sufficiency",
            "exact_principal_selector_from_centered_eta_class",
            "exact_depressed_cubic_intrinsic_spectrum",
            "mass_eigenstate_row_policy_nu1_nu2_nu3",
            "positive_load_balanced_least_distortion_midpoint_selector",
            "shape_closed_scale_invariant_left_basis",
            "pmns_from_shared_charged_and_intrinsic_bases",
        ],
        "current_isotropic_builder_facing_class": {
            "artifact": "oph_neutrino_only_edge_constant_centered_eta_class",
            "status": "exact_builder_facing_class_only",
            "edge_order": ["psi12", "psi23", "psi31"],
            "equivalence_class": "Any q_e = c > 0 gives the same centered class eta = 0.",
            "eta_e": {"psi12": 0.0, "psi23": 0.0, "psi31": 0.0},
            "mu_e_normalized": {"psi12": 1.0, "psi23": 1.0, "psi31": 1.0},
            "current_intrinsic_masses_gev": list(forward.get("masses_gev_sorted") or []),
            "current_delta_m21_sq_gev2": float(forward.get("delta_m21_sq_gev2") or 0.0),
            "current_delta_m31_sq_gev2": float(forward.get("delta_m31_sq_gev2") or 0.0),
            "ordering": forward.get("ordering_phase_certified"),
            "why_not_proof_facing": (
                "The isotropic neutrino-only bundle determines only the zero centered class. "
                "It does not determine the proof-facing overlap/gap certificate."
            ),
        },
        "demo_proof_facing_certificate_summary": {
            "eta_e": eta_payload,
        },
        "demo_intrinsic_result_summary": {
            "masses_gev_sorted": list(intrinsic_validation.get("collective_vector_actual_aligned") and [
                2.3929601069646055e-12,
                2.4048200109774875e-12,
                2.589606227283229e-12,
            ] or []),
            "delta_m21_sq_gev2": float(intrinsic_validation.get("solar_split_actual_gev2") or 0.0),
            "delta_m31_sq_gev2": float(intrinsic_validation.get("delta_m31_actual_gev2") or 0.0),
            "ordering": "normal_like_collective_dominance",
        },
        "exact_blocker_counts": {
            "same_label_proof_facing_continuous_dof_mod_common_scale": 5,
            "same_label_builder_facing_centered_eta_dof": 2,
            "charged_left_basis_artifact_dof_before_phase_quotients": 9,
        },
        "exact_blockers": exact_blockers,
        "neutrino_only_isotropy_obstruction": {
            "closed": True,
            "statement": (
                "The current forward neutrino bundle is exactly S_3-isotropic, so neutrino-only same-label readback stays edge-constant and cannot open the solar 1-2 split."
            ),
            "first_honest_solar_mover": "realized_arrow_pullback_from_flavor_gap_and_defect_certificates",
        },
        "no_hidden_discrete_branch": {
            "status": "closed" if repair_shape_closed else "not_yet_established",
            "open_discrete_blockers": [] if repair_shape_closed else ["branch_repair_or_discrete_selector_unknown"],
            "closed_discrete_witnesses": (
                [
                    "shape_closed_scale_invariant_left_basis",
                    "pmns_from_shared_charged_and_intrinsic_bases",
                    "cycle_basis_order_fixed_(f3,f1,f2)",
                    "holonomy_orientation_fixed_021",
                ]
                if repair_shape_closed
                else []
            ),
            "statement": (
                "The repaired weighted-cycle branch leaves only positive rescaling freedom; no unresolved discrete branch remains."
                if repair_shape_closed
                else "No hidden-discrete-branch theorem is not yet available before the repaired weighted-cycle branch closes."
            ),
        },
        "remaining_positive_scale_orbit": {
            "status": "open" if repair_shape_closed and absolute_normalization_open else "not_applicable",
            "group": "R_{>0}" if repair_shape_closed and absolute_normalization_open else None,
            "family_parameter": "lambda_nu > 0" if repair_shape_closed and absolute_normalization_open else None,
            "proof_obstruction": "positive_rescaling_nonidentifiability" if repair_shape_closed and absolute_normalization_open else None,
        },
        "live_continuation_branch_status": {
            "same_label_scalar_certificate_present": same_label_present,
            "shared_charged_left_basis_present": charged_basis_present,
            "pmns_present": pmns_present,
            "repair_artifact_present": repair_present,
            "status": (
                "physically_repaired_up_to_one_reduced_bridge_correction_invariant"
                if repair_shape_closed and absolute_normalization_open
                else (
                    "numerically_closed_but_quantitatively_wrong_branch"
                    if branch_repair_required
                    else "waiting_on_missing_upstream_artifacts"
                )
            ),
            "physical_branch_closed": physical_branch_closed,
            "no_hidden_discrete_branch": {
                "status": "closed" if repair_shape_closed else "not_yet_established",
                "open_discrete_blockers": [] if repair_shape_closed else ["branch_repair_or_discrete_selector_unknown"],
                "statement": (
                    "The repaired weighted-cycle branch already fixes the shared-basis discrete data; no unresolved discrete neutrino branch remains."
                    if repair_shape_closed
                    else "The current surface has not yet reduced the neutrino lane to a pure positive-scale orbit."
                ),
            },
            "remaining_positive_scale_orbit": {
                "status": "open" if repair_shape_closed and absolute_normalization_open else "not_applicable",
                "group": "R_{>0}" if repair_shape_closed and absolute_normalization_open else None,
                "family_parameter": "lambda_nu > 0" if repair_shape_closed and absolute_normalization_open else None,
            },
            "current_pmns_parameters": dict(
                (repair.get("pmns_observables") if repair_shape_closed else pmns.get("standard_pmns_parameters")) or {}
            ),
            "current_mass_splittings_gev2": (
                {
                    "status": "not_emitted_without_absolute_anchor",
                    "delta_m21_sq_gev2": None,
                    "delta_m31_sq_gev2": None,
                    "reason": (
                        "The repaired weighted-cycle branch closes only the dimensionless splitting pattern until "
                        "one positive reduced bridge-correction invariant is emitted."
                    ),
                }
                if repair_shape_closed and absolute_normalization_open
                else {
                    "status": "builder_intrinsic_snapshot",
                    **intrinsic_builder_mass_splittings,
                }
            ),
            "intrinsic_builder_mass_splittings_gev2": intrinsic_builder_mass_splittings,
            "repaired_branch_dimensionless_dm2": dict(repair.get("dimensionless_dm2") or {}),
            "compare_only_atmospheric_anchor": dict(repair.get("compare_only_atmospheric_anchor") or {}),
            "absolute_scale_no_go": (
                {
                    "status": "closed",
                    "theorem": "neutrino_weighted_cycle_absolute_scale_no_go",
                    "proof_obstruction": "positive_rescaling_nonidentifiability",
                    "statement": (
                        "The repaired weighted-cycle branch fixes only the dimensionless hierarchy and PMNS pattern. "
                        "For every lambda_nu > 0, the family m_i(lambda_nu) = lambda_nu * mhat_i and "
                        "Delta m^2_ij(lambda_nu) = lambda_nu^2 * Delta_hat_ij has the same dimensionless observables, "
                        "so no unique absolute scale is emitted without one further OPH theorem fixing lambda_nu."
                    ),
                    "absolute_family_parameter": "lambda_nu > 0",
                    "scale_free_mass_normal_form": repair_scale_free_mass,
                    "scale_free_dm2_normal_form_eV2": repair_scale_free_dm2,
                    "dimensionless_mass_family": list(repair_symbolic_family.get("absolute_masses") or []),
                    "dimensionless_dm2_family_eV2": dict(repair_symbolic_family.get("absolute_dm2") or {}),
                    "external_anchor_disallowed": {
                        "name": "atmospheric_delta_m32_sq",
                        "value_eV2": 2.438e-3,
                        "reason": "compare_only_external_oscillation_anchor",
                    },
                    "hard_separated_compare_only_adapter": {
                        "allowed_formula": "lambda_nu_cmp = sqrt(Delta m32^2_anchor / Delta_hat_32)",
                        "forbidden_feedback": "compare_only_anchor_must_not_feed_back_into_theorem_state_or_lambda_nu_emission",
                    },
                    "symbolic_absolute_outputs": {
                        "masses": [
                            "m1 = lambda_nu * mhat_1",
                            "m2 = lambda_nu * mhat_2",
                            "m3 = lambda_nu * mhat_3",
                        ],
                        "splittings": {
                            "21": "Delta m21^2 = lambda_nu^2 * Delta_hat_21",
                            "31": "Delta m31^2 = lambda_nu^2 * Delta_hat_31",
                            "32": "Delta m32^2 = lambda_nu^2 * Delta_hat_32",
                        },
                    },
                    "minimal_missing_object": "neutrino_weighted_cycle_absolute_attachment",
                    "sharper_attachment_object": (
                        bridge_candidate.get("exact_next_theorem_object") if bridge_candidate else None
                    ),
                    "smallest_exact_missing_object": reduced_bridge_object or None,
                    "current_attached_stack_irreducibility_theorem": (
                        {
                            "artifact": irreducibility.get("artifact"),
                            "status": irreducibility.get("status"),
                            "statement": irreducibility.get("theorem", {}).get("statement"),
                            "remaining_object": irreducibility.get("remaining_object"),
                            "reduced_remaining_object": irreducibility.get("reduced_remaining_object"),
                        }
                        if irreducibility
                        else None
                    ),
                    "immediate_theorem_gate": (
                        bridge_candidate.get("strictly_smaller_missing_clause") if bridge_candidate else None
                    ),
                    "compare_only_bridge_scalar_corridor": (
                        {
                            "artifact": bridge_scalar_corridor.get("artifact"),
                            "status": bridge_scalar_corridor.get("status"),
                            "primary_cross_route_corridor": bridge_scalar_corridor.get("primary_cross_route_corridor"),
                            "strongest_target_containing_bridge_scalar_corridor": bridge_scalar_corridor.get(
                                "strongest_target_containing_bridge_scalar_corridor"
                            ),
                            "shortlist_route_consensus_window": bridge_scalar_corridor.get("shortlist_route_consensus_window"),
                            "exact_reduced_correction_scalar": bridge_scalar_corridor.get("exact_reduced_correction_scalar"),
                            "bridge_correction_candidate_audit": bridge_scalar_corridor.get("bridge_correction_candidate_audit"),
                        }
                        if bridge_scalar_corridor
                        else None
                    ),
                    "absolute_amplitude_bridge_summary": (
                        {
                            "status": amplitude_bridge.get("status"),
                            "remaining_object": amplitude_bridge.get("remaining_object"),
                            "remaining_object_kind": amplitude_bridge.get("remaining_object_kind"),
                        }
                        if amplitude_bridge
                        else None
                    ),
                }
                if repair_shape_closed and absolute_normalization_open
                else None
            ),
        },
        "current_snapshot_scan": {
            "live_same_label_artifact_found": same_label_present,
            "live_charged_left_artifact_found": charged_basis_present,
            "live_pmns_artifact_found": pmns_present,
            "live_repair_artifact_found": repair_present,
        },
    }

    summary_payload = {
        "artifact": "oph_current_snapshot_blocker_summary_v8",
        "generated_utc": _timestamp(),
        "exact_remaining_blockers": [item["name"] for item in exact_blockers],
        "live_same_label_scalar_certificate_present": same_label_present,
        "shared_charged_left_basis_present": charged_basis_present,
        "pmns_present": pmns_present,
        "repair_artifact_present": repair_present,
        "physical_branch_closed": physical_branch_closed,
        "same_label_proof_facing_continuous_dof_mod_common_scale": 5,
        "same_label_builder_facing_centered_eta_dof": 2,
        "charged_left_basis_artifact_dof_before_phase_quotients": 9,
        "strongest_compare_only_bridge_scalar_corridor": (
            {
                "artifact": bridge_scalar_corridor.get("artifact"),
                "primary_cross_route_corridor": bridge_scalar_corridor.get("primary_cross_route_corridor"),
                "strongest_target_containing_bridge_scalar_corridor": bridge_scalar_corridor.get(
                    "strongest_target_containing_bridge_scalar_corridor"
                ),
                "shortlist_route_consensus_window": bridge_scalar_corridor.get("shortlist_route_consensus_window"),
            }
            if bridge_scalar_corridor
            else None
        ),
    }
    return exact_payload, summary_payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact neutrino blocker audit.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--certificate", default=str(CERTIFICATE_JSON))
    parser.add_argument("--pmns", default=str(PMNS_JSON))
    parser.add_argument("--charged-left", default=str(CHARGED_LEFT_JSON))
    parser.add_argument("--eta-demo", default=str(ETA_DEMO_JSON))
    parser.add_argument("--intrinsic-validation", default=str(INTRINSIC_VALIDATION_JSON))
    parser.add_argument("--repair", default=str(REPAIR_JSON))
    parser.add_argument("--amplitude-bridge", default=str(AMPLITUDE_BRIDGE_JSON))
    parser.add_argument("--bridge-candidate", default=str(BRIDGE_CANDIDATE_JSON))
    parser.add_argument("--irreducibility", default=str(IRREDUCIBILITY_JSON))
    parser.add_argument("--bridge-scalar-corridor", default=str(BRIDGE_SCALAR_CORRIDOR_JSON))
    parser.add_argument("--exact-output", default=str(DEFAULT_EXACT_OUT))
    parser.add_argument("--summary-output", default=str(DEFAULT_SUMMARY_OUT))
    args = parser.parse_args()

    exact_payload, summary_payload = build_exact_blockers(
        _load_json(Path(args.forward)),
        _load_json(Path(args.certificate)),
        _load_json(Path(args.pmns)),
        _load_json(Path(args.charged_left)),
        _load_json(Path(args.eta_demo)),
        _load_json(Path(args.intrinsic_validation)),
        _load_json(Path(args.repair)) if Path(args.repair).exists() else {},
        _load_json(Path(args.amplitude_bridge)) if Path(args.amplitude_bridge).exists() else None,
        _load_json(Path(args.bridge_candidate)) if Path(args.bridge_candidate).exists() else None,
        _load_json(Path(args.irreducibility)) if Path(args.irreducibility).exists() else None,
        _load_json(Path(args.bridge_scalar_corridor)) if Path(args.bridge_scalar_corridor).exists() else None,
    )

    exact_out = Path(args.exact_output)
    exact_out.parent.mkdir(parents=True, exist_ok=True)
    exact_out.write_text(json.dumps(exact_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary_out = Path(args.summary_output)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text(json.dumps(summary_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {exact_out}")
    print(f"saved: {summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
