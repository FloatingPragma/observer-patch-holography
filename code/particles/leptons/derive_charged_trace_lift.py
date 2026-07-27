#!/usr/bin/env python3
"""Audit and, only with a complete source certificate, certify the charged trace lift.

The current repository contains a centered same-family checksum, a numerical
same-label arrow readback, and conditional determinant-line algebra.  None of
those objects fixes the affine origin of the charged determinant line.  This
builder makes that boundary machine-readable and deliberately fails closed.

An optional future source certificate can close the gate, but only when it
supplies every source-side object named in ``SOURCE_CERTIFICATE_CONTRACT`` and
certifies both leakage and the attachment residual on the singleton interval
``[0, 0]``.  Merely having an interval which contains zero is not sufficient.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
RUNS = ROOT / "particles" / "runs"

EXACT_READOUT_JSON = RUNS / "leptons" / "lepton_current_family_exact_readout.json"
SAME_LABEL_READBACK_JSON = RUNS / "neutrino" / "realized_same_label_gap_defect_readback.json"
DETERMINANT_FRONTIER_JSON = RUNS / "leptons" / "charged_determinant_character_frontier.json"
D10_REDUCTION_JSON = RUNS / "leptons" / "charged_p_to_affine_anchor_reduction.json"
NORMALIZATION_NO_GO_JSON = RUNS / "leptons" / "charged_determinant_trace_normalization_no_go.json"
CHARGED_BUDGET_JSON = RUNS / "flavor" / "charged_budget_transport.json"
DEFAULT_SOURCE_CERTIFICATE = RUNS / "leptons" / "charged_trace_lift_source_certificate.json"
DEFAULT_OUT = RUNS / "leptons" / "charged_trace_lift_theorem.json"

EDGE_ORDER = ("psi12", "psi23", "psi31")
MASS_ORDER = ("electron", "muon", "tau")
FORBIDDEN_TARGET_ANCESTORS = (
    "particle_reference_values",
    "lepton_current_family_exact_readout",
    "charged_lepton_reference_targets",
    "pdg",
    "codata",
    "m_e",
    "m_mu",
    "m_tau",
    "target",
    "compare_only",
)

SOURCE_CERTIFICATE_CONTRACT = {
    "artifact": "oph_charged_trace_lift_source_certificate",
    "required_flags": ["theorem_grade", "source_only", "no_target_leak"],
    "required_objects": [
        "source-emitted charged central projector/idempotent",
        "zero-leakage D10 determinant factorization over charged and rest sectors",
        "numeric sector-isolated multiplicity vector M_ch on psi12/psi23/psi31",
        "source-closed, stage-indexed same-label q_psi(r) with a physical label map",
        "D10-fixed reference stage r_0",
        "normalized charged determinant-line basepoint/trivialization at r_0",
        "certified P interval and mass-space affine-anchor interval",
        "named electroweak-scale conversion into GeV",
        "no-target-leak ancestry receipt",
    ],
    "required_zero_intervals": [
        "factorization_lemma.leakage_bound.interval",
        "attachment_identity_residual.interval",
    ],
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _decimal(value: Any) -> Decimal | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError):
        return None
    return result if result.is_finite() else None


def _interval(value: Any) -> tuple[Decimal, Decimal] | None:
    if isinstance(value, Mapping):
        value = value.get("interval", [value.get("lower"), value.get("upper")])
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        return None
    lower, upper = _decimal(value[0]), _decimal(value[1])
    if lower is None or upper is None or lower > upper:
        return None
    return lower, upper


def _singleton_zero_interval(value: Any) -> bool:
    parsed = _interval(value)
    return parsed == (Decimal(0), Decimal(0))


def _all_numeric_multiplicities(value: Any) -> bool:
    if not isinstance(value, Mapping) or set(value) != set(EDGE_ORDER):
        return False
    return all(_decimal(value[edge]) is not None for edge in EDGE_ORDER)


def _ancestry_is_target_free(value: Any) -> bool:
    if not isinstance(value, list) or not value:
        return False
    lowered = [str(item).lower() for item in value]
    return not any(token in ancestor for ancestor in lowered for token in FORBIDDEN_TARGET_ANCESTORS)


def _source_certificate_checks(source: Mapping[str, Any] | None) -> dict[str, bool]:
    source = source or {}
    factorization = source.get("factorization_lemma") or {}
    constant = source.get("uncentered_lift_constant") or {}
    residual = source.get("attachment_identity_residual") or {}
    label_map = source.get("physical_label_map") or {}
    return {
        "certificate_present": bool(source),
        "artifact_type": source.get("artifact") == SOURCE_CERTIFICATE_CONTRACT["artifact"],
        "theorem_grade": source.get("theorem_grade") is True,
        "source_only": source.get("source_only") is True,
        "no_target_leak_flag": source.get("no_target_leak") is True,
        "no_target_leak_ancestry": _ancestry_is_target_free(source.get("ancestry")),
        "numeric_sector_isolated_M_ch": _all_numeric_multiplicities(
            source.get("multiplicity_vector_M_ch")
        ),
        "source_closed_stage_indexed_q": source.get("source_closed_stage_indexed_q") is True,
        "physical_label_map": (
            isinstance(label_map, Mapping)
            and set(label_map) == set(EDGE_ORDER)
            and set(label_map.values()) == set(MASS_ORDER)
        ),
        "charged_central_projector": bool(source.get("charged_central_projector")),
        "factorization_certified": factorization.get("status") == "certified",
        "zero_leakage": _singleton_zero_interval(factorization.get("leakage_bound")),
        "reference_stage_named": bool(constant.get("reference_stage")),
        "lift_constant_source_named": bool(constant.get("source_object_name")),
        "lift_constant_numeric": _decimal(constant.get("value")) is not None,
        "lift_constant_interval": _interval(constant.get("interval")) is not None,
        "determinant_scalar_interval": _interval(source.get("determinant_scalar_interval")) is not None,
        "P_interval": _interval(source.get("P_interval")) is not None,
        "mass_space_affine_anchor_interval": (
            _interval(source.get("mass_space_affine_anchor_log_gev_interval")) is not None
        ),
        "mass_space_anchor_source_named": bool(source.get("mass_space_anchor_source_object")),
        "residual_certified": residual.get("certified") is True,
        "residual_singleton_zero": _singleton_zero_interval(residual.get("interval")),
    }


def classify_source_certificate(source: Mapping[str, Any] | None) -> tuple[str, dict[str, bool]]:
    """Return the claim label and every promotion check for a candidate certificate."""

    checks = _source_certificate_checks(source)
    if all(checks.values()):
        return "trace_lift_certified", checks
    if not checks["certificate_present"]:
        return "no_go_confirmed_new_source_needed", checks
    return "open", checks


def _ratio_regression(exact_readout: Mapping[str, Any]) -> dict[str, Any]:
    centered = [float(value) for value in exact_readout["centered_log_shape_exact"]]
    witness = [float(value) for value in exact_readout["predicted_singular_values_abs"]]
    pairs = ((0, 1, "mu_over_e"), (1, 2, "tau_over_mu"), (0, 2, "tau_over_e"))
    rows: dict[str, Any] = {}
    max_relative_residual = 0.0
    for left, right, name in pairs:
        witness_ratio = witness[right] / witness[left]
        reconstructed_ratio = math.exp(centered[right] - centered[left])
        relative_residual = abs(reconstructed_ratio / witness_ratio - 1.0)
        max_relative_residual = max(max_relative_residual, relative_residual)
        rows[name] = {
            "witness": witness_ratio,
            "from_centered_family_functional": reconstructed_ratio,
            "relative_residual_abs": relative_residual,
        }
    passed = max_relative_residual <= 5.0e-15
    return {
        "status": "passed" if passed else "failed",
        "theorem_scope": exact_readout.get("theorem_scope"),
        "role": "target-anchored checksum only; never an ancestor of the trace-lift claim",
        "centered_log_shape": centered,
        "witness_masses_gev": witness,
        "ratios": rows,
        "max_relative_residual_abs": max_relative_residual,
    }


def _conditional_mass_rows(
    exact_readout: Mapping[str, Any], source: Mapping[str, Any] | None
) -> list[dict[str, Any]]:
    source = source or {}
    p_interval = _interval(source.get("P_interval"))
    anchor_interval = _interval(source.get("mass_space_affine_anchor_log_gev_interval"))
    if p_interval is None or anchor_interval is None:
        return []
    centered = [float(value) for value in exact_readout["centered_log_shape_exact"]]
    p_bounds = [str(p_interval[0]), str(p_interval[1])]
    rows = []
    for particle, centered_log in zip(MASS_ORDER, centered, strict=True):
        lower = math.exp(float(anchor_interval[0]) + centered_log)
        upper = math.exp(float(anchor_interval[1]) + centered_log)
        rows.append(
            {
                "particle": particle,
                "status": "conditional_on_P",
                "unit": "GeV",
                "P_interval": p_bounds,
                "mass_interval": [lower, upper],
                "formula": "m_i(P)=exp(A_ch^GeV(P)+ell_i_centered(P)) GeV",
                "centered_log_checksum": centered_log,
                "mass_space_anchor_source_object": source.get("mass_space_anchor_source_object"),
            }
        )
    return rows


def _current_source_audit(
    same_label: Mapping[str, Any],
    frontier: Mapping[str, Any],
    d10_reduction: Mapping[str, Any],
    charged_budget: Mapping[str, Any],
) -> dict[str, Any]:
    sector_witnesses = charged_budget.get("sector_isolation_witness_by_sector") or {}
    electron_sector = sector_witnesses.get("e") or {}
    q_values = same_label.get("q_e") or same_label.get("derived_q_e") or {}
    return {
        "same_label_readback": {
            "artifact": same_label.get("artifact"),
            "artifact_ref": _artifact_ref(SAME_LABEL_READBACK_JSON),
            "labels": list(q_values),
            "physical_label_assignment_certified": False,
            "stage_indexed_q_e_of_r": False,
            "proof_status": same_label.get("proof_status"),
            "source_only_physical_input_eligible": same_label.get(
                "source_only_physical_input_eligible"
            ),
            "missing_source_objects": (same_label.get("source_closure_status") or {}).get(
                "missing_objects", []
            ),
        },
        "multiplicity_vector": {
            "numeric_M_ch_present": False,
            "status": "formal_symbol_only",
            "evidence": frontier.get("determinant_character_candidate", {}).get("formula"),
        },
        "d10_determinant_landing": {
            "numeric_or_interval_s_det_present": False,
            "status": "conditional_reduction_only",
            "proof_status": d10_reduction.get("proof_status"),
            "forbidden_claims": d10_reduction.get("do_not_claim", []),
        },
        "charged_sector_isolation": {
            "d9_representation_role_projector_constructible": True,
            "source_emitted_d10_charged_projector_present": False,
            "d10_determinant_compatibility_present": False,
            "electron_budget_witness": {
                "sector_local_closed": electron_sector.get("sector_local_closed"),
                "status": electron_sector.get("status"),
                "value": electron_sector.get("value"),
            },
        },
        "reference_stage": {
            "r_0_named": False,
            "stage_indexed_readback_present": False,
            "d10_repair_square_present": False,
            "basepoint_determinant_trivialization_present": False,
        },
    }


def _representation_sector_factorization() -> dict[str, Any]:
    """Return the strongest sector lemma supported by the declared D9 carrier.

    This lemma is deliberately narrower than the issue-level D10 attachment.
    It isolates the gauge-equivariant charged Yukawa *channel* once a Yukawa
    tensor is supplied.  It neither supplies that tensor nor normalizes its
    determinant line.
    """

    return {
        "status": "certified_representation_channel_only_d10_attachment_open",
        "theorem_scope": (
            "declared D9 Standard-Model G-module, a chosen neutral Higgs vacuum, and a supplied "
            "gauge-equivariant charged Yukawa tensor"
        ),
        "joint_charge_spectrum": {
            "coordinates": ["C2_SU3", "C2_SU2", "Y"],
            "Q": ["4/3", "3/4", "1/6"],
            "u^c": ["4/3", "0", "-2/3"],
            "d^c": ["4/3", "0", "1/3"],
            "L": ["0", "3/4", "-1/2"],
            "e^c": ["0", "0", "1"],
            "H": ["0", "3/4", "1/2"],
        },
        "projector_construction": (
            "Choose a generic Z=a C2_SU3+b C2_SU2+c Y with distinct eigenvalues z_R on the "
            "six declared roles and set p_R=prod_{S!=R}(Z-z_S I)/(z_R-z_S). The charged "
            "Yukawa role is selected by p_L, p_e^c and the neutral-Higgs invariant tensor."
        ),
        "charged_hom_space": (
            "Hom_G(L tensor H^dagger,e^c) tensor Hom(F_L,F_e) "
            "isomorphic to C tensor Mat_3(C)"
        ),
        "determinant_line": "Det_ch = det(F_e) tensor det(F_L)^*",
        "direct_sum_naturality": (
            "Det(D_ch direct_sum D_rest) isomorphic to Det(D_ch) tensor Det(D_rest)"
        ),
        "leakage_bound": {
            "interval": ["0", "0"],
            "kind": "exact_gauge_selection_rule_quark_leakage_on_supplied_Y_e",
            "certified_zero": True,
        },
        "does_not_certify": [
            "a fermionic D10 determinant landing",
            "rank-one e/mu/tau projectors inside the generation multiplicity space",
            "a map from psi12/psi23/psi31 to physical charged labels",
            "a nonzero normalized section of Det_ch",
            "the issue-level attachment residual",
        ],
    }


def _existing_axiom_independence_theorem() -> dict[str, Any]:
    return {
        "status": "proved_by_one_parameter_countermodel",
        "family": "Y_e^(kappa)(r) = exp(kappa) Y_e(r), for every real kappa and every r",
        "preserved_existing_antecedents": [
            "D9 gauge group, matter representations, Higgs role, hypercharge lattice, N_c and N_g",
            "anomaly cancellation and the existence of the allowed L H^dagger e^c Yukawa monomial",
            "Yukawa-completability, rank, and the declared structural packet",
            "D10 gauge-coupling, weak-scale, and electroweak pole outputs",
            "D9 representation-role projectors and exact quark-channel leakage zero",
            "all refinement squares when the same constant rescaling is used at every stage",
            "the centered charged log spectrum and every charged mass ratio",
            "the current q_psi arrow readback, because no existing theorem attaches it to Y_e",
        ],
        "changed_quantities": {
            "log_abs_det_Y_e": "log_abs_det_Y_e + 3 kappa",
            "charged_masses": "m_i -> exp(kappa) m_i",
            "normalization_defect": "N_det -> N_det + 3 kappa",
        },
        "generation_center_obstruction": {
            "commutant": "End_G(W_L tensor C^3) = I_W_L tensor Mat_3(C)",
            "center": "Z(End_G(W_L tensor C^3)) = C I",
            "consequence": (
                "The gauge/collar center selects the full L and e^c isotypic roles but cannot "
                "select three physical rank-one family lines or their determinant normalization."
            ),
        },
        "even_if_total_determinant_is_granted": (
            "For an invertible n_rest-dimensional rest block, simultaneously take "
            "D_rest -> exp(-3 kappa/n_rest) D_rest. The charged and rest log determinants "
            "shift by +3 kappa and -3 kappa, so the total and exact factorization remain fixed."
        ),
        "logical_conclusion": (
            "The normalized attachment is not a consequence of the currently declared antecedents. "
            "A proof using only those antecedents would take the same value on this family and is "
            "therefore contradicted by the displayed determinant shift."
        ),
    }


def _candidate_escape_audit(same_label: Mapping[str, Any]) -> dict[str, Any]:
    q_values = same_label.get("q_e") or same_label.get("derived_q_e") or {}
    q_logs = {
        edge: math.log(float(q_values[edge]))
        for edge in EDGE_ORDER
        if edge in q_values and float(q_values[edge]) > 0.0
    }
    symmetric_character = sum(q_logs.values()) if len(q_logs) == len(EDGE_ORDER) else None
    stage5_log_det_yukawa = 0.5 * math.log(2.0) - 14.0 * math.log(6.0)
    return {
        "d10_determinant_premise": {
            "status": "contradicted_by_canonical_artifacts",
            "finding": (
                "The live D10 package emits gauge couplings, v(P), and electroweak pole/readout "
                "quantities. charged_p_to_affine_anchor_reduction.json is an if-landing reduction "
                "and explicitly forbids claiming physical Y_e(P), s_det(P), A_ch(P), or masses."
            ),
            "numeric_or_interval_matter_s_det_present": False,
        },
        "reference_stage_q_candidate": {
            "status": "fails",
            "finding": (
                "The current q_psi values are template-dependent arrows from refinement level 0 "
                "to 1, not q_psi(r_0). On an identity/reference arrow the declared defect "
                "1-overlap^2 is zero, hence q=sqrt(gap*defect)=0 and log q is not a finite "
                "determinant-line trivialization."
            ),
        },
        "stage5_continuation": {
            "status": "diagnostic_only_not_a_trace_lift_certificate",
            "assumptions": [
                "uniform Z6 ensemble epsilon=1/6",
                "balanced circulant/Koide carrier",
                "phase delta=(N_c+1)/(2 N_c N_g)=2/9",
                "hand-built exponent prescription n=(7,4,3), including a residual scan for n_e",
                "extra factor 2^(1/6) and a sorted-root-to-e/mu/tau assignment",
            ],
            "algebraic_consequence_inside_that_ansatz": {
                "det_M_e_GeV_cubed": "v(P)^3 / (2 * 6^14)",
                "log_abs_det_Y_e": "(1/2) log 2 - 14 log 6",
                "log_abs_det_Y_e_numeric": stage5_log_det_yukawa,
            },
            "current_q_diagnostic": {
                "q_psi": {edge: q_values.get(edge) for edge in EDGE_ORDER},
                "natural_symmetric_M_candidate": {edge: 1 for edge in EDGE_ORDER},
                "S_(1,1,1)": symmetric_character,
                "residual_against_stage5_log_det_Y_e": (
                    stage5_log_det_yukawa - symmetric_character
                    if symmetric_character is not None
                    else None
                ),
            },
            "why_not_promotable": (
                "The determinant formula is algebra inside a quarantined phenomenological "
                "continuation; it is not emitted by D9/D10 or attached to q_psi."
            ),
        },
        "thomson_inversion": {
            "status": "circular",
            "finding": (
                "paper_math first constructs the Stage-5 charged masses and then evaluates "
                "Delta_lep from those masses with mass_source=internal_stage5_continuation. "
                "Inverting that transport only recovers its own mass input."
            ),
        },
    }


def build_artifact(
    exact_readout: Mapping[str, Any],
    same_label: Mapping[str, Any],
    frontier: Mapping[str, Any],
    d10_reduction: Mapping[str, Any],
    normalization_no_go: Mapping[str, Any],
    charged_budget: Mapping[str, Any],
    source_certificate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    claim_label, certificate_checks = classify_source_certificate(source_certificate)
    certified = claim_label == "trace_lift_certified"
    candidate_factorization = (source_certificate or {}).get("factorization_lemma") or {}
    candidate_constant = (source_certificate or {}).get("uncentered_lift_constant") or {}
    candidate_residual = (source_certificate or {}).get("attachment_identity_residual") or {}

    if certified:
        factorization_lemma = {
            **candidate_factorization,
            "status": "certified",
            "leakage_bound": {
                "interval": ["0", "0"],
                "certified_zero": True,
            },
        }
        uncentered_lift_constant = {
            **candidate_constant,
            "status": "fixed_by_named_source_object",
        }
        attachment_residual = {
            **candidate_residual,
            "name": "N_det(P)",
            "formula": "s_det(P) - sum_e M_e^ch log q_e(P)",
            "interval": ["0", "0"],
            "certified_zero": True,
            "computable": True,
        }
    else:
        factorization_lemma = _representation_sector_factorization()
        uncentered_lift_constant = {
            "status": "absent_from_declared_source_objects",
            "source_object_name": None,
            "required_source_object_class": "normalized_charged_determinant_line_basepoint_trivialization",
            "value": None,
            "interval": None,
            "reference_stage": None,
        }
        attachment_residual = {
            "name": "N_det(P)",
            "formula": "s_det(P) - sum_e M_e^ch log q_e(P)",
            "value": None,
            "interval": [None, None],
            "interval_kind": "full_affine_orbit_unbounded",
            "certified_zero": False,
            "computable": False,
            "blocking_reasons": [
                "no numeric sector-isolated M_ch",
                "q_psi is source-open, arrow-indexed, and lacks a certified e/mu/tau map",
                "no theorem-grade numeric or interval D10 matter landing to s_det(P)",
                "the D9 zero-leakage channel lemma does not attach the D10 package or normalize its line",
                "no D10-fixed reference stage and basepoint normalization",
            ],
        }

    ratio_regression = _ratio_regression(exact_readout)
    conditional_rows = _conditional_mass_rows(exact_readout, source_certificate) if certified else []
    conditional_on_p = (
        certified and ratio_regression["status"] == "passed" and len(conditional_rows) == 3
    )
    return {
        "artifact": "oph_charged_trace_lift_theorem",
        "generated_utc": _timestamp(),
        "issue": 546,
        "claim_label": claim_label,
        "proof_status": claim_label,
        "source_only": certified,
        "factorization_lemma": factorization_lemma,
        "uncentered_lift_constant": uncentered_lift_constant,
        "attachment_identity_residual": attachment_residual,
        "exact_ratio_regression": ratio_regression,
        "declared_source_audit": _current_source_audit(
            same_label, frontier, d10_reduction, charged_budget
        ),
        "existing_axiom_independence_theorem": _existing_axiom_independence_theorem(),
        "candidate_escape_audit": _candidate_escape_audit(same_label),
        "source_certificate": {
            "artifact_ref": (
                _artifact_ref(DEFAULT_SOURCE_CERTIFICATE) if source_certificate is not None else None
            ),
            "contract": SOURCE_CERTIFICATE_CONTRACT,
            "checks": certificate_checks,
        },
        "strengthened_no_go": {
            "applies": not certified,
            "theorem": (
                "Even granting an exact total D10 determinant and a zero-leakage additive split "
                "s_det = s_det^ch + s_det^rest, the declared source objects do not select the charged "
                "affine origin. For every kappa, s_det^ch -> s_det^ch + 3 kappa, "
                "s_det^rest -> s_det^rest - 3 kappa, and mu -> mu + kappa preserve the total D10 "
                "scalar, central factorization, centered family functional, and q readback, while "
                "N_det -> N_det + 3 kappa. A q value at an unnamed/current arrow does not break this "
                "symmetry; a normalized physical basepoint attachment is required."
            ),
            "preserved_declared_objects": [
                "total D10 determinant scalar, even if granted",
                "zero-leakage additive sector factorization, even if granted",
                "centered charged family functional and all charged ratios",
                "same-label q_psi arrow readback and every formal S_M",
                "physical identity-mode equalizer on each fixed fiber",
            ],
            "changed_object": "N_det -> N_det + 3 kappa",
            "independent_current_blockers": [
                "numeric M_ch is not source-emitted",
                "q_psi is not source-closed or stage-indexed",
                "D10 -> matter s_det(P) is conditional rather than emitted",
                "the exact D9 charged-channel projector is not a D10 determinant attachment",
                "reference-stage determinant-line trivialization is absent",
            ],
            "prior_no_go_artifact": normalization_no_go.get("artifact"),
            "prior_no_go_artifact_ref": _artifact_ref(NORMALIZATION_NO_GO_JSON),
        },
        "missing_source_object_class": {
            "name": "theorem_grade_sector_isolated_normalized_charged_determinant_line_attachment",
            "must_carry": SOURCE_CERTIFICATE_CONTRACT["required_objects"],
            "normalization_equation": "s_det^ch(P,r_0) = sum_psi M_ch[psi] log q_psi(P,r_0)",
            "refinement_equation": "3 mu(r) = sum_psi M_ch[psi] log q_psi(r)",
        },
        "conditional_mass_rows": conditional_rows,
        "promotion": {
            "conditional_on_P_allowed": conditional_on_p,
            "public_promotion_allowed": (
                conditional_on_p
                and (source_certificate or {}).get("anchor_bridge_verdict")
                == "scheme_bridge_certified"
            ),
            "reason": (
                "attachment residual is not certified zero"
                if not certified
                else "trace lift certified; public rows still require the #545 anchor-bridge verdict"
            ),
        },
        "claim_boundary": (
            "No absolute charged-lepton mass is emitted. The current exact e/mu/tau triple remains "
            "a target-anchored ratio checksum. Under the requested no-new-axiom constraint, the "
            "kappa countermodel proves that the normalized attachment and absolute masses cannot "
            "be inferred from the currently declared antecedents."
            if not certified
            else "The trace lift is certified conditionally on P; public promotion remains separately gated."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged sector-isolated trace-lift theorem.")
    parser.add_argument("--exact-readout", default=str(EXACT_READOUT_JSON))
    parser.add_argument("--same-label-readback", default=str(SAME_LABEL_READBACK_JSON))
    parser.add_argument("--determinant-frontier", default=str(DETERMINANT_FRONTIER_JSON))
    parser.add_argument("--d10-reduction", default=str(D10_REDUCTION_JSON))
    parser.add_argument("--normalization-no-go", default=str(NORMALIZATION_NO_GO_JSON))
    parser.add_argument("--charged-budget", default=str(CHARGED_BUDGET_JSON))
    parser.add_argument("--source-certificate", default=str(DEFAULT_SOURCE_CERTIFICATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_path = Path(args.source_certificate)
    source_certificate = _load_json(source_path) if source_path.exists() else None
    artifact = build_artifact(
        _load_json(Path(args.exact_readout)),
        _load_json(Path(args.same_label_readback)),
        _load_json(Path(args.determinant_frontier)),
        _load_json(Path(args.d10_reduction)),
        _load_json(Path(args.normalization_no_go)),
        _load_json(Path(args.charged_budget)),
        source_certificate,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
