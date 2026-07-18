#!/usr/bin/env python3
"""Build the particle-derivation gap ledger.

The ledger separates the new compressed P-trunk from the remaining theorem,
codepath, and execution gaps in the particle program. It is intentionally
static and explicit: changing a status should be a conscious edit, not an
accidental side effect of a numeric rebuild.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_TRUNK = ROOT / "P_derivation" / "runtime" / "p_closure_trunk_current.json"
HIERARCHY_ROOT = PARTICLES_ROOT / "hierarchy"
HIERARCHY_WITNESS = HIERARCHY_ROOT / "computations" / "hierarchy_numeric_witness.json"
HIERARCHY_KRAWCZYK = HIERARCHY_ROOT / "certificates" / "R_U_krawczyk_certificate.json"
HIERARCHY_DAG = HIERARCHY_ROOT / "certificates" / "DAG_U.json"
HIERARCHY_RESONANCE = HIERARCHY_ROOT / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
HIERARCHY_EW_CAPACITY = HIERARCHY_ROOT / "certificates" / "R_EW_global_capacity_certificate.json"
HIERARCHY_NATURALITY = HIERARCHY_ROOT / "issue_332_rg_naturality_certificate.json"
CHARGED_TRACE_LIFT_THEOREM = (
    PARTICLES_ROOT / "runs" / "leptons" / "charged_trace_lift_theorem.json"
)
QUARK_SIGMA_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_sigma_source_nonidentifiability_obstruction.json"
)
QUARK_AXIOM_LEVEL_YUKAWA_OBSTRUCTION = (
    PARTICLES_ROOT
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
)
QUARK_SCHEME_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_running_mass_scheme_convention_obstruction.json"
)
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "particle_derivation_gap_ledger.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "DERIVATION_GAP_LEDGER.md"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_p_trunk_summary() -> dict[str, Any]:
    if not P_TRUNK.exists():
        return {
            "artifact_path": P_TRUNK.relative_to(ROOT).as_posix(),
            "exists": False,
            "claim_status": "not_emitted",
            "may_feed_live_particle_predictions": False,
        }
    payload = json.loads(P_TRUNK.read_text(encoding="utf-8"))
    return {
        "artifact_path": P_TRUNK.relative_to(ROOT).as_posix(),
        "exists": True,
        "claim_status": payload.get("claim_status"),
        "P": payload.get("fixed_point_candidate", {}).get("P"),
        "alpha_inv": payload.get("fixed_point_candidate", {}).get("alpha_inv"),
        "source_report_mode": payload.get("source_report_mode"),
        "may_feed_live_particle_predictions": payload.get("consumer_policy", {}).get(
            "may_feed_live_particle_predictions",
            False,
        ),
    }


def _load_hierarchy_summary() -> dict[str, Any]:
    required = (
        HIERARCHY_WITNESS,
        HIERARCHY_KRAWCZYK,
        HIERARCHY_DAG,
        HIERARCHY_RESONANCE,
        HIERARCHY_EW_CAPACITY,
        HIERARCHY_NATURALITY,
    )
    if not all(path.exists() for path in required):
        return {
            "artifact_path": HIERARCHY_ROOT.relative_to(ROOT).as_posix(),
            "exists": False,
            "claim_status": "not_emitted",
            "may_feed_local_hierarchy_claim": False,
        }

    witness = json.loads(HIERARCHY_WITNESS.read_text(encoding="utf-8"))
    krawczyk = json.loads(HIERARCHY_KRAWCZYK.read_text(encoding="utf-8"))
    dag = json.loads(HIERARCHY_DAG.read_text(encoding="utf-8"))
    resonance = json.loads(HIERARCHY_RESONANCE.read_text(encoding="utf-8"))
    ew_capacity = json.loads(HIERARCHY_EW_CAPACITY.read_text(encoding="utf-8"))
    naturality = json.loads(HIERARCHY_NATURALITY.read_text(encoding="utf-8"))
    public = witness["public_endpoint_branch"]
    source_audit = witness["source_audit_branch"]
    exact_capacity = ew_capacity["exact_capacity_fixed_point"]
    return {
        "artifact_path": HIERARCHY_ROOT.relative_to(ROOT).as_posix(),
        "exists": True,
        "claim_status": "exact_conditional_local_global_hierarchy_and_closed_naturality_certificate",
        "may_feed_local_hierarchy_claim": False,
        "may_feed_conditional_local_hierarchy_claim": True,
        "may_feed_naturality_claim": True,
        "local_global_resonance_status": resonance["status"],
        "full_theorem_grade_resonance_promoted": resonance["full_theorem_grade_resonance_promoted"],
        "work_in_progress_receipts": resonance["work_in_progress_receipts"],
        "ew_capacity_status": ew_capacity["status"],
        "N_CRC_EW": exact_capacity["N_CRC_EW"],
        "bridge_residual": exact_capacity["bridge_residual"],
        "fixed_point_residual_x": exact_capacity["fixed_point_residual_x"],
        "v_over_E_cell_source": exact_capacity["v_over_E_cell_source"],
        "epsilon_H": naturality["epsilon_H"],
        "epsilon_H_interval": naturality["epsilon_H_interval"],
        "public_endpoint_branch": {
            "P": public["P_C"],
            "alpha_U": public["alpha_U_display"],
            "v_over_E_star": public["v_over_E_star"],
            "log10_E_star_over_v": public["log10_E_star_over_v"],
        },
        "source_audit_branch": {
            "P": source_audit["P_cand"],
            "alpha_U": source_audit["alpha_U"],
            "v_over_E_star": source_audit["v_over_E_star"],
        },
        "interval": krawczyk["I_U"],
        "krawczyk_image": krawczyk["K_I"],
        "derivative_interval": krawczyk["derivative_interval"],
        "krawczyk_interior": krawczyk["inclusion"]["K_I_subset_interior_I_U"],
        "dag_forbidden_paths": dag["validation_result"]["forbidden_paths_to_protected_targets"],
        "boundary": (
            "The selected local P -> alpha_U -> v/E_star lane and the Higgs naturality defect "
            "epsilon_H=0 are exact on their declared branches. The local/global resonance bridge "
            "is exact under the screen premises listed in work_in_progress_receipts. Separate "
            "non-promoted gates are the public Thomson endpoint transport, theorem-grade "
            "W/Z promotion, charged-lepton absolute masses, source-only hadron masses, Strong CP, "
            "and the full no-G clock stack for SI gravity."
        ),
    }


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited")


def _is_singleton_zero_interval(value: Any) -> bool:
    if isinstance(value, dict):
        value = value.get("interval", [value.get("lower"), value.get("upper")])
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        return False
    try:
        lower = Decimal(str(value[0]))
        upper = Decimal(str(value[1]))
        return lower.is_finite() and upper.is_finite() and lower == 0 and upper == 0
    except (InvalidOperation, TypeError, ValueError):
        return False


def _artifact_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def _charged_trace_lift_gate() -> dict[str, Any]:
    """Classify #546 without trusting its top-level claim label alone."""

    payload: dict[str, Any] = {}
    if CHARGED_TRACE_LIFT_THEOREM.exists():
        try:
            payload = json.loads(CHARGED_TRACE_LIFT_THEOREM.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            payload = {}

    factorization = payload.get("factorization_lemma") or {}
    leakage = factorization.get("leakage_bound") or {}
    constant = payload.get("uncentered_lift_constant") or {}
    residual = payload.get("attachment_identity_residual") or {}
    source_checks = (payload.get("source_certificate") or {}).get("checks") or {}
    all_source_checks_pass = bool(source_checks) and all(value is True for value in source_checks.values())
    certified = all(
        (
            payload.get("claim_label") == "trace_lift_certified",
            payload.get("source_only") is True,
            factorization.get("status") == "certified",
            leakage.get("certified_zero") is True,
            _is_singleton_zero_interval(leakage),
            bool(constant.get("source_object_name")),
            constant.get("value") is not None,
            residual.get("computable") is True,
            residual.get("certified_zero") is True,
            _is_singleton_zero_interval(residual.get("interval")),
            all_source_checks_pass,
        )
    )

    if certified:
        status = "trace_lift_certified_conditional_on_P"
        boundary = (
            "The #546 artifact certifies zero charged-sector leakage, fixes the uncentered lift "
            "constant from a named source object, and proves the singleton residual interval "
            "N_det(P)=[0,0]. Charged mass rows may enter the conditional-on-P surface; public "
            "promotion remains separately gated by the #545 anchor-bridge verdict."
        )
        next_action = (
            "Propagate the certified conditional-on-P mass intervals and require the #545 "
            "anchor-bridge verdict before public promotion."
        )
    else:
        status = "closed_current_corpus_charged_end_to_end_no_go"
        boundary = (
            "The #546 audit proves exact representation-level isolation of a supplied D9 charged "
            "Yukawa channel, with zero quark leakage, but no D10 determinant attachment. Its "
            "existing-axiom countermodel Y_e -> exp(kappa)Y_e preserves every declared source "
            "antecedent and all ratios while shifting log|det Y_e| by 3 kappa. The corpus also "
            "lacks a numeric source-emitted M_ch, source-closed stage-indexed q, a numeric D10 "
            "matter s_det landing, and a finite normalized reference-stage attachment."
        )
        next_action = (
            "Keep charged masses suppressed. Under the no-new-axiom rule, reopen only if an "
            "declared determinant-sensitive source observable is exhibited that breaks "
            "the kappa countermodel; the gate flips only when N_det has the certified "
            "singleton interval [0,0]."
        )

    return {
        "id": "charged.determinant-normalization-transport",
        "lane": "Charged leptons",
        "status": status,
        "github_issue": 546,
        "closed_issue_refs": [201],
        "title": "Certify or refute the sector-isolated charged trace lift",
        "current_boundary": boundary,
        "next_action": next_action,
        "trace_lift_claim_label": payload.get("claim_label", "artifact_missing_or_invalid"),
        "trace_lift_artifact": _artifact_path(CHARGED_TRACE_LIFT_THEOREM),
        "attachment_identity_residual": residual,
        "target_surfaces": [
            "code/particles/leptons/derive_charged_trace_lift.py",
            "code/particles/runs/leptons/charged_trace_lift_theorem.json",
        ],
    }


def _load_artifact(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _quark_sigma_obstruction_gate() -> dict[str, Any]:
    payload = _load_artifact(QUARK_SIGMA_OBSTRUCTION)
    axiom_payload = _load_artifact(QUARK_AXIOM_LEVEL_YUKAWA_OBSTRUCTION)
    axiom_certified = all(
        (
            axiom_payload.get("proof_status") == "closed_axiom_level_nondefinability_theorem",
            axiom_payload.get("additional_axioms_used") is False,
            (axiom_payload.get("reference_data_policy") or {}).get(
                "no_target_leak_by_construction"
            )
            is True,
            (axiom_payload.get("counterfamily") or {}).get("parameter_space")
            == "(lambda_u,lambda_d) in (R_{>0})^2",
        )
    )
    certified = all(
        (
            payload.get("proof_status") == "closed_exact_current_corpus_obstruction",
            payload.get("theorem_grade_obstruction") is True,
            payload.get("issue_377_acceptance_met_as_obstruction") is True,
            payload.get("issue_379_acceptance_met_as_obstruction") is True,
            payload.get("issue_380_acceptance_met_as_obstruction") is True,
            (payload.get("dependency_audit") or {}).get("no_target_leak") is True,
            (payload.get("exact_ray_classification") or {}).get("fiber") == "(R_{>0})^2",
            axiom_certified,
        )
    )
    return {
        "id": "quark.source-spread-identifiability",
        "lane": "Quarks",
        "status": (
            "closed_current_corpus_two_modulus_nonidentifiability_obstruction"
            if certified
            else "open_source_spread_obstruction_artifact_invalid"
        ),
        "github_issue": 377,
        "related_github_issues": [379, 380],
        "title": "Classify the source-only up/down quark spread fiber",
        "current_boundary": (
            "The target-free source equations fix the two ordered profile rays but leave their up/down positive "
            "spans independent. The compatible fiber is (R_{>0})^2, and its free rescaling action changes the "
            "affine sector means and mass textures. The edge path begins at a hand-written transport template and "
            "does not select either modulus. The stronger Axioms-1--5 theorem constructs physically inequivalent "
            "Yukawa packages with the same MAR score, so the stated MAR order does not break this action."
            if certified
            else "The target-free spread-obstruction artifact is absent or fails its no-target certificate."
        ),
        "next_action": (
            "Keep all numeric quark rows suppressed. With additional axioms excluded, reopen only if a theorem from "
            "the existing MaxEnt/refinement data refutes the equal-MAR-score counterfamily and emits both spreads."
        ),
        "obstruction_artifact": _artifact_path(QUARK_SIGMA_OBSTRUCTION),
        "axiom_level_obstruction_artifact": _artifact_path(QUARK_AXIOM_LEVEL_YUKAWA_OBSTRUCTION),
        "axiom_level_obstruction_certified": axiom_certified,
        "target_surfaces": ["code/particles/flavor", "particle paper quark section"],
    }


def _quark_scheme_obstruction_gate() -> dict[str, Any]:
    payload = _load_artifact(QUARK_SCHEME_OBSTRUCTION)
    acceptance = payload.get("issue_acceptance") or {}
    certified = all(
        (
            payload.get("proof_status")
            == "closed_structural_finite_renormalization_nonidentifiability_obstruction",
            (acceptance.get("381") or {}).get("acceptance_met_as_sharper_obstruction") is True,
            (acceptance.get("382") or {}).get("acceptance_met_as_sharper_obstruction") is True,
            (payload.get("reference_data_policy") or {}).get("no_target_leak_by_construction") is True,
            (payload.get("stored_matrix_dimensional_audit") or {}).get(
                "certified_physical_yukawa_matrices"
            )
            is False,
        )
    )
    return {
        "id": "quark.running-scheme-and-yukawa-normalization",
        "lane": "Quarks",
        "status": (
            "closed_structural_scheme_nonidentifiability_obstruction"
            if certified
            else "open_scheme_obstruction_artifact_invalid"
        ),
        "github_issue": 381,
        "related_github_issues": [382],
        "title": "Separate source physics from running-mass and Yukawa coordinate conventions",
        "current_boundary": (
            "Finite renormalizations and scale changes preserve physical amplitudes while changing running-mass "
            "coordinates. Declaring an external comparison chart is allowed, but the current source signature also "
            "fails to emit the underlying RGI mass vector or six-flavor Yukawa trajectory. "
            "The stored six-row packet mixes light MSbar coordinates at 2 GeV, charm and bottom self-scale "
            "coordinates, and a separate top pole extraction coordinate. Its GeV-valued matrices are mass textures, "
            "not dimensionless physical Yukawa matrices."
            if certified
            else "The scheme-obstruction artifact is absent or fails its value-free acceptance checks."
        ),
        "next_action": (
            "Emit an RG-covariant mass trajectory or invariant. Apply a declared comparison chart afterward, and "
            "require one common scale, threshold transport, top conversion, running v(mu), and y=sqrt(2)m/v before "
            "using the phrase physical Yukawa matrix."
        ),
        "obstruction_artifact": _artifact_path(QUARK_SCHEME_OBSTRUCTION),
        "target_surfaces": ["code/particles/flavor", "particle paper quark section"],
    }


def build_gap_rows() -> list[dict[str, Any]]:
    return [
        {
            "id": "pclosure.compressed-trunk-artifact",
            "lane": "P closure",
            "status": "closed_canonical_guarded_candidate_artifact",
            "github_issue": 224,
            "closed_issue_refs": [224],
            "title": "Make the five-equation P trunk the canonical audit artifact",
            "current_boundary": (
                "The compressed trunk is the canonical guarded audit and compare artifact. It remains "
                "outside the certified particle-root surface until the source spectral measure payload "
                "and interval-certificate gates close."
            ),
            "next_action": "Keep emitting p_closure_trunk_current.json and block live prediction promotion.",
            "target_surfaces": ["code/P_derivation", "code/particles"],
        },
        {
            "id": "d10.ward-projected-thomson-endpoint",
            "lane": "D10 electromagnetic endpoint",
            "status": "closed_blocker_isolated_endpoint_package",
            "github_issue": 223,
            "successor_github_issue": 235,
            "title": "Close the Ward-projected U(1)_Q endpoint package",
            "current_boundary": (
                "The structured-running law is a continuation. The endpoint package computes the "
                "residual inverse-alpha transport packet and isolates the first non-internalized "
                "object. The source-residual non-identifiability boundary is closed in issue #235; "
                "the source-emitted QCD spectral map, scheme remainder, and interval error control "
                "are stage-gated."
            ),
            "next_action": (
                "Keep the package as the closed blocker-isolation artifact for issue #223."
            ),
            "target_surfaces": [
                "code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md",
                "code/P_derivation/runtime/thomson_endpoint_package_current.json",
                "code/particles/calibration",
            ],
        },
        {
            "id": "d10.source-residual-map-and-interval-certificate",
            "lane": "D10 electromagnetic endpoint",
            "status": "closed_blocker_isolated_source_residual_no_go",
            "github_issue": 235,
            "closed_issue_refs": [223],
            "title": "Populate the source spectral measure payload and interval certificate",
            "current_boundary": (
                "The endpoint package fixes the target residual and the current corpus proves "
                "non-identifiability of R_Q(P) from the existing D10 invariant packet. No OPH "
                "source payload emits the Ward-projected hadronic spectral measure or same-scheme "
                "electroweak remainder. The source-spectral reduction theorem is emitted, and the "
                "screening-invariant no-go rejects fitted c_Q and detuning-only shortcuts."
            ),
            "next_action": (
                "Populate the Ward-projected source spectral measure payload, including rho_had(s;P) "
                "or an equivalent spectral primitive, matching remainder, certified quadrature bounds, "
                "and the interval certificate for the final map."
            ),
            "target_surfaces": [
                "code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md",
                "code/P_derivation/SOURCE_SPECTRAL_THEOREM.md",
                "code/P_derivation/runtime/source_spectral_theorem_current.json",
                "code/P_derivation/runtime/thomson_endpoint_contract_current.json",
                "code/P_derivation/runtime/screening_invariant_no_go_current.json",
                "code/P_derivation/runtime/fine_structure_interval_certificate_current.json",
                "code/P_derivation/runtime/r_q_residual_contract_current.json",
                "code/particles/hadron/ward_projected_spectral_measure.schema.json",
            ],
        },
        {
            "id": "d10.rg-matching-threshold-scheme",
            "lane": "D10 running and matching",
            "status": "closed_declared_convention_contract",
            "github_issue": 32,
            "title": "Internalize RG matching, threshold placement, and scheme conversion",
            "current_boundary": (
                "The D10 branch uses declared running/matching conventions. Issue #32 is closed as "
                "a declared-convention contract, not as an OPH derivation of every coefficient, "
                "threshold, and conversion."
            ),
            "next_action": (
                "Keep the declared-convention status visible in prediction surfaces and require "
                "a separate theorem before treating those conventions as OPH-derived."
            ),
            "target_surfaces": ["paper compact D10 section", "code/P_derivation", "code/particles/calibration"],
        },
        {
            "id": "pclosure.certified-codepath-adoption",
            "lane": "P closure",
            "status": "closed_guarded_codepath_adoption",
            "github_issue": 224,
            "closed_issue_refs": [224],
            "title": "Replace particle P consumers with the certified P root",
            "current_boundary": (
                "Particle status surfaces consume the canonical guarded trunk for audit and compare "
                "surfaces. Live particle-root promotion remains blocked while the trunk remains "
                "candidate metadata."
            ),
            "next_action": (
                "Switch live particle builders only after the source spectral measure payload emits "
                "R_Q(P), the interval certificate proves the full map, and the compressed trunk is "
                "promoted beyond guarded candidate metadata."
            ),
            "target_surfaces": ["code/particles/scripts", "code/particles/runs/status", "WebProjects OPH summaries"],
        },
        _charged_trace_lift_gate(),
        _quark_sigma_obstruction_gate(),
        _quark_scheme_obstruction_gate(),
        {
            "id": "quark.selected-class-vs-global-classification",
            "lane": "Quarks",
            "status": "selected_class_descent_closed_global_classification_no_go",
            "github_issue": 199,
            "related_github_issues": [378],
            "title": "Keep selected-fiber descent distinct from global frame classification",
            "current_boundary": (
                "Representative independence is proved on the selected bridge fiber. It neither selects the two "
                "spread moduli nor certifies physical Yukawa matrices. The stronger class-uniform/global "
                "classification lane is a corpus-limited no-go because no source-emitted ambient public-frame "
                "classifier or quotient-intrinsic spread law exists."
            ),
            "next_action": (
                "Keep selected-fiber descent scoped to representative independence. Global classification may reopen "
                "only for a source-emitted ambient public-frame classifier; numeric mass promotion is governed by the "
                "separate spread and scheme obstructions."
            ),
            "target_surfaces": ["code/particles/flavor", "particle paper quark section"],
        },
        {
            "id": "neutrino.pmns-status-and-absolute-rows",
            "lane": "Neutrinos",
            "status": "rejected_candidate_source_basis_and_kernel_open",
            "github_issue": 513,
            "related_github_issues": [117, 508, 513],
            "title": "Rebuild the neutrino lane from source-closed operator and basis data",
            "current_boundary": (
                "The weighted-cycle point is a target-informed candidate, not a physical PMNS derivation. "
                "Its shared-basis construction cancels the inserted charged basis algebraically; the stored "
                "charged artifact is open and nearly degenerate; the same-label scalar inputs inherit an "
                "underived family kernel and candidate line lift; and the point fails the official NuFIT 6.1 "
                "correlated theta23/delta_CP profile. Absolute-splitting rows are fitted diagnostics."
            ),
            "next_action": (
                "Derive a source-closed neutrino operator, stable physical charged-lepton left basis, and "
                "mass-eigenstate label/order rule without oscillation-target feedback; freeze that construction "
                "before evaluating a later likelihood. Keep all present weighted-cycle, bridge, and exact-adapter "
                "numbers rejected or compare-only."
            ),
            "target_surfaces": ["code/particles/neutrino", "RESULTS_STATUS.md"],
        },
        {
            "id": "qcd.strong-cp-angle",
            "lane": "Strong CP",
            "status": "open_theta_qcd_bar_theta_vanishing_gap",
            "github_issue": 155,
            "title": "Keep the strong-CP branch explicit until the physical invariant is emitted",
            "current_boundary": (
                "The selected-class quark audit wrapper carries target-anchored mass textures on the public class "
                "f_P. The two spread moduli are non-identifiable from the source corpus, and the dimensionful "
                "mixed-scheme matrices are not certified physical Yukawas. "
                "The available corpus does not derive theta_QCD, does not emit the physical anomaly-invariant "
                "bar(theta), and does not prove that the physical strong-CP phase vanishes."
            ),
            "next_action": (
                "Keep strong CP explicit as an open branch. First emit a source-only quark mass matrix at one "
                "declared scale with physical determinant-line phase data. Then fix the topological-angle "
                "contribution and prove that the anomaly-invariant strong-CP phase vanishes."
            ),
            "target_surfaces": ["paper particle discussion", "README.md", "code/particles status surfaces"],
        },
        {
            "id": "calibration.direct-top-bridge",
            "lane": "D11/top codomain",
            "status": "closed_current_corpus_codomain_no_go",
            "github_issue": 207,
            "title": "Bridge the cross-section target-audit coordinate to the auxiliary direct-top PDG row",
            "current_boundary": (
                "The target-audit top coordinate uses the PDG cross-section codomain Q007TP4. The auxiliary "
                "direct-top entry Q007TP is a separate extraction codomain and remains compare-only; "
                "the available corpus emits no extraction-response map or uncertainty-propagation certificate."
            ),
            "next_action": (
                "Keep both top coordinates compare-only. Reopen only "
                "for a concrete source-side extraction-response kernel."
            ),
            "target_surfaces": ["code/particles/calibration", "code/particles/runs/status"],
        },
        {
            "id": "hadron.production-backend-systematics",
            "lane": "Hadrons",
            "status": "source_backend_absent_empirical_policy_emitted",
            "github_issue": 153,
            "related_github_issues": [153, 157],
            "title": "Execute the production hadron backend and publish systematics",
            "current_boundary": (
                "Source-only hadron prediction requires a working production hadron backend. "
                "Issues #153 and #157 are source-backend boundaries because there is no production "
                "hadron backend in the local environment. A credible backend is gated on "
                "OPH hardware such as GLORB/Echosahedron, outside local Python and Chrome workers. "
                "The empirical closure surface uses a separate e+e- payload class."
            ),
            "next_action": (
                "Keep source-only hadron rows suppressed. Use empirical hadron closure rows only "
                "through the documented e+e- spectral payload. Promote source-only hadron rows "
                "only after a working OPH hadron backend emits production hadron output, "
                "Ward-projected spectral data, and systematics."
            ),
            "target_surfaces": ["code/particles/hadron", "code/particles/qcd"],
        },
        _empirical_ee_gate(),
        {
            "id": "d10.repair-tuple-selection",
            "lane": "D10/D11 electroweak masses",
            "status": "conditional_selection_theorem_axioms_named",
            "github_issue": 521,
            "title": "Repair-tuple selection beneath the free quadratic family",
            "current_boundary": (
                "The unconditional theorem leaves the free family tau2 = -c eta^2, "
                "delta_n = d (1 - beta) eta^2. The selection-theorem artifact "
                "(runs/calibration/d10_repair_tuple_selection_theorem.json) proves (c, d) = "
                "(sqrt(3)/2, 3/2) uniquely under the named color-balanced descent axioms and "
                "records the discrimination table: the spread between the two on-disk candidates "
                "is below current experimental resolution in every observable, so the selection "
                "is a theorem task. Conditional H/top/W/Z rows with the selection-and-P envelope "
                "live in runs/calibration/conditional_ew_predictions_current.json."
            ),
            "next_action": (
                "Derive the descent axioms A2 and A3 from the realized carrier package (the "
                "color-singlet projection weight and the coherent neutral color sum), which "
                "closes the selection unconditionally and unblocks the D10 promotion review."
            ),
            "target_surfaces": [
                "code/particles/runs/calibration/d10_repair_tuple_selection_theorem.json",
                "code/particles/runs/calibration/conditional_ew_predictions_current.json",
                "code/particles/calibration",
            ],
        },
        {
            "id": "d10.same-scheme-anchor-bridge",
            "lane": "P closure / D10 electromagnetic endpoint",
            "status": "structure_resolved_reduces_to_source_hadron_backend",
            "github_issue": 545,
            "title": "Source-side electroweak scheme bridge for the anchor a0(P)",
            "current_boundary": (
                "The anchor scheme-bridge analysis (runtime/anchor_scheme_bridge_current.json) "
                "identifies the certified anchor gap as the hadronic and higher-order running "
                "deficit of the OPH one-loop unification anchor: the standard reference value "
                "alpha^-1(m_Z) = 128.939 (5-flavor on-shell) minus the OPH one-loop anchor 128.308 "
                "is 0.631, at the lower edge of the certified interval [0.649, 0.855]. The "
                "anchor-exactness no-go (route B) is false: the anchor is a one-loop value with an "
                "understood deficit, not a source-chain failure. A same-scheme bridge (route A) that "
                "fills the deficit needs the measured hadronic running (empirical class, "
                "carried by the payload) or the OPH source hadronic spectral measure."
            ),
            "next_action": (
                "A source-only anchor bridge reduces to the OPH hadronic spectral measure, blocked on "
                "the hadron backend (#425). No source-only scheme-bridge theorem exists on the current "
                "corpus; #545 stays open as that reduction. The empirical-class bridge is the built "
                "payload/endpoint pair."
            ),
            "target_surfaces": [
                "code/P_derivation/runtime/anchor_scheme_bridge_current.json",
                "code/P_derivation/runtime/empirical_thomson_endpoint_current.json",
                "code/P_derivation/EMPIRICAL_HADRON_SCHEME_BRIDGE.md",
            ],
        },
    ]


def _empirical_ee_gate() -> dict[str, Any]:
    payload = ROOT / "particles" / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"
    endpoint = ROOT / "P_derivation" / "runtime" / "empirical_thomson_endpoint_current.json"
    populated = payload.exists() and endpoint.exists()
    if populated:
        return {
            "id": "hadron.empirical-ee-spectral-closure",
            "lane": "Hadrons",
            "status": "payload_populated_endpoint_evaluated_gap_anchor_localized",
            "github_issue": None,
            "title": "Empirical e+e- -> hadrons payload and Thomson endpoint evaluation",
            "current_boundary": (
                "The empirical dispersion payload and the empirical Thomson endpoint artifact are on "
                "disk. The payload interval for the hadronic transport excludes the value required to "
                "reach the measured endpoint with the frozen source anchor; the certified discrepancy "
                "is the same-scheme anchor gap recorded in the endpoint artifact. Row class stays "
                "oph_plus_empirical_hadron_closure; nothing here is a source-only theorem."
            ),
            "next_action": (
                "Emit the source-side electroweak scheme bridge for a0(P) that produces the certified "
                "anchor-gap interval; refine the payload with experiment-level tables when a finer "
                "compilation is ingested."
            ),
            "target_surfaces": [
                "code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json",
                "code/P_derivation/runtime/empirical_thomson_endpoint_current.json",
                "code/P_derivation/empirical_thomson_endpoint.py",
            ],
        }
    return {
        "id": "hadron.empirical-ee-spectral-closure",
        "lane": "Hadrons",
        "status": "policy_scaffold_emitted_dataset_absent",
        "github_issue": None,
        "title": "Populate the empirical e+e- -> hadrons payload for closure rows",
        "current_boundary": (
            "The empirical output class is declared in docs/HADRON.md. The source registry and schema "
            "exist, while the integrated e+e- spectral dataset and dispersion artifact are absent."
        ),
        "next_action": (
            "Populate oph_empirical_ee_hadronic_spectral_measure from PDG, HEPData, alphaQED, "
            "or an equivalent documented compilation, then feed the empirical Thomson endpoint builder."
        ),
        "target_surfaces": [
            "docs/HADRON.md",
            "code/particles/hadron/empirical_ee_hadrons_sources.yaml",
            "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json",
            "code/P_derivation",
        ],
    }


def build_bundles() -> list[dict[str, Any]]:
    return [
        {
            "id": "electroweak-root-closure-bundle",
            "status": "endpoint_package_closed_source_measure_payload_absent",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.source-residual-map-and-interval-certificate",
                "d10.rg-matching-threshold-scheme",
                "d10.same-scheme-anchor-bridge",
                "pclosure.certified-codepath-adoption",
            ],
            "promotion_question": (
                "Can one source-emitted map Delta_Th(P), with declared matching and interval bounds, "
                "certify the compressed P trunk as the particle root without importing alpha(0)?"
            ),
            "result": (
                "Constructive result. The admissible endpoint object is explicit and the endpoint "
                "package computes the residual inverse-alpha packet. The source-spectral reduction "
                "theorem is emitted. Delta_Th(P) must split into "
                "source lepton transport, a Ward-projected hadronic spectral density rho_had(s;P), "
                "a certified electroweak/scheme remainder, RG/matching certificates, quadrature bounds, "
                "and an interval-level fixed-point certificate. The local implementation targets are "
                "P_derivation/runtime/thomson_endpoint_contract_current.json and "
                "P_derivation/runtime/source_spectral_theorem_current.json."
            ),
        },
        {
            "id": "spectrum-source-bundle",
            "status": "closed_current_corpus_source_boundaries_emitted",
            "gap_ids": [
                "charged.determinant-normalization-transport",
                "quark.source-spread-identifiability",
                "quark.running-scheme-and-yukawa-normalization",
                "quark.selected-class-vs-global-classification",
                "neutrino.pmns-status-and-absolute-rows",
            ],
            "promotion_question": (
                "Is there one OPH excitation dictionary and sector-isolated trace-lift theorem that "
                "explains the charged affine anchor and breaks the quark two-modulus spread action while also deriving "
                "a source-closed neutrino operator, charged basis, and mass-label rule without target fitting?"
            ),
            "result": (
                "No promotion. Charged leptons are closed as a corpus-limited no-go by the end-to-end "
                "impossibility theorem: the same-family witness and conditional algebraic readout remain, "
                "but no theorem-grade A_ch(P) is emitted. The quark source fiber retains two free positive spread "
                "moduli, its running-coordinate chart is conventional, and its GeV matrices are mass textures. "
                "Selected-fiber descent and the global-classification no-go do not change those facts. The "
                "weighted-cycle neutrino point is a "
                "rejected target-informed candidate; no physical PMNS, ordering, Majorana, or absolute-mass "
                "row is emitted."
            ),
        },
        {
            "id": "strong-cp-closure-bundle",
            "status": "open_physical_invariant_gap",
            "gap_ids": [
                "qcd.strong-cp-angle",
            ],
            "promotion_question": (
                "Can a source-only common-scale quark mass matrix be extended to the physical strong-CP invariant, "
                "including theta_QCD, bar(theta), and a vanishing theorem on the realized branch?"
            ),
            "result": (
                "No promotion. The selected-class quark wrapper carries only target-anchored mixed-scheme mass "
                "textures. The source spread pair is non-identifiable and the physical dimensionless Yukawa "
                "normalization is absent. "
                "The available corpus does not emit the determinant-line phase contribution, the bare theta_QCD "
                "term, or a theorem forcing the physical strong-CP phase to vanish."
            ),
        },
        {
            "id": "qcd-thomson-backend-bundle",
            "status": "source_backend_boundary_empirical_policy_emitted",
            "gap_ids": [
                "d10.ward-projected-thomson-endpoint",
                "hadron.production-backend-systematics",
                "hadron.empirical-ee-spectral-closure",
            ],
            "promotion_question": (
                "Can the source-only backend or the empirical e+e- spectral payload supply the "
                "Ward-projected hadronic term needed by the Thomson endpoint with row class recorded?"
            ),
            "result": (
                "Constructive result with two surfaces. The source-only primitive remains "
                "production_ward_projected_hadronic_spectral_measure_export and requires a real OPH "
                "hadron backend. The empirical surface uses a separate e+e- payload class and "
                "cannot promote the source-only theorem."
            ),
        },
        {
            "id": "top-codomain-bridge-bundle",
            "status": "closed_current_corpus_codomain_no_go",
            "gap_ids": [
                "calibration.direct-top-bridge",
            ],
            "promotion_question": (
                "Can the cross-section target-audit top coordinate be mapped into the auxiliary direct-top extraction codomain "
                "without using Q007TP as a calibration input?"
            ),
            "result": (
                "No-go result. Q007TP4 remains a target-audit coordinate. The auxiliary direct-top row "
                "Q007TP is compare-only because the available corpus emits no source-side extraction-response "
                "kernel into that codomain."
            ),
        },
        {
            "id": "particle-root-integration-gate",
            "status": "keep_candidate_with_constructive_next_artifacts",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.rg-matching-threshold-scheme",
                "pclosure.certified-codepath-adoption",
                "charged.determinant-normalization-transport",
                "quark.source-spread-identifiability",
                "quark.running-scheme-and-yukawa-normalization",
                "quark.selected-class-vs-global-classification",
                "neutrino.pmns-status-and-absolute-rows",
                "calibration.direct-top-bridge",
                "hadron.production-backend-systematics",
                "hadron.empirical-ee-spectral-closure",
            ],
            "promotion_question": (
                "Do the returned packets jointly close the endpoint, matching, interval, and source-object "
                "requirements strongly enough to promote the compressed trunk into particle builders?"
            ),
            "result": (
                "No promotion. The first wave emits constructive next artifacts, so the compressed P trunk "
                "remains candidate/audit metadata until those artifacts are populated and certified."
            ),
        },
    ]


def build_ledger() -> dict[str, Any]:
    rows = build_gap_rows()
    return {
        "artifact": "oph_particle_derivation_gap_ledger",
        "generated_utc": _now_utc(),
        "purpose": "Systematic claim-safe queue after the five-equation P-trunk simplification.",
        "p_trunk": _load_p_trunk_summary(),
        "electroweak_hierarchy": _load_hierarchy_summary(),
        "bundles": build_bundles(),
        "rows": rows,
        "promotion_policy": {
            "compressed_p_trunk_is_certified_prediction_root": False,
            "p_trunk_issue_closed_as_guarded_candidate_adoption": True,
            "reason": (
                "The canonical guarded trunk is emitted for audit and compare surfaces. Live "
                "prediction promotion waits on the populated source spectral measure payload, "
                "same-scheme remainder, and interval certificate."
            ),
            "hadron_backend_in_current_local_scope": False,
            "hadron_backend_scope_reason": (
                "Source-only production hadrons require a real OPH hardware backend. Issues #153/#157 "
                "are source-backend boundaries; local surrogate output is non-promoting and empirical "
                "hadron closure uses a separate e+e- payload class."
            ),
            "empirical_hadron_closure_class_declared": True,
            "empirical_hadron_closure_source_only": False,
            "torus_mode_language_allowed_in_pipeline": False,
            "address_remaining_blockers_one_by_one": False,
            "obstruction_only_worker_result_allowed": True,
        },
    }


def render_markdown(ledger: dict[str, Any]) -> str:
    p_trunk = ledger["p_trunk"]
    hierarchy = ledger["electroweak_hierarchy"]
    lines = [
        "# Particle Derivation Gap Ledger",
        "",
        f"Generated: `{ledger['generated_utc']}`",
        "",
        ledger["purpose"],
        "",
        "## P-Trunk Claim Boundary",
        "",
        f"- Artifact: `{p_trunk['artifact_path']}`",
        f"- Exists: `{p_trunk['exists']}`",
        f"- Claim label: `{p_trunk['claim_status']}`",
        f"- May feed promoted particle predictions: `{p_trunk['may_feed_live_particle_predictions']}`",
    ]
    if p_trunk.get("P") is not None:
        lines.extend(
            [
                f"- Candidate P: `{p_trunk['P']}`",
                f"- Candidate alpha^-1: `{p_trunk['alpha_inv']}`",
                f"- Source report mode: `{p_trunk['source_report_mode']}`",
            ]
        )
    lines.extend(
        [
            "",
        "## Electroweak Hierarchy Certificate",
        "",
        f"- Artifact: `{hierarchy['artifact_path']}`",
        f"- Exists: `{hierarchy['exists']}`",
        f"- Claim label: `{hierarchy['claim_status']}`",
            f"- May feed local hierarchy claim: `{hierarchy['may_feed_local_hierarchy_claim']}`",
            f"- May feed Higgs naturality claim: `{hierarchy.get('may_feed_naturality_claim', False)}`",
        ]
    )
    if hierarchy.get("exists"):
        public = hierarchy["public_endpoint_branch"]
        source_audit = hierarchy["source_audit_branch"]
        interval = hierarchy["interval"]
        k_image = hierarchy["krawczyk_image"]
        derivative = hierarchy["derivative_interval"]
        epsilon_interval = hierarchy["epsilon_H_interval"]
        lines.extend(
            [
                f"- Public endpoint P: `{public['P']}`",
                f"- Public endpoint alpha_U: `{public['alpha_U']}`",
                f"- Public endpoint v/E_star: `{public['v_over_E_star']}`",
                f"- Public endpoint log10(E_star/v): `{public['log10_E_star_over_v']}`",
                f"- Source-audit P: `{source_audit['P']}`",
                f"- Source-audit alpha_U: `{source_audit['alpha_U']}`",
                f"- Source-audit v/E_star: `{source_audit['v_over_E_star']}`",
                f"- R_U interval: `[{interval['lower']}, {interval['upper']}]`",
                f"- Krawczyk image: `[{k_image['lower']}, {k_image['upper']}]`",
                f"- Derivative interval: `[{derivative['lower']}, {derivative['upper']}]`",
                f"- Krawczyk interior inclusion: `{hierarchy['krawczyk_interior']}`",
                f"- Forbidden DAG paths into protected targets: `{hierarchy['dag_forbidden_paths']}`",
                f"- Local/global resonance status: `{hierarchy['local_global_resonance_status']}`",
                f"- Full theorem-grade resonance promoted: `{hierarchy['full_theorem_grade_resonance_promoted']}`",
                f"- Work-in-progress receipts: `{hierarchy['work_in_progress_receipts']}`",
                f"- Conditional EW bridge capacity (modulo F and CP-1 to CP-3): "
                f"`{hierarchy['N_CRC_EW']}`",
                f"- Bridge residual: `{hierarchy['bridge_residual']}`",
                f"- Fixed-point residual in log capacity: `{hierarchy['fixed_point_residual_x']}`",
                f"- Source v/E_cell: `{hierarchy['v_over_E_cell_source']}`",
                f"- Higgs naturality defect: `epsilon_H={hierarchy['epsilon_H']}`",
                f"- Higgs naturality interval: `[{epsilon_interval[0]}, {epsilon_interval[1]}]`",
                f"- Boundary: {hierarchy['boundary']}",
                "- Physical capacity comparison remains open until F and CP-1 to "
                "CP-3 close and the joint cosmological posterior is propagated.",
            ]
        )
    lines.extend(
        [
            "",
            "## Bundle Claim Gates",
            "",
            "Claim gates are grouped into coupled closure packets rather than a one-blocker-at-a-time queue.",
            "",
            "| Bundle | Claim label | Gaps | Promotion question |",
            "| --- | --- | --- | --- |",
        ]
    )
    for bundle in ledger["bundles"]:
        gaps = ", ".join(f"`{gap}`" for gap in bundle["gap_ids"])
        lines.append(
            f"| `{bundle['id']}` | `{_display_status(bundle['status'])}` | {gaps} | {bundle['promotion_question']} |"
        )
    lines.extend(
        [
            "",
            "## Bundle Packet Results",
            "",
        ]
    )
    for bundle in ledger["bundles"]:
        lines.append(f"- `{bundle['id']}`: `{_display_status(bundle['status'])}`. {bundle['result']}")
    lines.extend(
        [
            "",
            "## Claim Gates",
            "",
            "| ID | Lane | Claim label | Gate |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in ledger["rows"]:
        lines.append(
            f"| `{row['id']}` | {row['lane']} | `{_display_status(row['status'])}` | {row['next_action']} |"
        )
    lines.extend(
        [
            "",
            "## Claim Policy",
            "",
            "- The compressed P trunk is an audit/candidate artifact until the endpoint and certificate gates close.",
            "- Claim gates are handled as coupled bundles rather than isolated one-off fixes.",
            "- The particle pipeline must keep compare-only, continuation, selected-class, and theorem-grade rows mechanically distinct.",
            "- Golden-ratio torus or resonance language is not a derivation input unless a separate representation-to-spectrum theorem is supplied.",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the particle derivation gap ledger.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ledger = build_ledger()
    json_text = json.dumps(ledger, indent=2, sort_keys=True) + "\n"

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.write_text(render_markdown(ledger) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
