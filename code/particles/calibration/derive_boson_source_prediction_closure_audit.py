#!/usr/bin/env python3
"""Compose the fail-closed W/Z/H source-prediction closure audit.

This module does not derive a new mass formula.  It binds the existing pixel,
D10, D11, hierarchy, RG, adapter, and provenance records into one receipt that
answers the stronger question: does the current repository contain a complete,
single-branch, source-only prediction of the physical W, Z, and Higgs complex
poles?  The default answer is deliberately fail-closed.

The receipt keeps four pixel branches and four mass surfaces distinct, records
the difference between the raw neutral color trace N_c and the normalized
repair-chart coefficient d=N_c/2, and refuses promotion while any absolute-
scale, D10 uniqueness, concrete RG, branch-rigidity, pole, uncertainty, or
full source-ancestry gate remains open.
"""

from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any


CODE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = CODE_ROOT.parent
PARTICLES = CODE_ROOT / "particles"
CAL_RUNS = PARTICLES / "runs" / "calibration"
HIERARCHY = PARTICLES / "hierarchy"
P_RUNTIME = CODE_ROOT / "P_derivation" / "runtime"

LEGACY_BASIS_SOURCE = (
    PARTICLES / "calibration" / "legacy_d10" / "particle_masses_paper_d10_d11.py"
)
QT_CERTIFICATE_SCHEMA = PARTICLES / "calibration" / "d10_ew_quotient_path_certificate.schema.json"
DEFAULT_OUT = CAL_RUNS / "boson_source_prediction_closure_audit.json"
MANIFEST_SUPPORT_PATHS = [
    Path(__file__).resolve(),
    PARTICLES / "calibration" / "boson_source_prediction_closure_audit.schema.json",
    PARTICLES / "calibration" / "derive_d10_ew_target_free_repair_value_law.py",
    PARTICLES / "calibration" / "derive_d10_ew_quotient_transport_receipt.py",
    PARTICLES / "calibration" / "derive_d11_live_exact_split_pair_theorem.py",
    PARTICLES / "calibration" / "derive_color_amplitude_loop_split.py",
    QT_CERTIFICATE_SCHEMA,
]

INPUT_PATHS: dict[str, Path] = {
    "p_closure_trunk": P_RUNTIME / "p_closure_trunk_current.json",
    "measured_endpoint_calibration": P_RUNTIME / "measured_endpoint_calibration_current.json",
    "rg_matching_threshold_contract": P_RUNTIME / "rg_matching_threshold_contract_current.json",
    "d10_observable_family": CAL_RUNS / "d10_ew_observable_family.json",
    "d10_source_transport_pair": CAL_RUNS / "d10_ew_source_transport_pair.json",
    "d10_target_free_repair_candidate": CAL_RUNS / "d10_ew_target_free_repair_value_law.json",
    "d10_quotient_transport_receipt": CAL_RUNS / "d10_ew_quotient_transport_receipt.json",
    "d10_minimal_conditional_theorem": CAL_RUNS / "d10_ew_minimal_conditional_theorem.json",
    "d10_repair_tuple_selection": CAL_RUNS / "d10_repair_tuple_selection_theorem.json",
    "color_amplitude_loop_split": CAL_RUNS / "color_amplitude_loop_split.json",
    "d10_freeze_once_adapter": (
        CAL_RUNS / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
    ),
    "d11_declared_surface": CAL_RUNS / "d11_declared_calibration_surface.json",
    "d11_conditional_split": CAL_RUNS / "d11_live_exact_split_pair_theorem.json",
    "d11_reference_adapter": CAL_RUNS / "d11_reference_exact_adapter.json",
    "hierarchy_numeric_witness": HIERARCHY / "computations" / "hierarchy_numeric_witness.json",
    "hierarchy_declared_dag": HIERARCHY / "certificates" / "DAG_U.json",
    "public_pixel_certificate": HIERARCHY / "certificates" / "R_P_public_pixel_certificate.json",
    "source_audit_pixel_certificate": (
        HIERARCHY / "certificates" / "R_P_source_audit_pixel_certificate.json"
    ),
    "wz_boundary_certificate": HIERARCHY / "certificates" / "R_WZ_boundary_certificate.json",
    "ht_declared_surface_certificate": (
        HIERARCHY / "certificates" / "R_HT_declared_surface_certificate.json"
    ),
    "no_g_clock_certificate": HIERARCHY / "certificates" / "R_gamma_noG_DAG_certificate.json",
}


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _repo_relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _file_receipt(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    return {
        "path": _repo_relative(path),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
    }


def _literal_constants(path: Path, names: set[str]) -> dict[str, Any]:
    """Read literal module constants without importing or executing legacy code."""

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    values: dict[str, Any] = {}
    for node in tree.body:
        name: str | None = None
        value_node: ast.expr | None = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                name = target.id
                value_node = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name = node.target.id
            value_node = node.value
        if name in names and value_node is not None:
            values[name] = ast.literal_eval(value_node)
    missing = names - values.keys()
    if missing:
        raise ValueError(f"missing literal constants in {path}: {sorted(missing)}")
    return values


def load_inputs() -> tuple[
    dict[str, dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
    list[dict[str, Any]],
]:
    payloads = {name: _load_json(path) for name, path in INPUT_PATHS.items()}
    manifest_paths = sorted(
        {*INPUT_PATHS.values(), LEGACY_BASIS_SOURCE, *MANIFEST_SUPPORT_PATHS},
        key=lambda path: _repo_relative(path),
    )
    manifest = sorted((_file_receipt(path) for path in manifest_paths), key=lambda item: item["path"])
    constants = _literal_constants(LEGACY_BASIS_SOURCE, {"P_DEFAULT", "E_PLANCK_GEV"})
    # The reviewed correspondence has been integrated into the repository
    # receipt semantics. It remains intentionally absent from runtime inputs so
    # a standalone clone can rebuild this audit without workspace-only files.
    return payloads, manifest, constants, []


def _gate(
    theorem_id: str,
    title: str,
    *,
    theorem_status: str,
    closure_requirement_status: str,
    gate_passed: bool,
    evidence: list[str],
    open_requirements: list[str],
    gate_role: str = "required",
) -> dict[str, Any]:
    return {
        "theorem_id": theorem_id,
        "title": title,
        "theorem_status": theorem_status,
        "closure_requirement_status": closure_requirement_status,
        "gate_role": gate_role,
        "gate_passed_for_full_source_prediction": gate_passed,
        "evidence": evidence,
        "open_requirements": open_requirements,
    }


def build_artifact(
    payloads: dict[str, dict[str, Any]],
    input_manifest: list[dict[str, Any]],
    legacy_constants: dict[str, Any],
    *,
    correspondence_manifest: list[dict[str, Any]] | None = None,
    generated_utc: str | None = None,
) -> dict[str, Any]:
    p_trunk = payloads["p_closure_trunk"]
    endpoint = payloads["measured_endpoint_calibration"]
    rg_contract = payloads["rg_matching_threshold_contract"]
    observable = payloads["d10_observable_family"]
    source_pair = payloads["d10_source_transport_pair"]
    repair = payloads["d10_target_free_repair_candidate"]
    qt_receipt = payloads["d10_quotient_transport_receipt"]
    minimal = payloads["d10_minimal_conditional_theorem"]
    selection = payloads["d10_repair_tuple_selection"]
    color_split = payloads["color_amplitude_loop_split"]
    freeze_adapter = payloads["d10_freeze_once_adapter"]
    d11_surface = payloads["d11_declared_surface"]
    d11_split = payloads["d11_conditional_split"]
    d11_adapter = payloads["d11_reference_adapter"]
    hierarchy = payloads["hierarchy_numeric_witness"]
    hierarchy_dag = payloads["hierarchy_declared_dag"]
    public_pixel = payloads["public_pixel_certificate"]
    source_pixel = payloads["source_audit_pixel_certificate"]
    wz_boundary = payloads["wz_boundary_certificate"]
    ht_boundary = payloads["ht_declared_surface_certificate"]
    no_g_clock = payloads["no_g_clock_certificate"]

    public_branch = hierarchy["public_endpoint_branch"]
    source_audit_branch = hierarchy["source_audit_branch"]
    core = observable["core_source"]
    trunk_candidate = p_trunk["fixed_point_candidate"]
    endpoint_public = endpoint["calibrated_values"]

    source_branches = {
        "legacy_d10_P_1p63094": {
            "branch_id": "legacy_d10_P_1p63094",
            "P": str(core["p"]),
            "alpha_U": str(core["alpha_u"]),
            "v_over_E_star": None,
            "source_class": "legacy_rounded_calibration_carrier_with_external_dimensionful_scale",
            "source_artifact": _repo_relative(INPUT_PATHS["d10_observable_family"]),
            "eligible_for_source_only_mass_prediction": False,
            "reason": "The legacy branch uses rounded P and the externally normalized E_PLANCK_GEV constant.",
        },
        "public_endpoint_P_C": {
            "branch_id": "public_endpoint_P_C",
            "P": public_branch["P_C"],
            "P_full_precision_endpoint_record": endpoint_public["P_from_outer_equation"],
            "alpha_U": public_branch["alpha_U_display"],
            "v_over_E_star": public_branch["v_over_E_star"],
            "source_class": "public_endpoint_conditioned_branch",
            "source_artifact": _repo_relative(INPUT_PATHS["public_pixel_certificate"]),
            "eligible_for_source_only_mass_prediction": False,
            "reason": public_pixel["status"],
        },
        "source_audit_P_cand": {
            "branch_id": "source_audit_P_cand",
            "P": source_audit_branch["P_cand"],
            "alpha_U": source_audit_branch["alpha_U"],
            "v_over_E_star": source_audit_branch["v_over_E_star"],
            "source_class": "source_audit_witness_not_full_endpoint_proof",
            "source_artifact": _repo_relative(INPUT_PATHS["source_audit_pixel_certificate"]),
            "eligible_for_source_only_mass_prediction": False,
            "reason": source_pixel["status"],
        },
        "compressed_p_trunk_candidate": {
            "branch_id": "compressed_p_trunk_candidate",
            "P": trunk_candidate["P"],
            "alpha_U": None,
            "v_over_E_star": None,
            "source_class": p_trunk["claim_status"],
            "source_artifact": _repo_relative(INPUT_PATHS["p_closure_trunk"]),
            "eligible_for_source_only_mass_prediction": bool(
                p_trunk["consumer_policy"]["may_feed_live_particle_predictions"]
            ),
            "reason": "The compressed candidate trunk explicitly forbids live-prediction consumption.",
        },
    }

    carrier_quintet = source_pair["compact_hypercharge_only_mass_slice"]["coherent_output_quintet"]
    repair_quintet = repair["coherent_emitted_quintet"]
    freeze_quintet = freeze_adapter["coherent_repaired_quintet"]
    split_pair = d11_split["exact_split_pair"]
    boundary_values = wz_boundary["values"]
    mass_surfaces = [
        {
            "surface_id": "d10_selected_current_carrier",
            "branch_id": "legacy_d10_P_1p63094",
            "classification": "selected_carrier_tree_or_running_mass_coordinate",
            "values": {
                "W_GeV": carrier_quintet["MW_pole"],
                "Z_GeV": carrier_quintet["MZ_pole"],
            },
            "legacy_value_keys": ["MW_pole", "MZ_pole"],
            "physical_mass_semantics": "tree_level_algebraic_coordinate_not_certified_complex_pole",
            "source_unique": False,
            "physical_pole_certified": False,
            "promotion_allowed": False,
            "source_artifact": _repo_relative(INPUT_PATHS["d10_source_transport_pair"]),
        },
        {
            "surface_id": "d10_running_tree_repair_candidate",
            "branch_id": "legacy_d10_P_1p63094",
            "classification": "candidate_only_retrospective_selection_not_source_unique",
            "values": {
                "W_GeV": repair_quintet["MW_pole"],
                "Z_GeV": repair_quintet["MZ_pole"],
            },
            "boundary_adapter_alias_values": {
                "W_GeV": boundary_values["M_W_GeV"],
                "Z_GeV": boundary_values["M_Z_GeV"],
            },
            "legacy_value_keys": ["MW_pole", "MZ_pole"],
            "physical_mass_semantics": "tree_level_algebraic_coordinate_not_certified_complex_pole",
            "source_unique": False,
            "physical_pole_certified": False,
            "promotion_allowed": bool(repair["promotion_allowed"]),
            "source_artifact": _repo_relative(INPUT_PATHS["d10_target_free_repair_candidate"]),
        },
        {
            "surface_id": "d10_freeze_once_reference_adapter",
            "branch_id": "legacy_d10_P_1p63094",
            "classification": "inverse_target_adapter_compare_only",
            "values": {
                "W_GeV": freeze_quintet["MW_pole"],
                "Z_GeV": freeze_quintet["MZ_pole"],
            },
            "legacy_value_keys": ["MW_pole", "MZ_pole"],
            "physical_mass_semantics": "reference_frozen_algebraic_coordinate_not_prediction",
            "source_unique": False,
            "physical_pole_certified": False,
            "promotion_allowed": bool(freeze_adapter["prediction_promotion_allowed"]),
            "source_artifact": _repo_relative(INPUT_PATHS["d10_freeze_once_adapter"]),
        },
        {
            "surface_id": "d11_declared_surface_conditional_split",
            "branch_id": "legacy_d10_P_1p63094",
            "classification": "conditional_declared_surface_coordinate",
            "values": {
                "H_GeV": split_pair["mH_gev"],
                "top_companion_GeV": split_pair["mt_pole_gev"],
            },
            "legacy_value_keys": ["mH_gev", "mt_pole_gev"],
            "physical_mass_semantics": "declared_surface_linearized_readout_not_certified_complex_pole",
            "source_unique": False,
            "physical_pole_certified": False,
            "promotion_allowed": bool(d11_split["prediction_promotion_allowed"]),
            "source_artifact": _repo_relative(INPUT_PATHS["d11_conditional_split"]),
        },
    ]

    n_c = int(minimal["realized_color_count"])
    alpha_y = float(minimal["basis"]["alphaY_mz"])
    alpha_2 = float(minimal["basis"]["alpha2_mz"])
    beta_ew = float(minimal["basis"]["beta_EW"])
    eta = float(minimal["basis"]["eta_source"])
    alpha_sum = alpha_2 + alpha_y
    raw_trace_weight = float(n_c)
    normalized_d = raw_trace_weight / 2.0
    delta_n = normalized_d * (1.0 - beta_ew) * eta * eta
    uplift_from_chart = alpha_sum * delta_n
    uplift_from_raw_trace = raw_trace_weight * alpha_y * eta * eta
    normalization_lhs = alpha_sum * (1.0 - beta_ew)
    normalization_rhs = 2.0 * alpha_y
    normalization = {
        "N_c": n_c,
        "charged_channel": {
            "raw_amplitude_multiplicity": math.sqrt(n_c),
            "generator_half_weight": 0.5,
            "normalized_chart_coefficient_c": math.sqrt(n_c) / 2.0,
            "formula": "c=sqrt(N_c)/2",
        },
        "neutral_channel": {
            "raw_hypercharge_trace_weight": raw_trace_weight,
            "normalized_chart_coefficient_d": normalized_d,
            "chart_formula": "delta_n=d*(1-beta_EW)*eta_source^2",
            "physical_uplift_formula": "delta_alphaY_perp=N_c*alphaY_mz*eta_source^2",
            "normalization_identity": "(alpha2_mz+alphaY_mz)*(1-beta_EW)=2*alphaY_mz",
            "normalization_identity_lhs": normalization_lhs,
            "normalization_identity_rhs": normalization_rhs,
            "normalization_identity_residual": normalization_lhs - normalization_rhs,
            "delta_n_value": delta_n,
            "uplift_from_chart": uplift_from_chart,
            "uplift_from_raw_trace": uplift_from_raw_trace,
            "uplift_identity_residual": uplift_from_chart - uplift_from_raw_trace,
        },
        "conditionality": {
            "applies_to_surface": "historical_color_balanced_quadratic_descent_only",
            "does_not_derive_surface": "d10_running_tree_repair_candidate",
            "arithmetic_normalization_closed": True,
            "transport_channel_identification_closed": False,
            "source_only_theorem_emitted": bool(
                color_split["guards"]["source_only_theorem_emitted"]
            ),
            "promotion_allowed": False,
        },
    }

    # The dedicated QT receipt owns the canonical-helper evaluation, fibre
    # identities, and Jacobian checks.  This audit only composes its status with
    # the wider scale/RG/pole/DAG gates and binds the reviewed correspondence.
    d10_quotient_transport_review = dict(qt_receipt)
    d10_quotient_transport_review["reviewed_material"] = correspondence_manifest or []

    legacy_basis_path = _repo_relative(LEGACY_BASIS_SOURCE)
    absolute_scale = {
        "status": "open_independent_dimensionful_scale_theorem",
        "source_closed": False,
        "legacy_dimensionful_constant": {
            "name": "E_PLANCK_GEV",
            "value": legacy_constants["E_PLANCK_GEV"],
            "source_class": "externally_normalized_legacy_constant",
            "source_path": legacy_basis_path,
        },
        "legacy_P_default": legacy_constants["P_DEFAULT"],
        "no_g_clock_status": no_g_clock["status"],
        "missing_clock_components": [
            key
            for key, status in no_g_clock["component_status"].items()
            if "not supplied" in status
        ],
        "promotion_allowed": False,
    }

    constructive_ids = [item["id"] for item in rg_contract["constructive_objects"]]
    concrete_rg_receipt_present = all(
        key in rg_contract
        for key in ("scheme_lock", "threshold_map", "beta_provenance_table", "matching_interval_composition_certificate")
    )
    rg_matching = {
        "status": rg_contract["status"],
        "declared_convention_contract_present": True,
        "concrete_transport_receipt_present": concrete_rg_receipt_present,
        "required_objects": constructive_ids,
        "loop_order_thresholds_matching_and_scheme_fully_frozen": False,
        "interval_composition_certificate_present": False,
        "promotion_allowed": False,
        "source_artifact": _repo_relative(INPUT_PATHS["rg_matching_threshold_contract"]),
    }

    physical_pole = {
        "status": "open_no_wzh_complex_pole_certificate",
        "required_mass_convention": "s_B=(M_B-i*Gamma_B/2)^2",
        "analytic_continuation_and_riemann_sheet_fixed": False,
        "self_energy_functions_present": False,
        "contour_uniqueness_or_fixed_order_pole_certificate_present": False,
        "widths_present": False,
        "truncation_refinement_uncertainty_bounds_present": False,
        "BRST_complete_mixing_blocks_present": False,
        "Ward_Slavnov_Taylor_and_Nielsen_receipts_present": False,
        "physical_residue_certificate_present": False,
        "analytic_field_redefinition_invariance_checked": False,
        "legacy_fields_named_pole_are_physical_pole_certificates": False,
        "promotion_allowed": False,
        "required_objects": [
            "BRST_complete_physical_inverse_two_point_blocks_Gamma_W_Gamma_Z_Gamma_H",
            "Ward_Slavnov_Taylor_projectors_and_Nielsen_identities",
            "complex_determinant_zero_on_declared_Riemann_sheets",
            "uniqueness_and_stability_certificate",
            "nonzero_physical_residues_widths_and_uncertainty_budget",
            "analytic_continuation_and_physical_riemann_sheet_selection",
        ],
    }

    dag_nodes = [
        {"id": "measured_thomson_endpoint", "class": "external_measurement"},
        {"id": "public_endpoint_P_C", "class": "endpoint_conditioned_branch"},
        {"id": "source_audit_P_cand", "class": "source_audit_branch"},
        {"id": "compressed_p_trunk_candidate", "class": "candidate_branch"},
        {"id": "legacy_d10_P_1p63094", "class": "legacy_calibration_branch"},
        {"id": "legacy_E_PLANCK_GEV", "class": "external_dimensionful_constant"},
        {"id": "d10_observable_family", "class": "declared_running_surface"},
        {"id": "d10_source_transport_pair", "class": "carrier_builder"},
        {"id": "d10_selected_current_carrier", "class": "conditional_mass_surface"},
        {"id": "d10_running_tree_repair_candidate", "class": "candidate_mass_surface"},
        {"id": "measured_WZ_targets", "class": "external_measurement"},
        {"id": "d10_freeze_once_reference_adapter", "class": "inverse_adapter"},
        {"id": "d11_declared_surface", "class": "declared_constants_and_jacobian"},
        {"id": "d11_declared_surface_conditional_split", "class": "conditional_mass_surface"},
        {"id": "measured_HT_targets", "class": "external_measurement"},
        {"id": "d11_reference_adapter", "class": "inverse_adapter"},
    ]
    dag_edges = [
        {"from": "measured_thomson_endpoint", "to": "public_endpoint_P_C"},
        {"from": "legacy_d10_P_1p63094", "to": "d10_observable_family"},
        {"from": "legacy_E_PLANCK_GEV", "to": "d10_observable_family"},
        {"from": "d10_observable_family", "to": "d10_source_transport_pair"},
        {"from": "d10_source_transport_pair", "to": "d10_selected_current_carrier"},
        {"from": "d10_source_transport_pair", "to": "d10_running_tree_repair_candidate"},
        {"from": "measured_WZ_targets", "to": "d10_freeze_once_reference_adapter"},
        {"from": "d10_running_tree_repair_candidate", "to": "d11_declared_surface_conditional_split"},
        {"from": "d11_declared_surface", "to": "d11_declared_surface_conditional_split"},
        {"from": "measured_HT_targets", "to": "d11_reference_adapter"},
    ]
    source_dag = {
        "status": "declared_artifact_graph_acyclic_full_runtime_and_human_ancestry_open",
        "nodes": dag_nodes,
        "edges": dag_edges,
        "protected_outputs": [
            "d10_selected_current_carrier",
            "d10_running_tree_repair_candidate",
            "d11_declared_surface_conditional_split",
        ],
        "forbidden_target_sources": [
            "measured_thomson_endpoint",
            "measured_WZ_targets",
            "measured_HT_targets",
        ],
        "declared_hierarchy_dag_status": hierarchy_dag["status"],
        "declared_hierarchy_dag_pass": bool(hierarchy_dag["validation_result"]["pass"]),
        "single_end_to_end_source_branch_selected": False,
        "mixed_branch_consumption_allowed": False,
        "full_file_level_runtime_ancestry_certified": False,
        "human_formula_selection_ancestry_certified": False,
        "prospective_pre_reference_manifest_present": False,
        "promotion_allowed": False,
    }

    paths = {name: _repo_relative(path) for name, path in INPUT_PATHS.items()}
    theorem_gates = [
        _gate(
            "T01_structural_electroweak_underdetermination",
            "Structural electroweak underdetermination",
            theorem_status="closed_generic_no_go",
            closure_requirement_status="satisfied_as_claim_boundary",
            gate_passed=True,
            evidence=["continuous (g2,gY,mu2,lambdaH) counterfamily encoded by this audit"],
            open_requirements=[],
        ),
        _gate(
            "T02_absolute_scale_no_go",
            "Independent absolute-scale theorem",
            theorem_status="closed_generic_no_go",
            closure_requirement_status="open_independent_scale_or_clock",
            gate_passed=False,
            evidence=[legacy_basis_path, paths["no_g_clock_certificate"]],
            open_requirements=["source-normalized E_star or independent clock/length/energy theorem"],
        ),
        _gate(
            "T03_electroweak_running_mass_readout",
            "Electroweak running-mass chart",
            theorem_status="closed_algebraic_chart",
            closure_requirement_status="satisfied_for_tree_or_running_coordinates_only",
            gate_passed=True,
            evidence=[paths["d10_source_transport_pair"]],
            open_requirements=["physical-pole promotion remains T10"],
        ),
        _gate(
            "T04_inverse_adapters_are_not_predictions",
            "Inverse adapters are comparison surfaces",
            theorem_status="closed_classification_theorem",
            closure_requirement_status="satisfied_fail_closed",
            gate_passed=True,
            evidence=[paths["d10_freeze_once_adapter"], paths["d11_reference_adapter"]],
            open_requirements=[],
        ),
        _gate(
            "T05_d10_source_uniqueness",
            "D10 source-uniqueness criterion",
            theorem_status="exact_conditional_QT1_QT5_implication_plus_closed_underdetermination_boundary",
            closure_requirement_status="open_QT1_QT5_source_entailment_and_same_branch_packet",
            gate_passed=False,
            evidence=[
                paths["d10_quotient_transport_receipt"],
                paths["d10_minimal_conditional_theorem"],
                paths["d10_repair_tuple_selection"],
            ],
            open_requirements=[
                "emit QT3 to fix the displayed quotient-path coefficients",
                "emit QT5 to exclude additional output-changing deformations",
                "emit the remaining QT1, QT2, and QT4 finite-carrier premises",
                "bind the theorem to one strict P source branch",
            ],
        ),
        _gate(
            "T06_conditional_color_weights",
            "Conditional color-amplitude and color-trace weights",
            theorem_status="conditional_arithmetic_closed_channel_identification_open",
            closure_requirement_status="open_and_distinct_from_complete_QT_value_law",
            gate_passed=False,
            evidence=[paths["color_amplitude_loop_split"], paths["d10_repair_tuple_selection"]],
            open_requirements=["prove amplitude/trace channel factorization and scalar-channel exhaustion"],
            gate_role="alternative_model",
        ),
        _gate(
            "T07_rg_threshold_matching_uniqueness",
            "Frozen RG, threshold, matching, and scheme transport",
            theorem_status="closed_generic_ode_theorem",
            closure_requirement_status="declared_contract_only_concrete_receipt_open",
            gate_passed=False,
            evidence=[paths["rg_matching_threshold_contract"]],
            open_requirements=constructive_ids,
        ),
        _gate(
            "T08_conditional_d11_higgs_closure",
            "Conditional D11 Higgs/top split",
            theorem_status=d11_split["proof_status"],
            closure_requirement_status="conditional_declared_surface_upstream_open",
            gate_passed=False,
            evidence=[paths["d11_declared_surface"], paths["d11_conditional_split"], paths["ht_declared_surface_certificate"]],
            open_requirements=[
                "promotable D10 source repair",
                "DS1--DS5 finite D11 split-character and rigidity certificate",
                "source provenance for D11 core/Jacobian",
                "BRST-complete pole map",
            ],
        ),
        _gate(
            "T09_branch_rigidity_mar",
            "Branch rigidity or target-independent MAR selection",
            theorem_status="generic_criterion_recorded",
            closure_requirement_status="open_admissible_class_enumeration_and_unique_argmin",
            gate_passed=False,
            evidence=[paths["d10_quotient_transport_receipt"], paths["d10_repair_tuple_selection"]],
            open_requirements=["enumerate admissible D10/D11 deformations", "prove unique target-free selector"],
        ),
        _gate(
            "T10_complex_pole_promotion",
            "Physical complex-pole promotion and stability",
            theorem_status="generic_rouche_criterion_recorded",
            closure_requirement_status="open_no_self_energy_or_pole_certificate",
            gate_passed=False,
            evidence=[],
            open_requirements=physical_pole["required_objects"],
        ),
        _gate(
            "T11_deterministic_source_separation",
            "Deterministic source DAG and prospective provenance",
            theorem_status="declared_partial_dag_only",
            closure_requirement_status="open_full_chain_and_human_selection_ancestry",
            gate_passed=False,
            evidence=[paths["hierarchy_declared_dag"]],
            open_requirements=["full file-level DAG", "branch-coherent root", "pre-reference manifest/signature"],
        ),
        _gate(
            "T12_full_wzh_source_prediction",
            "Capstone W/Z/H source prediction",
            theorem_status="not_closed",
            closure_requirement_status="open",
            gate_passed=False,
            evidence=[paths["wz_boundary_certificate"], paths["ht_declared_surface_certificate"]],
            open_requirements=[
                "T02",
                "T05 via one selected D10 model (QT path law or color-balanced alternative)",
                "T07",
                "T08",
                "T09",
                "T10",
                "T11",
            ],
        ),
    ]

    failed_gate_ids = [
        gate["theorem_id"]
        for gate in theorem_gates
        if gate["gate_role"] == "required"
        and not gate["gate_passed_for_full_source_prediction"]
    ]
    promotion = {
        "overall_status": "source_incomplete_no_end_to_end_wzh_pole_prediction",
        "prediction_promotion_allowed": False,
        "all_required_gates_passed": False,
        "failed_gate_ids": failed_gate_ids,
        "shared_blockers": [
            "no_single_promotable_P_branch_feeds_the_current_D10_D11_mass_chain",
            "independent_dimensionful_scale_or_clock_theorem_open",
            "D10_source_model_selection_open_QT_path_law_or_alternative",
            "concrete_RG_threshold_matching_scheme_receipt_open",
            "D11_DS1_DS5_split_character_core_Jacobian_provenance_and_rigidity_open",
            "complex_pole_and_uncertainty_certificate_open",
            "full_target_free_dependency_DAG_and_prospective_manifest_open",
        ],
        "alternative_model_policy": (
            "Close one source-emitted D10 repair model: either the complete QT1--QT5 quotient-path "
            "law or a separately derived color-balanced alternative. The two are alternatives, "
            "not cumulative premises."
        ),
        "strongest_clean_numeric_statement": "dimensionless_electroweak_hierarchy_v_over_E_star_on_declared_branch",
        "wz_status": "source_incomplete",
        "higgs_status": "conditional_declared_surface_candidate",
    }

    return {
        "artifact": "oph_boson_source_prediction_closure_audit",
        "schema_version": "1.0.0",
        "generated_utc": generated_utc or _now_utc(),
        "claim": {
            "observables": ["W", "Z", "H"],
            "required_mass_convention": "complex_propagator_pole",
            "overall_status": promotion["overall_status"],
            "prediction_promotion_allowed": False,
        },
        "source_branches": source_branches,
        "branch_consistency": {
            "required_rule": "one branch_id must be used from P through the physical pole leaves",
            "single_end_to_end_branch_selected": False,
            "mixed_branch_consumption_allowed": False,
            "branch_count": len(source_branches),
            "all_branch_ids_distinct": len({item["branch_id"] for item in source_branches.values()})
            == len(source_branches),
        },
        "mass_surfaces": mass_surfaces,
        "repair_coefficient_normalization": normalization,
        "d10_quotient_transport_review": d10_quotient_transport_review,
        "absolute_scale_contract": absolute_scale,
        "rg_matching_contract": rg_matching,
        "physical_pole_contract": physical_pole,
        "source_dependency_dag": source_dag,
        "theorem_gates": theorem_gates,
        "promotion_decision": promotion,
        "comparison_adapters": {
            "d11_reference_adapter_status": d11_adapter["proof_status"],
            "d11_reference_adapter_promotable": bool(d11_adapter["promotable"]),
            "wz_boundary_status": wz_boundary["status"],
            "ht_boundary_status": ht_boundary["status"],
            "selection_status": selection["status"],
            "d11_declared_surface_promotion_allowed": bool(d11_surface["predictive_promotion_allowed"]),
        },
        "input_manifest": input_manifest,
        "reviewed_correspondence_manifest": correspondence_manifest or [],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the fail-closed W/Z/H source-closure audit.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    payloads, manifest, constants, correspondence_manifest = load_inputs()
    artifact = build_artifact(
        payloads,
        manifest,
        constants,
        correspondence_manifest=correspondence_manifest,
    )
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered, encoding="utf-8")
    if args.print_json:
        print(rendered, end="")
    else:
        print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
