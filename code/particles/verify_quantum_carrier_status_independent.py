#!/usr/bin/env python3
"""Independent exact reconstruction of the issue-552 carrier status packet.

The verifier imports no producer code. It reconstructs the projection from the
full classical parents, validates the bounded source-frontier artifacts, and
rebuilds every receipt field. Unknown fields and altered prose fail equality.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


REPO_ROOT = Path(__file__).resolve().parents[2]
CLASSICAL_PATH = REPO_ROOT / "code/particles/runs/status/carrier_mode_acceptance.json"
COLOR_PATH = REPO_ROOT / "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
LOCAL_DOMAIN_PATH = REPO_ROOT / "code/a5_closure/manifests/stage4_receipt.json"
EW_FRONTIER_PATH = (
    REPO_ROOT
    / "code/particles/hierarchy/higgs_yukawa_source_frontier/outputs"
    / "higgs_yukawa_source_frontier.json"
)
EINSTEIN_PATH = REPO_ROOT / "code/geometry/runs/realized_branch_receipt_report.json"
QUANTUM_EFT_PATH = (
    REPO_ROOT
    / "code/particles/calibration/wz_native_source_packet/outputs"
    / "source_parent_inventory.json"
)
QCD_RESOURCE_PATH = (
    REPO_ROOT
    / "code/particles/runs/qcd/hadron_source_backend/qcd_ensemble"
    / "solver_on_standby.json"
)
PROJECTION_PATH = (
    REPO_ROOT / "code/particles/manifests/quantum_carrier_status_source_projection.json"
)
RECEIPT_PATH = REPO_ROOT / "code/particles/runs/status/quantum_carrier_status.json"

ORDER = ("photon", "gluon", "graviton")
DIMENSION = 4
QUANTUM_EVIDENCE = (
    "quantization_constructed_from_oph",
    "positive_physical_hilbert_space",
    "physical_two_point_or_hamiltonian_spectrum_constructed",
    "positive_residue_massless_pole",
    "asymptotic_or_deconfined_particle_state",
)

CLASSICAL_ROWS = {
    "photon": {
        "label": "Electromagnetic carrier",
        "branch": "unbroken_deconfined_Maxwell_branch_with_Z_A_positive",
        "physical_degrees_of_freedom": (
            "two transverse classical polarizations in four dimensions"
        ),
        "multiplicity_factor": 1,
        "multiplicity_role": "u1_lie_algebra_generator",
    },
    "gluon": {
        "label": "Color gauge carrier",
        "branch": (
            "pure_Yang-Mills_quadratic_expansion_before_nonperturbative_confinement"
        ),
        "physical_degrees_of_freedom": (
            "two perturbative transverse modes per Lie-algebra generator"
        ),
        "multiplicity_factor": 8,
        "multiplicity_role": "su3_adjoint_generator",
    },
    "graviton": {
        "label": "Einstein tensor carrier",
        "branch": (
            "pure_Einstein-Hilbert_linearization_about_a_suitable_Ricci-flat_background"
        ),
        "physical_degrees_of_freedom": (
            "two transverse-traceless classical tensor polarizations in four dimensions"
        ),
        "multiplicity_factor": 1,
        "multiplicity_role": "symmetric_metric_tensor_field",
    },
}

ROW_CONFIG = {
    "photon": {
        "verdict": "NOT_EVALUABLE_NO_SOURCE_SELECTED_MAXWELL_QUANTUM_SECTOR",
        "blocking_frontier": [
            "source_selected_unbroken_u1_quantum_maxwell_sector",
            "finite_source_to_lorentzian_quantum_eft_construction",
        ],
        "state_space": "a positive physical photon state space after the gauge quotient",
        "observable_algebra": (
            "a source-derived gauge-invariant electromagnetic observable algebra"
        ),
        "classical_observable_algebra": (
            "gauge-invariant classical Maxwell field-strength observables"
        ),
        "gauge_quotient": (
            "a quantum BRST or equivalent physical quotient with positivity"
        ),
        "vacuum": (
            "a source-selected unbroken U(1) quantum vacuum with a deconfined photon sector"
        ),
        "spectral_object": (
            "a gauge-invariant physical two-point or Hamiltonian spectral object"
        ),
        "residue": (
            "a nonzero positive residue against a physical electromagnetic current"
        ),
        "phase_status": "DECLARED_CLASSICAL_BRANCH_ONLY",
        "open_interfaces": [
            {
                "issue": 630,
                "requirement": (
                    "source selection of the electroweak action, unbroken electromagnetic phase, and physical current"
                ),
            },
            {
                "issue": 635,
                "requirement": "finite-source to Lorentzian quantum-EFT construction",
            },
        ],
        "resource_boundaries": [],
        "strongest_supported_statement": (
            "On the declared unbroken deconfined Maxwell branch, the hard quadratic mass parameter is zero and the classical reduced system has two transverse modes. This is not a photon pole or rest-mass prediction."
        ),
    },
    "gluon": {
        "verdict": "NOT_EVALUABLE_NO_QCD",
        "blocking_frontier": [
            "finite_source_to_lorentzian_quantum_eft_construction",
            "source_derived_qcd_physical_spectral_sector",
        ],
        "state_space": (
            "a positive gauge-invariant QCD physical state space in the selected phase"
        ),
        "observable_algebra": "a source-derived gauge-invariant QCD observable algebra",
        "classical_observable_algebra": (
            "classical gauge-invariant Yang-Mills polynomial observables"
        ),
        "gauge_quotient": "a nonperturbative physical gauge quotient with positivity",
        "vacuum": "a source-selected QCD vacuum and phase",
        "spectral_object": (
            "a gauge-invariant QCD Hamiltonian or two-point spectral object"
        ),
        "residue": (
            "a nonzero physical residue in an admitted deconfined carrier channel"
        ),
        "phase_status": "NONPERTURBATIVE_QCD_SECTOR_ABSENT_FROM_DECLARED_CORPUS",
        "open_interfaces": [
            {
                "issue": 635,
                "requirement": "finite-source to Lorentzian quantum-EFT construction",
            }
        ],
        "resource_boundaries": [
            {
                "issue": 425,
                "requirement": (
                    "source-only QCD physical spectral production beyond the recorded resource boundary"
                ),
            },
            {
                "issue": 294,
                "requirement": (
                    "continuum reconstruction and nontriviality beyond the recorded bounded exit"
                ),
            },
        ],
        "strongest_supported_statement": (
            "On the declared perturbative Yang-Mills quadratic branch, each of the eight color generators has two transverse modes and zero hard quadratic mass parameter, giving sixteen conditional classical modes. This supplies neither an asymptotic colored gluon nor a QCD spectral prediction."
        ),
    },
    "graviton": {
        "verdict": "NOT_EVALUABLE_NO_INHABITED_EINSTEIN_QUANTUM_CARRIER",
        "blocking_frontier": [
            "inhabited_source_derived_einstein_tower",
            "finite_source_to_lorentzian_linearized_quantum_carrier",
        ],
        "state_space": (
            "a positive physical state space for linearized Einstein perturbations"
        ),
        "observable_algebra": (
            "a source-derived linearized diffeomorphism-invariant observable algebra"
        ),
        "classical_observable_algebra": (
            "linearized curvature observables on the declared background"
        ),
        "gauge_quotient": (
            "a quantum quotient of linearized diffeomorphisms with a positive physical sector"
        ),
        "vacuum": (
            "an inhabited source-derived Einstein tower with a selected quantum background or vacuum"
        ),
        "spectral_object": (
            "a physical linearized spin-two Hamiltonian or two-point spectral object"
        ),
        "residue": (
            "a nonzero positive spin-two residue against a physical conserved stress tensor"
        ),
        "phase_status": "DECLARED_CLASSICAL_BACKGROUND_ONLY",
        "open_interfaces": [
            {
                "issue": 503,
                "requirement": "one inhabited source-derived typed Einstein tower",
            },
            {
                "issue": 635,
                "requirement": "finite-source to Lorentzian quantum-EFT construction",
            },
        ],
        "resource_boundaries": [],
        "strongest_supported_statement": (
            "On the declared pure Einstein-Hilbert linearization about a suitable Ricci-flat background, the hard quadratic mass parameter is zero and the transverse-traceless classical system has two modes. This is not a graviton Hilbert-space, pole, or rest-mass prediction."
        ),
    },
}


class VerificationError(ValueError):
    pass


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise VerificationError(f"object required at {path}")
    return value


def _digest(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _pin(path: Path, role: str) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "bytes": len(raw),
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "role": role,
        "sha256": _digest(raw),
    }


def _receipt_hash(receipt: Mapping[str, Any]) -> str:
    core = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    raw = json.dumps(
        core,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return _digest(raw)


def _boundary_evidence() -> list[dict[str, Any]]:
    local_domain = _load(LOCAL_DOMAIN_PATH)
    if (
        local_domain.get("schema") != "oph.local-domain-stage4.v1"
        or local_domain.get("issue") != 634
        or local_domain.get("verdict") != "ATTAINED"
        or local_domain.get("physical_promotion_allowed") is not False
        or local_domain.get("clause_verdicts", {}).get(
            "continuum_promotion_blocked_and_recorded"
        )
        is not True
    ):
        raise VerificationError("finite local-domain boundary changed")

    ew = _load(EW_FRONTIER_PATH)
    if (
        ew.get("schema") != "oph.higgs_yukawa_source_frontier.v1"
        or ew.get("issue") != 630
        or ew.get("status")
        != "BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN"
        or ew.get("physical_source_action_emitted") is not False
        or ew.get("promotion_allowed") is not False
    ):
        raise VerificationError("electroweak source frontier changed")

    einstein = _load(EINSTEIN_PATH)
    if (
        einstein.get("artifact") != "einstein_branch_realized_receipt_evaluation"
        or einstein.get("issue") != 503
        or einstein.get("realized_geometric_branch_certified_nonempty") is not False
    ):
        raise VerificationError("inhabited Einstein boundary changed")

    quantum_eft = _load(QUANTUM_EFT_PATH)
    interfaces = quantum_eft.get("required_interfaces")
    if not isinstance(interfaces, list):
        raise VerificationError("quantum-EFT required-interface list missing")
    selected = [
        row
        for row in interfaces
        if isinstance(row, dict)
        and row.get("gate_id") == "finite_to_lorentzian_quantum_eft_transfer"
    ]
    if (
        quantum_eft.get("schema") != "oph.wz.source_parent_inventory.v1"
        or quantum_eft.get("provenance_issue") != 594
        or quantum_eft.get("promotion_allowed") is not False
        or len(selected) != 1
        or selected[0].get("provenance_issues") != [635]
        or selected[0].get("classification") != "not_supplied"
        or selected[0].get("supplied_evidence") != []
        or selected[0].get("terminal_for_dimensionless_output") is not True
    ):
        raise VerificationError("quantum-EFT source frontier changed")

    qcd = _load(QCD_RESOURCE_PATH)
    if (
        qcd.get("artifact") != "qcd_solver_on_standby_receipt"
        or qcd.get("status")
        != "SOLVER_COMPILED_AND_SMOKE_BLOCKED_INVOCATION_GATED_ON_SOURCE_PARAMETERS"
        or qcd.get("promotion_allowed") is not False
        or qcd.get("external_targets_used") != []
    ):
        raise VerificationError("QCD resource boundary changed")

    return [
        {
            "boundary_id": "finite_local_domain",
            "evidence_class": "bounded_local_receipt",
            "issue_context": [634, 635],
            "pin": _pin(
                LOCAL_DOMAIN_PATH,
                "finite domain with continuum promotion blocked",
            ),
            "status": "FINITE_DOMAIN_ATTAINED__CONTINUUM_PROMOTION_BLOCKED",
        },
        {
            "boundary_id": "electroweak_source_action",
            "evidence_class": "bounded_local_frontier",
            "issue_context": [630],
            "pin": _pin(
                EW_FRONTIER_PATH,
                "nonpromoting scalar and Yukawa source frontier",
            ),
            "status": "POSITIVE_SOURCE_ACTION_OPEN",
        },
        {
            "boundary_id": "inhabited_einstein_tower",
            "evidence_class": "bounded_local_frontier",
            "issue_context": [503],
            "pin": _pin(
                EINSTEIN_PATH,
                "partial Einstein branch receipts without full inhabited antecedent",
            ),
            "status": "FULL_INHABITED_ANTECEDENT_NOT_CERTIFIED",
        },
        {
            "boundary_id": "lorentzian_quantum_eft_transfer",
            "evidence_class": "required_interface_not_supplied",
            "issue_context": [635],
            "pin": _pin(
                QUANTUM_EFT_PATH,
                "source inventory with the quantum-EFT transfer classified required and not supplied",
            ),
            "status": "REQUIRED_INTERFACE_NOT_SUPPLIED_ON_DECLARED_SOURCE_INVENTORY",
        },
        {
            "boundary_id": "qcd_spectral_resource",
            "evidence_class": "bounded_local_resource_receipt",
            "issue_context": [425, 294],
            "pin": _pin(
                QCD_RESOURCE_PATH,
                "standby solver with no source spectral output",
            ),
            "status": "NO_PROMOTED_QCD_SPECTRAL_OUTPUT_ON_DECLARED_LOCAL_PATH",
        },
    ]


def _projection() -> dict[str, Any]:
    classical = _load(CLASSICAL_PATH)
    color = _load(COLOR_PATH)
    if (
        classical.get("artifact") != "oph_massless_carrier_mode_acceptance"
        or classical.get("schema") != "oph_carrier_mode_quantum_particle_gate_v1"
        or classical.get("status")
        != "classical_action_branch_modes_recorded_quantum_particle_gate_open"
    ):
        raise VerificationError("classical carrier parent identity changed")
    try:
        color_dimension = color["closure"]["derived_block_dimensions"]["even_block_su3"]
        dimensions_verified = color["port_to_generator_map"]["block_dimensions_verified"]
    except (KeyError, TypeError) as exc:
        raise VerificationError("color parent is incomplete") from exc
    if color_dimension != 8 or dimensions_verified is not True:
        raise VerificationError("color adjoint is not fixed at dimension eight")

    rows = classical.get("carriers")
    if not isinstance(rows, list):
        raise VerificationError("classical carrier rows missing")
    by_id = {row.get("carrier_id"): row for row in rows if isinstance(row, dict)}
    if set(by_id) != set(ORDER):
        raise VerificationError("classical carrier set changed")

    projected: list[dict[str, Any]] = []
    counts: list[int] = []
    for carrier_id in ORDER:
        row = by_id[carrier_id]
        expected = CLASSICAL_ROWS[carrier_id]
        for key in ("label", "branch", "physical_degrees_of_freedom"):
            if row.get(key) != expected[key]:
                raise VerificationError(f"{carrier_id}: {key} changed")
        if (
            row.get("branch_is_additional_input_not_group_output") is not True
            or row.get("abstract_symmetry_group_alone_sufficient") is not False
            or row.get("particle_promotion_allowed") is not False
        ):
            raise VerificationError(f"{carrier_id}: branch boundary changed")
        classical_gate = row.get("classical_carrier_gate", {})
        quantum_gate = row.get("quantum_particle_gate", {})
        evidence = row.get("evidence", {})
        if (
            classical_gate.get("passed") is not True
            or classical_gate.get("missing") != []
            or classical_gate.get("status")
            != "conditional_pass_on_declared_action_phase_branch"
        ):
            raise VerificationError(f"{carrier_id}: classical gate changed")
        if (
            quantum_gate.get("passed") is not False
            or quantum_gate.get("missing") != list(QUANTUM_EVIDENCE)
        ):
            raise VerificationError(f"{carrier_id}: quantum gate changed")
        if any(evidence.get(key) is not False for key in QUANTUM_EVIDENCE):
            raise VerificationError(f"{carrier_id}: unsupported quantum evidence appeared")
        if row.get("hard_quadratic_mass_parameter_squared") != 0:
            raise VerificationError(f"{carrier_id}: hard mass parameter changed")

        factor = expected["multiplicity_factor"]
        if carrier_id == "gluon" and factor != color_dimension:
            raise VerificationError("color multiplicity differs from the adjoint dimension")
        total = 2 * factor
        counts.append(total)
        projected.append(
            {
                "branch": row["branch"],
                "branch_is_additional_input_not_group_output": True,
                "carrier_id": carrier_id,
                "classical_gate_passed": True,
                "classical_gate_status": classical_gate["status"],
                "continuum_spacetime_dimension": DIMENSION,
                "hard_quadratic_mass_parameter_squared": 0,
                "label": row["label"],
                "modes_per_carrier_component": 2,
                "multiplicity_factor": factor,
                "multiplicity_role": expected["multiplicity_role"],
                "quantum_gate_passed": False,
                "quantum_missing_receipts": list(QUANTUM_EVIDENCE),
                "total_classical_modes": total,
            }
        )

    return {
        "classical_carriers": projected,
        "classical_mode_vector": counts,
        "classical_mode_vector_order": list(ORDER),
        "comparison_policy": {
            "blind_prediction_eligible": False,
            "comparison_data_present": False,
            "comparison_values_consumed": False,
            "laboratory_values_present": False,
            "status_packet_only": True,
            "target_named_rows_present": True,
        },
        "continuum_spacetime_dimension": DIMENSION,
        "declared_boundary_evidence": _boundary_evidence(),
        "schema": "oph.quantum_carrier_status_source_projection.v2",
        "scope": (
            "Comparison-value-free projection of three target-named conditional classical carrier rows and the exact color-adjoint dimension. The mode vector contains differently typed four-dimensional propagating-mode counts. It is neither a particle-count vector nor a uniform gauge-algebra-dimension vector. The boundary evidence is limited to the pinned declared corpus."
        ),
        "source_pins": [
            _pin(
                CLASSICAL_PATH,
                "conditional classical carrier modes and unpassed quantum evidence",
            ),
            _pin(
                COLOR_PATH,
                "exact color-adjoint dimension used only to total the perturbative gluon modes",
            ),
        ],
    }


def _capabilities(carrier_id: str) -> dict[str, Any]:
    config = ROW_CONFIG[carrier_id]
    return {
        "state_space": {
            "available_classical_object": "conditional reduced classical transverse phase space",
            "classical_object_available": True,
            "physical_quantum_object_available": False,
            "required_quantum_object": config["state_space"],
            "status": "CLASSICAL_ONLY",
        },
        "observable_algebra": {
            "available_classical_object": config["classical_observable_algebra"],
            "classical_object_available": True,
            "physical_quantum_object_available": False,
            "required_quantum_object": config["observable_algebra"],
            "status": "CLASSICAL_ONLY",
        },
        "gauge_quotient": {
            "available_classical_object": (
                "conditional classical gauge fixing and reduced mode quotient"
            ),
            "classical_object_available": True,
            "physical_quantum_object_available": False,
            "required_quantum_object": config["gauge_quotient"],
            "status": "CLASSICAL_ONLY",
        },
        "vacuum": {
            "available_declared_object": (
                "the conditional action, phase, or background branch named by the classical parent"
            ),
            "branch_declared": True,
            "required_quantum_object": config["vacuum"],
            "source_selected_quantum_object_available": False,
            "status": "DECLARED_BRANCH_ONLY",
        },
        "spectral_object": {
            "available_classical_object": (
                "conditional classical characteristic or reduced Hamiltonian spectrum"
            ),
            "classical_object_available": True,
            "positive_physical_quantum_object_available": False,
            "required_quantum_object": config["spectral_object"],
            "status": "CLASSICAL_ONLY",
        },
        "physical_current_residue": {
            "nonzero_positive_residue_available": False,
            "required_quantum_object": config["residue"],
            "status": "MISSING_SOURCE_CONSTRUCTION",
        },
        "refinement_control": {
            "available": False,
            "required_object": (
                "a source-derived refinement family preserving the physical pole classification"
            ),
            "status": "MISSING_SOURCE_CONSTRUCTION",
        },
        "phase_or_asymptotic_sector": {
            "declared_classical_branch_available": True,
            "physical_quantum_sector_available": False,
            "status": config["phase_status"],
        },
    }


def _row(source: Mapping[str, Any]) -> dict[str, Any]:
    carrier_id = str(source["carrier_id"])
    config = ROW_CONFIG[carrier_id]
    return {
        "blocking_frontier": config["blocking_frontier"],
        "capabilities": _capabilities(carrier_id),
        "carrier_id": carrier_id,
        "classical_baseline": {
            "branch": source["branch"],
            "branch_is_additional_input_not_group_output": source[
                "branch_is_additional_input_not_group_output"
            ],
            "classical_gate_status": source["classical_gate_status"],
            "continuum_spacetime_dimension": source["continuum_spacetime_dimension"],
            "exact_total_mode_count": source["total_classical_modes"],
            "hard_quadratic_mass_parameter_squared": source[
                "hard_quadratic_mass_parameter_squared"
            ],
            "modes_per_carrier_component": source["modes_per_carrier_component"],
            "multiplicity_factor": source["multiplicity_factor"],
            "multiplicity_role": source["multiplicity_role"],
            "particle_claim": False,
        },
        "label": source["label"],
        "open_interfaces": config["open_interfaces"],
        "particle_promotion_allowed": False,
        "resource_boundaries": config["resource_boundaries"],
        "strongest_supported_statement": config["strongest_supported_statement"],
        "verdict": config["verdict"],
        "verdict_class": "EXPLICIT_NOT_EVALUABLE",
    }


def _expected_receipt(projection: Mapping[str, Any]) -> dict[str, Any]:
    projection_raw = PROJECTION_PATH.read_bytes()
    rows = [_row(row) for row in projection["classical_carriers"]]
    receipt: dict[str, Any] = {
        "allowed_exit_policy": (
            "Each carrier row exits only through a source-positive pole receipt, a rigorous negative theorem, or an explicit NOT_EVALUABLE verdict. No classical zero or conditional propagator can promote a quantum particle."
        ),
        "all_rows_at_allowed_exit": all(
            row["verdict_class"] == "EXPLICIT_NOT_EVALUABLE"
            and row["particle_promotion_allowed"] is False
            for row in rows
        ),
        "artifact": "oph_quantum_carrier_status",
        "blind_prediction_eligible": False,
        "classical_mode_vector": projection["classical_mode_vector"],
        "classical_mode_vector_order": projection["classical_mode_vector_order"],
        "classical_mode_vector_scope": (
            "Exact conditional four-dimensional propagating-mode totals on three differently typed declared branches: two Maxwell modes from one U(1) generator, sixteen perturbative color modes from eight SU(3) adjoint generators, and two Einstein transverse-traceless modes from one metric tensor field. This is neither a particle-count vector nor a uniform gauge-algebra-dimension vector."
        ),
        "comparison_values_consumed": False,
        "continuum_spacetime_dimension": DIMENSION,
        "github_issue": 552,
        "parent_dependency_context": [
            {
                "issue": 634,
                "role": (
                    "finite local source domain; it does not supply a Lorentzian physical quantum state space or pole"
                ),
            },
            {
                "issue": 635,
                "role": "finite-source to Lorentzian quantum-EFT bridge",
            },
            {
                "issue": 630,
                "role": "electroweak source integration and physical action selection",
            },
            {
                "issue": 503,
                "role": "inhabited source-derived typed Einstein tower",
            },
        ],
        "producer_scope": (
            "Bounded classification of the pinned declared corpus. It proves no exhaustive absence claim outside that projection."
        ),
        "rows": rows,
        "schema": "oph.quantum_carrier_status.v2",
        "source_projection_pin": {
            "bytes": len(projection_raw),
            "path": PROJECTION_PATH.relative_to(REPO_ROOT).as_posix(),
            "sha256": _digest(projection_raw),
        },
        "status": "THREE_ROW_EXPLICIT_NOT_EVALUABLE",
        "target_named_status_rows": True,
    }
    receipt["receipt_sha256"] = _receipt_hash(receipt)
    return receipt


def verify() -> dict[str, Any]:
    expected_projection = _projection()
    projection = _load(PROJECTION_PATH)
    if projection != expected_projection:
        raise VerificationError("source projection differs from exact reconstruction")
    if projection["classical_mode_vector"] != [2, 16, 2]:
        raise VerificationError("conditional classical mode vector changed")

    expected_receipt = _expected_receipt(projection)
    receipt = _load(RECEIPT_PATH)
    if receipt.get("receipt_sha256") != _receipt_hash(receipt):
        raise VerificationError("receipt hash mismatch")
    if receipt != expected_receipt:
        raise VerificationError("receipt differs from exact independent reconstruction")
    return receipt


def main() -> int:
    try:
        receipt = verify()
    except VerificationError as exc:
        print(f"QUANTUM_CARRIER_STATUS_INDEPENDENT_INVALID: {exc}")
        return 1
    print(
        "QUANTUM_CARRIER_STATUS_INDEPENDENT_VALID "
        f"{receipt['receipt_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
