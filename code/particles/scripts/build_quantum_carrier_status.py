#!/usr/bin/env python3
"""Build the comparison-value-free issue-552 quantum-carrier status packet.

The parent carrier receipt proves only conditional classical mode statements.
This producer preserves those exact statements and classifies the quantum pole
interfaces row by row. It consumes no laboratory mass, pole, width, or other
comparison value. The named carrier rows make this a status audit, rather than
a blind prediction packet.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


CODE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = CODE_ROOT.parent
PARTICLES_ROOT = CODE_ROOT / "particles"
PROJECTION_PATH = (
    PARTICLES_ROOT / "manifests" / "quantum_carrier_status_source_projection.json"
)
CLASSICAL_PARENT_PATH = (
    PARTICLES_ROOT / "runs" / "status" / "carrier_mode_acceptance.json"
)
COLOR_PARENT_PATH = (
    CODE_ROOT / "a5_closure" / "receipts" / "port_current_inner_reference.receipt.json"
)
LOCAL_DOMAIN_BOUNDARY_PATH = (
    CODE_ROOT / "a5_closure" / "manifests" / "stage4_receipt.json"
)
EW_SOURCE_FRONTIER_PATH = (
    PARTICLES_ROOT
    / "hierarchy"
    / "higgs_yukawa_source_frontier"
    / "outputs"
    / "higgs_yukawa_source_frontier.json"
)
EINSTEIN_BRANCH_PATH = (
    CODE_ROOT / "geometry" / "runs" / "realized_branch_receipt_report.json"
)
QUANTUM_EFT_FRONTIER_PATH = (
    PARTICLES_ROOT
    / "calibration"
    / "wz_native_source_packet"
    / "outputs"
    / "source_parent_inventory.json"
)
QCD_RESOURCE_BOUNDARY_PATH = (
    PARTICLES_ROOT
    / "runs"
    / "qcd"
    / "hadron_source_backend"
    / "qcd_ensemble"
    / "solver_on_standby.json"
)
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "quantum_carrier_status.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "QUANTUM_CARRIER_STATUS.md"

SCHEMA = "oph.quantum_carrier_status.v2"
ARTIFACT = "oph_quantum_carrier_status"
STATUS = "THREE_ROW_EXPLICIT_NOT_EVALUABLE"
ROW_ORDER = ("photon", "gluon", "graviton")
CONTINUUM_SPACETIME_DIMENSION = 4
QUANTUM_MISSING = (
    "quantization_constructed_from_oph",
    "positive_physical_hilbert_space",
    "physical_two_point_or_hamiltonian_spectrum_constructed",
    "positive_residue_massless_pole",
    "asymptotic_or_deconfined_particle_state",
)

EXPECTED_CLASSICAL_FIELDS = {
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
        "state_space": (
            "a positive physical photon state space after the gauge quotient"
        ),
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
                "requirement": (
                    "finite-source to Lorentzian quantum-EFT construction"
                ),
            },
        ],
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
        "observable_algebra": (
            "a source-derived gauge-invariant QCD observable algebra"
        ),
        "classical_observable_algebra": (
            "classical gauge-invariant Yang-Mills polynomial observables"
        ),
        "gauge_quotient": (
            "a nonperturbative physical gauge quotient with positivity"
        ),
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
                "requirement": (
                    "finite-source to Lorentzian quantum-EFT construction"
                ),
            },
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
                "requirement": (
                    "one inhabited source-derived typed Einstein tower"
                ),
            },
            {
                "issue": 635,
                "requirement": (
                    "finite-source to Lorentzian quantum-EFT construction"
                ),
            },
        ],
        "strongest_supported_statement": (
            "On the declared pure Einstein-Hilbert linearization about a suitable Ricci-flat background, the hard quadratic mass parameter is zero and the transverse-traceless classical system has two modes. This is not a graviton Hilbert-space, pole, or rest-mass prediction."
        ),
    },
}


class CertificateError(ValueError):
    """Fail-closed status-packet error with a stable mutation-test code."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise CertificateError("MISSING", f"required input is absent: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("JSON", f"cannot read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CertificateError("TYPE", f"top-level JSON object required: {path}")
    return value


def _sha256_bytes(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _file_pin(path: Path, role: str) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "bytes": len(raw),
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "role": role,
        "sha256": _sha256_bytes(raw),
    }


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _pretty_bytes(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )


def _receipt_hash(value: Mapping[str, Any]) -> str:
    core = {key: item for key, item in value.items() if key != "receipt_sha256"}
    return _sha256_bytes(_canonical_bytes(core))


def _declared_boundary_evidence() -> list[dict[str, Any]]:
    """Validate and pin the bounded local evidence used by the status audit."""

    local_domain = _load_json(LOCAL_DOMAIN_BOUNDARY_PATH)
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
        raise CertificateError(
            "LOCAL_DOMAIN_BOUNDARY",
            "the finite local-domain boundary left its declared scope",
        )

    ew_frontier = _load_json(EW_SOURCE_FRONTIER_PATH)
    if (
        ew_frontier.get("schema") != "oph.higgs_yukawa_source_frontier.v1"
        or ew_frontier.get("issue") != 630
        or ew_frontier.get("status")
        != "BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN"
        or ew_frontier.get("physical_source_action_emitted") is not False
        or ew_frontier.get("promotion_allowed") is not False
    ):
        raise CertificateError(
            "EW_SOURCE_BOUNDARY",
            "the electroweak source frontier left its nonpromoting scope",
        )

    einstein = _load_json(EINSTEIN_BRANCH_PATH)
    if (
        einstein.get("artifact") != "einstein_branch_realized_receipt_evaluation"
        or einstein.get("issue") != 503
        or einstein.get("realized_geometric_branch_certified_nonempty") is not False
    ):
        raise CertificateError(
            "EINSTEIN_BOUNDARY",
            "the full inhabited Einstein antecedent was promoted",
        )

    quantum_eft = _load_json(QUANTUM_EFT_FRONTIER_PATH)
    required_interfaces = quantum_eft.get("required_interfaces")
    if not isinstance(required_interfaces, list):
        raise CertificateError(
            "QUANTUM_EFT_BOUNDARY", "required-interface list is absent"
        )
    quantum_rows = [
        row
        for row in required_interfaces
        if isinstance(row, dict)
        and row.get("gate_id") == "finite_to_lorentzian_quantum_eft_transfer"
    ]
    if (
        quantum_eft.get("schema") != "oph.wz.source_parent_inventory.v1"
        or quantum_eft.get("provenance_issue") != 594
        or quantum_eft.get("promotion_allowed") is not False
        or len(quantum_rows) != 1
        or quantum_rows[0].get("provenance_issues") != [635]
        or quantum_rows[0].get("classification") != "not_supplied"
        or quantum_rows[0].get("supplied_evidence") != []
    ):
        raise CertificateError(
            "QUANTUM_EFT_BOUNDARY",
            "the finite-to-Lorentzian quantum-EFT frontier changed",
        )

    qcd = _load_json(QCD_RESOURCE_BOUNDARY_PATH)
    if (
        qcd.get("artifact") != "qcd_solver_on_standby_receipt"
        or qcd.get("status")
        != "SOLVER_COMPILED_AND_SMOKE_BLOCKED_INVOCATION_GATED_ON_SOURCE_PARAMETERS"
        or qcd.get("promotion_allowed") is not False
        or qcd.get("external_targets_used") != []
    ):
        raise CertificateError(
            "QCD_RESOURCE_BOUNDARY",
            "the local QCD standby receipt emitted a promoted source result",
        )

    return [
        {
            "boundary_id": "finite_local_domain",
            "evidence_class": "bounded_local_receipt",
            "issue_context": [634, 635],
            "pin": _file_pin(
                LOCAL_DOMAIN_BOUNDARY_PATH,
                "finite domain with continuum promotion blocked",
            ),
            "status": "FINITE_DOMAIN_ATTAINED__CONTINUUM_PROMOTION_BLOCKED",
        },
        {
            "boundary_id": "electroweak_source_action",
            "evidence_class": "bounded_local_frontier",
            "issue_context": [630],
            "pin": _file_pin(
                EW_SOURCE_FRONTIER_PATH,
                "nonpromoting scalar and Yukawa source frontier",
            ),
            "status": "POSITIVE_SOURCE_ACTION_OPEN",
        },
        {
            "boundary_id": "inhabited_einstein_tower",
            "evidence_class": "bounded_local_frontier",
            "issue_context": [503],
            "pin": _file_pin(
                EINSTEIN_BRANCH_PATH,
                "partial Einstein branch receipts without full inhabited antecedent",
            ),
            "status": "FULL_INHABITED_ANTECEDENT_NOT_CERTIFIED",
        },
        {
            "boundary_id": "lorentzian_quantum_eft_transfer",
            "evidence_class": "declared_open_interface",
            "issue_context": [635],
            "pin": _file_pin(
                QUANTUM_EFT_FRONTIER_PATH,
                "source inventory with the quantum-EFT transfer typed open",
            ),
            "status": "OPEN_ON_DECLARED_SOURCE_INVENTORY",
        },
        {
            "boundary_id": "qcd_spectral_resource",
            "evidence_class": "bounded_local_resource_receipt",
            "issue_context": [425, 294],
            "pin": _file_pin(
                QCD_RESOURCE_BOUNDARY_PATH,
                "standby solver with no source spectral output",
            ),
            "status": "NO_PROMOTED_QCD_SPECTRAL_OUTPUT_ON_DECLARED_LOCAL_PATH",
        },
    ]


def reconstruct_source_projection() -> dict[str, Any]:
    """Reconstruct the admitted comparison-value-free fields and boundaries."""

    classical = _load_json(CLASSICAL_PARENT_PATH)
    color = _load_json(COLOR_PARENT_PATH)
    if (
        classical.get("artifact") != "oph_massless_carrier_mode_acceptance"
        or classical.get("schema") != "oph_carrier_mode_quantum_particle_gate_v1"
        or classical.get("status")
        != "classical_action_branch_modes_recorded_quantum_particle_gate_open"
    ):
        raise CertificateError("CLASSICAL_PARENT", "unexpected carrier parent identity")

    try:
        color_dimension = color["closure"]["derived_block_dimensions"]["even_block_su3"]
        color_dimensions_verified = color["port_to_generator_map"][
            "block_dimensions_verified"
        ]
    except (KeyError, TypeError) as exc:
        raise CertificateError("COLOR_PARENT", "missing exact color-dimension fields") from exc
    if color_dimension != 8 or color_dimensions_verified is not True:
        raise CertificateError(
            "COLOR_PARENT", "the color parent does not certify an eight-dimensional adjoint"
        )

    rows = classical.get("carriers")
    if not isinstance(rows, list):
        raise CertificateError("CLASSICAL_ROWS", "carrier list is absent")
    by_id = {row.get("carrier_id"): row for row in rows if isinstance(row, dict)}
    if set(by_id) != set(ROW_ORDER):
        raise CertificateError("CLASSICAL_ROWS", "expected exactly photon, gluon, graviton")

    projected_rows: list[dict[str, Any]] = []
    mode_vector: list[int] = []
    for carrier_id in ROW_ORDER:
        row = by_id[carrier_id]
        expected = EXPECTED_CLASSICAL_FIELDS[carrier_id]
        for key in ("label", "branch", "physical_degrees_of_freedom"):
            expected_value = expected[key]
            if row.get(key) != expected_value:
                raise CertificateError(
                    "CLASSICAL_FIELD", f"{carrier_id}.{key} left the frozen branch"
                )
        if (
            row.get("branch_is_additional_input_not_group_output") is not True
            or row.get("abstract_symmetry_group_alone_sufficient") is not False
            or row.get("particle_promotion_allowed") is not False
        ):
            raise CertificateError(
                "CLASSICAL_BOUNDARY",
                f"{carrier_id} lost its declared-input or nonpromotion boundary",
            )
        classical_gate = row.get("classical_carrier_gate")
        quantum_gate = row.get("quantum_particle_gate")
        if not isinstance(classical_gate, dict) or not isinstance(quantum_gate, dict):
            raise CertificateError("GATE_TYPE", f"{carrier_id} gate object missing")
        if (
            classical_gate.get("passed") is not True
            or classical_gate.get("missing") != []
            or classical_gate.get("status")
            != "conditional_pass_on_declared_action_phase_branch"
        ):
            raise CertificateError("CLASSICAL_GATE", f"{carrier_id} classical gate failed")
        if quantum_gate.get("passed") is not False:
            raise CertificateError("QUANTUM_GATE", f"{carrier_id} quantum gate promoted")
        missing = quantum_gate.get("missing")
        if missing != list(QUANTUM_MISSING):
            raise CertificateError(
                "QUANTUM_GATE", f"{carrier_id} missing-receipt vector changed"
            )
        evidence = row.get("evidence")
        if not isinstance(evidence, dict):
            raise CertificateError("EVIDENCE", f"{carrier_id} evidence object missing")
        if any(evidence.get(key) is not False for key in QUANTUM_MISSING):
            raise CertificateError("EVIDENCE", f"{carrier_id} has unsupported quantum evidence")
        if row.get("hard_quadratic_mass_parameter_squared") != 0:
            raise CertificateError("CLASSICAL_MASS", f"{carrier_id} mass parameter changed")

        multiplicity_factor = expected["multiplicity_factor"]
        if carrier_id == "gluon" and multiplicity_factor != color_dimension:
            raise CertificateError(
                "COLOR_PARENT",
                "the configured color multiplicity differs from the exact adjoint dimension",
            )
        total_modes = 2 * multiplicity_factor
        mode_vector.append(total_modes)
        projected_rows.append(
            {
                "branch": row["branch"],
                "branch_is_additional_input_not_group_output": True,
                "carrier_id": carrier_id,
                "classical_gate_passed": True,
                "classical_gate_status": classical_gate["status"],
                "continuum_spacetime_dimension": CONTINUUM_SPACETIME_DIMENSION,
                "hard_quadratic_mass_parameter_squared": 0,
                "label": row["label"],
                "modes_per_carrier_component": 2,
                "multiplicity_factor": multiplicity_factor,
                "multiplicity_role": expected["multiplicity_role"],
                "quantum_gate_passed": False,
                "quantum_missing_receipts": list(missing),
                "total_classical_modes": total_modes,
            }
        )

    return {
        "classical_carriers": projected_rows,
        "classical_mode_vector": mode_vector,
        "classical_mode_vector_order": list(ROW_ORDER),
        "continuum_spacetime_dimension": CONTINUUM_SPACETIME_DIMENSION,
        "declared_boundary_evidence": _declared_boundary_evidence(),
        "schema": "oph.quantum_carrier_status_source_projection.v2",
        "scope": (
            "Comparison-value-free projection of three target-named conditional classical carrier rows and the exact color-adjoint dimension. The mode vector contains differently typed four-dimensional propagating-mode counts. It is neither a particle-count vector nor a uniform gauge-algebra-dimension vector. The boundary evidence is limited to the pinned declared corpus."
        ),
        "source_pins": [
            _file_pin(
                CLASSICAL_PARENT_PATH,
                "conditional classical carrier modes and unpassed quantum evidence",
            ),
            _file_pin(
                COLOR_PARENT_PATH,
                "exact color-adjoint dimension used only to total the perturbative gluon modes",
            ),
        ],
        "comparison_policy": {
            "blind_prediction_eligible": False,
            "comparison_data_present": False,
            "comparison_values_consumed": False,
            "laboratory_values_present": False,
            "status_packet_only": True,
            "target_named_rows_present": True,
        },
    }


def load_and_validate_projection() -> dict[str, Any]:
    projection = _load_json(PROJECTION_PATH)
    expected = reconstruct_source_projection()
    if projection != expected:
        raise CertificateError(
            "PROJECTION",
            "committed comparison-value-free projection differs from reconstruction",
        )
    return projection


def _capability_inventory(carrier_id: str) -> dict[str, Any]:
    config = ROW_CONFIG[carrier_id]
    return {
        "state_space": {
            "available_classical_object": (
                "conditional reduced classical transverse phase space"
            ),
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


def _build_row(source_row: Mapping[str, Any]) -> dict[str, Any]:
    carrier_id = str(source_row["carrier_id"])
    config = ROW_CONFIG[carrier_id]
    return {
        "capabilities": _capability_inventory(carrier_id),
        "carrier_id": carrier_id,
        "classical_baseline": {
            "branch": source_row["branch"],
            "branch_is_additional_input_not_group_output": source_row[
                "branch_is_additional_input_not_group_output"
            ],
            "classical_gate_status": source_row["classical_gate_status"],
            "continuum_spacetime_dimension": source_row[
                "continuum_spacetime_dimension"
            ],
            "exact_total_mode_count": source_row["total_classical_modes"],
            "hard_quadratic_mass_parameter_squared": source_row[
                "hard_quadratic_mass_parameter_squared"
            ],
            "modes_per_carrier_component": source_row["modes_per_carrier_component"],
            "multiplicity_factor": source_row["multiplicity_factor"],
            "multiplicity_role": source_row["multiplicity_role"],
            "particle_claim": False,
        },
        "blocking_frontier": config["blocking_frontier"],
        "label": source_row["label"],
        "open_interfaces": config["open_interfaces"],
        "particle_promotion_allowed": False,
        "resource_boundaries": config.get("resource_boundaries", []),
        "strongest_supported_statement": config["strongest_supported_statement"],
        "verdict": config["verdict"],
        "verdict_class": "EXPLICIT_NOT_EVALUABLE",
    }


def build_payload() -> dict[str, Any]:
    projection = load_and_validate_projection()
    payload = build_payload_without_validation(projection)
    payload["receipt_sha256"] = _receipt_hash(payload)
    validate_payload(payload)
    return payload


def validate_payload(payload: Mapping[str, Any]) -> None:
    if not isinstance(payload, Mapping):
        raise CertificateError("PAYLOAD_TYPE", "receipt must be a JSON object")
    recorded_hash = payload.get("receipt_sha256")
    if not isinstance(recorded_hash, str) or recorded_hash != _receipt_hash(payload):
        raise CertificateError("RECEIPT_HASH", "receipt hash mismatch")

    projection = load_and_validate_projection()
    if payload.get("schema") != SCHEMA or payload.get("artifact") != ARTIFACT:
        raise CertificateError("IDENTITY", "unexpected receipt identity")
    if payload.get("github_issue") != 552 or payload.get("status") != STATUS:
        raise CertificateError("STATUS", "unexpected issue or aggregate status")
    if payload.get("comparison_values_consumed") is not False:
        raise CertificateError("COMPARISON", "laboratory comparison values are forbidden")
    if payload.get("blind_prediction_eligible") is not False:
        raise CertificateError("BLINDNESS", "a target-named status packet is not blind")
    if payload.get("all_rows_at_allowed_exit") is not True:
        raise CertificateError("EXIT", "not every row reached an allowed exit")
    if payload.get("classical_mode_vector_order") != list(ROW_ORDER):
        raise CertificateError("MODE_ORDER", "classical mode order changed")
    if payload.get("classical_mode_vector") != [2, 16, 2]:
        raise CertificateError("MODE_VECTOR", "classical mode vector changed")
    if payload.get("continuum_spacetime_dimension") != 4:
        raise CertificateError("DIMENSION", "the declared four-dimensional premise changed")

    projection_raw = PROJECTION_PATH.read_bytes()
    expected_projection_pin = {
        "bytes": len(projection_raw),
        "path": PROJECTION_PATH.relative_to(REPO_ROOT).as_posix(),
        "sha256": _sha256_bytes(projection_raw),
    }
    if payload.get("source_projection_pin") != expected_projection_pin:
        raise CertificateError("PROJECTION_PIN", "source projection pin mismatch")

    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != 3:
        raise CertificateError("ROWS", "exactly three carrier rows required")
    expected_rows = [_build_row(row) for row in projection["classical_carriers"]]
    if rows != expected_rows:
        raise CertificateError("ROWS", "typed carrier rows differ from reconstruction")

    expected = {
        **build_payload_without_validation(projection),
        "receipt_sha256": recorded_hash,
    }
    if dict(payload) != expected:
        raise CertificateError("PAYLOAD", "unexpected field or value in receipt")


def build_payload_without_validation(projection: Mapping[str, Any]) -> dict[str, Any]:
    """Construct the receipt core without recursion through ``validate_payload``."""

    projection_raw = PROJECTION_PATH.read_bytes()
    rows = [_build_row(row) for row in projection["classical_carriers"]]
    all_rows_at_allowed_exit = all(
        row["verdict_class"] == "EXPLICIT_NOT_EVALUABLE"
        and row["particle_promotion_allowed"] is False
        for row in rows
    )
    return {
        "allowed_exit_policy": (
            "Each carrier row exits only through a source-positive pole receipt, a rigorous negative theorem, or an explicit NOT_EVALUABLE verdict. No classical zero or conditional propagator can promote a quantum particle."
        ),
        "all_rows_at_allowed_exit": all_rows_at_allowed_exit,
        "artifact": ARTIFACT,
        "blind_prediction_eligible": False,
        "classical_mode_vector": projection["classical_mode_vector"],
        "classical_mode_vector_order": projection["classical_mode_vector_order"],
        "classical_mode_vector_scope": (
            "Exact conditional four-dimensional propagating-mode totals on three differently typed declared branches: two Maxwell modes from one U(1) generator, sixteen perturbative color modes from eight SU(3) adjoint generators, and two Einstein transverse-traceless modes from one metric tensor field. This is neither a particle-count vector nor a uniform gauge-algebra-dimension vector."
        ),
        "comparison_values_consumed": False,
        "continuum_spacetime_dimension": CONTINUUM_SPACETIME_DIMENSION,
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
        "schema": SCHEMA,
        "source_projection_pin": {
            "bytes": len(projection_raw),
            "path": PROJECTION_PATH.relative_to(REPO_ROOT).as_posix(),
            "sha256": _sha256_bytes(projection_raw),
        },
        "status": STATUS,
        "target_named_status_rows": True,
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Quantum Carrier Status",
        "",
        "This packet classifies the photon, gluon, and graviton rows at the quantum-particle gate. It consumes no laboratory comparison values. The named rows make it ineligible as a blind prediction.",
        "",
        "The exact conditional classical mode vector is `(2, 16, 2)` in four dimensions: two Maxwell modes from one U(1) generator, sixteen perturbative color modes from eight SU(3) adjoint generators, and two Einstein transverse-traceless modes from one metric tensor field. The entries are differently typed propagating-mode totals. They are neither particle counts nor a uniform gauge-algebra-dimension vector.",
        "",
        "| Carrier | Exact conditional classical modes | Quantum verdict | Blocking frontier |",
        "| --- | ---: | --- | --- |",
    ]
    for row in payload["rows"]:
        baseline = row["classical_baseline"]
        blockers = ", ".join(f"`{item}`" for item in row["blocking_frontier"])
        lines.append(
            f"| {row['label']} | `{baseline['exact_total_mode_count']}` | "
            f"`{row['verdict']}` | {blockers} |"
        )
    lines.extend(
        [
            "",
            "Each row separately records the required state space, observable algebra, gauge quotient, vacuum, spectral object, physical-current residue, refinement control, and phase or asymptotic sector. The available objects stop at the conditional classical layer. No row has a positive physical quantum spectrum, a nonzero physical-current residue, or a source-derived refinement family that preserves a physical pole classification.",
            "",
            "The gluon row does not promote a colored gauge-potential pole to an asymptotic particle in confining QCD. The photon and graviton rows do not turn a zero hard quadratic action parameter into a measured zero rest mass.",
            "",
            "The not-evaluable verdicts are bounded to the pinned declared corpus. The resource rows record missing local production paths and do not prove that every possible completion fails.",
            "",
            f"Receipt: `{payload['receipt_sha256']}`.",
            "",
        ]
    )
    return "\n".join(lines)


def validate_committed(
    json_path: Path = DEFAULT_JSON_OUT,
    markdown_path: Path = DEFAULT_MD_OUT,
) -> dict[str, Any]:
    expected = build_payload()
    actual = _load_json(json_path)
    validate_payload(actual)
    if json_path.read_bytes() != _pretty_bytes(expected):
        raise CertificateError("JSON_BYTES", "committed receipt is not byte-exact")
    expected_markdown = render_markdown(expected).encode("utf-8")
    if markdown_path.read_bytes() != expected_markdown:
        raise CertificateError("MARKDOWN_BYTES", "generated status page is stale")
    return actual


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build or validate the comparison-value-free issue-552 status packet."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="write the canonical outputs")
    mode.add_argument(
        "--validate-only",
        action="store_true",
        help="validate the committed outputs without writing",
    )
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    parser.add_argument("--markdown-out", type=Path, default=DEFAULT_MD_OUT)
    args = parser.parse_args()

    try:
        if args.write:
            projection = reconstruct_source_projection()
            PROJECTION_PATH.parent.mkdir(parents=True, exist_ok=True)
            PROJECTION_PATH.write_bytes(_pretty_bytes(projection))
            payload = build_payload()
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_bytes(_pretty_bytes(payload))
            args.markdown_out.write_text(render_markdown(payload), encoding="utf-8")
            print(f"wrote {PROJECTION_PATH.relative_to(REPO_ROOT)}")
            print(f"wrote {args.json_out.relative_to(REPO_ROOT)}")
            print(f"wrote {args.markdown_out.relative_to(REPO_ROOT)}")
        else:
            validate_committed(args.json_out, args.markdown_out)
        print("QUANTUM_CARRIER_STATUS_VALID")
    except (CertificateError, OSError) as exc:
        print(f"QUANTUM_CARRIER_STATUS_INVALID: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
