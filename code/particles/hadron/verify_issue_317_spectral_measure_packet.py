#!/usr/bin/env python3
"""Verifier for OPH issue #317: conditional Ward-projected spectral-measure theorem and contract packet.

Scope (re-scoped per the PR #597 review): this packet does NOT construct a
physical source-only QCD ensemble and does NOT numerically produce the
physical hadronic spectral measure. Those production objects belong to the
production backend lane (issue #425). This packet defines and verifies the
conditional mathematical construction and the machine-readable contract
that the production export must satisfy:

    SpectralMeasure_Q : ValidQCDSourcePacket -> WardProjectedSpectralMeasure

The theorem is genuinely conditional: ten explicitly typed premises (gauge
quotient ensemble, conserved current, Ward certificate, reflection-positive
transfer certificate, correlators with joint covariance, normalization
conventions, level/residue records, threshold typing, subtraction
convention, and no-measured-HVP provenance) and eight conclusions, each
tied to the premises it uses. The premises are specified and validated in
representation by the strict contract validator; their physical QCD values
are not produced here, and the certificate says so.

Executed implementation witnesses (they test the implementation; they do
not replace the reflection-positivity, transfer, continuum, or
physical-source premises):

1. Gauge witness: both the local-local transverse correlator and the
   conserved(sink)-local(source) transverse correlator - the latter
   containing the actual link-inserted point-split conserved current - are
   recomputed on a random SU(3) gauge rotation of a rough interacting
   background and must match, together with the plaquette.
2. Ward witness: the exact backward-difference divergence of the conserved
   point-split current vanishes off the source at solver tolerance.
3. Normalization witness: the conserved-local estimator
   Z_V^eff = C_CL/(2 kappa C_LL) equals 1 on the free field (convention
   anchor; not an interacting plateau determination).
4. Diagnostic demonstrator: the constructor runs end-to-end on the
   committed real quenched ensemble and emits a machine-readable
   single-atom object with correlator covariance and the joint jackknife
   covariance of the extracted observables. Its correlator sign check is
   explicitly labeled "no statistically significant sign violation
   detected"; it is NOT a proof of spectral positivity and is not used as
   one anywhere in the certificate.
5. Acceptance gate: the strict fail-closed validator
   (ward_projected_spectral_measure_validator.py) accepts a conformant
   synthetic payload and rejects an adversarial battery covering every
   declared semantic requirement of the contract.
6. Downstream compatibility: every gate-approved payload variant embeds as
   source_measure into the downstream source-transport validator
   (code/P_derivation/thomson_spectral_transport.py) with zero reasons, so
   the gate defines a contract the stated consumer actually accepts.
7. Typing checks: empirical comparison data confined to the non-source,
   non-promoting row class; higher-point and radiative objects separately
   typed (live artifact reads).
8. Fail-closed probes: the physical-availability predicate is executed on
   synthetic unknown/renamed/missing/blocked status combinations and must
   report unavailable on every one.

The certificate separates three verdicts, computed independently:

- contract_certified: the theorem-and-contract layer above, all
  machine-checked here.
- physical_source_payload_available: computed fail-closed from explicit
  success-state allowlists over the live production-lane artifacts plus
  the requirement that an actual gate-approved production payload exists;
  currently false.
- physical_promotion_allowed: requires the payload to be available and a
  production source-transport payload accepted by the full downstream
  validator; currently false.

Run:
    python3 code/particles/hadron/verify_issue_317_spectral_measure_packet.py \
        --check --output \
        code/particles/runs/hadron/ward_projected_spectral_measure_proof_packet.json
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lattice_backend.core import (  # noqa: E402
    average_plaquette,
    cold_start,
    gauge_transform,
    random_su3_near_identity,
    sweep,
)
from lattice_backend.dirac import WilsonClover, point_propagator  # noqa: E402
from lattice_backend.conserved_vector import (  # noqa: E402
    conserved_local_correlator,
    transverse_vector_correlator,
    ward_divergence_offsource_max,
    zv_effective,
)
from lattice_backend.vector_correlator import fold_correlator  # noqa: E402

import ward_projected_spectral_measure_validator as strict_validator  # noqa: E402

ROOT = HERE.parents[1]
P_DERIVATION = ROOT / "P_derivation"
if str(P_DERIVATION) not in sys.path:
    sys.path.insert(0, str(P_DERIVATION))

from thomson_spectral_transport import (  # noqa: E402
    _validate_source_measure,
    source_payload_forbidden_keys,
    validate_source_transport_payload,
)

SCHEMA_PATH = HERE / "ward_projected_spectral_measure.schema.json"
EMPIRICAL_SCHEMA_PATH = HERE / "empirical_ward_projected_spectral_measure.schema.json"
ENSEMBLE_NPZ = ROOT / "particles" / "runs" / "hadron" / "hybrid_ir_ensembleA_2026-07-16.npz"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_proof_packet.json"
BACKEND_DIR = ROOT / "particles" / "runs" / "qcd" / "hadron_source_backend"
COMPARISON_MANIFEST = BACKEND_DIR / "controls" / "comparison_data_manifest.json"
Q4_RECEIPT = BACKEND_DIR / "higher_point" / "Q4_HLbL_receipt.json"
XI_LEDGER = BACKEND_DIR / "endpoint" / "Xi_same_scheme.json"
BASE_MEASURE = BACKEND_DIR / "qcd_ensemble" / "base_measure.json"
WARD_CURRENT_DEFINITION = BACKEND_DIR / "currents" / "ward_current_definition.json"
READINESS_REPORT = ROOT / "particles" / "runs" / "hadron" / "hadron_production_readiness_report.json"

# Declared production artifact locations. The physical availability verdict
# requires an actual gate-approved payload at these paths; absence fails
# closed.
PRODUCTION_PAYLOAD_PATH = (
    ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure.production.json"
)
PRODUCTION_TRANSPORT_PAYLOAD_PATH = (
    ROOT / "P_derivation" / "runtime" / "source_transport_payload.production.json"
)

FORBIDDEN_TARGETS = list(strict_validator.FORBIDDEN_TARGETS)

GAUGE_INVARIANCE_TOL = 1e-8
WARD_RELATIVE_TOL = 1e-9
ZV_FREE_FIELD_TOL = 5e-3
EFFMASS_WINDOW_START = 3

# Explicit success-state allowlists for the physical production lane. Any
# status not literally in its allowlist - including unknown, renamed,
# missing, None, or blocked values - fails closed.
STATUS_SUCCESS_ALLOWLISTS: dict[str, tuple[str, ...]] = {
    "export_bundle_status": ("complete",),
    "closure_grade": ("execution_complete",),
    "base_measure_status": ("POPULATED",),
    "ward_current_status": ("SOURCE_CERTIFICATE_VERIFIED",),
}


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rel(path: Path) -> str:
    return path.relative_to(ROOT.parent).as_posix()


# ---------------------------------------------------------------------------
# Conditional construction theorem: typed premises -> typed conclusions.
# ---------------------------------------------------------------------------


def build_theorem() -> dict[str, Any]:
    """The conditional theorem SpectralMeasure_Q with typed premises and conclusions."""
    premises = [
        {
            "id": "P1",
            "name": "gauge_invariant_quotient_ensemble_and_action",
            "statement": (
                "a declared gauge-invariant quotient ensemble with a positive Gibbs "
                "source law at zero density and theta = 0 (a normal form is not a "
                "probability law; source-law non-substitution theorem)"
            ),
            "contract_fields": ["provenance.source_inputs", "branch.flavors", "backend"],
            "physical_population_status": "REQUIRED_NOT_POPULATED (qcd_ensemble/base_measure.json)",
        },
        {
            "id": "P2",
            "name": "conserved_current_or_certified_renormalization_map",
            "statement": (
                "a conserved electromagnetic hadronic current, or a local current "
                "with a certified renormalization map to the conserved one"
            ),
            "contract_fields": ["ward_projected_residues[].current_normalization"],
            "physical_population_status": (
                "MISSING_SOURCE_CERTIFICATE (currents/ward_current_definition.json)"
            ),
        },
        {
            "id": "P3",
            "name": "ward_identity_certificate",
            "statement": "a Ward-identity certificate for the declared current",
            "contract_fields": ["projection.ward_projected"],
            "physical_population_status": (
                "MISSING_SOURCE_CERTIFICATE (currents/ward_current_definition.json)"
            ),
        },
        {
            "id": "P4",
            "name": "reflection_positive_state_and_transfer_certificate",
            "statement": (
                "a reflection-positive Euclidean state and transfer-operator "
                "certificate; reflection positivity is a theorem for the unimproved "
                "Wilson action (Osterwalder-Seiler) and is not available at finite "
                "spacing for the clover-improved engine, so this premise is never "
                "supplied by the diagnostic demonstrator"
            ),
            "contract_fields": ["rho_had_or_measure.positivity_status"],
            "physical_population_status": "not populated (vacuum/ transfer certificates absent)",
        },
        {
            "id": "P5",
            "name": "finite_volume_correlators_with_joint_covariance",
            "statement": (
                "finite-volume current correlators and their joint covariance from "
                "one common resampling pipeline"
            ),
            "contract_fields": ["covariance.row_basis", "covariance.matrix"],
            "physical_population_status": "not populated (production correlator arrays absent)",
        },
        {
            "id": "P6",
            "name": "current_normalization_conventions",
            "statement": (
                "current-normalization conventions: Z_V = 1 for the exactly "
                "conserved current, the declared matching certificate for a local one"
            ),
            "contract_fields": ["ward_projected_residues[].current_normalization"],
            "physical_population_status": "convention declared; production matching not populated",
        },
        {
            "id": "P7",
            "name": "finite_volume_level_and_residue_records",
            "statement": "finite-volume level and residue records with identifier integrity",
            "contract_fields": ["finite_volume_levels", "ward_projected_residues"],
            "physical_population_status": "not populated (production levels absent)",
        },
        {
            "id": "P8",
            "name": "threshold_and_channel_typing",
            "statement": "threshold and channel typing for the declared vector channel",
            "contract_fields": ["finite_volume_levels[].channel"],
            "physical_population_status": "typed; production spectrum not populated",
        },
        {
            "id": "P9",
            "name": "declared_subtraction_and_transport_convention",
            "statement": (
                "the declared once-subtracted transport kernel "
                "mZ(P)^2/(3 pi s (s + mZ(P)^2)) with the same-subtraction lock to a0(P)"
            ),
            "contract_fields": ["transport_moment_certificate.kernel"],
            "physical_population_status": "convention declared; production certificate not populated",
        },
        {
            "id": "P10",
            "name": "no_measured_hvp_or_target_provenance",
            "statement": (
                "provenance proving that no measured HVP or hadronic target entered "
                "the source construction"
            ),
            "contract_fields": ["provenance", "external_targets_used"],
            "physical_population_status": "enforced by the strict validator on any payload",
        },
    ]
    for premise in premises:
        premise["premise_type"] = "declared_premise_validated_in_representation_only"
        premise["physical_value_produced_here"] = False

    conclusions = [
        {
            "id": "C1",
            "name": "current_correlator_gauge_invariant",
            "statement": (
                "each hop term psibar(x+mu)(1 +- gamma_mu) U_mu(x)^dag psi(x) of the "
                "point-split current is a color trace of gauge-covariant bilinears "
                "with the compensating link insertion, hence the current correlator "
                "is exactly gauge invariant"
            ),
            "uses": ["P1", "P2"],
            "implementation_witness": "machine_witnesses.gauge_and_ward.conserved_current_gauge_invariance",
        },
        {
            "id": "C2",
            "name": "ward_projection_transverse",
            "statement": (
                "the current is the exact Noether current of the U(1) phase rotation "
                "of the Wilson-clover action (the clover term is phase-invariant), so "
                "q_mu Pi_Q^{mu nu} = 0 and the two-current correlator has the "
                "transverse form Pi_Q^{mu nu} = (q^mu q^nu - q^2 eta^{mu nu}) Pi_Q(q^2)"
            ),
            "uses": ["P2", "P3"],
            "implementation_witness": "machine_witnesses.gauge_and_ward.ward_identity",
        },
        {
            "id": "C3",
            "name": "positive_locally_finite_spectral_measure",
            "statement": (
                "given the reflection-positive transfer certificate, the spectral "
                "decomposition gives G(t) = sum_n w_n (e^{-E_n t} + e^{-E_n (T-t)}) "
                "with w_n proportional to |<0|J_Q|n>|^2 >= 0, so rho_Q(s) = "
                "sum_n Z_n delta(s - E_n^2) is a positive locally finite measure; "
                "without P4 this conclusion is not asserted, and no finite-sample "
                "sign check substitutes for it"
            ),
            "uses": ["P4", "P5", "P7"],
            "implementation_witness": None,
        },
        {
            "id": "C4",
            "name": "threshold_support",
            "statement": (
                "the support of rho_Q begins at the lowest state in the declared "
                "vector channel above the vacuum (two-pion channel in infinite volume)"
            ),
            "uses": ["P4", "P7", "P8"],
            "implementation_witness": None,
        },
        {
            "id": "C5",
            "name": "normalization",
            "statement": (
                "the exactly conserved current has Z_V = 1; a local current is "
                "renormalized through the declared matching certificate "
                "Z_V^eff = C_CL/(2 kappa C_LL), anchored at exactly 1 on the free field"
            ),
            "uses": ["P2", "P6"],
            "implementation_witness": "machine_witnesses.zv_free_field_anchor",
        },
        {
            "id": "C6",
            "name": "scheme_confinement",
            "statement": (
                "the spectral measure is scheme-independent spectral data; the "
                "endpoint subtraction convention does not change the measure, and "
                "scheme dependence is confined to the declared transport kernel and "
                "matching remainder (the kernel identity Delta_had = 4 pi Pihat(mZ^2) "
                "is exact; lattice_backend/vector_correlator.py)"
            ),
            "uses": ["P9"],
            "implementation_witness": None,
        },
        {
            "id": "C7",
            "name": "covariance_propagation",
            "statement": (
                "linear spectral moments and nonlinear downstream quantities carry "
                "covariance propagated from one common resampling pipeline; the "
                "contract requires the covariance block and the seven-row budget "
                "ledger with finite ordered bound intervals"
            ),
            "uses": ["P5"],
            "implementation_witness": "machine_witnesses.demonstrator_measure.joint_extraction_covariance",
        },
        {
            "id": "C8",
            "name": "higher_point_and_radiative_separately_typed",
            "statement": (
                "higher-point (four-current, HLbL-class) and radiative (QED/EW) "
                "objects remain separately typed artifacts and cannot be substituted "
                "by the two-current measure"
            ),
            "uses": ["P8", "P10"],
            "implementation_witness": "machine_witnesses.higher_point_typing",
        },
    ]

    return {
        "theorem_id": "SpectralMeasureQ_ConditionalConstruction",
        "form": "SpectralMeasure_Q : ValidQCDSourcePacket -> WardProjectedSpectralMeasure",
        "statement": (
            "Given a source packet satisfying the declared gauge, current, "
            "reflection-positivity, transfer, normalization, threshold, and "
            "covariance premises (P1-P10), the corresponding Ward-projected "
            "positive spectral measure exists and has the stated structural "
            "properties (C1-C8). This issue specifies the premises and validates "
            "their representation; it does not produce their physical QCD values."
        ),
        "premises": premises,
        "conclusions": conclusions,
        "conditionality": (
            "every conclusion holds only under the premises it lists; the theorem "
            "remains explicitly conditional whenever reflection positivity (P4), "
            "continuum existence, the scattering pushforward, or any other "
            "physical premise is absent. Diagnostic or quenched data exercise the "
            "code paths but cannot satisfy the physical premises or promote a "
            "hadronic number."
        ),
    }


def theorem_structure_check(theorem: dict[str, Any]) -> dict[str, Any]:
    """Machine check that the theorem is stated with typed premises and conclusions."""
    premises = theorem.get("premises", [])
    conclusions = theorem.get("conclusions", [])
    premise_ids = [p.get("id") for p in premises]
    checks = {
        "ten_premises_typed": len(premises) == 10
        and all(
            p.get("premise_type") == "declared_premise_validated_in_representation_only"
            and p.get("physical_value_produced_here") is False
            and p.get("statement")
            and p.get("contract_fields")
            for p in premises
        ),
        "eight_conclusions_typed": len(conclusions) == 8
        and all(c.get("statement") and c.get("uses") for c in conclusions),
        "conclusion_uses_reference_declared_premises": all(
            set(c.get("uses", [])) <= set(premise_ids) and c.get("uses")
            for c in conclusions
        ),
        "positivity_conclusion_requires_reflection_positivity_premise": any(
            c["id"] == "C3" and "P4" in c["uses"] for c in conclusions
        ),
        "threshold_conclusion_requires_reflection_positivity_premise": any(
            c["id"] == "C4" and "P4" in c["uses"] for c in conclusions
        ),
        "explicit_conditionality_statement": "conditional" in theorem.get("conditionality", ""),
    }
    return {"checks": checks, "passed": bool(all(checks.values()))}


# ---------------------------------------------------------------------------
# Live artifact checks: typing criteria.
# ---------------------------------------------------------------------------


def empirical_typing_check(gate: dict[str, Any]) -> dict[str, Any]:
    """Criterion: empirical comparison data confined to a non-source, non-promoting row class.

    Machine-checked in three layers against live artifacts:

    1. Blinding of the source construction: no empirical dataset is attached
       to the source lane (comparison manifest empty), and the source export
       cannot admit empirical input (the strict gate fails closed on every
       measured-HVP path; executed in the acceptance-gate witness).
    2. Comparison channel typing: the only surface where empirical data may
       ever attach to this lane is the declared-empirical companion row
       class, pinned non-source by schema constants.
    3. Disclosure: the empirical companion surface allows public final
       values (usable_for_public_final_values const true). That use is
       confined to the disclosed oph_plus_empirical_hadron_closure row
       class, which by the same constants can never be promoted as an OPH
       source theorem and never satisfies the production constructive
       artifact. The check records the disclosure explicitly.
    """
    manifest = _load_json(COMPARISON_MANIFEST)
    schema = _load_json(EMPIRICAL_SCHEMA_PATH)
    production_schema = _load_json(SCHEMA_PATH)
    guards = schema["properties"]["guards"]["properties"]
    production_guards = production_schema["properties"]["guards"]["properties"]
    gate_rejections = {
        row["control_id"]: row["rejected"] for row in gate["negative_controls"]
    }
    checks = {
        "comparison_manifest_empty": manifest.get("comparison_data") == [],
        "comparison_manifest_status_no_data": manifest.get("status") == "NO_COMPARISON_DATA_ATTACHED",
        "production_schema_rejects_compare_only_const_false": (
            production_guards["compare_only_external_endpoint"].get("const") is False
        ),
        "production_schema_rejects_surrogate_const_false": (
            production_guards["surrogate_hadron_artifact"].get("const") is False
        ),
        "gate_rejects_measured_hvp_comparison_endpoint": bool(
            gate_rejections.get("measured_hvp_comparison_endpoint")
        ),
        "gate_rejects_ee_to_hadrons_target_leak": bool(
            gate_rejections.get("forbidden_target_leak_ee_to_hadrons")
        ),
        "gate_rejects_nested_measured_hvp_provenance": bool(
            gate_rejections.get("nested_measured_hvp_provenance")
        ),
        "row_class_const_empirical_closure": (
            schema["properties"]["row_class"].get("const") == "oph_plus_empirical_hadron_closure"
        ),
        "empirical_surface_not_source_only_const": guards["source_only"].get("const") is False,
        "empirical_use_disclosed_const_true": (
            guards["external_cross_section_data_used"].get("const") is True
        ),
        "promotable_as_oph_source_theorem_const_false": (
            guards["promotable_as_oph_source_theorem"].get("const") is False
        ),
        "satisfies_production_constructive_next_artifact_const_false": (
            guards["satisfies_production_constructive_next_artifact"].get("const") is False
        ),
    }
    public_use_disclosure = {
        "usable_for_public_final_values_const": guards["usable_for_public_final_values"].get("const"),
        "confinement": (
            "public use is a property of the disclosed empirical-closure row class "
            "only; the same schema constants forbid promoting that surface as an "
            "OPH source theorem or substituting it for the production artifact, "
            "so with respect to the source spectral measure the empirical data "
            "remain comparison-only"
        ),
    }
    blinding_status = {
        "comparison_event_performed": bool(manifest.get("comparison_data")),
        "statement": (
            "no empirical dataset is attached to the source lane, so no comparison "
            "has occurred and there is nothing unblinded; the comparison event is "
            "production-scope and must append to the comparison manifest after the "
            "source payload is emitted"
        ),
    }
    return {
        "checked_artifacts": {
            "comparison_manifest": _rel(COMPARISON_MANIFEST),
            "empirical_companion_schema": _rel(EMPIRICAL_SCHEMA_PATH),
            "production_export_schema": _rel(SCHEMA_PATH),
        },
        "checks": checks,
        "public_use_disclosure": public_use_disclosure,
        "blinding_status": blinding_status,
        "passed": bool(all(checks.values())),
    }


def higher_point_typing_check() -> dict[str, Any]:
    """Criterion: higher-point and radiative corrections stay separately typed.

    Machine-checked for both halves against live artifacts:

    - higher-point half: the four-current receipt records the two-point
      measure as insufficient for HLbL-class objects, non-promoting, with no
      external targets;
    - radiative half: the same-scheme remainder ledger Xi_same_scheme is a
      distinct typed artifact whose required term rows include the QED and
      EW corrections, non-promoting, with no external targets; the
      production export schema additionally excludes QED from the two-point
      branch (branch.qed limited to off/explicitly_superseded).

    The Xi ledger's source certificate is itself still missing
    (MISSING_SOURCE_CERTIFICATE); that is a production-level content gap
    recorded explicitly. The typing separation checked here is exactly that
    the radiative rows live in that distinct ledger rather than inside the
    two-point measure.
    """
    receipt = _load_json(Q4_RECEIPT)
    xi = _load_json(XI_LEDGER)
    production_schema = _load_json(SCHEMA_PATH)
    qed_enum = (
        production_schema["properties"]["branch"]["properties"]["qed"].get("enum") or []
    )
    checks = {
        "q4_status_two_point_measure_insufficient": (
            receipt.get("status") == "TWO_POINT_MEASURE_INSUFFICIENT"
        ),
        "q4_promotion_not_allowed": receipt.get("promotion_allowed") is False,
        "q4_no_external_targets": receipt.get("external_targets_used") == [],
        "xi_ledger_distinct_artifact": xi.get("artifact") == "Xi_same_scheme",
        "xi_required_terms_include_qed_and_ew": (
            {"QED", "EW"} <= set(xi.get("required_terms", []))
        ),
        "xi_promotion_not_allowed": xi.get("promotion_allowed") is False,
        "xi_no_external_targets": xi.get("external_targets_used") == [],
        "production_branch_excludes_qed": "on" not in qed_enum and "off" in qed_enum,
    }
    return {
        "checked_artifacts": {
            "q4_hlbl_receipt": _rel(Q4_RECEIPT),
            "xi_same_scheme_ledger": _rel(XI_LEDGER),
            "production_export_schema": _rel(SCHEMA_PATH),
        },
        "checks": checks,
        "xi_ledger_source_certificate_status": xi.get("status"),
        "passed": bool(all(checks.values())),
    }


# ---------------------------------------------------------------------------
# Physical production lane: fail-closed availability and promotion verdicts.
# ---------------------------------------------------------------------------


def read_live_production_state() -> dict[str, Any]:
    """Read the live production-lane statuses and payload locations."""
    readiness = _load_json(READINESS_REPORT)
    base_measure = _load_json(BASE_MEASURE)
    ward_current = _load_json(WARD_CURRENT_DEFINITION)
    remaining = readiness.get("exact_remaining_runtime_object", {})
    closure = readiness.get("closure_status", {})
    return {
        "export_bundle_status": remaining.get("status"),
        "closure_grade": closure.get("closure_grade"),
        "public_unsuppression_ready": closure.get("public_unsuppression_ready"),
        "base_measure_status": base_measure.get("status"),
        "ward_current_status": ward_current.get("status"),
        "production_payload_path": PRODUCTION_PAYLOAD_PATH,
        "checked_artifacts": {
            "production_backend_export_bundle": _rel(READINESS_REPORT),
            "qcd_base_measure": _rel(BASE_MEASURE),
            "ward_current_source_certificate": _rel(WARD_CURRENT_DEFINITION),
        },
    }


def physical_source_payload_verdict(state: dict[str, Any] | None = None) -> dict[str, Any]:
    """Fail-closed availability decision for the physical production payload.

    Success is decided ONLY by explicit allowlists and by the existence of a
    gate-approved production payload; unknown, renamed, missing, malformed,
    blocked, or contradictory statuses all fail closed with typed reasons.
    """
    state = state if state is not None else read_live_production_state()
    reasons: list[str] = []
    for key, allowlist in STATUS_SUCCESS_ALLOWLISTS.items():
        value = state.get(key)
        if value not in allowlist:
            reasons.append(f"status_not_in_success_allowlist:{key}={value!r}")
    if state.get("public_unsuppression_ready") is not True:
        reasons.append(
            "public_unsuppression_ready_not_true:"
            f"{state.get('public_unsuppression_ready')!r}"
        )
    payload_path = state.get("production_payload_path")
    if not isinstance(payload_path, Path) or not payload_path.exists():
        reasons.append(f"production_payload_absent:{Path(str(payload_path)).name}")
    else:
        try:
            payload = _load_json(payload_path)
        except (json.JSONDecodeError, OSError):
            payload = None
        if payload is None:
            reasons.append("production_payload_unreadable")
        else:
            result = strict_validator.validate_production_payload(payload)
            if not result.accepted:
                reasons.append(
                    "production_payload_rejected_by_strict_gate:"
                    + ";".join(result.reasons[:5])
                )
    return {
        "available": not reasons,
        "reasons": reasons,
        "success_allowlists": {k: list(v) for k, v in STATUS_SUCCESS_ALLOWLISTS.items()},
        "additional_positive_requirements": [
            "public_unsuppression_ready is exactly true",
            "a production payload exists at "
            + PRODUCTION_PAYLOAD_PATH.relative_to(ROOT.parent).as_posix()
            + " and passes the strict production gate",
        ],
        "live_status_snapshot": {
            key: state.get(key)
            for key in (
                "export_bundle_status",
                "closure_grade",
                "public_unsuppression_ready",
                "base_measure_status",
                "ward_current_status",
            )
        },
        "checked_artifacts": state.get("checked_artifacts", {}),
    }


def physical_promotion_verdict(payload_verdict: dict[str, Any]) -> dict[str, Any]:
    """Fail-closed promotion decision; requires availability plus the downstream certificate."""
    reasons: list[str] = []
    if not payload_verdict["available"]:
        reasons.append("physical_source_payload_not_available")
    if not PRODUCTION_TRANSPORT_PAYLOAD_PATH.exists():
        reasons.append(
            "production_source_transport_payload_absent:"
            + PRODUCTION_TRANSPORT_PAYLOAD_PATH.name
        )
    else:
        transport = _load_json(PRODUCTION_TRANSPORT_PAYLOAD_PATH)
        validation = validate_source_transport_payload(transport)
        if not validation.promotion_allowed:
            reasons.append(
                "downstream_source_transport_validator_blocked:"
                + ";".join(validation.reasons[:5])
            )
    return {
        "allowed": not reasons,
        "reasons": reasons,
        "positive_requirements": [
            "physical_source_payload_available is true",
            "a production source-transport payload exists at "
            + PRODUCTION_TRANSPORT_PAYLOAD_PATH.relative_to(ROOT.parent).as_posix()
            + " and is accepted by validate_source_transport_payload "
            "(code/P_derivation/thomson_spectral_transport.py)",
        ],
    }


def fail_closed_probe_battery() -> dict[str, Any]:
    """Execute the availability predicate on adversarial status combinations.

    Every probe uses a state that is NOT in the success allowlists (unknown,
    renamed, missing, None, blocked, or contradictory); the predicate must
    report unavailable on all of them.
    """
    all_success = {
        "export_bundle_status": "complete",
        "closure_grade": "execution_complete",
        "public_unsuppression_ready": True,
        "base_measure_status": "POPULATED",
        "ward_current_status": "SOURCE_CERTIFICATE_VERIFIED",
        "production_payload_path": PRODUCTION_PAYLOAD_PATH,
    }
    probes: list[tuple[str, dict[str, Any]]] = [
        ("status_blocked", {**all_success, "export_bundle_status": "blocked"}),
        ("status_missing_string", {**all_success, "base_measure_status": "MISSING"}),
        ("status_none", {**all_success, "ward_current_status": None}),
        ("status_renamed", {**all_success, "export_bundle_status": "done"}),
        ("closure_grade_unknown", {**all_success, "closure_grade": "complete_pending_review"}),
        ("public_unsuppression_not_ready", {**all_success, "public_unsuppression_ready": False}),
        ("public_unsuppression_stringly_true", {**all_success, "public_unsuppression_ready": "true"}),
        ("all_keys_missing", {}),
        (
            "all_statuses_success_but_payload_absent",
            dict(all_success),
        ),
    ]
    rows = []
    all_failed_closed = True
    for probe_id, state in probes:
        verdict = physical_source_payload_verdict(state)
        rows.append(
            {
                "probe_id": probe_id,
                "reported_unavailable": not verdict["available"],
                "reasons": verdict["reasons"][:3],
            }
        )
        all_failed_closed = all_failed_closed and not verdict["available"]
    return {
        "statement": (
            "the availability predicate is executed on unknown, renamed, missing, "
            "None, blocked, and contradictory status combinations and must report "
            "unavailable on every one; success requires literal allowlist membership "
            "plus an existing gate-approved payload"
        ),
        "probes": rows,
        "passed": bool(all_failed_closed),
    }


# ---------------------------------------------------------------------------
# Executed lattice witnesses (small, seeded, deterministic; non-promoting).
# ---------------------------------------------------------------------------


def gauge_and_ward_witness(
    shape: tuple[int, int, int, int] = (8, 2, 2, 2),
    beta: float = 5.5,
    kappa: float = 0.12,
    n_sweeps: int = 6,
    cg_tol: float = 1e-12,
    seed: int = 3170,
) -> dict[str, Any]:
    """Gauge invariance (local-local AND conserved-local) and the exact Ward identity.

    The conserved-local comparison exercises the actual link-inserted
    point-split conserved current whose gauge invariance theorem part C1
    claims; the local-local comparison covers the source-side local current.
    """
    rng = np.random.default_rng(seed)
    u = cold_start(shape)
    for _ in range(n_sweeps):
        sweep(rng, u, beta, n_or=1)
    op = WilsonClover(u, kappa=kappa, c_sw=1.0)
    prop, _ = point_propagator(op, shape, tol=cg_tol)
    g_ll = transverse_vector_correlator(prop)
    g_cl = conserved_local_correlator(prop, op.ubc, kappa)
    plaq = average_plaquette(u)

    g_field = random_su3_near_identity(rng, shape, eps=1.0)
    u_rot = gauge_transform(u, g_field)
    op_rot = WilsonClover(u_rot, kappa=kappa, c_sw=1.0)
    prop_rot, _ = point_propagator(op_rot, shape, tol=cg_tol)
    g_ll_rot = transverse_vector_correlator(prop_rot)
    g_cl_rot = conserved_local_correlator(prop_rot, op_rot.ubc, kappa)
    plaq_rot = average_plaquette(u_rot)

    ll_scale = float(np.max(np.abs(g_ll)))
    ll_defect = float(np.max(np.abs(g_ll_rot - g_ll))) / ll_scale
    cl_scale = float(np.max(np.abs(g_cl)))
    cl_defect = float(np.max(np.abs(g_cl_rot - g_cl))) / cl_scale
    plaq_defect = abs(plaq_rot - plaq)

    ward_defect, ward_scale = ward_divergence_offsource_max(prop, op.ubc, kappa, nu=1)
    ward_relative = ward_defect / ward_scale

    return {
        "non_promoting_diagnostic": True,
        "background": {
            "lattice_shape_TXYZ": list(shape),
            "beta": beta,
            "kappa": kappa,
            "n_sweeps": n_sweeps,
            "cg_tol": cg_tol,
            "rng_seed": seed,
            "plaquette": plaq,
        },
        "local_current_gauge_invariance": {
            "observable": (
                "local-local transverse correlator "
                "sum_x <psibar gamma_k psi(x) psibar gamma_k psi(0)>, k = 1..3"
            ),
            "statement": (
                "the local-local transverse U(1) vector correlator and the plaquette "
                "are invariant under U_mu(x) -> g(x) U_mu(x) g(x+mu)^dag for random SU(3) g"
            ),
            "correlator_relative_defect": ll_defect,
            "plaquette_absolute_defect": plaq_defect,
            "tolerance": GAUGE_INVARIANCE_TOL,
            "passed": bool(ll_defect < GAUGE_INVARIANCE_TOL and plaq_defect < GAUGE_INVARIANCE_TOL),
        },
        "conserved_current_gauge_invariance": {
            "observable": (
                "conserved(sink)-local(source) transverse correlator C_CL(t): the sink "
                "operator is the actual link-inserted point-split conserved current "
                "V_mu(x) = kappa [psibar(x+mu)(1+gamma_mu) U_mu(x)^dag psi(x) - "
                "psibar(x)(1-gamma_mu) U_mu(x) psi(x+mu)], mu = 1..3"
            ),
            "statement": (
                "the correlator containing the link-inserted point-split conserved "
                "current is invariant under the same random SU(3) gauge rotation; this "
                "is the observable whose gauge invariance conclusion C1 asserts"
            ),
            "correlator_relative_defect": cl_defect,
            "tolerance": GAUGE_INVARIANCE_TOL,
            "passed": bool(cl_defect < GAUGE_INVARIANCE_TOL),
        },
        "ward_identity": {
            "observable": (
                "backward-difference divergence of <V_mu^cons(x) J_nu^loc(0)> over all "
                "off-source sites"
            ),
            "statement": (
                "the conserved point-split current V_mu satisfies the exact backward-difference "
                "Ward identity sum_mu [V_mu(x) - V_mu(x-mu)] = 0 off the contact point"
            ),
            "max_offsource_divergence": ward_defect,
            "correlator_scale": ward_scale,
            "relative_defect": ward_relative,
            "tolerance": WARD_RELATIVE_TOL,
            "passed": bool(ward_relative < WARD_RELATIVE_TOL),
        },
    }


def zv_free_field_witness(
    shape: tuple[int, int, int, int] = (16, 2, 2, 2),
    kappa: float = 0.13,
    cg_tol: float = 1e-11,
) -> dict[str, Any]:
    """Free-field anchor: Z_V^eff = C_CL/(2 kappa C_LL) equals 1 on the free field.

    This anchors the sign and normalization convention (conclusion C5). It
    is a free-field convention check, not an interacting plateau
    determination; the interacting normalization certificate is premise P6
    content for the production lane.
    """
    u = cold_start(shape)
    op = WilsonClover(u, kappa=kappa, c_sw=1.0)
    prop, _ = point_propagator(op, shape, tol=cg_tol)
    g_ll = transverse_vector_correlator(prop)
    g_cl = conserved_local_correlator(prop, op.ubc, kappa)
    zv = zv_effective(g_cl, g_ll, kappa)
    plateau = zv[4:8]
    max_dev = float(np.max(np.abs(plateau - 1.0)))
    return {
        "non_promoting_diagnostic": True,
        "lattice_shape_TXYZ": list(shape),
        "kappa": kappa,
        "statement": (
            "on the free field the local-current renormalization estimator "
            "Z_V^eff(t) = C_CL(t)/(2 kappa C_LL(t)) equals 1, anchoring sign and "
            "normalization convention; not an interacting plateau determination"
        ),
        "window_t": [4, 7],
        "window_values": [float(v) for v in plateau],
        "max_abs_deviation_from_one": max_dev,
        "tolerance": ZV_FREE_FIELD_TOL,
        "passed": bool(max_dev < ZV_FREE_FIELD_TOL),
    }


# ---------------------------------------------------------------------------
# Diagnostic demonstrator from the committed real quenched ensemble.
# ---------------------------------------------------------------------------


def _cosh_effective_mass(t: int, ratio: float, half: int) -> float:
    """Solve cosh(m*(T/2 - t))/cosh(m*(T/2 - t - 1)) = ratio by bisection."""
    if ratio <= 1.0:
        return float("nan")
    lo, hi = 1e-8, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        val = math.cosh(mid * (half - t)) / math.cosh(mid * (half - t - 1))
        if val > ratio:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _effmass_window(folded: np.ndarray, start: int = EFFMASS_WINDOW_START) -> list[int]:
    half = len(folded) - 1
    ts: list[int] = []
    for t in range(start, half):
        if folded[t] > 0.0 and folded[t + 1] > 0.0 and folded[t] > folded[t + 1]:
            ts.append(t)
        elif ts:
            break
    return ts


def _windowed_mass(folded: np.ndarray, start: int = EFFMASS_WINDOW_START) -> float:
    half = len(folded) - 1
    vals = [
        _cosh_effective_mass(t, folded[t] / folded[t + 1], half)
        for t in _effmass_window(folded, start)
    ]
    vals = [v for v in vals if math.isfinite(v)]
    return float(np.mean(vals)) if vals else float("nan")


def _windowed_zv(folded_cl: np.ndarray, folded_ll: np.ndarray, kappa: float) -> float:
    """Windowed average of Z_V^eff(t); a diagnostic, not a validated plateau."""
    half = len(folded_ll) - 1
    vals = [
        float(folded_cl[t] / (2.0 * kappa * folded_ll[t]))
        for t in range(EFFMASS_WINDOW_START, half)
        if folded_ll[t] != 0.0
    ]
    return float(np.mean(vals))


def _atom_pipeline(
    ll_mean: np.ndarray, cl_mean: np.ndarray, kappa: float, t_extent: int
) -> dict[str, float]:
    """Single-atom constructor: level, Z_V, Ward-projected residue weight.

    Single-state model for the folded correlator, f(t) = A (e^{-m t} +
    e^{-m (T-t)}).  For an R-normalized spectral atom rho_R = w delta(s - m^2)
    the amplitude is A = w m / (24 pi^2)
    (``vector_correlator.synthetic_atom_correlator``), so w_hop = 24 pi^2 A / m.
    The correlator is built from hopping-normalized propagators; the physical
    correlator is (2 kappa)^2 times the hopping one (psi_phys =
    sqrt(2 kappa) psi_hop; ``conserved_vector.py`` docstring: "the analysis
    layer applies that factor"), and the local current renormalizes with
    Z_V^2: w_phys = Z_V^2 (2 kappa)^2 w_hop.
    """
    folded_ll = fold_correlator(ll_mean)
    folded_cl = fold_correlator(cl_mean)
    m_v = _windowed_mass(folded_ll)
    z_v = _windowed_zv(folded_cl, folded_ll, kappa)
    window = _effmass_window(folded_ll)
    amps = [
        folded_ll[t] / (math.exp(-m_v * t) + math.exp(-m_v * (t_extent - t)))
        for t in window
    ]
    amp = float(np.mean(amps))
    w_hop = 24.0 * math.pi**2 * amp / m_v
    w_phys = z_v * z_v * (2.0 * kappa) ** 2 * w_hop
    return {"am_v": m_v, "z_v": z_v, "amplitude_hop": amp, "weight_phys": w_phys}


def demonstrator_measure(npz_path: Path = ENSEMBLE_NPZ) -> dict[str, Any]:
    """Machine-readable diagnostic single-atom object from the committed real ensemble.

    This is an executed implementation witness only: it demonstrates that
    the constructor pipeline runs end-to-end on real source-only lattice
    data and emits a machine-readable object with covariance. It is
    diagnostic scale, non-promoting, quenched, and is rejected by the
    production gate by construction. Its correlator sign check is not a
    proof of spectral positivity and is not used as one.
    """
    d = np.load(npz_path)
    kappa = float(d["kappa"])
    t_extent = int(d["shape"][0])
    g_ll = np.asarray(d["g_ll"])
    g_cl = np.asarray(d["g_cl"])
    n_cfg = g_ll.shape[0]

    central = _atom_pipeline(g_ll.mean(axis=0), g_cl.mean(axis=0), kappa, t_extent)

    extraction_keys = list(central)
    samples: dict[str, list[float]] = {k: [] for k in extraction_keys}
    folded_samples = []
    for i in range(n_cfg):
        ll_i = np.delete(g_ll, i, axis=0).mean(axis=0)
        cl_i = np.delete(g_cl, i, axis=0).mean(axis=0)
        res = _atom_pipeline(ll_i, cl_i, kappa, t_extent)
        for key, value in res.items():
            samples[key].append(value)
        folded_samples.append(fold_correlator(ll_i))
    errors = {}
    for key, vals in samples.items():
        arr = np.array(vals)
        center = arr.mean()
        errors[key] = float(np.sqrt((n_cfg - 1) / n_cfg * np.sum((arr - center) ** 2)))

    # Joint jackknife covariance of ALL extracted (nonlinear) observables
    # from the one common leave-one-out pipeline.
    sample_matrix = np.array([samples[key] for key in extraction_keys])
    centers = sample_matrix.mean(axis=1, keepdims=True)
    dev_extract = sample_matrix - centers
    joint_cov = (n_cfg - 1) / n_cfg * dev_extract @ dev_extract.T
    joint_eigvals = np.linalg.eigvalsh(0.5 * (joint_cov + joint_cov.T))
    joint_scale = float(np.max(np.abs(joint_cov)))
    joint_psd = bool(float(joint_eigvals.min()) >= -1e-12 * max(1.0, joint_scale))

    folded_arr = np.array(folded_samples)
    folded_center = folded_arr.mean(axis=0)
    dev = folded_arr - folded_center
    covariance = (n_cfg - 1) / n_cfg * dev.T @ dev
    eigvals = np.linalg.eigvalsh(covariance)
    cov_scale = float(np.max(np.abs(covariance)))
    cov_psd = bool(float(eigvals.min()) >= -1e-12 * cov_scale)

    # Diagnostic sign check ONLY. Pointwise positivity of a finite noisy
    # correlator does not imply a positive spectral measure (e.g.
    # G(t) = e^{-t} - 0.1 e^{-2t} > 0 for t >= 0 with a negative spectral
    # coefficient). The positive-measure conclusion C3 depends on the
    # reflection-positivity/transfer premise P4, which the finite-spacing
    # clover demonstrator does not supply.
    folded_mean = fold_correlator(g_ll.mean(axis=0))
    folded_err = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    significant = np.abs(folded_mean) > 2.0 * folded_err
    positive_where_significant = bool(np.all(folded_mean[significant] > 0.0))
    negative_only_within_noise = bool(
        np.all((folded_mean >= 0.0) | (np.abs(folded_mean) <= 2.0 * folded_err))
    )
    sign_check_passed = positive_where_significant and negative_only_within_noise
    weight_positive = bool(central["weight_phys"] > 0.0)
    sign_check_report = {
        "statement": (
            "no statistically significant sign violation was detected in the "
            "diagnostic correlator"
        ),
        "not_a_proof_of_spectral_positivity": True,
        "why_not_a_proof": (
            "pointwise positivity of G(t) does not imply a positive spectral "
            "measure (counterexample: G(t) = e^{-t} - 0.1 e^{-2t} > 0 for t >= 0 "
            "with a negative coefficient); spectral positivity is conclusion C3, "
            "conditional on the reflection-positive transfer premise P4, which is "
            "not available for the finite-spacing clover engine"
        ),
        "rule": (
            "G(t) > 0 on every statistically significant timeslice (|G| > 2 sigma) "
            "and |G| <= 2 sigma wherever the central value is negative"
        ),
        "timeslices": [
            {
                "t": int(t),
                "value": float(folded_mean[t]),
                "jackknife_error": float(folded_err[t]),
                "significant": bool(significant[t]),
            }
            for t in range(len(folded_mean))
        ],
        "positive_where_significant": positive_where_significant,
        "negative_values_consistent_with_zero_at_2_sigma": negative_only_within_noise,
        "noise_note": (
            "the vector channel enters noise at t >= 6 on this volume (declared in "
            "the hybrid IR bracket note); negative central values there are "
            "statistical fluctuations around zero"
        ),
    }

    # quadrature/window budget: shift the extraction window start by +1
    window_shift_mass = abs(
        _windowed_mass(folded_mean, EFFMASS_WINDOW_START + 1) - central["am_v"]
    )

    unquantified = {
        "status": "not_quantified_at_diagnostic_scale",
        "production_scope": "production backend export budgets (dependency_note)",
    }
    am_v = central["am_v"]

    return {
        "artifact": "oph_diagnostic_ward_projected_spectral_measure_demonstrator",
        "non_promoting_diagnostic": True,
        "role": (
            "executed implementation witness that the constructor pipeline runs "
            "end-to-end on source-only lattice data and emits a machine-readable "
            "object; diagnostic scale, non-promoting, not a production payload, "
            "rejected by the production gate by construction"
        ),
        "ensemble": {
            "npz": npz_path.relative_to(ROOT).as_posix(),
            "ensemble_id": str(d["ensemble_id"]),
            "lattice_shape_TXYZ": [int(v) for v in d["shape"]],
            "beta": float(d["beta"]),
            "kappa": kappa,
            "n_configs": n_cfg,
            "rng_seed": int(d["seed"]),
            "flavors": "quenched_nf0",
            "external_inputs_used": False,
            "input_note": (
                "bare parameters (beta, kappa) only; the extracted level and residue "
                "are whatever the lattice dynamics produce; no measured HVP, R(s), "
                "alpha, or hadron-mass input enters the constructor"
            ),
        },
        "units": "lattice units a = 1; no scale setting at diagnostic scale (declared)",
        "charge_convention": (
            "single flavor, unit charge (R-normalization N_c Q^2 = 3); the U(1)_Q "
            "flavor charge weighting is a declared production-side factor"
        ),
        "finite_volume_levels": [
            {
                "ensemble_id": str(d["ensemble_id"]),
                "channel": "vector_transverse_unit_charge",
                "levels": [
                    {
                        "level_id": "vector_ground_A",
                        "s": am_v * am_v,
                        "energy": am_v,
                        "energy_jackknife_error": errors["am_v"],
                        "weight": central["weight_phys"],
                        "weight_jackknife_error": errors["weight_phys"],
                    }
                ],
            }
        ],
        "ward_projected_residues": [
            {
                "level_id": "vector_ground_A",
                "residue": central["weight_phys"],
                "current_normalization": (
                    "Z_V^eff = C_CL/(2 kappa C_LL) from the exactly conserved "
                    f"point-split current; Z_V = {central['z_v']:.6f} "
                    f"+- {errors['z_v']:.6f} (windowed average, diagnostic)"
                ),
            }
        ],
        "rho_had_or_measure": {
            "representation": "primitive_spectral_measure",
            "support_variable": "s",
            "pushforward_rule": (
                "rho(s) = w delta(s - (a m_V)^2) from the single-state model of the "
                "folded transverse correlator; the production pushforward to "
                "rho_Q(s;P) is the declared Luscher-class finite-volume map "
                "(rho_scattering/ scaffolding, production backend export lane)"
            ),
            "positivity_status": (
                "diagnostic_sign_check_passed_not_a_positivity_certificate"
                if sign_check_passed and weight_positive
                else "DIAGNOSTIC_SIGN_CHECK_FAILED"
            ),
            "sign_check_report": sign_check_report,
        },
        "covariance": {
            "object": (
                "jackknife covariance of the folded transverse correlator (this "
                "block is exactly the correlator covariance and nothing more)"
            ),
            "dimension": int(covariance.shape[0]),
            "matrix": [[float(x) for x in row] for row in covariance],
            "min_eigenvalue": float(eigvals.min()),
            "max_eigenvalue": float(eigvals.max()),
            "positive_semidefinite": cov_psd,
        },
        "joint_extraction_covariance": {
            "object": (
                "joint leave-one-out jackknife covariance of the nonlinear "
                "extracted observables from the one common resampling pipeline"
            ),
            "row_basis": extraction_keys,
            "dimension": len(extraction_keys),
            "matrix": [[float(x) for x in row] for row in joint_cov],
            "min_eigenvalue": float(joint_eigvals.min()),
            "max_eigenvalue": float(joint_eigvals.max()),
            "positive_semidefinite": joint_psd,
        },
        "extraction": {
            "am_vector": am_v,
            "am_vector_jackknife_error": errors["am_v"],
            "z_v": central["z_v"],
            "z_v_jackknife_error": errors["z_v"],
            "amplitude_hop": central["amplitude_hop"],
            "weight_physical_units": central["weight_phys"],
            "weight_jackknife_error": errors["weight_phys"],
            "window_rule": (
                "cosh effective mass over the positive decaying window from t = 3; "
                "Z_V from the windowed average of Z_V^eff(t) over the same range "
                "(not a validated plateau)"
            ),
        },
        "diagnostic_limitations": [
            "no correlated fit and no goodness-of-fit statistic",
            "no excited-state stability test and no model-selection test",
            "no model-truncation uncertainty budget",
            "Z_V is a windowed average, not a statistically validated plateau",
            "single-state model at toy volume; values carry O(50%) jackknife errors",
            (
                "the emitted correlator covariance is the covariance of the folded "
                "correlator only; the joint covariance of the extracted observables "
                "is reported separately in joint_extraction_covariance"
            ),
        ],
        "systematics": {
            "statistical_budget": {
                "method": "leave-one-out jackknife through the full constructor",
                "relative_error_weight": errors["weight_phys"] / central["weight_phys"],
                "relative_error_energy": errors["am_v"] / central["am_v"],
            },
            "quadrature_budget": {
                "method": "extraction-window start shifted by +1",
                "abs_mass_shift": window_shift_mass,
            },
            "continuum_budget": unquantified,
            "finite_volume_budget": unquantified,
            "chiral_budget": unquantified,
            "current_matching_budget": unquantified,
            "endpoint_remainder_budget": unquantified,
        },
        "guards": {
            "execution_class": "real_lattice_diagnostic_toy_scale",
            "promotion_allowed": False,
            "production_schema_conformant": False,
            "production_schema_nonconformance_reason": (
                "quenched branch, no typed provenance manifest, no transport-moment "
                "certificate, unquantified budget rows: the strict production gate "
                "rejects this demonstrator on multiple independent grounds, by "
                "construction"
            ),
            "surrogate_hadron_artifact": False,
            "target_anchored": False,
            "external_inputs_used": False,
        },
        "checks_passed": bool(
            sign_check_passed
            and weight_positive
            and cov_psd
            and joint_psd
            and math.isfinite(am_v)
            and am_v > 0.0
        ),
    }


# ---------------------------------------------------------------------------
# Acceptance gate: strict validator with a full adversarial battery.
# ---------------------------------------------------------------------------


def build_conformant_payload() -> dict[str, Any]:
    """Synthetic production payload conforming to the strict export contract."""
    return {
        "artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
        "format_version": 2,
        "profile_id": "issue_317_acceptance_gate_synthetic_conformant_payload",
        "branch": {"flavors": "2+1", "qed": "off"},
        "projection": {
            "lane": "U(1)_Q",
            "ward_projected": True,
            "external_kinematics": "zero-momentum rest-frame finite-volume levels",
            "zero_momentum_endpoint_compatible": True,
        },
        "backend": {
            "family": "synthetic_gate_check",
            "name": "acceptance_gate",
            "version": "1",
            "git_commit": "0000000000",
            "run_id": "gate",
            "build_id": "gate",
            "machine": "gate",
        },
        "provenance": {
            "source_inputs": [
                {"kind": "source_ensemble", "identifier": "gate_ens"},
                {"kind": "declared_convention", "identifier": "Z_V ledger, Q = T3 + Y"},
            ],
            "measured_hvp_input_present": False,
            "target_calibration_present": False,
        },
        "external_targets_used": [],
        "finite_volume_levels": [
            {
                "ensemble_id": "gate_ens",
                "channel": "U(1)_Q_vector",
                "levels": [
                    {"level_id": "L0", "s": 1.0, "energy": 1.0, "weight": 0.5}
                ],
            }
        ],
        "ward_projected_residues": [
            {"level_id": "L0", "residue": 0.5, "current_normalization": "Z_V ledger"}
        ],
        "rho_had_or_measure": {
            "representation": "rho_had",
            "support_variable": "s",
            "pushforward_rule": "declared Luscher-class finite-volume to continuum map",
            "positivity_status": "certified_positive",
            "positivity_certificate": (
                "reflection-positive transfer certificate reference (premise P4 "
                "content, supplied by the production packet)"
            ),
        },
        "covariance": {
            "row_basis": ["L0_energy", "L0_weight"],
            "dimension": 2,
            "matrix": [[1.0e-4, 0.0], [0.0, 1.0e-4]],
        },
        "transport_moment_certificate": {
            "kernel": "mZ(P)^2/(3*pi*s*(s+mZ(P)^2))",
            "Delta_had_image": {"lo": 4.30, "hi": 4.40},
            "quadrature_error_bound": 1.0e-6,
            "tail_bound": 1.0e-6,
        },
        "systematics": {
            "statistical_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "continuum_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "finite_volume_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "chiral_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "current_matching_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "quadrature_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
            "endpoint_remainder_budget": {"bound_interval": {"lo": 0.0, "hi": 0.01}},
        },
        "guards": {
            "stable_channel_only": False,
            "surrogate_hadron_artifact": False,
            "compare_only_external_endpoint": False,
        },
    }


def build_negative_controls() -> list[dict[str, Any]]:
    """Adversarial payload mutations covering every declared semantic requirement."""
    controls: list[dict[str, Any]] = []

    def add(control_id: str, requirements: list[str], mutate: Callable[[dict[str, Any]], None]) -> None:
        payload = build_conformant_payload()
        mutate(payload)
        controls.append(
            {
                "control_id": control_id,
                "semantic_requirements": requirements,
                "payload": payload,
            }
        )

    def level0(p: dict[str, Any]) -> dict[str, Any]:
        return p["finite_volume_levels"][0]["levels"][0]

    add(
        "measured_hvp_comparison_endpoint",
        ["guard_flags_exactly_false"],
        lambda p: p["guards"].__setitem__("compare_only_external_endpoint", True),
    )
    add(
        "forbidden_target_leak_ee_to_hadrons",
        ["no_measured_hvp_or_target_inputs"],
        lambda p: p.__setitem__("external_targets_used", ["EE_TO_HADRONS"]),
    )
    add(
        "nested_measured_hvp_provenance",
        ["no_measured_hvp_or_target_inputs", "closed_object_boundaries"],
        lambda p: p.__setitem__("measured_hvp_input", {"source": "EE_TO_HADRONS"}),
    )
    add(
        "external_targets_used_supplied_as_string",
        ["no_measured_hvp_or_target_inputs"],
        lambda p: p.__setitem__("external_targets_used", "EE_TO_HADRONS"),
    )
    add(
        "s_not_energy_squared",
        ["s_equals_energy_squared"],
        lambda p: level0(p).__setitem__("s", 999.0),
    )
    add(
        "dangling_residue_reference",
        ["residue_level_reference_integrity"],
        lambda p: p["ward_projected_residues"][0].__setitem__("level_id", "NO_SUCH_LEVEL"),
    )
    add(
        "duplicate_level_id",
        ["unique_level_identifiers"],
        lambda p: p["finite_volume_levels"][0]["levels"].append(
            {"level_id": "L0", "s": 4.0, "energy": 2.0, "weight": 0.25}
        ),
    )
    add(
        "nan_energy",
        ["finite_numeric_values"],
        lambda p: level0(p).__setitem__("energy", float("nan")),
    )
    add(
        "infinite_weight",
        ["finite_numeric_values"],
        lambda p: level0(p).__setitem__("weight", float("inf")),
    )
    add(
        "nan_residue",
        ["finite_numeric_values"],
        lambda p: p["ward_projected_residues"][0].__setitem__("residue", float("nan")),
    )
    add(
        "all_budgets_status_not_quantified",
        ["finite_ordered_bound_intervals"],
        lambda p: p.__setitem__(
            "systematics",
            {name: {"status": "not_quantified"} for name in p["systematics"]},
        ),
    )
    add(
        "missing_chiral_budget",
        ["complete_budget_rows"],
        lambda p: p["systematics"].pop("chiral_budget"),
    )
    add(
        "unordered_bound_interval",
        ["finite_ordered_bound_intervals"],
        lambda p: p["systematics"]["statistical_budget"].__setitem__(
            "bound_interval", {"lo": 0.02, "hi": 0.01}
        ),
    )
    add(
        "nonfinite_bound_interval",
        ["finite_ordered_bound_intervals"],
        lambda p: p["systematics"]["continuum_budget"].__setitem__(
            "bound_interval", {"lo": 0.0, "hi": "inf"}
        ),
    )
    add(
        "missing_transport_moment_certificate",
        ["transport_moment_certificate_complete"],
        lambda p: p.pop("transport_moment_certificate"),
    )
    add(
        "covariance_asymmetric",
        ["covariance_dimension_and_symmetry"],
        lambda p: p["covariance"].__setitem__("matrix", [[1.0e-4, 5.0e-5], [0.0, 1.0e-4]]),
    )
    add(
        "covariance_not_positive_semidefinite",
        ["covariance_positive_semidefinite"],
        lambda p: p["covariance"].__setitem__("matrix", [[1.0e-4, 2.0e-4], [2.0e-4, 1.0e-4]]),
    )
    add(
        "covariance_dimension_mismatch",
        ["covariance_dimension_and_symmetry"],
        lambda p: p["covariance"].__setitem__("dimension", 3),
    )
    add(
        "weight_residue_inconsistent",
        ["weight_residue_consistency"],
        lambda p: p["ward_projected_residues"][0].__setitem__("residue", 0.4),
    )
    add(
        "negative_residue",
        ["nonnegative_weights_and_residues"],
        lambda p: (
            p["ward_projected_residues"][0].__setitem__("residue", -0.5),
            level0(p).__setitem__("weight", -0.5),
        )
        and None,
    )
    add(
        "negative_level_weight",
        ["nonnegative_weights_and_residues"],
        lambda p: level0(p).__setitem__("weight", -0.1),
    )
    add(
        "nonpositive_energy",
        ["positive_energies"],
        lambda p: (level0(p).__setitem__("energy", 0.0), level0(p).__setitem__("s", 0.0)) and None,
    )
    add(
        "empty_level_support",
        ["nonempty_level_support"],
        lambda p: p.__setitem__("finite_volume_levels", []),
    )
    add(
        "quenched_branch",
        ["strict_branch_typing"],
        lambda p: p["branch"].__setitem__("flavors", "quenched"),
    )
    add(
        "surrogate_hadron_artifact",
        ["guard_flags_exactly_false"],
        lambda p: p["guards"].__setitem__("surrogate_hadron_artifact", True),
    )
    add(
        "stable_channel_only_export",
        ["guard_flags_exactly_false"],
        lambda p: p["guards"].__setitem__("stable_channel_only", True),
    )
    add(
        "guard_stringly_false",
        ["guard_flags_exactly_false"],
        lambda p: p["guards"].__setitem__("compare_only_external_endpoint", "false"),
    )
    add(
        "ward_projection_dropped",
        ["ward_projection_required"],
        lambda p: p["projection"].__setitem__("ward_projected", False),
    )
    add(
        "provenance_manifest_missing",
        ["complete_typed_provenance"],
        lambda p: p.pop("provenance"),
    )
    add(
        "provenance_unknown_input_kind",
        ["complete_typed_provenance"],
        lambda p: p["provenance"]["source_inputs"].append(
            {"kind": "external_dataset", "identifier": "unknown"}
        ),
    )
    add(
        "provenance_hvp_flag_true",
        ["complete_typed_provenance"],
        lambda p: p["provenance"].__setitem__("measured_hvp_input_present", True),
    )
    add(
        "provenance_empty_source_inputs",
        ["complete_typed_provenance"],
        lambda p: p["provenance"].__setitem__("source_inputs", []),
    )
    add(
        "unknown_root_key",
        ["closed_object_boundaries"],
        lambda p: p.__setitem__("undeclared_extra_block", {"note": "unknown"}),
    )
    add(
        "unknown_guard_key",
        ["closed_object_boundaries"],
        lambda p: p["guards"].__setitem__("undeclared_guard", False),
    )
    add(
        "positivity_status_unknown_string",
        ["positivity_status_allowlisted"],
        lambda p: p["rho_had_or_measure"].__setitem__("positivity_status", "checked_numerically"),
    )
    add(
        "positivity_status_sign_check_substitution",
        ["positivity_status_allowlisted"],
        lambda p: p["rho_had_or_measure"].__setitem__(
            "positivity_status", "diagnostic_sign_check_passed"
        ),
    )
    return controls


def run_acceptance_gate() -> dict[str, Any]:
    """Strict-gate witness: conformant accepted, full adversarial battery rejected."""
    schema = strict_validator.load_schema()
    conformant = strict_validator.validate_production_payload(build_conformant_payload(), schema)
    control_rows = []
    all_rejected = True
    coverage: dict[str, list[str]] = {name: [] for name in strict_validator.SEMANTIC_REQUIREMENTS}
    for control in build_negative_controls():
        result = strict_validator.validate_production_payload(control["payload"], schema)
        rejected = not result.accepted
        control_rows.append(
            {
                "control_id": control["control_id"],
                "semantic_requirements": control["semantic_requirements"],
                "rejected": rejected,
                "rejection_reasons": list(result.reasons[:4]),
            }
        )
        all_rejected = all_rejected and rejected
        if rejected:
            for requirement in control["semantic_requirements"]:
                coverage.setdefault(requirement, []).append(control["control_id"])
    uncovered = sorted(
        name for name in strict_validator.SEMANTIC_REQUIREMENTS if not coverage.get(name)
    )
    return {
        "schema": SCHEMA_PATH.relative_to(ROOT).as_posix(),
        "validator": "code/particles/hadron/ward_projected_spectral_measure_validator.py",
        "semantic_requirements": list(strict_validator.SEMANTIC_REQUIREMENTS),
        "declared_tolerances": {
            "s_equals_energy_squared_rel_tol": strict_validator.S_EQUALS_ENERGY_SQUARED_REL_TOL,
            "weight_residue_rel_tol": strict_validator.WEIGHT_RESIDUE_REL_TOL,
            "covariance_symmetry_rel_tol": strict_validator.COVARIANCE_SYMMETRY_REL_TOL,
            "covariance_psd_rel_tol": strict_validator.COVARIANCE_PSD_REL_TOL,
        },
        "conformant_payload_accepted": bool(conformant.accepted),
        "conformant_payload_reasons": list(conformant.reasons),
        "negative_controls": control_rows,
        "all_negative_controls_rejected": bool(all_rejected),
        "adversarial_coverage": {
            requirement: controls for requirement, controls in coverage.items()
        },
        "uncovered_semantic_requirements": uncovered,
        "adversarial_coverage_complete": not uncovered,
        "passed": bool(conformant.accepted and all_rejected and not uncovered),
    }


def build_gate_approved_variants() -> list[tuple[str, dict[str, Any]]]:
    """Gate-approved payload variants for the downstream implication."""
    base = build_conformant_payload()

    multi = copy.deepcopy(base)
    multi["profile_id"] = "issue_317_gate_multi_level_variant"
    multi["finite_volume_levels"][0]["levels"].append(
        {
            "level_id": "L1",
            "s": 2.25,
            "energy": 1.5,
            "weight": 0.25,
            "energy_jackknife_error": 0.01,
            "weight_jackknife_error": 0.005,
        }
    )
    multi["ward_projected_residues"].append(
        {
            "level_id": "L1",
            "residue": 0.25,
            "current_normalization": "Z_V ledger",
            "residue_jackknife_error": 0.005,
        }
    )
    multi["covariance"] = {
        "row_basis": ["L0_energy", "L0_weight", "L1_energy", "L1_weight"],
        "dimension": 4,
        "matrix": [
            [1.0e-4, 0.0, 0.0, 0.0],
            [0.0, 1.0e-4, 0.0, 0.0],
            [0.0, 0.0, 1.0e-4, 0.0],
            [0.0, 0.0, 0.0, 1.0e-4],
        ],
    }

    two_ensembles = copy.deepcopy(base)
    two_ensembles["profile_id"] = "issue_317_gate_two_ensemble_variant"
    two_ensembles["finite_volume_levels"].append(
        {
            "ensemble_id": "gate_ens_B",
            "channel": "U(1)_Q_vector",
            "levels": [{"level_id": "B0", "s": 1.21, "energy": 1.1, "weight": 0.3}],
        }
    )
    two_ensembles["ward_projected_residues"].append(
        {"level_id": "B0", "residue": 0.3, "current_normalization": "Z_V ledger"}
    )

    string_bounds = copy.deepcopy(base)
    string_bounds["profile_id"] = "issue_317_gate_string_decimal_bounds_variant"
    string_bounds["transport_moment_certificate"] = {
        "kernel": "mZ(P)^2/(3*pi*s*(s+mZ(P)^2))",
        "Delta_had_image": {"lo": "4.30", "hi": "4.40"},
        "quadrature_error_bound": "1e-40",
        "tail_bound": "1e-40",
    }
    for budget in string_bounds["systematics"].values():
        budget["bound_interval"] = {"lo": "0", "hi": "1e-30"}

    return [
        ("conformant_baseline", base),
        ("multi_level", multi),
        ("two_ensembles", two_ensembles),
        ("string_decimal_bounds", string_bounds),
    ]


def downstream_compatibility_witness() -> dict[str, Any]:
    """Executed implication: gate accepts payload => downstream validator accepts source_measure.

    Every gate-approved payload variant is embedded as ``source_measure``
    into the downstream source-transport contract
    (code/P_derivation/thomson_spectral_transport.py). The witness requires
    zero downstream source-measure reasons and zero downstream forbidden
    keys for every variant, so the contract certified here is one the
    stated consumer actually accepts.
    """
    schema = strict_validator.load_schema()
    rows = []
    all_ok = True
    for variant_id, payload in build_gate_approved_variants():
        gate_result = strict_validator.validate_production_payload(payload, schema)
        downstream_reasons: list[str] = []
        _validate_source_measure({"source_measure": payload}, downstream_reasons)
        forbidden = sorted(source_payload_forbidden_keys({"source_measure": payload}))
        ok = bool(gate_result.accepted and not downstream_reasons and not forbidden)
        rows.append(
            {
                "variant_id": variant_id,
                "gate_accepted": bool(gate_result.accepted),
                "gate_reasons": list(gate_result.reasons[:4]),
                "downstream_source_measure_reasons": downstream_reasons,
                "downstream_forbidden_keys": forbidden,
                "implication_holds": ok,
            }
        )
        all_ok = all_ok and ok
    return {
        "implication": (
            "packet gate accepts payload => downstream source-transport validator "
            "accepts the resulting source_measure"
        ),
        "downstream_validator": (
            "code/P_derivation/thomson_spectral_transport.py "
            "(_validate_source_measure and source_payload_forbidden_keys; the same "
            "module whose validate_source_transport_payload gates the Thomson "
            "endpoint promotion)"
        ),
        "variants": rows,
        "passed": bool(all_ok),
    }


# ---------------------------------------------------------------------------
# Derivation chain.
# ---------------------------------------------------------------------------


def build_derivation_chain() -> list[dict[str, Any]]:
    return [
        {
            "step": 1,
            "premise": "P1: declared gauge-invariant quotient ensemble and action",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "receipts": "particles/runs/qcd/hadron_source_backend/qcd_ensemble/",
            "conclusion": (
                "the source law mu_r^QCD is a declared positive Gibbs measure on the "
                "gauge quotient at zero density and theta = 0, with refinement maps and "
                "no-target-leak controls; a normal form is not a probability law "
                "(source-law non-substitution theorem); its physical population is "
                "production-lane content (REQUIRED_NOT_POPULATED)"
            ),
        },
        {
            "step": 2,
            "premise": "P4: reflection positivity and transfer operator (declared, not supplied here)",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "receipts": "particles/runs/qcd/hadron_source_backend/vacuum/",
            "conclusion": (
                "given a reflection-positive Euclidean slab, the transfer operator's "
                "spectral decomposition supplies the Hilbert-space representation used "
                "for positivity; reflection positivity is a theorem for the unimproved "
                "Wilson action (Osterwalder-Seiler) and is not available at finite "
                "spacing for the clover-improved engine, so conclusion C3 stays "
                "conditional on this premise and no finite-sample sign check "
                "substitutes for it"
            ),
        },
        {
            "step": 3,
            "premise": "P2, P3: conserved current and Ward certificate",
            "source_artifact": "code/particles/hadron/lattice_backend/conserved_vector.py",
            "conclusion": (
                "the point-split conserved U(1)_Q current V_mu(x) = kappa [psibar(x+mu) "
                "(1+gamma_mu) U_mu(x)^dag psi(x) - psibar(x)(1-gamma_mu) U_mu(x) psi(x+mu)] "
                "is gauge invariant and exactly conserved (C1, C2; executed witnesses "
                "on the actual conserved-current observable); the local current carries "
                "Z_V^eff = C_CL/(2 kappa C_LL) with the free-field anchor exactly 1 (C5)"
            ),
        },
        {
            "step": 4,
            "premise": "Ward projection",
            "source_artifact": "code/P_derivation/SOURCE_SPECTRAL_THEOREM.md",
            "conclusion": (
                "current conservation gives q_mu Pi_Q^{mu nu} = 0, so the two-current "
                "correlator has the transverse scalar form "
                "Pi_Q^{mu nu} = (q^mu q^nu - q^2 eta^{mu nu}) Pi_Q(q^2) (C2)"
            ),
        },
        {
            "step": 5,
            "premise": "P4, P5, P7: positive spectral measure with thresholds (conditional)",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "conclusion": (
                "under the reflection-positivity/transfer premise the smeared-current "
                "quadratic spectral measure is positive and locally finite; in finite "
                "volume it is the atomic measure rho_Q(s) = sum_n Z_n delta(s - E_n^2), "
                "Z_n >= 0, supported above the vacuum from the lowest vector-channel "
                "level (C3, C4)"
            ),
        },
        {
            "step": 6,
            "premise": "P8: resonances and continuum typing",
            "source_artifact": "code/particles/hadron/rho_scattering/",
            "conclusion": (
                "two-body channels map to infinite volume through the declared "
                "Luscher-class pushforward; unstable resonances and thresholds are "
                "treated explicitly at that stage (production backend export lane)"
            ),
        },
        {
            "step": 7,
            "premise": "P9: scheme lock and transport kernel",
            "source_artifact": "code/particles/hadron/lattice_backend/vector_correlator.py",
            "conclusion": (
                "the once-subtracted spacelike kernel identity "
                "Delta_had = 4 pi Pihat(mZ^2) is exact; the measure itself is "
                "scheme-independent and the same-subtraction lock ties the transport "
                "to a0(P) in the declared convention (C6)"
            ),
        },
        {
            "step": 8,
            "premise": "P5: covariance and uncertainty propagation",
            "source_artifact": "code/particles/hadron/verify_issue_317_spectral_measure_packet.py",
            "conclusion": (
                "one common leave-one-out jackknife pipeline propagates the correlator "
                "covariance and the joint covariance of the extracted observables "
                "(executed witness; both blocks emitted and labeled separately); "
                "production payloads must carry the covariance block and the seven-row "
                "budget ledger with finite ordered bound intervals (C7)"
            ),
        },
        {
            "step": 9,
            "premise": "machine-readable emission with a strict fail-closed gate",
            "source_artifact": "code/particles/hadron/ward_projected_spectral_measure_validator.py",
            "conclusion": (
                "the strict typed validator positively enforces finiteness, identifier "
                "integrity, s = E^2, residue/weight consistency, typed provenance, "
                "covariance integrity, ordered budget intervals, the transport-moment "
                "certificate, branch typing, and literal guard flags; unknown, missing, "
                "malformed, contradictory, or non-finite values fail closed (executed "
                "adversarial battery with complete requirement coverage)"
            ),
        },
        {
            "step": 10,
            "premise": "downstream consumer compatibility",
            "source_artifact": "code/P_derivation/thomson_spectral_transport.py",
            "conclusion": (
                "every gate-approved payload embeds as source_measure into the "
                "downstream source-transport validator with zero reasons (executed "
                "implication witness), so the certified contract is the one the "
                "stated consumer accepts"
            ),
        },
        {
            "step": 11,
            "premise": "P10 and typing: blinded comparison, higher-point, radiative",
            "source_artifact": "code/particles/runs/qcd/hadron_source_backend/",
            "conclusion": (
                "empirical compilations are typed oph_plus_empirical_hadron_closure with "
                "promotable_as_oph_source_theorem = false and enter only through the "
                "comparison manifest; the two-current measure is a marginal of the "
                "hadronic precision functor, with four-current and transition objects "
                "separately typed (Q4_HLbL_receipt: TWO_POINT_MEASURE_INSUFFICIENT) and "
                "QED/EW radiative corrections confined to the Xi_Q remainder ledger (C8)"
            ),
        },
    ]


# ---------------------------------------------------------------------------
# Packet assembly.
# ---------------------------------------------------------------------------


def _witness_labels_non_promoting(
    witnesses: dict[str, Any], demonstrator: dict[str, Any]
) -> dict[str, Any]:
    checks = {
        "gauge_and_ward_labeled": witnesses.get("non_promoting_diagnostic") is True,
        "zv_anchor_labeled": True,  # filled by caller with the zv witness
        "demonstrator_labeled": (
            demonstrator.get("non_promoting_diagnostic") is True
            and demonstrator["guards"]["promotion_allowed"] is False
        ),
        "sign_check_labeled_not_a_proof": (
            demonstrator["rho_had_or_measure"]["sign_check_report"][
                "not_a_proof_of_spectral_positivity"
            ]
            is True
        ),
    }
    return checks


def build_packet() -> dict[str, Any]:
    theorem = build_theorem()
    theorem_check = theorem_structure_check(theorem)
    witnesses = gauge_and_ward_witness()
    zv_witness = zv_free_field_witness()
    demonstrator = demonstrator_measure()
    gate = run_acceptance_gate()
    downstream = downstream_compatibility_witness()
    empirical_typing = empirical_typing_check(gate)
    higher_point_typing = higher_point_typing_check()
    probes = fail_closed_probe_battery()
    payload_verdict = physical_source_payload_verdict()
    promotion_verdict = physical_promotion_verdict(payload_verdict)

    non_promoting_labels = _witness_labels_non_promoting(witnesses, demonstrator)
    non_promoting_labels["zv_anchor_labeled"] = zv_witness.get("non_promoting_diagnostic") is True

    # Deterministic re-emission agreement on the cheap acceptance-relevant
    # blocks (theorem, gate, downstream implication); the full stored-vs-fresh
    # agreement including the seeded numeric witnesses is enforced by
    # test_issue_317_spectral_measure_packet.py.
    theorem_again = build_theorem()
    gate_again = run_acceptance_gate()
    downstream_again = downstream_compatibility_witness()
    reemission_agree = bool(
        theorem_again == theorem
        and gate_again == gate
        and downstream_again == downstream
    )

    contract_checks = {
        "theorem_typed_conditional": theorem_check["passed"],
        "conserved_current_gauge_witness": witnesses["conserved_current_gauge_invariance"]["passed"],
        "local_current_gauge_witness": witnesses["local_current_gauge_invariance"]["passed"],
        "ward_identity_witness": witnesses["ward_identity"]["passed"],
        "zv_free_field_anchor": zv_witness["passed"],
        "demonstrator_pipeline": demonstrator["checks_passed"],
        "strict_gate_with_full_adversarial_coverage": gate["passed"],
        "downstream_consumer_compatibility": downstream["passed"],
        "empirical_typing": empirical_typing["passed"],
        "higher_point_typing": higher_point_typing["passed"],
        "fail_closed_status_probes": probes["passed"],
        "diagnostic_witnesses_labeled_non_promoting": all(non_promoting_labels.values()),
        "deterministic_reemission_agreement": reemission_agree,
    }
    contract_certified = bool(all(contract_checks.values()))

    verdicts = {
        "contract_certified": {
            "value": contract_certified,
            "meaning": (
                "the conditional construction theorem and the machine-readable "
                "contract are stated, machine-checked, fail-closed, and "
                "downstream-compatible; no physical claim is included"
            ),
            "computed_from": contract_checks,
        },
        "physical_source_payload_available": {
            "value": bool(payload_verdict["available"]),
            "meaning": (
                "a gate-approved physical production payload exists and the "
                "production-lane artifacts are in their explicit success states "
                "(fail-closed allowlists); supplied by the production backend lane"
            ),
            "computed_from": payload_verdict,
        },
        "physical_promotion_allowed": {
            "value": bool(promotion_verdict["allowed"]),
            "meaning": (
                "the physical payload is available AND the downstream "
                "source-transport certificate accepts the production transport "
                "payload; only then may any hadronic number be promoted"
            ),
            "computed_from": promotion_verdict,
        },
    }

    criteria = {
        "conditional_theorem_stated_with_typed_premises_and_conclusions": {
            "packet_level_passed": bool(theorem_check["passed"]),
            "machine_checks": "theorem_structure_check (theorem.premises P1-P10, theorem.conclusions C1-C8)",
        },
        "gauge_invariance_and_ward_identity_proved_for_point_split_current": {
            "packet_level_passed": bool(
                witnesses["conserved_current_gauge_invariance"]["passed"]
                and witnesses["local_current_gauge_invariance"]["passed"]
                and witnesses["ward_identity"]["passed"]
            ),
            "machine_checks": (
                "machine_witnesses.gauge_and_ward: the gauge witness now includes the "
                "conserved(sink)-local(source) correlator containing the actual "
                "link-inserted point-split conserved current, plus the exact "
                "backward-difference Ward identity"
            ),
        },
        "positivity_and_threshold_conditional_without_sign_check_substitution": {
            "packet_level_passed": bool(
                theorem_check["checks"][
                    "positivity_conclusion_requires_reflection_positivity_premise"
                ]
                and theorem_check["checks"][
                    "threshold_conclusion_requires_reflection_positivity_premise"
                ]
                and demonstrator["rho_had_or_measure"]["sign_check_report"][
                    "not_a_proof_of_spectral_positivity"
                ]
                is True
                and gate["adversarial_coverage"].get("positivity_status_allowlisted")
            ),
            "machine_checks": (
                "conclusions C3/C4 are typed conditional on premise P4 (reflection "
                "positivity + transfer); the demonstrator's finite-sample check is "
                "labeled a diagnostic sign check that is not a positivity proof; the "
                "production gate rejects sign-check positivity statuses "
                "(positivity_status_sign_check_substitution control)"
            ),
        },
        "normalization_and_scheme_dependence_stated_conditionally": {
            "packet_level_passed": bool(
                zv_witness["passed"] and theorem_check["passed"]
            ),
            "machine_checks": (
                "conclusions C5 (Z_V = 1 conserved / declared matching certificate) and "
                "C6 (scheme confinement to the declared kernel) typed against premises "
                "P6/P9; free-field anchor executed as a convention witness, explicitly "
                "not an interacting plateau determination"
            ),
        },
        "validator_fails_closed_on_all_invalid_payload_classes": {
            "packet_level_passed": bool(
                gate["all_negative_controls_rejected"] and gate["conformant_payload_accepted"]
            ),
            "machine_checks": (
                "strict validator battery: malformed, non-finite, inconsistent, "
                "target-leaking (incl. nested provenance and stringly-typed target "
                "lists), quenched, surrogate, and unquantified payloads all rejected "
                "with typed reasons"
            ),
        },
        "gate_validates_consistency_relations": {
            "packet_level_passed": bool(
                all(
                    gate["adversarial_coverage"].get(name)
                    for name in (
                        "s_equals_energy_squared",
                        "unique_level_identifiers",
                        "residue_level_reference_integrity",
                        "weight_residue_consistency",
                        "covariance_dimension_and_symmetry",
                        "covariance_positive_semidefinite",
                        "finite_ordered_bound_intervals",
                    )
                )
                and gate["conformant_payload_accepted"]
            ),
            "machine_checks": (
                "executed controls for s = E^2, identifier integrity, residue/weight "
                "consistency, covariance integrity, and complete uncertainty intervals"
            ),
        },
        "gate_approved_payloads_accepted_by_downstream_validator": {
            "packet_level_passed": bool(downstream["passed"]),
            "machine_checks": (
                "machine_witnesses.downstream_compatibility: every gate-approved "
                "variant embeds as source_measure with zero downstream reasons"
            ),
        },
        "empirical_data_confined_to_nonsource_nonpromoting_row_class": {
            "packet_level_passed": bool(empirical_typing["passed"]),
            "machine_checks": "machine_witnesses.empirical_typing (live artifact reads)",
        },
        "higher_point_and_radiative_corrections_separately_typed": {
            "packet_level_passed": bool(higher_point_typing["passed"]),
            "machine_checks": "machine_witnesses.higher_point_typing (live artifact reads)",
        },
        "diagnostic_witnesses_explicitly_labeled_non_promoting": {
            "packet_level_passed": bool(all(non_promoting_labels.values())),
            "machine_checks": str(non_promoting_labels),
        },
        "certificate_distinguishes_three_verdicts": {
            "packet_level_passed": bool(
                set(verdicts)
                == {
                    "contract_certified",
                    "physical_source_payload_available",
                    "physical_promotion_allowed",
                }
                and all(isinstance(v["value"], bool) for v in verdicts.values())
                and verdicts["contract_certified"]["value"] is contract_certified
                and (
                    verdicts["physical_source_payload_available"]["value"]
                    or bool(
                        verdicts["physical_source_payload_available"]["computed_from"]["reasons"]
                    )
                )
                and (
                    verdicts["physical_promotion_allowed"]["value"]
                    or bool(verdicts["physical_promotion_allowed"]["computed_from"]["reasons"])
                )
            ),
            "machine_checks": (
                "verdicts.contract_certified / physical_source_payload_available / "
                "physical_promotion_allowed are computed independently; a false "
                "physical verdict must carry typed reasons"
            ),
        },
        "unknown_or_missing_status_values_fail_closed": {
            "packet_level_passed": bool(probes["passed"]),
            "machine_checks": (
                "machine_witnesses.fail_closed_probes: unknown, renamed, missing, "
                "None, blocked, contradictory, and payload-absent states all report "
                "unavailable"
            ),
        },
        "adversarial_tests_cover_all_semantic_requirements": {
            "packet_level_passed": bool(gate["adversarial_coverage_complete"]),
            "machine_checks": (
                "acceptance_gate.adversarial_coverage maps every declared semantic "
                "requirement to at least one rejected control; "
                "uncovered_semantic_requirements is empty"
            ),
        },
        "stored_and_fresh_certificates_agree": {
            "packet_level_passed": bool(reemission_agree),
            "machine_checks": (
                "deterministic re-emission agreement executed at build time on the "
                "theorem, gate, and downstream blocks; full stored-vs-fresh agreement "
                "including the seeded numeric witnesses is enforced by "
                "test_issue_317_spectral_measure_packet.py::"
                "test_stored_packet_matches_fresh_emission"
            ),
        },
    }
    accepted = bool(
        contract_certified and all(row["packet_level_passed"] for row in criteria.values())
    )

    return {
        "issue": 317,
        "certificate_id": "issue-317-hadronic-spectral-measure-proof-packet-v2",
        "artifact": "oph_ward_projected_spectral_measure_proof_packet",
        "status": (
            "conditional_theorem_and_contract_certified_physical_source_payload_unavailable"
        ),
        "accepted": accepted,
        "acceptance_scope": (
            "'accepted' certifies the re-scoped theorem-and-contract layer only: the "
            "conditional construction theorem with typed premises and conclusions, "
            "the strict fail-closed production contract, downstream-consumer "
            "compatibility, typing confinement, and the fail-closed separation of "
            "the physical verdicts. It asserts no physical source payload, no "
            "production positivity, and no promotable hadronic number; those are "
            "the physical_source_payload_available and physical_promotion_allowed "
            "verdicts, both currently false with typed reasons."
        ),
        "theorem": theorem,
        "theorem_structure_check": theorem_check,
        "derivation_chain": build_derivation_chain(),
        "machine_witnesses": {
            "gauge_and_ward": witnesses,
            "zv_free_field_anchor": zv_witness,
            "demonstrator_measure": demonstrator,
            "acceptance_gate": gate,
            "downstream_compatibility": downstream,
            "empirical_typing": empirical_typing,
            "higher_point_typing": higher_point_typing,
            "fail_closed_probes": probes,
        },
        "verdicts": verdicts,
        "acceptance_criteria_status": criteria,
        "issue_closure_condition": {
            "re_scoped_closure_rule": (
                "per the maintainer's re-scoping on the PR review (2026-07-22), the "
                "issue closes when the theorem-and-contract acceptance list passes; "
                "closure of the contract issue must not be interpreted as completion "
                "of the production backend lane or as production of a physical "
                "hadronic spectral measure"
            ),
            "theorem_and_contract_acceptance_list_passes": accepted,
            "original_scope_production_condition": {
                "author_condition": (
                    "backend export bundle, unquenching execution, runtime receipt, "
                    "and uncertainty ledger present"
                ),
                "maps_to": "verdicts.physical_source_payload_available (fail-closed allowlists)",
                "met_locally": bool(payload_verdict["available"]),
            },
            "closing_keyword_policy": (
                "the PR references this issue with a non-closing reference until the "
                "issue body itself adopts the re-scoped acceptance list"
            ),
        },
        "claim_boundary": {
            "closed_here": (
                "the conditional construction theorem SpectralMeasure_Q with ten "
                "typed premises and eight typed conclusions, the strict fail-closed "
                "machine-readable production contract (schema plus semantic "
                "validator) with complete adversarial coverage, the executed "
                "implication that gate-approved payloads are accepted by the "
                "downstream source-transport validator, executed implementation "
                "witnesses on the actual conserved point-split current, and the "
                "typing confinement of empirical, higher-point, and radiative objects"
            ),
            "not_closed_here": [
                "a populated physical 2+1-flavor QCD ensemble or any production payload "
                "(production backend execution, dependency_note)",
                "production spectral positivity: conclusion C3 remains conditional on "
                "the reflection-positive transfer certificate (premise P4)",
                "production-scale resonance extraction through the Luscher-class map",
                "continuum-limit and infinite-volume numerical extrapolation",
                "the source QCD parameter map P -> (g_3, quark masses, scheme)",
                "a source-derived hadronic numerical value or physical promotion of "
                "any hadronic number into the fine-structure endpoint",
            ],
            "scope": (
                "This packet is the conditional proof and contract layer of the "
                "hadronic spectral measure: it fixes what the production export must "
                "satisfy, proves the structural claims conditional on the typed "
                "premises, and emits no physical hadronic prediction. The diagnostic "
                "demonstrator is non-promoting and is rejected by the production "
                "gate it demonstrates."
            ),
        },
        "dependency_note": {
            "contract_design_edge": {
                "direction": "#317 contract specification -> #425 implementation",
                "meaning": (
                    "this packet fixes the theorem and the strict export contract "
                    "that the production backend must satisfy; it consumes no "
                    "production result"
                ),
            },
            "execution_dependency_edge": {
                "direction": (
                    "#425 physical payload -> physical_source_payload_available -> "
                    "physical_promotion_allowed"
                ),
                "meaning": (
                    "the production backend's physical payload is required before "
                    "the physical availability verdict can pass, and the downstream "
                    "source-transport certificate is additionally required before "
                    "any promotion; the production lane is therefore upstream of "
                    "issue completion in the execution sense, while this packet is "
                    "upstream in the specification sense"
                ),
            },
            "gauge_carrier_selection": (
                "the physical gauge/matter carrier selection (#590) gates the "
                "production lane's physical promotion, not this packet"
            ),
        },
        "dependency_artifacts": {
            "export_schema": "code/particles/hadron/ward_projected_spectral_measure.schema.json",
            "strict_validator": "code/particles/hadron/ward_projected_spectral_measure_validator.py",
            "constructive_contract": "code/particles/runs/hadron/ward_projected_spectral_measure_contract.json",
            "nonidentifiability_obstruction": "code/particles/hadron/derive_ward_projected_spectral_measure_obstruction.py",
            "endpoint_reduction_theorem": "code/P_derivation/source_spectral_theorem.py (WardProjectedHadronicSpectralEmission_Q)",
            "downstream_source_transport_validator": "code/P_derivation/thomson_spectral_transport.py",
            "backend_spec_note": "physics-problems/hadronic_precision_endpoint.md",
            "receipt_bundle": "code/particles/runs/qcd/hadron_source_backend/ (claim: SOURCE_PROTOTYPE_NOT_PROMOTED)",
            "lattice_engine": "code/particles/hadron/lattice_backend/",
            "demonstrator_ensemble": "code/particles/runs/hadron/hybrid_ir_ensembleA_2026-07-16.npz",
        },
        "consumer_artifacts": {
            "production_backend_export": (
                "production backend export bundle conforming to the strict schema "
                "and validator certified here (dependency_note); expected at "
                + PRODUCTION_PAYLOAD_PATH.relative_to(ROOT.parent).as_posix()
            ),
            "thomson_endpoint_lane": (
                "code/P_derivation/source_spectral_theorem.py consumes a conforming "
                "payload for the source-only Thomson endpoint interval certificate"
            ),
        },
        "forbidden_inputs": FORBIDDEN_TARGETS,
        "verifier_command": (
            "python3 code/particles/hadron/verify_issue_317_spectral_measure_packet.py "
            "--check --output "
            "code/particles/runs/hadron/ward_projected_spectral_measure_proof_packet.json"
        ),
        "volatile": {
            "generated_utc": _now_utc(),
            "note": "timestamp only; excluded from acceptance logic",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify the issue #317 conditional spectral-measure theorem and contract packet."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--check", action="store_true", help="exit nonzero unless the packet is accepted")
    args = parser.parse_args()

    packet = build_packet()
    text = json.dumps(packet, indent=1) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")

    witnesses = packet["machine_witnesses"]
    gw = witnesses["gauge_and_ward"]
    print(
        "gauge invariance defects: local-local "
        f"{gw['local_current_gauge_invariance']['correlator_relative_defect']:.2e}, "
        f"conserved-local {gw['conserved_current_gauge_invariance']['correlator_relative_defect']:.2e}"
    )
    print(f"ward identity relative defect: {gw['ward_identity']['relative_defect']:.2e}")
    print(f"Z_V free-field max deviation: {witnesses['zv_free_field_anchor']['max_abs_deviation_from_one']:.2e}")
    demo = witnesses["demonstrator_measure"]["extraction"]
    print(
        f"demonstrator (diagnostic, non-promoting): a*m_V = {demo['am_vector']:.4f} "
        f"+- {demo['am_vector_jackknife_error']:.4f}, Z_V = {demo['z_v']:.4f}"
    )
    gate = witnesses["acceptance_gate"]
    print(
        f"strict gate: conformant accepted = {gate['conformant_payload_accepted']}, "
        f"controls rejected = {gate['all_negative_controls_rejected']} "
        f"({len(gate['negative_controls'])} controls), "
        f"coverage complete = {gate['adversarial_coverage_complete']}"
    )
    print(f"downstream compatibility: {witnesses['downstream_compatibility']['passed']}")
    print(f"fail-closed probes: {witnesses['fail_closed_probes']['passed']}")
    verdicts = packet["verdicts"]
    print(f"contract_certified: {verdicts['contract_certified']['value']}")
    print(
        "physical_source_payload_available: "
        f"{verdicts['physical_source_payload_available']['value']}"
    )
    print(f"physical_promotion_allowed: {verdicts['physical_promotion_allowed']['value']}")
    print(f"packet accepted (theorem-and-contract layer): {packet['accepted']}")
    print(f"wrote {out_path}")
    if args.check and not packet["accepted"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
