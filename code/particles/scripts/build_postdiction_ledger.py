#!/usr/bin/env python3
"""Aggregate the certified postdiction rows into one ledger.

The ledger is a deterministic aggregator.  Numeric values and measured
references are read live from their parent artifacts.  Structural rows are
derived from validated structured parents, with any direct algebraic
corollary identified as such.  A missing or inconsistent parent is a hard
failure, not a silently absent row.

Section one records the forced-structure layer: machine-checked finite
theorems and executable certificates that precede or constrain numeric lanes,
including the icosahedral gauge packet and generic observer-law boundaries.
Lean-backed rows record module paths and exact declaration names; executable
rows record their structured artifacts.  The builder rejects a missing
receipt and records the declared hypothesis boundaries of each owning paper.

The numeric sections carry the per-lane claim discipline of their parents:
interval rows report containment of the compare-only witness, conditional
rows carry their declared premises, chart coordinates keep their
NOT_EVALUABLE physical-comparison status, and the quark absolute-mass row
is an obstruction theorem rather than a number.

Run:
    python3 code/particles/scripts/build_postdiction_ledger.py
writes code/particles/runs/status/postdiction_ledger.json and
docs/POSTDICTION_LEDGER.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parent
PARTICLES = SCRIPTS.parent
CODE = PARTICLES.parent
REPO = CODE.parent
RUNS = PARTICLES / "runs"
RUNTIME = CODE / "P_derivation" / "runtime"
LEAN_SCREEN = REPO / "Lean" / "Screen"

PARENTS = {
    "mass_surface": RUNS / "status" / "source_only_mass_prediction_surface.json",
    "conditional_ew": RUNS / "calibration" / "conditional_ew_predictions_current.json",
    "endpoint": RUNTIME / "empirical_thomson_endpoint_current.json",
    "anchor_bridge": RUNTIME / "anchor_scheme_bridge_current.json",
    "kappa_rectangle": RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json",
    "kappa_coherent": RUNS / "leptons" / "charged_kappa_interval_coherent_closure.json",
    "koide_balance": RUNS / "leptons" / "koide_balance_comparison.json",
    "clebsch_lane": RUNS / "flavor" / "down_type_register_clebsch_lane.json",
    "clebsch_selection": RUNS / "flavor" / "clebsch_register_pairing_selection.json",
    "fiber_obstruction": RUNS / "flavor" / "quark_spread_fiber_structure_transport_obstruction.json",
    "matter_receipt": CODE / "a5_closure" / "receipts" / "super_tannakian_matter_reference.receipt.json",
    "matter_menu": CODE / "a5_closure" / "manifests" / "matter_menu_spectral_ledger_reference.json",
    "port_current": CODE / "a5_closure" / "receipts" / "port_current_inner_reference.receipt.json",
    "axis_center_descent": CODE / "a5_closure" / "receipts" / "axis_center_descent_reference.receipt.json",
    "carrier_modes": RUNS / "status" / "carrier_mode_acceptance.json",
    "quantum_carrier_status": RUNS / "status" / "quantum_carrier_status.json",
    "alpha_hvp_verdict": PARTICLES / "alpha_hvp_audit" / "outputs" / "alpha_hvp_class_verdict.json",
    "hadron_payload": RUNS / "hadron" / "empirical_ee_hadronic_spectral_measure.json",
    "solver_standby": RUNS / "qcd" / "hadron_source_backend" / "qcd_ensemble" / "solver_on_standby.json",
    "carrier_class_dispersion": CODE / "a5_fingerprint" / "runtime"
    / "carrier_class_dispersion_receipt.json",
    "carrier_frequency_speed": CODE / "a5_fingerprint" / "runtime"
    / "carrier_frequency_speed_receipt.json",
    "gauge_kinetic_invariant_forms": CODE / "e9_kinetic"
    / "gauge_kinetic_invariant_forms.certificate.json",
    "oriented_face_bracket_selector": CODE / "b14_jacobi"
    / "oriented_face_bracket_selector.certificate.json",
    "invariant_metric_phase": CODE / "b14_jacobi"
    / "invariant_metric_phase.certificate.json",
}

LEAN_RECEIPTS = {
    "A2HolonomyBridge": LEAN_SCREEN / "A2HolonomyBridge.lean",
    "A5OPH": LEAN_SCREEN / "A5OPH.lean",
    "A5CharacterField": LEAN_SCREEN / "A5CharacterField.lean",
    "A5SixAxes": LEAN_SCREEN / "A5SixAxes.lean",
    "Z6Exact": LEAN_SCREEN / "Z6Exact.lean",
    "Z6Descent": LEAN_SCREEN / "Z6Descent.lean",
    "A5CouplingSymmetry": LEAN_SCREEN / "A5CouplingSymmetry.lean",
    "A5PortAction": LEAN_SCREEN / "A5PortAction.lean",
    "PortFrameGram": LEAN_SCREEN / "PortFrameGram.lean",
    "ExteriorSelection": LEAN_SCREEN / "ExteriorSelection.lean",
    "TimeOrderLedger": REPO / "Lean" / "Time" / "TimeOrderLedger.lean",
    "ObserverHistory": REPO / "Lean" / "Time" / "ObserverHistory.lean",
    "ClockReadout": REPO / "Lean" / "Time" / "ClockReadout.lean",
    "WorldlineRealization": REPO / "Lean" / "Time" / "WorldlineRealization.lean",
    "ProperTimeCalibration": REPO / "Lean" / "Time" / "ProperTimeCalibration.lean",
    "ClockComparison": REPO / "Lean" / "Time" / "ClockComparison.lean",
    "ConsensusTower": REPO / "Lean" / "Tower" / "ConsensusTower.lean",
    "PublicWorldQuotient": REPO / "Lean" / "Tower"
    / "PublicWorldQuotient.lean",
    "FixedPointEndpoint": REPO / "Lean" / "Tower"
    / "FixedPointEndpoint.lean",
    "CanonicalLorentzModule": REPO / "Lean" / "Geometry"
    / "CanonicalLorentzModule.lean",
    "CelestialNullCone": REPO / "Lean" / "Geometry"
    / "CelestialNullCone.lean",
    "ObserverFrameHyperboloid": REPO / "Lean" / "Geometry"
    / "ObserverFrameHyperboloid.lean",
    "ObserverRestSpace": REPO / "Lean" / "Geometry"
    / "ObserverRestSpace.lean",
    "EinsteinTensorBridge": REPO / "Lean" / "Geometry"
    / "EinsteinTensorBridge.lean",
    "LorentzOverlapCocycle": REPO / "Lean" / "Geometry"
    / "LorentzOverlapCocycle.lean",
    "EventGermDisplacement": REPO / "Lean" / "Geometry"
    / "EventGermDisplacement.lean",
    "CelestialSoldering": REPO / "Lean" / "Geometry"
    / "CelestialSoldering.lean",
    "EventFrameSoldering": REPO / "Lean" / "Geometry"
    / "EventFrameSoldering.lean",
    "SpatialReadbackSoldering": REPO / "Lean" / "Geometry"
    / "SpatialReadbackSoldering.lean",
    "FiniteCausalObserverNet": REPO / "Lean" / "QFT"
    / "FiniteCausalObserverNet.lean",
    "ObserverNetDescent": REPO / "Lean" / "QFT"
    / "ObserverNetDescent.lean",
    "RichFibreWitness": REPO / "Lean" / "QFT" / "RichFibreWitness.lean",
    "RichFibreRegionalNet": REPO / "Lean" / "QFT"
    / "RichFibreRegionalNet.lean",
    "PublicRecordAlgebra": REPO / "Lean" / "EventAlgebra"
    / "PublicRecordAlgebra.lean",
    "NoBroadcastingAdapter": REPO / "Lean" / "EventAlgebra"
    / "NoBroadcastingAdapter.lean",
    "PartitionAverageCP": REPO / "Lean" / "EventAlgebra"
    / "PartitionAverageCP.lean",
    "TwoScalePublicRepair": REPO / "Lean" / "EventAlgebra"
    / "TwoScalePublicRepair.lean",
    "PoissonizedRepair": REPO / "Lean" / "Thermodynamics"
    / "PoissonizedRepair.lean",
    "PoissonizedRepairOperatorExp": REPO / "Lean" / "Thermodynamics"
    / "PoissonizedRepairOperatorExp.lean",
    "ConditionalExpectationGenerator": REPO / "Lean" / "Dynamics"
    / "ConditionalExpectationGenerator.lean",
    "ChoiCPTP": REPO / "Lean" / "Dynamics" / "ChoiCPTP.lean",
    "PublicMarkov": REPO / "Lean" / "Dynamics" / "PublicMarkov.lean",
    "PublicAutomorphism": REPO / "Lean" / "Dynamics"
    / "PublicAutomorphism.lean",
    "PrivateInner": REPO / "Lean" / "Dynamics" / "PrivateInner.lean",
    "FiniteBornFrame": REPO / "Lean" / "EventAlgebra"
    / "FiniteBornFrame.lean",
    "FiniteEffectClosureBoundary": REPO / "Lean" / "EventAlgebra"
    / "FiniteEffectClosureBoundary.lean",
    "Robertson": REPO / "Lean" / "EventAlgebra" / "Robertson.lean",
    "Superselection": REPO / "Lean" / "EventAlgebra"
    / "Superselection.lean",
    "ExteriorComponentBridge": LEAN_SCREEN / "ExteriorComponentBridge.lean",
    "QuantumMatterIntegration": LEAN_SCREEN / "QuantumMatterIntegration.lean",
    "B10EdgeCenterAction": LEAN_SCREEN / "B10EdgeCenterAction.lean",
    "HolonomyInterference": LEAN_SCREEN / "HolonomyInterference.lean",
    "FiniteConditionalRepair": REPO / "Lean" / "Thermodynamics"
    / "FiniteConditionalRepair.lean",
    "StationaryRealization": REPO / "Lean" / "Thermodynamics"
    / "StationaryRealization.lean",
    "FirstLawIdentity": REPO / "Lean" / "Thermodynamics"
    / "FirstLawIdentity.lean",
    "FluctuationTheorems": REPO / "Lean" / "Thermodynamics"
    / "FluctuationTheorems.lean",
    "CapFirstLaw": REPO / "Lean" / "Thermodynamics"
    / "CapFirstLaw.lean",
    "EinsteinPremiseLink": REPO / "Lean" / "Thermodynamics"
    / "EinsteinPremiseLink.lean",
    "PartitionPinchingCP": REPO / "Lean" / "EventAlgebra"
    / "PartitionPinchingCP.lean",
    "RegionalContinuity": LEAN_SCREEN / "RegionalContinuity.lean",
    "DiscreteGauss": LEAN_SCREEN / "DiscreteGauss.lean",
    "ProtectedCharge": REPO / "Lean" / "Dynamics" / "ProtectedCharge.lean",
    "WardLimitManifest": REPO / "Lean" / "Dynamics" / "WardLimitManifest.lean",
    "GreenKubo": REPO / "Lean" / "Thermodynamics" / "GreenKubo.lean",
    "GraphDiffusion": REPO / "Lean" / "Thermodynamics" / "GraphDiffusion.lean",
    "DependencyCone": REPO / "Lean" / "ObserverPatchHolography"
    / "Locality" / "DependencyCone.lean",
    "NoSignalling": REPO / "Lean" / "ObserverPatchHolography"
    / "Locality" / "NoSignalling.lean",
    "AdaptiveScheduler": REPO / "Lean" / "ObserverPatchHolography"
    / "Locality" / "AdaptiveScheduler.lean",
    "PathGibbs": REPO / "Lean" / "InformationProjection" / "PathGibbs.lean",
    "DiscreteEulerLagrange": REPO / "Lean" / "Variational"
    / "DiscreteEulerLagrange.lean",
    "DiscreteNoether": REPO / "Lean" / "Variational"
    / "DiscreteNoether.lean",
    "FiniteHistoryBridge": REPO / "Lean" / "Variational"
    / "FiniteHistoryBridge.lean",
    "RealizedHistoryLegendreNoGo": REPO / "Lean" / "Variational"
    / "RealizedHistoryLegendreNoGo.lean",
    "RecordMajorization": REPO / "Lean" / "EventAlgebra"
    / "RecordMajorization.lean",
    "SpectralEntropyBoundary": REPO / "Lean" / "EventAlgebra"
    / "SpectralEntropyBoundary.lean",
    "PortGramRepairBand": LEAN_SCREEN / "PortGramRepairBand.lean",
    "PortGramRepairCovariance": LEAN_SCREEN / "PortGramRepairCovariance.lean",
    "PrimitivePortFrameQuotient": LEAN_SCREEN / "PrimitivePortFrameQuotient.lean",
    "PortGramA5Isometry": LEAN_SCREEN / "PortGramA5Isometry.lean",
    "RepairWordCarrierReadout": LEAN_SCREEN / "RepairWordCarrierReadout.lean",
    "SeamCurrentCarrierQuotient": LEAN_SCREEN / "SeamCurrentCarrierQuotient.lean",
    "A5CarrierClassBand": LEAN_SCREEN / "A5CarrierClassBand.lean",
    "CarrierFrequencySpeed": LEAN_SCREEN / "CarrierFrequencySpeed.lean",
    "GaugeKineticInvariantForms": LEAN_SCREEN / "GaugeKineticInvariantForms.lean",
    "OrientedFaceBracketSelector": LEAN_SCREEN / "OrientedFaceBracketSelector.lean",
    "OrientedFaceInvariantMetric": LEAN_SCREEN / "OrientedFaceInvariantMetric.lean",
    "OperationalOverlapEvidence": REPO / "Lean" / "QFT"
    / "OperationalOverlapEvidence.lean",
    "CommonReferenceObstruction": REPO / "Lean" / "Thermodynamics"
    / "CommonReferenceObstruction.lean",
    "FiniteWebBornNoGo": REPO / "Lean" / "EventAlgebra"
    / "FiniteWebBornNoGo.lean",
    "SourceContextTomographyNoGo": REPO / "Lean" / "QFT"
    / "SourceContextTomographyNoGo.lean",
    "SourcePhaseLiftBridge": REPO / "Lean" / "QFT"
    / "SourcePhaseLiftBridge.lean",
    "ConjugationGauge": REPO / "Lean" / "QFT"
    / "ConjugationGauge.lean",
    "RepairCurrentOrientation": REPO / "Lean" / "Thermodynamics"
    / "RepairCurrentOrientation.lean",
    "SourceOrientedCompletion": REPO / "Lean" / "QFT"
    / "SourceOrientedCompletion.lean",
    "TwoFactorHistoryBinding": REPO / "Lean" / "QFT"
    / "TwoFactorHistoryBinding.lean",
}

DEFAULT_OUT = RUNS / "status" / "postdiction_ledger.json"
DEFAULT_MD = REPO / "docs" / "POSTDICTION_LEDGER.md"


def _load(key: str, override: Path | None = None) -> dict[str, Any]:
    path = override or PARENTS[key]
    if not path.exists():
        raise SystemExit(f"postdiction ledger parent missing: {key} at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _lean_receipt(
    *modules: str,
    declarations: dict[str, tuple[str, ...]] | None = None,
) -> list[str]:
    """Resolve Lean modules and fail closed on missing named declarations."""

    refs = []
    declarations = declarations or {}
    unknown = set(declarations) - set(modules)
    if unknown:
        raise SystemExit(
            "postdiction ledger declaration map names unbound modules: "
            + ", ".join(sorted(unknown))
        )
    for m in modules:
        path = LEAN_RECEIPTS[m]
        if not path.exists():
            raise SystemExit(f"postdiction ledger Lean receipt missing: {path}")
        source = path.read_text(encoding="utf-8")
        for declaration in declarations.get(m, ()):
            pattern = (
                rf"(?m)^\s*(?:theorem|lemma)\s+{re.escape(declaration)}"
                rf"(?:\s|:)"
            )
            if re.search(pattern, source) is None:
                raise SystemExit(
                    "postdiction ledger Lean declaration missing: "
                    f"{m}.{declaration} in {path}"
                )
        refs.append(path.relative_to(REPO).as_posix())
    return refs


def _rel(key: str) -> str:
    return PARENTS[key].relative_to(REPO).as_posix()


def _canonical_self_digest(payload: dict[str, Any]) -> str:
    body = {key: value for key, value in payload.items() if key != "receipt_sha256"}
    raw = (
        json.dumps(
            body,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _certificate_self_digest(payload: dict[str, Any]) -> str:
    body = {
        key: value
        for key, value in payload.items()
        if key != "certificate_sha256"
    }
    raw = json.dumps(
        body,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _quantum_carrier_status_row(packet: dict[str, Any]) -> dict[str, Any]:
    core = {
        key: value for key, value in packet.items() if key != "receipt_sha256"
    }
    digest = "sha256:" + hashlib.sha256(
        json.dumps(
            core,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    expected_verdicts = {
        "photon": "NOT_EVALUABLE_NO_SOURCE_SELECTED_MAXWELL_QUANTUM_SECTOR",
        "gluon": "NOT_EVALUABLE_NO_QCD",
        "graviton": "NOT_EVALUABLE_NO_INHABITED_EINSTEIN_QUANTUM_CARRIER",
    }
    rows = packet.get("rows")
    if not isinstance(rows, list):
        raise SystemExit("quantum-carrier status rows are absent")
    by_id = {
        row.get("carrier_id"): row for row in rows if isinstance(row, dict)
    }
    if (
        packet.get("schema") != "oph.quantum_carrier_status.v2"
        or packet.get("github_issue") != 552
        or packet.get("status") != "THREE_ROW_EXPLICIT_NOT_EVALUABLE"
        or packet.get("receipt_sha256") != digest
        or packet.get("comparison_values_consumed") is not False
        or packet.get("blind_prediction_eligible") is not False
        or packet.get("target_named_status_rows") is not True
        or packet.get("all_rows_at_allowed_exit") is not True
        or packet.get("continuum_spacetime_dimension") != 4
        or packet.get("classical_mode_vector_order")
        != ["photon", "gluon", "graviton"]
        or packet.get("classical_mode_vector") != [2, 16, 2]
        or set(by_id) != set(expected_verdicts)
    ):
        raise SystemExit("quantum-carrier status packet left its typed boundary")
    expected_multiplicities = {
        "photon": (1, "u1_lie_algebra_generator", 2),
        "gluon": (8, "su3_adjoint_generator", 16),
        "graviton": (1, "symmetric_metric_tensor_field", 2),
    }
    for carrier_id, verdict in expected_verdicts.items():
        row = by_id[carrier_id]
        capabilities = row.get("capabilities", {})
        baseline = row.get("classical_baseline", {})
        factor, role, total = expected_multiplicities[carrier_id]
        if (
            row.get("verdict") != verdict
            or row.get("verdict_class") != "EXPLICIT_NOT_EVALUABLE"
            or row.get("particle_promotion_allowed") is not False
            or not isinstance(row.get("blocking_frontier"), list)
            or not row["blocking_frontier"]
            or baseline.get("particle_claim") is not False
            or baseline.get("continuum_spacetime_dimension") != 4
            or baseline.get("multiplicity_factor") != factor
            or baseline.get("multiplicity_role") != role
            or baseline.get("exact_total_mode_count") != total
            or "gauge_algebra_dimension" in baseline
            or capabilities.get("state_space", {}).get(
                "physical_quantum_object_available"
            )
            is not False
            or capabilities.get("spectral_object", {}).get(
                "positive_physical_quantum_object_available"
            )
            is not False
            or capabilities.get("physical_current_residue", {}).get(
                "nonzero_positive_residue_available"
            )
            is not False
        ):
            raise SystemExit(
                f"quantum-carrier row {carrier_id} left its typed boundary"
            )
    return {
        "artifact_ref": _rel("quantum_carrier_status"),
        "classical_mode_vector": packet["classical_mode_vector"],
        "classical_mode_vector_order": packet["classical_mode_vector_order"],
        "receipt_sha256": packet["receipt_sha256"],
        "rows": [
            {
                "blocking_frontier": by_id[carrier_id]["blocking_frontier"],
                "carrier_id": carrier_id,
                "verdict": by_id[carrier_id]["verdict"],
            }
            for carrier_id in ("photon", "gluon", "graviton")
        ],
        "scope": (
            "The exact (2,16,2) vector contains differently typed conditional "
            "four-dimensional propagating-mode totals. Every quantum-particle "
            "row is explicitly not evaluable on the pinned declared corpus. "
            "The packet is target-named, comparison-value-free, and ineligible "
            "as a blind prediction."
        ),
    }


def _forced_structure(
    matter: dict[str, Any],
    matter_menu: dict[str, Any],
    port_current: dict[str, Any],
    axis_center_descent: dict[str, Any],
    carrier_modes: dict[str, Any],
    carrier_class: dict[str, Any],
    carrier_frequency: dict[str, Any],
    gauge_kinetic: dict[str, Any],
    oriented_face: dict[str, Any],
    invariant_metric: dict[str, Any],
) -> list[dict[str, Any]]:
    spectrum = matter["realized_package"]["charge_spectrum"]
    sm_spectrum = {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 6}
    scan = matter["matter_selection_scan"]
    classification = matter_menu["subset_classification"]
    scan_agreement = matter_menu["scan_agreement"]
    menu_verdicts = matter_menu["verdicts"]
    lean_cross_reference = classification["lean_cross_reference"]
    subsets_enumerated = classification["subsets_enumerated"]
    component_count = subsets_enumerated.bit_length() - 1
    survivors = classification["survivors"]
    survivor_dimensions = {row["dimension"] for row in survivors}
    if (
        2**component_count != subsets_enumerated
        or classification["survivor_count"] != len(survivors)
        or classification["survivor_count"] != scan["survivor_count"]
        or subsets_enumerated != scan["subsets_enumerated"]
        or not classification["survivors_are_conjugate_pair"]
        or not scan["survivors_are_conjugate_pair"]
        or survivor_dimensions != {matter["realized_package"]["dimension"]}
        or scan["derived_block_charges"]
        != matter_menu["declared_algebra"]["derived_block_charges"]
        or not scan_agreement["lean_masks_match_ledger"]
        or not scan_agreement["matter_lift_scan_matches_ledger"]
        or menu_verdicts["menu_completeness_inside_declared_algebra"] != "exact"
        or lean_cross_reference["agreement"] is not True
    ):
        raise SystemExit(
            "matter-menu, matter-lift, and Lean cross-reference parents disagree"
        )

    frequency_theorem = carrier_frequency.get("generic_exact_theorem", {})
    frequency_physical = carrier_frequency.get("physical_boundary", {})
    frequency_exposure = carrier_frequency.get("exposure_boundary", {})
    frequency_supports = carrier_frequency.get("exact_support_instantiations", {})
    if (
        carrier_frequency.get("schema") != "oph.carrier_frequency_speed.v1"
        or carrier_frequency.get("status")
        != "EXACT_POSITIVE_TIGHT_FRAME_FREQUENCY_CONTRACTION__FZ11_FZ12_INSTANTIATED__PHYSICAL_BRIDGES_OPEN"
        or carrier_frequency.get("receipt_sha256")
        != _canonical_self_digest(carrier_frequency)
        or frequency_theorem.get("certified_upper_constant") != "1"
        or frequency_theorem.get("feature_identity")
        != "Lambda_a(k)=||Phi_a(k)||^2"
        or frequency_theorem.get("frequency_bound")
        != "|Omega_a(k)-Omega_a(p)| <= |k-p|"
        or frequency_supports.get("vertex12", {}).get("tight_constant") != "4"
        or frequency_supports.get("edge30", {}).get("tight_constant") != "10"
        or any(value is not False for value in frequency_physical.values())
        or frequency_exposure.get("comparison_inputs") != []
        or any(
            frequency_exposure.get(key) is not False
            for key in (
                "comparison_data_read",
                "public_measurement_read",
                "score_emitted",
                "verdict_emitted",
            )
        )
        or carrier_frequency.get("branch_bindings", {}).get(
            "new_prediction_payload"
        )
        is not False
    ):
        raise SystemExit(
            "carrier-frequency parent has left its exact auxiliary boundary"
        )

    kinetic_families = gauge_kinetic.get("families", {})
    if (
        gauge_kinetic.get("schema")
        != "oph.e9.gauge_kinetic_invariant_forms.v1"
        or gauge_kinetic.get("issue") != 716
        or gauge_kinetic.get("certificate_sha256")
        != _certificate_self_digest(gauge_kinetic)
        or kinetic_families.get("F", {}).get("carrier_invariant_dimension") != 3
        or kinetic_families.get("F", {}).get("ad_invariant_dimension") != 2
        or kinetic_families.get("G", {}).get("carrier_invariant_dimension") != 3
        or kinetic_families.get("G", {}).get("ad_invariant_dimension") != 2
        or kinetic_families.get("P", {}).get("carrier_invariant_dimension") != 2
        or kinetic_families.get("P", {}).get("ad_invariant_dimension") != 2
        or gauge_kinetic.get("mirror_common_control", {}).get(
            "not_source_selected"
        )
        is not True
    ):
        raise SystemExit(
            "gauge-kinetic invariant-form parent left its exact bounded boundary"
        )

    face_source = oriented_face.get("source_face_bracket", {})
    face_jacobi = oriented_face.get("jacobi_failure", {})
    face_selector = oriented_face.get(
        "orthogonal_compact_locus_discriminator", {}
    )
    face_norms = oriented_face.get("endpoint_norm_robustness", {})
    if (
        oriented_face.get("schema")
        != "oph.b14.oriented_face_bracket_selector.v1"
        or oriented_face.get("issue") != 705
        or oriented_face.get("certificate_sha256")
        != _certificate_self_digest(oriented_face)
        or face_source.get("oriented_face_count") != 20
        or face_source.get("identity")
        != "B_face = 60 * R13 exactly in all 12*12*12 tensor coordinates"
        or face_jacobi.get("nonzero_count") != 240
        or face_jacobi.get("positive_count") != 120
        or face_jacobi.get("negative_count") != 120
        or face_selector.get("unique_nearest_family") != "G"
        or face_selector.get("metric_is_source_derived") is not False
        or face_selector.get(
            "minimum_hs_or_jacobi_repair_is_source_derived"
        )
        is not False
        or face_norms.get("repair_norm_or_minimization_rule_is_source_derived")
        is not False
        or face_norms.get("l1", {}).get(
            "unique_nearest_compact_family_by_minimum_or_infimum"
        )
        != "G"
        or face_norms.get("linfinity", {}).get(
            "unique_nearest_compact_family"
        )
        != "G"
        or "infimum" not in face_norms.get("l1", {}).get(
            "families", {}
        ).get("F", {}).get("compact_attainment", "")
    ):
        raise SystemExit(
            "oriented-face bracket parent left its conditional discriminator boundary"
        )
    exterior_declarations = tuple(lean_cross_reference["theorems"])
    exterior_path = LEAN_RECEIPTS["ExteriorSelection"].relative_to(REPO).as_posix()
    if lean_cross_reference["file"] != exterior_path:
        raise SystemExit("matter-menu Lean source path does not match the ledger binding")
    realized_fields = matter["realized_package"]["fields"]
    field_order = ("Q", "u_c", "d_c", "L", "e_c")
    field_summary = ", ".join(
        f"{name}: {realized_fields[name]['charge']} x"
        f"{realized_fields[name]['dimension']}"
        for name in field_order
    )

    current_map = port_current["port_to_generator_map"]
    current_closure = port_current["closure"]
    derived_dimensions = current_closure["derived_block_dimensions"]
    color_adjoint_dimension = derived_dimensions["even_block_su3"]
    weak_adjoint_dimension = derived_dimensions["kernel_block_so3"]
    abelian_dimension = current_closure["center_dimension"]
    if (
        port_current["claim_boundary"]["status"]
        != "proved_conditional_on_declared_response_representation"
        or port_current["physical_source_gate"]["passed"] is not False
        or port_current["physical_source_gate"][
            "target_blind_impulse_readback_recomputed"
        ]
        is not True
        or port_current["source_definedness"][
            "response_model_declared_as_branch_premise"
        ]
        is not True
        or port_current["source_definedness"][
            "physical_response_source_bound"
        ]
        is not False
        or current_map["compact_lie_type"]
        != "u(3) (+) so(3) = u(1) (+) su(3) (+) su(2)"
        or current_map["block_dimensions_verified"] is not True
        or current_map["injective"] is not True
        or color_adjoint_dimension + weak_adjoint_dimension
        != current_closure["derived_dimension"]
        or (
            color_adjoint_dimension
            + weak_adjoint_dimension
            + abelian_dimension
            != current_map["image_real_dimension"]
        )
    ):
        raise SystemExit("port-current parent does not realize the declared product algebra")
    adjoint_branching = {
        "color_adjoint": [color_adjoint_dimension, 1, 0],
        "weak_adjoint": [1, weak_adjoint_dimension, 0],
        "abelian": [1, 1, 0],
        "mixed_xy_bifundamental_dimension": 0,
    }

    descent_gate = axis_center_descent["physical_global_form_gate"]
    declared_loop = axis_center_descent["carrier_deck_and_declared_loop_system"]
    tensor_kernel = axis_center_descent["kernel_on_realized_tensors"]
    effective_image = axis_center_descent["maximal_effective_image"]
    z6_bridge = axis_center_descent["two_z6_constructions"]
    if (
        axis_center_descent["claim_boundary"]["status"]
        != "conditional_exact_arithmetic_with_physical_global_form_open"
        or descent_gate["passed"] is not False
        or descent_gate["laboratory_global_form_attachment"] is not False
        or descent_gate["theta_periodicity_derived"] is not False
        or descent_gate["axis_relation_lattice_source_selected"] is not False
        or descent_gate["complete_character_category_source_derived"] is not False
        or descent_gate["same_source_loop_to_tensor_kernel_identification"]
        is not False
        or tensor_kernel["kernel_order"] != 6
        or tensor_kernel["matches_emitted_kernel_data"] is not True
        or effective_image["group"] != "(SU(3) x SU(2) x U(1)) / Z6"
        or declared_loop["six_axis_class_group_order"] != 6
        or declared_loop["declared_coefficient_system_menu"] != list(range(6))
        or declared_loop["axis_relation_lattice_source_selected"] is not False
        or z6_bridge["physical_loop_intertwiner_derived"] is not False
        or z6_bridge["conditional_algebraic_intertwiner_verified"] is not True
    ):
        raise SystemExit(
            "axis-centre descent parent has left its conditional arithmetic boundary"
        )

    carrier_by_id = {
        row["carrier_id"]: row for row in carrier_modes["carriers"]
    }
    if set(carrier_by_id) != {"photon", "gluon", "graviton"}:
        raise SystemExit("carrier-mode packet does not contain the expected rows")
    if not all(
        row["classical_carrier_gate"]["passed"]
        and not row["quantum_particle_gate"]["passed"]
        and row["hard_quadratic_mass_parameter_squared"] == 0
        for row in carrier_by_id.values()
    ):
        raise SystemExit("carrier-mode packet has left its declared boundary")

    class_statements = carrier_class.get("class_statements", {})
    eighth_order = carrier_class.get("eighth_order", {})
    if (
        carrier_class.get("schema") != "oph.carrier_class_dispersion.v1"
        or class_statements.get("sign_law") != "C4 < 0 for every member"
        or "10/21" not in class_statements.get("isotropic_floor", "")
        or "[-16/135, 16/75]" not in class_statements.get("rank_six_band", "")
        or eighth_order.get("kernel_universal_constant") != "256/75"
        or eighth_order.get("cross_order_identity")
        != "5 D6 B0 = 12 B6 D0"
        or eighth_order.get("cross_order_ratio_identity")
        != "D6/D0 = (12/5)(B6/B0)"
        or eighth_order.get("multi_radius_negative_control", {}).get(
            "lock_residual"
        )
        != "-57344/10360225"
        or eighth_order.get("zero_mixture_control", {}).get("B6_over_B0") != "0"
        or eighth_order.get("zero_mixture_control", {}).get("D6_over_D0") != "0"
    ):
        raise SystemExit(
            "carrier-class parent does not carry the exact sign, floor, band, "
            "and division-free eighth-order lock"
        )

    rows = [
        {
            "id": "gauge_lie_algebra",
            "statement": (
                "The certified twelve-port carrier has module 1+3+3'+5 and "
                "one fixed line. Complete compact port response from A1 and "
                "endogenous proper-carrier transport from A2 force the abstract "
                "Lie type u(1)+su(2)+su(3). Target-blind impulse and readback "
                "separately determine R=-J. The charged-double-triplet matrices "
                "are an exact declared witness, while ordered source tomography "
                "and same-current holonomy remain open"
            ),
            "observed_counterpart": "Standard Model gauge Lie algebra su(3)+su(2)+u(1)",
            "match": "axiom-forced abstract Lie type; conditional matrix witness",
            "artifact_ref": _rel("port_current"),
            "machine_checked_steps": (
                "the carrier action and fixed-space dimension are exact; Lean "
                "checks the A2 holonomy-to-inner-action bridge, the centreless "
                "four-factor fixed-space exclusion, triviality of A5-actions "
                "on at most four objects, and unique "
                "partitions 11 = 8+3 and 12 = 3+3+3+3 over the compact-simple "
                "dimension list {3, 8, 10}; no compact semisimple algebra in "
                "dimensions 1, 2, 4, 5, 7; the characteristic-centre step; "
                "A5 not a subgroup of SU(2); the noncentrality witness "
                "[iS, iT] != 0; Galois stability of the two three-dimensional "
                "characters over Q(sqrt 5) with the centre-dimension list "
                "{0, 1, 5, 6, 7, 11, 12}; six-axis 2-transitivity and the "
                "dimension count of the dimension-six branch"
            ),
            "lean_declarations": {
                "A2HolonomyBridge": [
                    "internalImplementation_of_holonomy",
                    "four_factor_fixed_dimension_ne_one",
                    "compact_product_dimensions_of_fixed_space",
                ],
                "A5OPH": [
                    "sum_eq_eleven",
                    "sum_eq_twelve",
                    "action_trivial_of_card_le_four",
                    "sum_not_mem_excluded",
                    "quintet_noncentral",
                ],
                "A5CharacterField": [
                    "multiplicities_equal_of_galoisStable",
                    "centreDim_mem_trichotomy_list",
                ],
                "A5SixAxes": [
                    "two_transitive",
                    "V5_irreducible",
                    "no_three_plus_three_split",
                ],
            },
            "lean_receipts": _lean_receipt(
                "A2HolonomyBridge",
                "A5OPH",
                "A5CharacterField",
                "A5SixAxes",
                declarations={
                    "A2HolonomyBridge": (
                        "internalImplementation_of_holonomy",
                        "four_factor_fixed_dimension_ne_one",
                        "compact_product_dimensions_of_fixed_space",
                    ),
                    "A5OPH": (
                        "sum_eq_eleven",
                        "sum_eq_twelve",
                        "action_trivial_of_card_le_four",
                        "sum_not_mem_excluded",
                        "quintet_noncentral",
                    ),
                    "A5CharacterField": (
                        "multiplicities_equal_of_galoisStable",
                        "centreDim_mem_trichotomy_list",
                    ),
                    "A5SixAxes": (
                        "two_transitive",
                        "V5_irreducible",
                        "no_three_plus_three_split",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the abstract theorem uses the A1 faithful complete compact "
                "response and A2 endogenous holonomy clauses. The direct "
                "receipt verifies a declared charged-double-triplet witness "
                "and does not reconstruct it from ordered source histories or "
                "identify a laboratory current. Compact reductivity and the "
                "compact-simple classification are declared classical inputs"
            ),
            "paper_ref": "Standard Model gauge paper, Compact-Lie trichotomy section",
        },
        {
            "id": "oriented_face_nearest_compact_discriminator",
            "statement": (
                "The equal-weight cyclic bracket of the twenty pinned oriented "
                "faces is exactly 60 times Reynolds basis vector R13 and is "
                "A5-equivariant. It fails Jacobi in exactly 240 of the 2640 "
                "independent output/input-triple coordinates, split 120 at +1 "
                "and 120 at -1. Conditional on the displayed 792-coordinate "
                "upper-triangular convention and the certified compact locus, "
                "exact primal-dual certificates make G the winning family for "
                "three edit norms. Distances to G, F, P are respectively "
                "30(sqrt(5)-1), 60, 60 in L1; "
                "(615-123 sqrt(5))/22, (615+123 sqrt(5))/22, 45 in squared "
                "L2; and (5-sqrt(5))/10, sqrt(5)/5, 1/2 in Linfinity. "
                "The F L1 value is an unattained compact-family infimum"
            ),
            "observed_counterpart": (
                "a finite source-incidence discriminator among the three compact "
                "bracket families"
            ),
            "match": (
                "exact conditional finite discriminator; no source repair law"
            ),
            "artifact_ref": _rel("oriented_face_bracket_selector"),
            "receipt_sha256": oriented_face["certificate_sha256"],
            "lean_declarations": {
                "OrientedFaceBracketSelector": [
                    "face_bracket_eq_sixty_r13",
                    "jacobi_failure_witness",
                    "unique_nearest_G",
                    "three_norm_unique_nearest_G",
                ],
            },
            "lean_receipts": _lean_receipt(
                "OrientedFaceBracketSelector",
                declarations={
                    "OrientedFaceBracketSelector": (
                        "face_bracket_eq_sixty_r13",
                        "jacobi_failure_witness",
                        "unique_nearest_G",
                        "three_norm_unique_nearest_G",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the equal-weight oriented-face construction is a declared "
                "deterministic rule applied to the pinned incidence orientation, "
                "not a rule forced by the OPH axioms. The comparison adds "
                "three basis-dependent coordinate norms; neither a norm, "
                "minimum-distance repair, nor any Jacobi-repair dynamics is "
                "source-derived. Exact code proves the 792-coordinate "
                "optimization, while Lean checks the serialized radical "
                "values and order. The result compares only with the "
                "classified compact locus and does not close B14 or select a "
                "physical gauge bracket"
            ),
            "paper_ref": "Standard Model gauge paper, bracket-selection boundary",
        },
        {
            "id": "invariant_metric_phase_diagram",
            "statement": (
                "The carrier splits multiplicity-free into 1+3+3'+5 and the "
                "commutant of the port action is exactly four-dimensional, "
                "spanned by the four symmetric spectral projectors, so the "
                "positive sector-scale cone is the complete family of "
                "invariant carrier inner products. Every induced bracket "
                "metric is channel-diagonal, and the squared distances from "
                "the face bracket to the classified compact families are "
                "exact three-term Laurent forms in the sector scales with "
                "the fixed-sector scale absent. P is strictly excluded for "
                "every invariant metric; every sector-balanced metric and "
                "every metric with beta/delta in [1/50, 6] selects G "
                "uniquely for all gamma and delta; F occupies the nonempty "
                "side d_F^2<d_G^2 of the exact three-scale tie surface, "
                "with witness (8,1,1); "
                "and Galois conjugation with the sector swap maps d_G to "
                "d_F exactly. A non-carrier-induced channel reweighting "
                "reverses the balanced-point selection"
            ),
            "observed_counterpart": (
                "a metric-robust phase diagram for the finite source-incidence "
                "discriminator over the complete invariant carrier-metric cone"
            ),
            "match": (
                "exact conditional phase diagram; no source metric or repair law"
            ),
            "artifact_ref": _rel("invariant_metric_phase"),
            "receipt_sha256": invariant_metric["certificate_sha256"],
            "lean_declarations": {
                "OrientedFaceInvariantMetric": [
                    "reference_G",
                    "dG2_lt_dP2",
                    "dF2_lt_dP2",
                    "balanced_unique_nearest_G",
                    "box_unique_nearest_G",
                    "F_wins_at_witness",
                ],
            },
            "lean_receipts": _lean_receipt(
                "OrientedFaceInvariantMetric",
                declarations={
                    "OrientedFaceInvariantMetric": (
                        "reference_G",
                        "dG2_lt_dP2",
                        "dF2_lt_dP2",
                        "balanced_unique_nearest_G",
                        "box_unique_nearest_G",
                        "F_wins_at_witness",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the closed forms are certificate content derived from the "
                "pinned tensors by the independently replayed producer; Lean "
                "proves the phase consequences as quantified real theorems "
                "over those forms. The nearest-point repair rule and the "
                "restriction to carrier-induced metrics are declared "
                "discriminator choices, the comparison is conditional on the "
                "classified compact locus, and no metric, bracket, current, "
                "or physical gauge structure is source-selected"
            ),
            "paper_ref": "Standard Model gauge paper, bracket-selection boundary",
        },
        {
            "id": "gauge_kinetic_invariant_form_drop",
            "statement": (
                "For each certified compact bracket, exact ad-invariance of a "
                "carrier-projector quadratic form leaves one coefficient per "
                "simple factor. The F and G families reduce the three "
                "carrier-invariant weights to two: F imposes "
                "w(3+)=sqrt(5) w(5), while G imposes "
                "w(3-)=sqrt(5) w(5). The P control remains two-to-two. "
                "Simultaneous F/G invariance leaves one ray but is an extra "
                "mirror-common premise"
            ),
            "observed_counterpart": (
                "the finite invariant quadratic-form shape of candidate gauge "
                "kinetic terms"
            ),
            "match": "exact representation-level finite theorem",
            "artifact_ref": _rel("gauge_kinetic_invariant_forms"),
            "receipt_sha256": gauge_kinetic["certificate_sha256"],
            "lean_declarations": {
                "GaugeKineticInvariantForms": [
                    "f_exact_two_parameter",
                    "g_exact_two_parameter",
                    "p_exact_two_parameter",
                    "mirror_common_extra_premise_one_ray",
                ],
            },
            "lean_receipts": _lean_receipt(
                "GaugeKineticInvariantForms",
                declarations={
                    "GaugeKineticInvariantForms": (
                        "f_exact_two_parameter",
                        "g_exact_two_parameter",
                        "p_exact_two_parameter",
                        "mirror_common_extra_premise_one_ray",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the brackets and carrier projectors are supplied certified "
                "finite objects. The theorem selects neither bracket nor "
                "overall or relative coupling coefficients, and it constructs "
                "no source action, continuum field, or laboratory current. "
                "The one-ray intersection cannot be promoted without the "
                "additional simultaneous-invariance premise"
            ),
            "paper_ref": "Standard Model gauge paper, kinetic-action boundary",
        },
        {
            "id": "two_factor_constructed_history_binding",
            "statement": (
                "For two finite additive step groups, a Gibbs kernel constructed "
                "from the weighted sum a q1+b q2 factorizes exactly, and its "
                "history action is the corresponding sum of factor actions. "
                "If each factor cost has one nonconstant direction, equality of "
                "two such constructed transition kernels identifies both "
                "multiplier-weighted coefficients; with nonzero multipliers, "
                "only one common scaling remains. Exact instances bind both P "
                "factors on 3^6 steps and the complete 8+3 dimensional F family "
                "on 3^11 steps"
            ),
            "observed_counterpart": (
                "a finite multifactor gauge-history action with an identifiable "
                "relative kinetic coefficient"
            ),
            "match": (
                "exact constructed-kernel compatibility; no independent source law"
            ),
            "lean_declarations": {
                "TwoFactorHistoryBinding": [
                    "twoFactor_kernel_identifies_weighted_coefficients",
                    "twoFactor_kernel_only_common_multiplier_scaling",
                    "fullP_action_reproduces_law",
                    "fullP_kernel_relative_coefficients_identifiable",
                    "fFamily_action_reproduces_law",
                    "fFamily_kernel_relative_coefficients_identifiable",
                ],
            },
            "lean_receipts": _lean_receipt(
                "TwoFactorHistoryBinding",
                declarations={
                    "TwoFactorHistoryBinding": (
                        "twoFactor_kernel_identifies_weighted_coefficients",
                        "twoFactor_kernel_only_common_multiplier_scaling",
                        "fullP_action_reproduces_law",
                        "fullP_kernel_relative_coefficients_identifiable",
                        "fFamily_action_reproduces_law",
                        "fFamily_kernel_relative_coefficients_identifiable",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "both transition kernels are Gibbs kernels constructed from the "
                "displayed two-factor costs. The theorem symbolically handles "
                "the 729- and 177147-state groups without a large enumeration. "
                "It does not identify either kernel with an independently "
                "source-produced process, select either coefficient or their "
                "ratio, add the abelian sector, or supply physical units, a "
                "continuum field, laboratory current, or prediction"
            ),
            "paper_ref": "Standard Model gauge paper, multifactor history binding",
        },
        {
            "id": "global_form_z6",
            "statement": (
                "The common central kernel on every declared tensor is the "
                "order-six diagonal subgroup, so the maximal faithful image "
                "of that representation is (SU(3) x SU(2) x U(1))/Z6. "
                "The six-axis class has order six only after diagonal and "
                "zero-sum coefficient relations are declared. Source selection "
                "of those relations, a complete character category, and a "
                "same-source loop-to-kernel theorem remain open"
            ),
            "observed_counterpart": (
                "Standard Model global gauge-group form and its charge "
                "quantization pattern"
            ),
            "match": "exact conditional kernel and maximal faithful image",
            "artifact_refs": [
                _rel("matter_receipt"),
                _rel("axis_center_descent"),
            ],
            "lean_declarations": {
                "Z6Exact": [
                    "gauge_eq_kernel",
                    "residue_surjective",
                    "representative_formula",
                ],
                "Z6Descent": [
                    "kernel_on_realized_weights",
                    "four_admissible_global_forms",
                    "sixAxis_generator_maps_to_kernel_generator",
                    "sixAxisToKernel_intertwines_involutions",
                    "sixAxisToKernel_injective",
                    "sixAxisToKernel_range",
                ],
            },
            "lean_receipts": _lean_receipt(
                "Z6Exact",
                "Z6Descent",
                declarations={
                    "Z6Exact": (
                        "gauge_eq_kernel",
                        "residue_surjective",
                        "representative_formula",
                    ),
                    "Z6Descent": (
                        "kernel_on_realized_weights",
                        "four_admissible_global_forms",
                        "sixAxis_generator_maps_to_kernel_generator",
                        "sixAxisToKernel_intertwines_involutions",
                        "sixAxisToKernel_injective",
                        "sixAxisToKernel_range",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the result is exact for the declared matter table, central "
                "descent congruence, axis coefficient relations, and line "
                "lattice. The source current, physical matter action, relation-"
                "lattice selection, character completeness, loop-to-kernel "
                "identity, laboratory attachment, and continuum quantum field "
                "theory remain outside this result"
            ),
            "paper_ref": "Standard Model gauge paper, Z6 global-form section",
        },
        {
            "id": "hypercharge_spectrum",
            "statement": (
                f"Inside the declared {component_count}-component "
                "exterior-response algebra, an exhaustive scan of all "
                f"{subsets_enumerated} subsets selects exactly one "
                "unordered charge-conjugate pair of nonempty chiral "
                "anomaly-free rank-"
                f"{matter['realized_package']['dimension']} projectors. "
                "Primitive determinant balance "
                "fixes the block charges up to conjugation, and the selected "
                f"representative has multiset {{{field_summary}}}"
            ),
            "observed_counterpart": (
                "Standard Model one-generation hypercharge assignment"
            ),
            "realized_spectrum": spectrum,
            "match": "exact" if spectrum == sm_spectrum else "MISMATCH",
            "artifact_refs": [_rel("matter_receipt"), _rel("matter_menu")],
            "subset_count": subsets_enumerated,
            "survivor_count": classification["survivor_count"],
            "survivor_dimension": matter["realized_package"]["dimension"],
            "derived_block_charges": matter_menu["declared_algebra"][
                "derived_block_charges"
            ],
            "lean_declarations": {"ExteriorSelection": list(exterior_declarations)},
            "lean_receipts": _lean_receipt(
                "ExteriorSelection",
                declarations={"ExteriorSelection": exterior_declarations},
            ),
            "hypothesis_boundary": (
                "the selection is exhaustive inside the declared exterior "
                "algebra; completeness beyond that algebra, selection of one "
                "charge-conjugate representative, light-sector attachment, "
                "family multiplicity, scalar content, and laboratory "
                "identification remain separate"
            ),
            "paper_ref": "zoo paper, matter lift section",
        },
        {
            "id": "finite_quantum_limitation_suite",
            "statement": (
                "A supplied finite density state and two Hermitian observables "
                "obey the ordinary-commutator Robertson inequality, with exact "
                "noncommuting saturation and zero-variance controls. For a "
                "supplied projective partition, equality after block pinching "
                "is exactly equality of every trace statistic against the full "
                "sector-preserving commutant. Partition averaging is the "
                "distinct commutative projector-span readout, factors through "
                "pinching, is surjective onto the supplied partition's public "
                "projector-span algebra, and commutes with the partition "
                "commutant. An exact rank-two control keeps the maps distinct. "
                "The same pinching is the positive uniform random-unitary "
                "average over all 2^k independently signed block reflections. "
                "A separate orthogonal-density witness has failed support "
                "inclusion but zero raw relative entropy under the totalized "
                "matrix logarithm"
            ),
            "observed_counterpart": (
                "finite uncertainty, partition-relative superselection, and "
                "the support-aware spectral-information boundary"
            ),
            "match": (
                "exact bounded finite package; source and "
                "physical-instrument attachments remain open"
            ),
            "lean_declarations": {
                "Robertson": [
                    "finite_state_robertson_commutator",
                    "neg_I_mul_commutator_expectation_eq_readout",
                    "pauliX_pauliY_ne_pauliY_pauliX",
                    "pauli_xy_noncommuting_control",
                    "pauliZ_pauliX_ne_pauliX_pauliZ",
                    "pauli_z_zero_variance_control",
                ],
                "Superselection": [
                    "partitionOperationallyEquivalent_iff_pinching_eq",
                    "trace_mul_eq_zero_of_partitionOffDiagonal",
                    "partitionPinching_partitionCorner_eq_zero",
                    "partitionAverage_partitionCorner_eq_zero",
                    "trace_partitionCorner_mul_eq_zero_of_mem_span",
                ],
                "B10EdgeCenterAction": [
                    "partitionCenterAdaptor_after_blockReadout",
                    "partitionCenterAdaptor_surjective",
                    "partitionCenterAdaptor_commutes_with_block",
                    "rankTwo_pinching_ne_average",
                ],
                "RecordMajorization": [
                    "recordSignAverage_eq_partitionPinching",
                    "globalSignAverage_not_binary_pinching",
                ],
                "SpectralEntropyBoundary": [
                    "binary_orthogonal_density_receipt",
                    "totalizedRelativeEntropy_binary_orthogonal_eq_zero",
                    "supportAware_not_totalizedRelativeEntropy",
                ],
            },
            "lean_receipts": _lean_receipt(
                "Robertson",
                "Superselection",
                "B10EdgeCenterAction",
                "RecordMajorization",
                "SpectralEntropyBoundary",
                declarations={
                    "Robertson": (
                        "finite_state_robertson_commutator",
                        "neg_I_mul_commutator_expectation_eq_readout",
                        "pauliX_pauliY_ne_pauliY_pauliX",
                        "pauli_xy_noncommuting_control",
                        "pauliZ_pauliX_ne_pauliX_pauliZ",
                        "pauli_z_zero_variance_control",
                    ),
                    "Superselection": (
                        "partitionOperationallyEquivalent_iff_pinching_eq",
                        "trace_mul_eq_zero_of_partitionOffDiagonal",
                        "partitionPinching_partitionCorner_eq_zero",
                        "partitionAverage_partitionCorner_eq_zero",
                        "trace_partitionCorner_mul_eq_zero_of_mem_span",
                    ),
                    "B10EdgeCenterAction": (
                        "partitionCenterAdaptor_after_blockReadout",
                        "partitionCenterAdaptor_surjective",
                        "partitionCenterAdaptor_commutes_with_block",
                        "rankTwo_pinching_ne_average",
                    ),
                    "RecordMajorization": (
                        "recordSignAverage_eq_partitionPinching",
                        "globalSignAverage_not_binary_pinching",
                    ),
                    "SpectralEntropyBoundary": (
                        "binary_orthogonal_density_receipt",
                        "totalizedRelativeEntropy_binary_orthogonal_eq_zero",
                        "supportAware_not_totalizedRelativeEntropy",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the state, observables, and projective partition are supplied "
                "finite inputs. Pinching lands in the generally noncommutative "
                "commutant, while averaging lands in the commutative projector "
                "span. The adaptor is relative to that supplied partition. "
                "The sign average is an exact precursor, not a proof of "
                "spectral majorization. The totalized-log countermodel rejects "
                "only that naive architecture; a support-aware extended "
                "divergence, pinching Pythagoras, constrained maximum entropy, "
                "and the publicization information chain remain open on #685. "
                "No source rule selects that partition, a state, an observable, "
                "a detector algebra, or a public instrument, and no "
                "edge-to-partition identification is constructed"
            ),
            "paper_ref": "finite-event-algebra paper",
        },
        {
            "id": "finite_exterior_component_bridge",
            "statement": (
                "The Mathlib exterior basis on the declared five-mode carrier "
                "has 32 labels and binds the ten non-vacuum, non-top component "
                "rows to their dimensions, charges, parity, conjugation, "
                "square-zero creation, and anticommutation. One explicit typed "
                "map assigns those rows to supplied partition sectors and "
                "central weights. The supplied weights define a sixth-root "
                "character action on the component-labelled finite product of "
                "mapped projector ranges with the exact six-element tensor "
                "kernel. The supplied "
                "anomaly-free exterior-degree parity support is nontrivial, "
                "invariant, and detects the same kernel"
            ),
            "observed_counterpart": (
                "finite exclusion and one-generation central-weight structure"
            ),
            "match": (
                "exact bounded finite action; source selection and physical "
                "matter attachment open"
            ),
            "lean_declarations": {
                "ExteriorComponentBridge": [
                    "exterior_basis_label_count",
                    "bidegree_count_table",
                    "componentDegree_exact_nontrivial_menu",
                    "component_dimension_binding",
                    "component_charge_binding",
                    "component_parity_binding",
                    "component_conjugation_binding",
                    "creation_square_zero",
                    "creation_actions_anticommute",
                ],
                "QuantumMatterIntegration": [
                    "even_component_weights_eq_matterWeights",
                    "kernel_on_exterior_component_weights",
                    "fractional_singlet_mutation_collapses_component_kernel",
                    "coordinate_diagonal_not_partitionOffDiagonal",
                    "coordinate_nonzero_offDiagonal_control",
                    "declaredBlockReadout_eq_iff_operationallyEquivalent",
                    "kernel_on_mapped_component_weights",
                    "bridge_selection_is_parity_sector",
                ],
                "B10EdgeCenterAction": [
                    "mappedCentralAction_zero",
                    "mappedCentralAction_add",
                    "mappedCentralAction_neg_comp",
                    "mappedCentralAction_eq_id_iff_component_phases_zero",
                    "mappedCentralAction_eq_id_iff",
                    "mappedCentralAction_kernel_card",
                    "selectedMappedMatter_support_is_parity",
                    "selectedMappedMatter_nontrivial",
                    "mappedCentralAction_preserves_selected",
                    "kernel_on_selected_mapped_components",
                    "selectedMappedCentralAction_eq_id_iff",
                    "no_selected_central_parameter_realizes_parity_sign",
                ],
            },
            "lean_receipts": _lean_receipt(
                "ExteriorComponentBridge",
                "QuantumMatterIntegration",
                "B10EdgeCenterAction",
                declarations={
                    "ExteriorComponentBridge": (
                        "exterior_basis_label_count",
                        "bidegree_count_table",
                        "componentDegree_exact_nontrivial_menu",
                        "component_dimension_binding",
                        "component_charge_binding",
                        "component_parity_binding",
                        "component_conjugation_binding",
                        "creation_square_zero",
                        "creation_actions_anticommute",
                    ),
                    "QuantumMatterIntegration": (
                        "even_component_weights_eq_matterWeights",
                        "kernel_on_exterior_component_weights",
                        "fractional_singlet_mutation_collapses_component_kernel",
                        "coordinate_diagonal_not_partitionOffDiagonal",
                        "coordinate_nonzero_offDiagonal_control",
                        "declaredBlockReadout_eq_iff_operationallyEquivalent",
                        "kernel_on_mapped_component_weights",
                        "bridge_selection_is_parity_sector",
                    ),
                    "B10EdgeCenterAction": (
                        "mappedCentralAction_zero",
                        "mappedCentralAction_add",
                        "mappedCentralAction_neg_comp",
                        "mappedCentralAction_eq_id_iff_component_phases_zero",
                        "mappedCentralAction_eq_id_iff",
                        "mappedCentralAction_kernel_card",
                        "selectedMappedMatter_support_is_parity",
                        "selectedMappedMatter_nontrivial",
                        "mappedCentralAction_preserves_selected",
                        "kernel_on_selected_mapped_components",
                        "selectedMappedCentralAction_eq_id_iff",
                        "no_selected_central_parameter_realizes_parity_sign",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the exterior carrier, component-to-sector map, central-weight "
                "labels, and separate selection mask are supplied finite inputs. "
                "The character action is imposed through those weights rather "
                "than derived from ambient-projector conjugation, and projector "
                "ranks are not bound to exterior multiplicities. The universal "
                "minus-one result is only a nonconflation control and constructs "
                "no physical fermion parity. No source rule selects a physical "
                "matter action, and the package proves no continuum "
                "spin-statistics, particle spectrum, physical global form, or "
                "laboratory charge"
            ),
            "paper_ref": "zoo paper, finite exterior component bridge",
        },
        {
            "id": "coupling_universality",
            "statement": (
                "A5-invariant readouts have port-independent group-averaged cap "
                "sums, so the per-cap ratio of any two averaged readouts is "
                "universal with zero spread"
            ),
            "observed_counterpart": (
                "universality clause of the Einstein-branch coupling law"
            ),
            "match": "structural",
            "lean_declarations": {
                "A5CouplingSymmetry": [
                    "groupAverage_port_independent",
                    "coupling_ratio_universal",
                ],
                "A5PortAction": ["transitive_on_ports"],
                "PortFrameGram": ["degree_five", "gram_sq"],
            },
            "lean_receipts": _lean_receipt(
                "A5CouplingSymmetry",
                "A5PortAction",
                "PortFrameGram",
                declarations={
                    "A5CouplingSymmetry": (
                        "groupAverage_port_independent",
                        "coupling_ratio_universal",
                    ),
                    "A5PortAction": ("transitive_on_ports",),
                    "PortFrameGram": ("degree_five", "gram_sq"),
                },
            ),
            "hypothesis_boundary": (
                "reduces the universality clause to A5-equivariance of the "
                "implemented source law; no coupling value is implied"
            ),
            "paper_ref": "Standard Model gauge paper, coupling symmetry section",
        },
        {
            "id": "intrinsic_rank_three_response_completion",
            "statement": (
                "The declared twelve-port repair mean selects an intrinsic "
                "rank-three Gram quotient as its normalized infinite-response "
                "limit. The antipodal-odd integer load quotient is Z^6, the "
                "thirty-seam boundary image is the even-sum D6 sublattice, and "
                "both signed modules embed densely into the same abstract "
                "three-dimensional Euclidean completion. The sixty proper "
                "carrier maps act faithfully and isometrically on that "
                "completion"
            ),
            "observed_counterpart": (
                "a three-dimensional local spatial carrier with proper "
                "icosahedral frame changes"
            ),
            "match": (
                "exact intrinsic metric completion; physical position, scale, "
                "refinement, and gluing open"
            ),
            "lean_declarations": {
                "PortGramRepairBand": [
                    "portGram_unique_lowest_positive_galois_maximal",
                    "selected_family_band_is_port_gram",
                ],
                "PortGramRepairCovariance": [
                    "normalizedKernel_tendsto_portGram",
                    "portGram_antipodal_quotient",
                ],
                "PrimitivePortFrameQuotient": [
                    "frameQuotient_finrank",
                    "quotientEquivVec3_preserves_gram",
                    "pointEuclideanFrame_denseRange",
                ],
                "RepairWordCarrierReadout": [
                    "loadPosition_denseRange",
                    "universalPosition_isometry",
                ],
                "SeamCurrentCarrierQuotient": [
                    "exists_seamCurrent_iff_even",
                    "d6Position_denseRange",
                    "d6Position_isometry",
                ],
                "PortGramA5Isometry": [
                    "selected_band_action_faithful",
                    "carrierRotation_isometry",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PortGramRepairBand",
                "PortGramRepairCovariance",
                "PrimitivePortFrameQuotient",
                "RepairWordCarrierReadout",
                "SeamCurrentCarrierQuotient",
                "PortGramA5Isometry",
                declarations={
                    "PortGramRepairBand": (
                        "portGram_unique_lowest_positive_galois_maximal",
                        "selected_family_band_is_port_gram",
                    ),
                    "PortGramRepairCovariance": (
                        "normalizedKernel_tendsto_portGram",
                        "portGram_antipodal_quotient",
                    ),
                    "PrimitivePortFrameQuotient": (
                        "frameQuotient_finrank",
                        "quotientEquivVec3_preserves_gram",
                        "pointEuclideanFrame_denseRange",
                    ),
                    "RepairWordCarrierReadout": (
                        "loadPosition_denseRange",
                        "universalPosition_isometry",
                    ),
                    "SeamCurrentCarrierQuotient": (
                        "exists_seamCurrent_iff_even",
                        "d6Position_denseRange",
                        "d6Position_isometry",
                    ),
                    "PortGramA5Isometry": (
                        "selected_band_action_faithful",
                        "carrierRotation_isometry",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the finite carrier, scalar repair mean, complete centered probe "
                "census, cumulative signed-load readback, and response-Gram "
                "topology are declared inputs. The result constructs no "
                "pathwise physical position, operational scale, cofinal "
                "refinement, overlap gluing, global space, clock, or field"
            ),
            "paper_ref": (
                "spacetime-recovery paper, repair-response metric completion"
            ),
        },
        {
            "id": "carrier_class_dispersion_band",
            "statement": (
                "Every member of the declared positive-weight scalar cosine "
                "class, whose full spatial symbol is the orbit sum, has C4 < 0 "
                "and B0/C4^2 at least 10/21, "
                "with equality exactly on one-radius support. Its anisotropic "
                "ranks one through five vanish and B6/B0 lies in "
                "[-16/135, 16/75] on the unique rotated I6 line. At eighth "
                "order no new angular shape appears, and every one-radius "
                "member obeys 5 D6 B0 = 12 B6 D0, equivalently "
                "D6/D0 = (12/5)(B6/B0). The polynomial form includes the "
                "exact zero-anisotropy mixture; multi-radius members retain "
                "radial-moment dependence"
            ),
            "observed_counterpart": (
                "a linked isotropic and rank-six vacuum-dispersion surface not "
                "fixed by the Standard Model with General Relativity"
            ),
            "match": (
                "exact class theorem; physical sector, frame, finite scale, "
                "readout, and comparison remain open"
            ),
            "artifact_ref": _rel("carrier_class_dispersion"),
            "lean_declarations": {
                "A5CarrierClassBand": [
                    "band_endpoints",
                    "tuned_zero",
                    "gap_zero_iff_single_radius",
                    "general_member_in_band",
                    "cross_order_lock",
                    "cross_order_polynomial",
                    "multi_radius_negative_control",
                ]
            },
            "lean_receipts": _lean_receipt(
                "A5CarrierClassBand",
                declarations={
                    "A5CarrierClassBand": (
                        "band_endpoints",
                        "tuned_zero",
                        "gap_zero_iff_single_radius",
                        "general_member_in_band",
                        "cross_order_lock",
                        "cross_order_polynomial",
                        "multi_radius_negative_control",
                    )
                },
            ),
            "hypothesis_boundary": (
                "the theorem applies to positive-weight finite mixtures of "
                "proper-carrier direction orbits with the declared cosine hop "
                "symbol, quadratic normalization, and no independent isotropic "
                "counterterm. No theorem identifies this class with a physical "
                "field or fixes its scale, frame, detector response, or "
                "exclusivity"
            ),
            "paper_ref": "flagship paper, carrier-class dispersion theorem",
        },
        {
            "id": "positive_cosine_frequency_contraction",
            "statement": (
                "Every normalized complete positive tight-frame cosine symbol "
                "has an exact sine-feature realization. The feature map is a "
                "Euclidean contraction, so its nonnegative auxiliary frequency "
                "is globally 1-Lipschitz at all momenta. Exact support bindings "
                "give t=4 and prefactor 1/(2a^2) for the FZ-11 vertex support, "
                "and t=10 and prefactor 1/(5a^2) for the FZ-12 edge support"
            ),
            "observed_counterpart": (
                "an all-momentum upper bound on an auxiliary carrier dispersion"
            ),
            "match": (
                "exact bounded finite theorem; physical position, frequency, "
                "clock, field, signal front, frame, scale, readout, and "
                "comparison remain open"
            ),
            "artifact_ref": _rel("carrier_frequency_speed"),
            "receipt_sha256": carrier_frequency["receipt_sha256"],
            "lean_declarations": {
                "CarrierFrequencySpeed": carrier_frequency["lean"]["theorems"]
            },
            "lean_receipts": _lean_receipt(
                "CarrierFrequencySpeed",
                declarations={
                    "CarrierFrequencySpeed": tuple(
                        carrier_frequency["lean"]["theorems"]
                    )
                },
            ),
            "hypothesis_boundary": (
                "the unit constant is a certified upper bound for the auxiliary "
                "norm in the selected Euclidean carrier chart, not an "
                "optimality theorem or a physical signal-speed claim. The "
                "receipt reads no comparison data and changes no frozen bytes"
            ),
            "paper_ref": (
                "screen-microphysics paper, positive-cosine frequency contraction"
            ),
        },
        {
            "id": "time_order_type_ledger",
            "statement": (
                "Universe closure, repair execution order, observer record "
                "order, modular parameter, worldline realization, clock "
                "readout, proper time, and optional global time are distinct "
                "formal types. In the committed source environment, a "
                "canonical witness matrix rejects all 56 ordered transitive "
                "coercions between distinct layers, and explicit named maps are required. "
                "Positive affine clock regraduation preserves strict "
                "record monotonicity, and an inhabited record with nonzero "
                "offset proves clock-origin nonuniqueness"
            ),
            "observed_counterpart": (
                "typed separation between operational ordering and physical time"
            ),
            "match": "exact formal boundary; physical time realization open",
            "lean_declarations": {
                "TimeOrderLedger": [
                    "canonicalLedgerKinds_pairwise",
                    "offsetGauge_ne",
                ],
            },
            "lean_receipts": _lean_receipt(
                "TimeOrderLedger",
                declarations={
                    "TimeOrderLedger": (
                        "canonicalLedgerKinds_pairwise",
                        "offsetGauge_ne",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the ledger constructs no public-world endpoint, worldline, "
                "physical clock, proper-time calibration, global time "
                "function, or modular-to-time identity and emits no prediction"
            ),
            "paper_ref": "observers paper, time and order interface",
        },
        {
            "id": "bounded_observer_time_calibration",
            "statement": (
                "A record history with a supplied strictly increasing natural-number "
                "rank gives an order-compatible scalar readout, which admits every "
                "strictly increasing regrading; an exact three-tick cubic control is "
                "not affine. "
                "After a unit-timelike affine event law is supplied, every "
                "precedence is future timelike across overlapping charts and "
                "along that same supplied history the positive clock increment "
                "is additive with square equal to the invariant Lorentz "
                "quadratic interval. One shared event leaves an affine "
                "comparison nonunique; two ordered event pairs determine the "
                "unique positive-affine interpolation of their four supplied "
                "readings. At a third shared event, affine consistency is "
                "equivalent to a cross-product equation and gives a "
                "nondegenerate no-new-fit-parameter check when its event and "
                "both readings differ from the anchors. A held-out reading "
                "requires separate predesignation and custody"
            ),
            "observed_counterpart": (
                "record-order data and conditional operational clock comparison"
            ),
            "match": (
                "exact bounded conditional algebra with finite controls; source "
                "physical clock open"
            ),
            "lean_declarations": {
                "ObserverHistory": [
                    "threeRecord_control",
                    "discreteConstantClock_not_injective",
                ],
                "ClockReadout": [
                    "cubicThreeTickClock_not_affine",
                    "throughTwoPoints_unique",
                ],
                "WorldlineRealization": [
                    "displacement_futureTimelike_in_chart",
                    "threeRecord_twoChart_futureTimelike",
                ],
                "ProperTimeCalibration": [
                    "properTimeBetween_sq_eq_interval_in_chart",
                    "properTimeBetween_add",
                ],
                "ClockComparison": [
                    "onePoint_not_unique",
                    "calibration_unique",
                    "affineConsistent_iff_crossMultiplication",
                ],
            },
            "lean_receipts": _lean_receipt(
                "ObserverHistory",
                "ClockReadout",
                "WorldlineRealization",
                "ProperTimeCalibration",
                "ClockComparison",
                declarations={
                    "ObserverHistory": (
                        "threeRecord_control",
                        "discreteConstantClock_not_injective",
                    ),
                    "ClockReadout": (
                        "cubicThreeTickClock_not_affine",
                        "throughTwoPoints_unique",
                    ),
                    "WorldlineRealization": (
                        "displacement_futureTimelike_in_chart",
                        "threeRecord_twoChart_futureTimelike",
                    ),
                    "ProperTimeCalibration": (
                        "properTimeBetween_sq_eq_interval_in_chart",
                        "properTimeBetween_add",
                    ),
                    "ClockComparison": (
                        "onePoint_not_unique",
                        "calibration_unique",
                        "affineConsistent_iff_crossMultiplication",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the event atlas, visibility, event map, clock, future-unit "
                "direction, affine unit-speed law, and shared-event equalities "
                "are supplied. No source history, refinement transport, "
                "physical instrument, SI unit, global time, modular-time "
                "identity, observable, decision rule, or prediction follows"
            ),
            "paper_ref": "observers paper, bounded observer-time calibration",
        },
        {
            "id": "finite_public_record_algebra_and_sharp_no_cloning",
            "statement": (
                "The span of a finite projective partition is a commutative "
                "matrix star subalgebra contained in its commutant and is "
                "star-algebra equivalent to complex functions on the nonzero "
                "projector labels. A common linear isometry can sharply copy "
                "two distinct states from one normalized blank only when "
                "they are orthogonal"
            ),
            "observed_counterpart": (
                "classical public records and the sharp-state copying boundary"
            ),
            "match": "exact finite theorem package",
            "lean_declarations": {
                "PublicRecordAlgebra": [
                    "publicSubalgebra_mul_comm",
                    "publicSubalgebra_le_commutant",
                    "recordSynthesisStarAlgHom_bijective",
                    "publicRecordFunctionEquiv_apply",
                ],
                "NoBroadcastingAdapter": [
                    "SharpCloneWitness.overlap_zero_or_one",
                    "SharpCloneWitness.eq_of_overlap_one",
                    "SharpCloneWitness.orthogonal_of_ne",
                    "NoBroadcastingAdapter.objective_pair_compatible",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PublicRecordAlgebra",
                "NoBroadcastingAdapter",
                declarations={
                    "PublicRecordAlgebra": (
                        "publicSubalgebra_mul_comm",
                        "publicSubalgebra_le_commutant",
                        "recordSynthesisStarAlgHom_bijective",
                        "publicRecordFunctionEquiv_apply",
                    ),
                    "NoBroadcastingAdapter": (
                        "SharpCloneWitness.overlap_zero_or_one",
                        "SharpCloneWitness.eq_of_overlap_one",
                        "SharpCloneWitness.orthogonal_of_ne",
                        "NoBroadcastingAdapter.objective_pair_compatible",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "zero projectors are removed before coordinate equivalence; "
                "the mixed-state no-broadcasting implication remains an "
                "explicit adapter premise, and no physical record-selection "
                "or measurement theorem follows"
            ),
            "paper_ref": "observers paper, public event algebra",
        },
        {
            "id": "finite_consensus_tower_interface",
            "statement": (
                "One typed directed refinement object carries finite observer "
                "and record fibres, observer record orders, private matrix "
                "algebras, commutative public star subalgebras, record "
                "representatives, certified states, and linear generators. "
                "Its refinement laws preserve every layer, with states "
                "restricting contravariantly by exact trace pairing. A "
                "constant adaptor reuses an existing projective partition "
                "and density state with discrete order and zero generator"
            ),
            "observed_counterpart": (
                "one common finite refinement substrate for observer theories"
            ),
            "match": "exact structural interface; source realization open",
            "lean_declarations": {
                "ConsensusTower": [
                    "public_mem_refine",
                    "refine_recordElement",
                    "refine_precedes",
                    "refine_generator",
                    "refine_state_pairing",
                    "constantConsensusTower_public",
                    "constantConsensusTower_recordElement",
                ],
            },
            "lean_receipts": _lean_receipt(
                "ConsensusTower",
                declarations={
                    "ConsensusTower": (
                        "public_mem_refine",
                        "refine_recordElement",
                        "refine_precedes",
                        "refine_generator",
                        "refine_state_pairing",
                        "constantConsensusTower_public",
                        "constantConsensusTower_recordElement",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the constant adaptor proves packaging only. No nonconstant "
                "source tower, repair endpoint, causal net, geometry, clock, "
                "continuum limit, physical evolution, or prediction follows"
            ),
            "paper_ref": "observers paper, consensus-tower root interface",
        },
        {
            "id": "finite_public_world_endpoint",
            "statement": (
                "A finite inhabited raw presentation has a literal kernel "
                "quotient whose equality is exactly public-readback equality "
                "and whose points are the realized signatures. Termination "
                "constructs a finite completed schedule; confluence, semantic "
                "fixed-point completeness, and explicit repair-output plus "
                "enabledness congruence make the consistent public endpoint "
                "independent of completed schedule and raw representative. "
                "The existing OPH Repair descends to an idempotent public map, "
                "and typed OPH and A3-regulator adaptors retain every premise"
            ),
            "observed_counterpart": (
                "an observer-independent public normal-form endpoint"
            ),
            "match": (
                "exact bounded conditional endpoint; source and limit open"
            ),
            "lean_declarations": {
                "PublicWorldQuotient": [
                    "toPublicWorld_eq_iff",
                    "publicSignature_injective",
                    "hiddenBit_distinct_but_publicly_equal",
                ],
                "FixedPointEndpoint": [
                    "public_endpoint_exists_unique_on_public_class",
                    "publicRepair_idempotent",
                    "lr_public_endpoint_exists_unique_on_gauge_class",
                    "representative_no_descended_repair",
                    "primitiveLR_endpoint_exists_unique_on_gauge_class",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PublicWorldQuotient",
                "FixedPointEndpoint",
                declarations={
                    "PublicWorldQuotient": (
                        "toPublicWorld_eq_iff",
                        "publicSignature_injective",
                        "hiddenBit_distinct_but_publicly_equal",
                    ),
                    "FixedPointEndpoint": (
                        "public_endpoint_exists_unique_on_public_class",
                        "publicRepair_idempotent",
                        "lr_public_endpoint_exists_unique_on_gauge_class",
                        "representative_no_descended_repair",
                        "primitiveLR_endpoint_exists_unique_on_gauge_class",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "CompletedSchedule is terminal finite completion, not infinite "
                "scheduler fairness. The A3 adaptor requires a caller-supplied "
                "seed and injective readback encoding. No source-selected "
                "physical world, cross-regulator naturality, continuum limit, "
                "clock, observable, decision rule, or prediction follows"
            ),
            "paper_ref": "observers paper, finite public-world endpoint",
        },
        {
            "id": "canonical_intrinsic_lorentz_module",
            "statement": (
                "The real Pauli-coordinate module Herm2 is exactly the "
                "four-dimensional space of two-by-two complex Hermitian "
                "matrices, with determinant equal to the Lorentz quadratic "
                "form of constructive inertia (+---). Positive future-null "
                "rays are set-equivalent to the unit two-sphere, the algebraic "
                "future-unit hyperboloid has three-dimensional positive rest "
                "spaces, and an explicit linear chart matches the existing "
                "Einstein coordinates with exactly the required sign flip"
            ),
            "observed_counterpart": (
                "four-dimensional Lorentzian event and observer-frame geometry"
            ),
            "match": (
                "exact intrinsic geometry and coordinate bridge; physical "
                "soldering remains outside C1"
            ),
            "lean_declarations": {
                "CanonicalLorentzModule": [
                    "det_toMatrix",
                    "isHermitian_iff_existsUnique_toMatrix",
                    "finrank_Herm2",
                    "time_axis_positive",
                    "spatial_axis_negative",
                ],
                "CelestialNullCone": [
                    "rayToCelestial_celestialToRay",
                    "celestialToRay_rayToCelestial",
                ],
                "ObserverFrameHyperboloid": [
                    "frame_time_sq_eq_one_add_spatial",
                    "frame_time_ge_one",
                ],
                "ObserverRestSpace": [
                    "finrank_restSpace",
                    "restMetric_pos",
                    "restMetric_self_eq_zero_iff",
                ],
                "EinsteinTensorBridge": [
                    "lorentzQ_eq_neg_einsteinQuad",
                    "lorentzQ_eq_zero_iff_einsteinQuad_eq_zero",
                    "isFutureNull_iff_einstein",
                ],
            },
            "lean_receipts": _lean_receipt(
                "CanonicalLorentzModule",
                "CelestialNullCone",
                "ObserverFrameHyperboloid",
                "ObserverRestSpace",
                "EinsteinTensorBridge",
                declarations={
                    "CanonicalLorentzModule": (
                        "det_toMatrix",
                        "isHermitian_iff_existsUnique_toMatrix",
                        "finrank_Herm2",
                        "time_axis_positive",
                        "spatial_axis_negative",
                    ),
                    "CelestialNullCone": (
                        "rayToCelestial_celestialToRay",
                        "celestialToRay_rayToCelestial",
                    ),
                    "ObserverFrameHyperboloid": (
                        "frame_time_sq_eq_one_add_spatial",
                        "frame_time_ge_one",
                    ),
                    "ObserverRestSpace": (
                        "finrank_restSpace",
                        "restMetric_pos",
                        "restMetric_self_eq_zero_iff",
                    ),
                    "EinsteinTensorBridge": (
                        "lorentzQ_eq_neg_einsteinQuad",
                        "lorentzQ_eq_zero_iff_einsteinQuad_eq_zero",
                        "isFutureNull_iff_einstein",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the celestial equivalence is set-level and the frames and "
                "rest spaces are algebraic. C1 alone supplies no observer-patch "
                "selection or soldering. The separate bounded C2 contract "
                "supplies algebraic overlap covariance only from declared "
                "inputs; source event realization, rods, clocks, physical "
                "spacetime, continuum attachment, observable, decision rule, "
                "and prediction remain open"
            ),
            "paper_ref": "observers paper, canonical Lorentz module",
        },
        {
            "id": "algebraic_event_frame_soldering",
            "statement": (
                "Coincidence-invariant Herm2 readback descends uniquely through "
                "an actual event setoid. Separately, one supplied time-oriented "
                "affine Lorentz overlap cocycle and a chart-coordinate family "
                "satisfying its overlap law induce translation-free "
                "displacement covariance, invariant intervals, future-null "
                "celestial transport, compatible event frames, and isometric "
                "transport of their positive rest spaces. The rank-three source "
                "FrameQuotient is linearly and isometrically identified with "
                "the standard internal rest fiber as a candidate local readback. "
                "The formal stack does not identify the quotient-descended "
                "readback with the supplied atlas coordinates. "
                "Exact controls exhibit a nontrivial translated future-null "
                "atlas and show that reflexive symmetric pairwise overlap need "
                "not be transitive"
            ),
            "observed_counterpart": (
                "event-frame and local rest-space Lorentz covariance"
            ),
            "match": (
                "exact bounded algebraic contract; source and physical receipts open"
            ),
            "lean_declarations": {
                "LorentzOverlapCocycle": [
                    "LorentzOverlapCocycle.act_cocycle",
                    "LorentzOverlapCocycle.act_reverse_left",
                    "LorentzOverlapCocycle.act_reverse_right",
                    "LorentzOverlapCocycle.act_sub_act",
                ],
                "EventGermDisplacement": [
                    "coincidenceInvariant_iff_existsUnique_descendedReadback",
                    "overlapControl_not_transitive",
                    "EventGermAtlas.displacement_reverse",
                    "EventGermAtlas.displacement_chain",
                    "EventGermAtlas.displacement_overlap",
                    "EventGermAtlas.interval_overlap",
                ],
                "CelestialSoldering": [
                    "OrientedLorentzEquiv.mapFutureNullRay_trans",
                    "OrientedLorentzEquiv.celestialAction_trans",
                    "EventGermAtlas.futureNullDisplacementRay_overlap",
                    "EventGermAtlas.celestialSolder_overlap",
                ],
                "EventFrameSoldering": [
                    "EventFrameSoldering.algebraicConsequences",
                    "EventFrameSoldering.expanded_handoff_iff_residual",
                    "control_chart_translation_nonzero",
                    "control_displacement_nonzero",
                    "control_displacement_futureNull",
                ],
                "SpatialReadbackSoldering": [
                    "restProjection_decomposition",
                    "OrientedLorentzEquiv.restProjection_covariant",
                    "OrientedLorentzEquiv.restEquiv_preserves_metric",
                    "frameQuotientEquivStandardRest_preserves_metric",
                    "EventFrameSoldering.displacement_time_add_spatial",
                    "EventFrameSoldering.spatialReadback_overlap",
                ],
            },
            "lean_receipts": _lean_receipt(
                "LorentzOverlapCocycle",
                "EventGermDisplacement",
                "CelestialSoldering",
                "EventFrameSoldering",
                "SpatialReadbackSoldering",
                declarations={
                    "LorentzOverlapCocycle": (
                        "LorentzOverlapCocycle.act_cocycle",
                        "LorentzOverlapCocycle.act_reverse_left",
                        "LorentzOverlapCocycle.act_reverse_right",
                        "LorentzOverlapCocycle.act_sub_act",
                    ),
                    "EventGermDisplacement": (
                        "coincidenceInvariant_iff_existsUnique_descendedReadback",
                        "overlapControl_not_transitive",
                        "EventGermAtlas.displacement_reverse",
                        "EventGermAtlas.displacement_chain",
                        "EventGermAtlas.displacement_overlap",
                        "EventGermAtlas.interval_overlap",
                    ),
                    "CelestialSoldering": (
                        "OrientedLorentzEquiv.mapFutureNullRay_trans",
                        "OrientedLorentzEquiv.celestialAction_trans",
                        "EventGermAtlas.futureNullDisplacementRay_overlap",
                        "EventGermAtlas.celestialSolder_overlap",
                    ),
                    "EventFrameSoldering": (
                        "EventFrameSoldering.algebraicConsequences",
                        "EventFrameSoldering.expanded_handoff_iff_residual",
                        "control_chart_translation_nonzero",
                        "control_displacement_nonzero",
                        "control_displacement_futureNull",
                    ),
                    "SpatialReadbackSoldering": (
                        "restProjection_decomposition",
                        "OrientedLorentzEquiv.restProjection_covariant",
                        "OrientedLorentzEquiv.restEquiv_preserves_metric",
                        "frameQuotientEquivStandardRest_preserves_metric",
                        "EventFrameSoldering.displacement_time_add_spatial",
                        "EventFrameSoldering.spatialReadback_overlap",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the coincidence setoid, invariant readback, affine Lorentz "
                "cocycle, overlap-compatible chart-coordinate family, and base "
                "frame are supplied. Quotient descent does not construct or "
                "identify that atlas family. The exact handoff retains "
                "source atlas realization, event population, certified "
                "separation, open rank-four charts, physical cone attachment, "
                "refinement naturality, semantic causal reachability, and an "
                "operational clock. Issues #693, #694, and #703 own those "
                "residuals; no physical spacetime, Einstein dynamics, "
                "observable, decision rule, or prediction follows"
            ),
            "paper_ref": "spacetime paper, bounded algebraic event-frame soldering",
        },
        {
            "id": "finite_causal_observer_net_interface",
            "statement": (
                "A proof-carrying finite enrichment of one consensus tower "
                "has declared region posets, overlaps, disjointness, local star "
                "subalgebras, isotony, algebraic locality, covariant refinement, "
                "compatible expectations, and idempotent regional repairs that "
                "explicitly fix local and declared-disjoint observables. Its "
                "relaxations compose exactly and preserve remote expectations. "
                "Supplied overlap restrictions support a unique restriction-gluing "
                "interface on declared nonempty subregion families, controls "
                "isolate missing premises, and a partition-and-state-parameterized "
                "commutative model proves conditional consistency. A retained "
                "source packet supplies four disjoint windows with split-fibre "
                "labels; a separately declared adapter constructs noncommutative "
                "regional blocks, exact coverage and gluing, and nonunital "
                "two-by-two matrix corners at every window"
            ),
            "observed_counterpart": (
                "a causal local quantum-observable net with overlap descent"
            ),
            "match": (
                "substantial conditional finite interface and declared-adapter "
                "coverage; source-attached operators and product split open"
            ),
            "lean_declarations": {
                "FiniteCausalObserverNet": [
                    "commute_of_disjoint",
                    "regionalExpectation_refine",
                    "relaxedRepair_compose",
                    "relaxedRepair_fixes_disjoint",
                    "relaxedRepair_remote_expectation",
                    "kraus_remote_marginal_invariant",
                    "fullM2_distinct_regions_not_local",
                    "idempotence_does_not_force_remote_fix",
                    "partitionPublicCausalNet_has_disjoint_pair",
                ],
                "ObserverNetDescent": [
                    "jointly_injective_of_unique_descent",
                    "no_descent_of_indistinguishable_global_sections",
                    "glue_restrict",
                    "glue_unique",
                    "partitionPublicTwoRegionCover_hasUniqueDescent",
                ],
                "RichFibreWitness": [
                    "richSupport_pairwise_disjoint",
                    "richSplit_census",
                ],
                "RichFibreRegionalNet": [
                    "richRegional_noncommutative_all",
                    "richDesignatedFactor",
                    "richWindowCover_coverageLaw",
                    "richDropCover_not_coverageLaw",
                    "richWindowCover_reconstruction",
                ],
            },
            "lean_receipts": _lean_receipt(
                "FiniteCausalObserverNet",
                "ObserverNetDescent",
                "RichFibreWitness",
                "RichFibreRegionalNet",
                declarations={
                    "FiniteCausalObserverNet": (
                        "commute_of_disjoint",
                        "regionalExpectation_refine",
                        "relaxedRepair_compose",
                        "relaxedRepair_fixes_disjoint",
                        "relaxedRepair_remote_expectation",
                        "kraus_remote_marginal_invariant",
                        "fullM2_distinct_regions_not_local",
                        "idempotence_does_not_force_remote_fix",
                        "partitionPublicCausalNet_has_disjoint_pair",
                    ),
                    "ObserverNetDescent": (
                        "jointly_injective_of_unique_descent",
                        "no_descent_of_indistinguishable_global_sections",
                        "glue_restrict",
                        "glue_unique",
                        "partitionPublicTwoRegionCover_hasUniqueDescent",
                    ),
                    "RichFibreWitness": (
                        "richSupport_pairwise_disjoint",
                        "richSplit_census",
                    ),
                    "RichFibreRegionalNet": (
                        "richRegional_noncommutative_all",
                        "richDesignatedFactor",
                        "richWindowCover_coverageLaw",
                        "richDropCover_not_coverageLaw",
                        "richWindowCover_reconstruction",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "ordinary isotony does not supply the star-homomorphic "
                "restriction retractions or unique gluing. `FiniteCover` has no "
                "joint-coverage axiom and the B4 helper has no regional-factor "
                "attachment. The rich-fibre source payload supplies window/class "
                "labels only: its block algebra, base-point state, restrictions, "
                "and repair are declared postprocessors. Conditional coverage is "
                "attained inside that adapter, but its nonunital matrix corners are "
                "not tensor factors or a `TensorSplitReceipt`. Issue #692 gates "
                "source-attached operator generation and a justified region-product "
                "or local-channel adapter; "
                "no CP/CPTP channel, scheduler locality, spacetime causality, "
                "time-slice property, continuum QFT, observable, decision rule, "
                "or prediction is supplied"
            ),
            "paper_ref": "observers paper, finite causal observer-net interface",
        },
        {
            "id": "finite_publicization_dynamics",
            "statement": (
                "An idempotent linear publicization map has an exactly "
                "solvable public-residual relaxation with multiplicative "
                "composition and exponential semigroup law. Its Poissonized "
                "closed form satisfies the initial-value and generator-flow "
                "identities and equals the literal Banach-algebra operator "
                "exponential for bounded idempotent endomorphisms. Partition "
                "averaging has an explicit normalized "
                "Kraus family and is formally CPTP. Partition pinching is "
                "also CPTP, and its generator equals the displayed projector-rate matrix "
                "dissipator with fixed algebra equal to the commutant at "
                "nonzero rate"
            ),
            "observed_counterpart": (
                "finite public/private relaxation and stable pointer algebra"
            ),
            "match": "exact finite linear and matrix identities",
            "lean_declarations": {
                "PartitionAverageCP": [
                    "ProjectivePartition.partitionAverageKraus_complete",
                    "partitionAverage_kraus_form",
                ],
                "TwoScalePublicRepair": [
                    "publicRelax_compose",
                    "publicRelaxTime_add",
                    "publicRelaxTime_residual",
                ],
                "PoissonizedRepair": [
                    "poissonizedRepair_add",
                    "repairGenerator_eq_zero_iff",
                    "hasDerivAt_poissonizedRepair_eq_generator",
                ],
                "PoissonizedRepairOperatorExp": [
                    "normedSpace_exp_smul_idempotent",
                    "normedSpace_exp_continuousRepairGenerator",
                    "normedSpace_exp_continuousRepairGenerator_apply",
                ],
                "ConditionalExpectationGenerator": [
                    "conditionalExpectationGenerator_eq_projectorGKSL",
                    "conditionalExpectationGenerator_eq_zero_iff_mem_commutant",
                    "multiCollarGenerator_eq_zero_iff_stableIntersection",
                ],
                "ChoiCPTP": [
                    "partitionAverage_isCPTP",
                    "partitionPinching_isCPTP",
                    "relaxationChannel_isCPTP",
                    "transposeMap_positive_tracePreserving_not_CP",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PartitionAverageCP",
                "TwoScalePublicRepair",
                "PoissonizedRepair",
                "PoissonizedRepairOperatorExp",
                "ConditionalExpectationGenerator",
                "ChoiCPTP",
                declarations={
                    "PartitionAverageCP": (
                        "ProjectivePartition.partitionAverageKraus_complete",
                        "partitionAverage_kraus_form",
                    ),
                    "TwoScalePublicRepair": (
                        "publicRelax_compose",
                        "publicRelaxTime_add",
                        "publicRelaxTime_residual",
                    ),
                    "PoissonizedRepair": (
                        "poissonizedRepair_add",
                        "repairGenerator_eq_zero_iff",
                        "hasDerivAt_poissonizedRepair_eq_generator",
                    ),
                    "PoissonizedRepairOperatorExp": (
                        "normedSpace_exp_smul_idempotent",
                        "normedSpace_exp_continuousRepairGenerator",
                        "normedSpace_exp_continuousRepairGenerator_apply",
                    ),
                    "ConditionalExpectationGenerator": (
                        "conditionalExpectationGenerator_eq_projectorGKSL",
                        "conditionalExpectationGenerator_eq_zero_iff_mem_commutant",
                        "multiCollarGenerator_eq_zero_iff_stableIntersection",
                    ),
                    "ChoiCPTP": (
                        "partitionAverage_isCPTP",
                        "partitionPinching_isCPTP",
                        "relaxationChannel_isCPTP",
                        "transposeMap_positive_tracePreserving_not_CP",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the formal CP/CPTP predicate covers partition averaging, "
                "partition pinching, and the nonnegative-time relaxation "
                "channel; a positive trace-preserving transpose control is "
                "not completely positive. The operator-exponential theorem "
                "uses bounded endomorphisms of a complete real normed space. Poisson rate and "
                "forward-time interpretations require nonnegative parameters. "
                "No source-derived rate, physical clock, or prediction is supplied"
            ),
            "paper_ref": "observers paper, finite publicization dynamics",
        },
        {
            "id": "finite_public_private_dynamics",
            "statement": (
                "Positive unital complex-linear maps of the finite active-"
                "record function algebra are exactly real row-stochastic "
                "kernels under the declared coordinatewise cone. Every public "
                "star automorphism is uniquely pullback by a label permutation, "
                "so every pointwise-continuous real-parameter group of arbitrary "
                "public star automorphisms is trivial. Every star automorphism "
                "of one finite full private "
                "endomorphism block is unitarily inner, and a supplied "
                "self-adjoint Hamiltonian generates a unitary real-parameter "
                "von Neumann flow"
            ),
            "observed_counterpart": (
                "classical-stochastic public and quantum-unitary private dynamics"
            ),
            "match": "substantial exact finite packet; global converse open",
            "lean_declarations": {
                "PublicMarkov": [
                    "recordMapOfKernel_injective",
                    "positive_unital_iff_stochastic",
                    "activeRecord_positive_unital_iff_stochastic",
                    "toPerm_eq_refl",
                    "function_action_eq",
                ],
                "PublicAutomorphism": [
                    "publicStarAutomorphism_is_labelPermutation",
                    "publicStarAutomorphism_labelPermutation_unique",
                    "toAut_eq_refl",
                ],
                "PrivateInner": [
                    "finitePrivateStarAutomorphism_inner",
                    "hasDerivAt_realVonNeumannFlow",
                    "hamiltonianPropagator_mem_unitary",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PublicMarkov",
                "PublicAutomorphism",
                "PrivateInner",
                declarations={
                    "PublicMarkov": (
                        "recordMapOfKernel_injective",
                        "positive_unital_iff_stochastic",
                        "activeRecord_positive_unital_iff_stochastic",
                        "toPerm_eq_refl",
                        "function_action_eq",
                    ),
                    "PublicAutomorphism": (
                        "publicStarAutomorphism_is_labelPermutation",
                        "publicStarAutomorphism_labelPermutation_unique",
                        "toAut_eq_refl",
                    ),
                    "PrivateInner": (
                        "finitePrivateStarAutomorphism_inner",
                        "hasDerivAt_realVonNeumannFlow",
                        "hamiltonianPropagator_mem_unitary",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "public star automorphisms are classified exactly as unique "
                "label permutations, so arbitrary pointwise-continuous public "
                "automorphism groups are trivial. Private innerness covers "
                "only one full matrix block. Arbitrary finite central sums, a "
                "coherent continuous unitary lift and single generator for every "
                "private automorphism group, source dynamics, physical time, and "
                "predictions are outside the attained packet"
            ),
            "paper_ref": "observers paper, finite public/private dynamics",
        },
        {
            "id": "finite_holonomy_character_phase",
            "statement": (
                "For reversal-compatible group labels on endpoint-typed finite "
                "paths, the ordered transport ratio of two common-endpoint "
                "paths equals the holonomy of their closed ratio loop, and "
                "every unitary character maps it to the exact relative phase. "
                "Recharting conjugates based holonomy and leaves character "
                "phase invariant. An explicit four-vertex path/face control "
                "has two flat declared triangular faces and a separate "
                "undeclared, unfilled loop with nontrivial holonomy and "
                "two-arm phase. A supplied cyclic ZMod n sector "
                "forces the character phase to be an nth root of unity"
            ),
            "observed_counterpart": (
                "finite algebraic two-path character phase and cyclic-root structure"
            ),
            "match": (
                "exact bounded algebraic packet; physical attachment open"
            ),
            "lean_declarations": {
                "HolonomyInterference": [
                    "transportRatio_eq_closedLoopHolonomy",
                    "relativeCharacterPhase_eq_closedLoopPhase",
                    "holonomy_rechart_conjugate",
                    "holonomy_rechart_invariant",
                    "characterPhase_rechart_invariant",
                    "exists_localTriangleFlat_globalHolonomy_nontrivial",
                    "long_reference_relativeCharacterPhase",
                    "exists_localTriangleFlat_relativeCharacterPhase_nontrivial",
                    "zmodCharacter_phase_pow_order",
                    "cyclicSectorCharacter_phase_pow_order",
                    "loopCharacterPhase_pow_order_of_cyclicHolonomy",
                    "cyclicLoopPhase_pow_order",
                ],
            },
            "lean_receipts": _lean_receipt(
                "HolonomyInterference",
                declarations={
                    "HolonomyInterference": (
                        "transportRatio_eq_closedLoopHolonomy",
                        "relativeCharacterPhase_eq_closedLoopPhase",
                        "holonomy_rechart_conjugate",
                        "holonomy_rechart_invariant",
                        "characterPhase_rechart_invariant",
                        "exists_localTriangleFlat_globalHolonomy_nontrivial",
                        "long_reference_relativeCharacterPhase",
                        "exists_localTriangleFlat_relativeCharacterPhase_nontrivial",
                        "zmodCharacter_phase_pow_order",
                        "cyclicSectorCharacter_phase_pow_order",
                        "loopCharacterPhase_pow_order_of_cyclicHolonomy",
                        "cyclicLoopPhase_pow_order",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the edge labels, reversal law, declared four-vertex path/face "
                "control, unitary "
                "character, cyclic sector map, and factorization are supplied. "
                "Local flatness quantifies only over the two declared faces of "
                "the four-vertex control; no puncture or noncontractibility is "
                "proved. The Aharonov--Bohm terminology is only an algebraic "
                "two-path analogy. No observer-source connection, physical "
                "gauge field, spacetime loop, charge, clock, flux, detector, "
                "laboratory fringe, or prediction is constructed"
            ),
            "paper_ref": "gauge paper, finite holonomy/interference packet",
        },
        {
            "id": "fixed_regulator_operational_overlap_evidence",
            "statement": (
                "At one finite regulator, two distinct operational observers "
                "can be bound through the E6 access cut so that each owner "
                "algebra is the declared accessible algebra, every committed "
                "record and own readout agrees after typed restriction to a "
                "proper meet, and one common restriction is nonzero and "
                "accessible to both. In the exact witness, owner-region "
                "records differ before restriction and agree only at the "
                "shared corner; a one-corner mutation falsifies the receipt"
            ),
            "observed_counterpart": (
                "the seven-clause bounded operational-observer and shared-record "
                "interface"
            ),
            "match": "exact fixed-regulator witness and negative control",
            "lean_declarations": {
                "OperationalOverlapEvidence": [
                    "operationalObservers_share_visible_event_record",
                    "operationalOverlap_good_evidence",
                    "operationalOverlap_good_common_mem_both",
                    "operationalOverlap_bad_not_evidence",
                ],
            },
            "lean_receipts": _lean_receipt(
                "OperationalOverlapEvidence",
                declarations={
                    "OperationalOverlapEvidence": (
                        "operationalObservers_share_visible_event_record",
                        "operationalOverlap_good_evidence",
                        "operationalOverlap_good_common_mem_both",
                        "operationalOverlap_bad_not_evidence",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the access cut, restrictions, observer values, and record "
                "packet are supplied finite data. The result is not a "
                "source-production theorem and proves no cross-regulator "
                "naturality, higher-overlap cocycle, continuum observer, "
                "laboratory attachment, or prediction"
            ),
            "paper_ref": "consensus paper, operational observer receipt",
        },
        {
            "id": "finite_born_frame_rank_gap",
            "statement": (
                "The twelve declared central port atoms form one classical "
                "context with an eleven-dimensional normalized weight simplex. "
                "The separate declared qubit adapter has six disjoint antipodal "
                "binary contexts: additive weights have affine dimension six, "
                "while the trace-one Hermitian/Born slice has dimension three "
                "and is characterized by three exact golden-ratio relations. "
                "Tomography is unique when a representation exists, but exact "
                "unit-interval controls show both nonrepresentation and a unique "
                "Hermitian representation whose matrix is nonpositive. On the "
                "full celestial sphere, the continuous normalized binary weight "
                "F(n)=(1+n_z^3)/2 is exactly non-affine; after affinity is "
                "supplied, dense probability tests force the coefficient into "
                "the closed unit ball. The current finite unsharp battery also "
                "fails: F_y(n)=(1+n_y^3)/2 is normalized, probability-valued, "
                "noncontextual on the whole web, and non-affine. The source-attached "
                "real S3 algebraic contexts obtained by applying a declared representation to source-realized gauge labels are not complex tomographically "
                "complete: distinct pure Pauli-Y states agree on every declared "
                "outcome and a missing complex Y projector separates them. "
                "For two source-attached algebraic projector candidates P,Q, the algebraic "
                "phase lift I/2-(2sqrt(3)/3)i(QP-PQ) is exactly that Y "
                "projector and completes fixed-trace tomography; every effect "
                "in a generous real/Kraus closure remains Y-blind and cannot "
                "produce it. Entrywise conjugation exchanges the two phase "
                "candidates while simultaneous state-effect conjugation "
                "preserves the real Born weight. A post-hoc raw-count product-gap "
                "diagnostic from retained B12 data supplies an exact reversal-odd "
                "bit and normalized cycle products, but neither its statistic nor "
                "designation rule was preregistered; its pairing with the phase "
                "torsor is declared and emits no source selection or validation"
            ),
            "observed_counterpart": (
                "finite noncontextual weight representation and the Born rule"
            ),
            "match": (
                "exact bounded rank-gap and phase-free no-gos plus an algebraic "
                "complex tomography target; public quantum instrument missing"
            ),
            "lean_declarations": {
                "FiniteBornFrame": [
                    "contextAdditive_unique_parameterization",
                    "exists_frameCentered_iff_frameRelations",
                    "hermitianRepresentation_unique",
                    "densityRepresentation_unique",
                    "exists_admissible_not_densityRepresentable",
                ],
                "FiniteEffectClosureBoundary": [
                    "continuous_celestialBinaryWeight",
                    "nonlinearBinaryWeight_mem_Icc",
                    "nonlinearBinaryWeight_antipodal_sum",
                    "nonlinearBinaryWeight_not_affine",
                    "dense_affine_probability_tests_force_closed_unit_ball",
                ],
                "FiniteWebBornNoGo": [
                    "planarCubicAssignment_noncontextual",
                    "finiteBuschGleasonInterface_false",
                    "current_finite_web_born_no_go",
                ],
                "SourceContextTomographyNoGo": [
                    "current_source_context_web_not_tomographically_complete",
                    "pauliY_context_distinguishes_states",
                ],
                "SourcePhaseLiftBridge": [
                    "sourcePhaseLift_eq_rhoYPlus",
                    "sourcePhaseLift_mem_complexSourceAlgebra",
                    "realSourceEffectClosure_not_tomographically_complete",
                    "sourcePhaseTomography_injective_on_equalTrace",
                    "sourcePhaseLift_boundary_summary",
                ],
                "ConjugationGauge": [
                    "bornWeight_re_matrixConj",
                    "conj_invisible_on_fixed_effects",
                    "yStates_are_conj_orbit",
                ],
                "RepairCurrentOrientation": [
                    "designatedPair_lexLeast",
                    "designatedCycle_lexLeast",
                    "designatedCycle_normalized_products",
                    "reversal_flips_orientation",
                    "reversibleControl_no_orientation",
                ],
                "SourceOrientedCompletion": [
                    "orientationApplicable_holds",
                    "reversal_selects_conjugate",
                    "completionTomography_injective_on_states",
                    "oriented_born_capstone",
                ],
            },
            "lean_receipts": _lean_receipt(
                "FiniteBornFrame",
                "FiniteEffectClosureBoundary",
                "FiniteWebBornNoGo",
                "SourceContextTomographyNoGo",
                "SourcePhaseLiftBridge",
                "ConjugationGauge",
                "RepairCurrentOrientation",
                "SourceOrientedCompletion",
                declarations={
                    "FiniteBornFrame": (
                        "contextAdditive_unique_parameterization",
                        "exists_frameCentered_iff_frameRelations",
                        "hermitianRepresentation_unique",
                        "densityRepresentation_unique",
                        "exists_admissible_not_densityRepresentable",
                    ),
                    "FiniteEffectClosureBoundary": (
                        "continuous_celestialBinaryWeight",
                        "nonlinearBinaryWeight_mem_Icc",
                        "nonlinearBinaryWeight_antipodal_sum",
                        "nonlinearBinaryWeight_not_affine",
                        "dense_affine_probability_tests_force_closed_unit_ball",
                    ),
                    "FiniteWebBornNoGo": (
                        "planarCubicAssignment_noncontextual",
                        "finiteBuschGleasonInterface_false",
                        "current_finite_web_born_no_go",
                    ),
                    "SourceContextTomographyNoGo": (
                        "current_source_context_web_not_tomographically_complete",
                        "pauliY_context_distinguishes_states",
                    ),
                    "SourcePhaseLiftBridge": (
                        "sourcePhaseLift_eq_rhoYPlus",
                        "sourcePhaseLift_mem_complexSourceAlgebra",
                        "realSourceEffectClosure_not_tomographically_complete",
                        "sourcePhaseTomography_injective_on_equalTrace",
                        "sourcePhaseLift_boundary_summary",
                    ),
                    "ConjugationGauge": (
                        "bornWeight_re_matrixConj",
                        "conj_invisible_on_fixed_effects",
                        "yStates_are_conj_orbit",
                    ),
                    "RepairCurrentOrientation": (
                        "designatedPair_lexLeast",
                        "designatedCycle_lexLeast",
                        "designatedCycle_normalized_products",
                        "reversal_flips_orientation",
                        "reversibleControl_no_orientation",
                    ),
                    "SourceOrientedCompletion": (
                        "orientationApplicable_holds",
                        "reversal_selects_conjugate",
                        "completionTomography_injective_on_states",
                        "oriented_born_capstone",
                    ),
                },
            ),
            "artifact_refs": [
                "code/born_frame/runtime/finite_born_frame_certificate.json",
                "code/born_frame/verify_finite_born_frame_independent.py",
                "code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json",
                "code/born_context_phase_lift/README.md",
                "code/born_context_phase_lift/verify_source_phase_lift.py",
                "code/born_context_phase_lift/test_source_phase_lift.py",
                "code/thermodynamics/repair_current_orientation/repair_current_payload.v3.json",
                "code/thermodynamics/repair_current_orientation/verify_repair_current_orientation.py",
                "code/thermodynamics/repair_current_orientation/test_verify_repair_current_orientation.py",
            ],
            "hypothesis_boundary": (
                "the central atoms and spinor projectors are different objects. "
                "The projector family is a declared mathematical adapter obtained "
                "by applying a two-dimensional representation to source-realized "
                "gauge labels, not a source-produced public quantum "
                "instrument. The celestial countermodel proves that continuity "
                "and normalized antipodal binary contexts still do not derive "
                "affinity; the transverse cubic refutes the current finite "
                "Busch--Gleason interface, and the Pauli-Y pair identifies the "
                "missing complex tomography direction. The exact phase lift "
                "constructs that direction only inside the complex operator "
                "algebra; the source has no phase producer, rotated or phase "
                "outcome receipts, common-preparation validation, or operational "
                "effect composition. Conjugation identifies one two-candidate "
                "orbit but does not prove that this orbit exhausts all hidden "
                "phase data. The repair-count bit is a post-hoc diagnostic on a "
                "locally hash-pinned B12 run; its statistic and designation rule "
                "were not preregistered, and the phase pairing is an arbitrary "
                "typed convention. The full-effect theorem still applies only "
                "after full coexistent-effect additivity is supplied. No physical "
                "Born derivation, observable, or prediction is emitted. Issue "
                "#702 owns the source-earned phase instrument, operational "
                "additivity, and public readback"
            ),
            "paper_ref": "observers paper, finite Born-frame rank audit",
        },
        {
            "id": "thermodynamic_four_law_package",
            "statement": (
                "With a faithful common reference and repaired-visible fibre "
                "supplied, Axiom 3 instantiated on states selects the Gibbs "
                "exponential family by the exact information-projection "
                "Pythagorean identity, and instantiated on transition "
                "distributions over the repaired-visible fibre selects "
                "weighted conditional resampling from the same reference. "
                "The kernel is stochastic, idempotent, reversible, "
                "stationary, and fixes fibre-measurable charges; relative "
                "entropy to the reference contracts under it; the exact "
                "first-law split carries its bilinear cross term; the "
                "excited Gibbs mass obeys the finite gap bound with "
                "entropy limit log g0; partition pinching has an explicit "
                "normalized projector Kraus family and is formally CPTP"
                "; more generally, every stochastic kernel preserving a "
                "faithful stationary reference contracts relative entropy "
                "even without detailed balance, with an exact lazy directed "
                "three-cycle as the nonreversible separation witness. On the "
                "current source artifact, however, the transition action has a "
                "nonconstant eigenmode with eigenvalue 665437/726948, whereas "
                "the candidate state-side heat-bath action is idempotent; every "
                "intertwiner kills that mode. Its stationary mass 7155/61511 is "
                "also not a deterministic pushforward of the equally weighted "
                "16384-state empirical table"
            ),
            "observed_counterpart": (
                "the zeroth, first, second, and third laws of "
                "thermodynamics"
            ),
            "match": "finite theorem package under named receipts",
            "lean_declarations": {
                "FiniteConditionalRepair": [
                    "gibbs_pythagorean",
                    "gibbs_minimizer",
                    "heatBath_row_optimal",
                    "heatBath_secondLaw",
                    "heatBath_detailedBalance",
                    "heatBath_fixes_fiberObservable",
                    "kl_push_le",
                    "excitedMass_le",
                    "excitedMass_lt_of_beta_large",
                    "gibbs_beta_injective",
                    "clausius",
                    "landauer",
                    "heatBath_preserves_pos",
                    "mixture_stochastic",
                    "mixture_stationary",
                    "block_entropy_le",
                ],
                "StationaryRealization": [
                    "stationary_secondLaw",
                    "constantObservable_fixed",
                    "directedLazy3_stationary",
                    "directedLazy3_not_detailedBalance",
                    "directedLazy3_secondLaw",
                ],
                "FirstLawIdentity": ["firstLaw_split"],
                "FluctuationTheorems": [
                    "integral_fluctuation",
                    "crooks_pointwise",
                    "crooks_level_set",
                    "sigma_mean_eq_kl_descent",
                    "correlation_symm",
                    "heatBath_integral_fluctuation",
                    "heatBath_crooks",
                    "heatBath_correlation_symm",
                ],
                "CapFirstLaw": [
                    "cap_firstLaw_exact",
                    "cap_firstLaw_split",
                    "cap_clausius_of_central_conserved",
                    "push_heatBath_fixes_mean",
                    "heatBath_cap_clausius",
                ],
                "EinsteinPremiseLink": [
                    "shannon_diff_eq_pairing_sub_kl",
                    "thermoFirstLawData_passes",
                    "thermo_first_law_on_simplex_tangent",
                    "repair_variation_mem_massZero",
                    "repair_variation_central_pairing_zero",
                ],
                "PartitionPinchingCP": [
                    "ProjectivePartition.kraus_complete",
                    "partitionPinching_kraus_form",
                ],
                "ChoiCPTP": [
                    "partitionPinching_isCPTP",
                ],
                "CommonReferenceObstruction": [
                    "mixingMode_eigenpair",
                    "no_nondegenerate_current_common_object_intertwiner",
                    "no_empirical_deterministic_stationary_pushforward",
                    "current_common_reference_obstruction_summary",
                ],
            },
            "lean_receipts": _lean_receipt(
                "FiniteConditionalRepair",
                "StationaryRealization",
                "FirstLawIdentity",
                "FluctuationTheorems",
                "CapFirstLaw",
                "EinsteinPremiseLink",
                "PartitionPinchingCP",
                "ChoiCPTP",
                "CommonReferenceObstruction",
                declarations={
                    "FiniteConditionalRepair": (
                        "gibbs_pythagorean",
                        "gibbs_minimizer",
                        "heatBath_row_optimal",
                        "heatBath_secondLaw",
                        "heatBath_detailedBalance",
                        "heatBath_fixes_fiberObservable",
                        "kl_push_le",
                        "excitedMass_le",
                        "excitedMass_lt_of_beta_large",
                        "gibbs_beta_injective",
                        "clausius",
                        "landauer",
                        "heatBath_preserves_pos",
                        "mixture_stochastic",
                        "mixture_stationary",
                        "block_entropy_le",
                    ),
                    "StationaryRealization": (
                        "stationary_secondLaw",
                        "constantObservable_fixed",
                        "directedLazy3_stationary",
                        "directedLazy3_not_detailedBalance",
                        "directedLazy3_secondLaw",
                    ),
                    "FirstLawIdentity": ("firstLaw_split",),
                    "FluctuationTheorems": (
                        "integral_fluctuation",
                        "crooks_pointwise",
                        "crooks_level_set",
                        "sigma_mean_eq_kl_descent",
                        "correlation_symm",
                        "heatBath_integral_fluctuation",
                        "heatBath_crooks",
                        "heatBath_correlation_symm",
                    ),
                    "CapFirstLaw": (
                        "cap_firstLaw_exact",
                        "cap_firstLaw_split",
                        "cap_clausius_of_central_conserved",
                        "push_heatBath_fixes_mean",
                        "heatBath_cap_clausius",
                    ),
                    "EinsteinPremiseLink": (
                        "shannon_diff_eq_pairing_sub_kl",
                        "thermoFirstLawData_passes",
                        "thermo_first_law_on_simplex_tangent",
                        "repair_variation_mem_massZero",
                        "repair_variation_central_pairing_zero",
                    ),
                    "PartitionPinchingCP": (
                        "ProjectivePartition.kraus_complete",
                        "partitionPinching_kraus_form",
                    ),
                    "ChoiCPTP": (
                        "partitionPinching_isCPTP",
                    ),
                    "CommonReferenceObstruction": (
                        "mixingMode_eigenpair",
                        "no_nondegenerate_current_common_object_intertwiner",
                        "no_empirical_deterministic_stationary_pushforward",
                        "current_common_reference_obstruction_summary",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "The exact obstruction closes issue #688 only as a bounded "
                "negative result for two direct mechanisms on the current "
                "artifact: a mixing-mode-retaining linear intertwiner into the "
                "idempotent heat bath and a deterministic empirical pushforward. "
                "It does not exclude stochastic, nonlinear, reverse-direction, "
                "dilated, or enriched-source constructions. Issue #725 owns the "
                "replacement common reference, collar, objective, and genuinely "
                "varying refinement family. Issue #703 separately owns physical energy "
                "and clock calibration beyond the attained central-interface "
                "modular split. The "
                "pinned 20-state collar table has an audit of all "
                "15 nonempty field-subset projections: its repair-load quotient "
                "is an eight-state raw ergodic nonreversible H-theorem "
                "probe, but the declared record charge is constant and "
                "the common reference is unidentified; its only fine-chain "
                "recurrent restriction is singleton freezeout. The "
                "strict-descent normalizer carries no entropy inequality"
            ),
            "paper_ref": "observers paper, thermodynamics section",
        },
        {
            "id": "finite_green_kubo_graph_transport",
            "statement": (
                "A finite reversible Markov kernel with a linear Poisson "
                "solver has a symmetric positive-semidefinite Green--Kubo "
                "matrix and an exact finite correlation sum with propagated "
                "remainder. One full-fibre repair projector has constant "
                "positive-lag correlation and cannot supply a nonzero decaying "
                "memory tail with stabilizing sums. Typed finite-graph Fick "
                "and Fourier updates obey exact source balance and source-free "
                "conservation"
            ),
            "observed_counterpart": (
                "Onsager symmetry, Green--Kubo response, Fick diffusion, and "
                "Fourier heat transport"
            ),
            "match": (
                "exact finite conditional transport structure; physical "
                "generator and coefficients open"
            ),
            "lean_declarations": {
                "GreenKubo": [
                    "dissipation_eq_dirichlet",
                    "greenKuboPair_symm_of_poisson",
                    "greenKuboPair_finite_matrix_psd",
                    "greenKuboPair_eq_integratedCorrelation_add_remainder",
                    "heatBath_integratedCorrelation_not_stable",
                    "heatBath_integratedCorrelation_eq_equalTime",
                    "binary_heatBath_integratedCorrelation_eq_one",
                    "identityKernel_no_poisson_of_ne_zero",
                ],
                "GraphDiffusion": [
                    "summation_by_parts",
                    "fickParticleAmountStep_bridge",
                    "fickParticleAmountStep_total_conservation",
                    "fourierEnergyStep_bridge",
                    "fourierEnergyStep_total_conservation",
                    "fick_flux_gradient_power_nonpositive",
                    "fourier_flux_gradient_power_nonpositive",
                    "negative_conductance_counterexample",
                    "negative_thermal_conductance_counterexample",
                    "twoVertex_fick_closed_step",
                    "twoVertex_fourier_closed_step",
                ],
            },
            "lean_receipts": _lean_receipt(
                "GreenKubo",
                "GraphDiffusion",
                declarations={
                    "GreenKubo": (
                        "dissipation_eq_dirichlet",
                        "greenKuboPair_symm_of_poisson",
                        "greenKuboPair_finite_matrix_psd",
                        "greenKuboPair_eq_integratedCorrelation_add_remainder",
                        "heatBath_integratedCorrelation_not_stable",
                        "heatBath_integratedCorrelation_eq_equalTime",
                        "binary_heatBath_integratedCorrelation_eq_one",
                        "identityKernel_no_poisson_of_ne_zero",
                    ),
                    "GraphDiffusion": (
                        "summation_by_parts",
                        "fickParticleAmountStep_bridge",
                        "fickParticleAmountStep_total_conservation",
                        "fourierEnergyStep_bridge",
                        "fourierEnergyStep_total_conservation",
                        "fick_flux_gradient_power_nonpositive",
                        "fourier_flux_gradient_power_nonpositive",
                        "negative_conductance_counterexample",
                        "negative_thermal_conductance_counterexample",
                        "twoVertex_fick_closed_step",
                        "twoVertex_fourier_closed_step",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the reversible kernel, linear Poisson solver, graph, distance, "
                "clock increment, volumes, heat capacities, and conductances "
                "are declared finite inputs. No theorem identifies the "
                "Green--Kubo coefficient with graph conductance. Issues #725, "
                "#693, #694, and #703 own the source evolution, physical equilibrium "
                "reference and conserved quantity, source-realized geometry, "
                "clock, and calibration; bounded algebraic C2 issue #690 is "
                "closed and this row emits no prediction-ladder "
                "entry"
            ),
            "paper_ref": "observers paper, finite transport theorem",
        },
        {
            "id": "finite_conservation_ward_precursor",
            "statement": (
                "On the exact twelve-port, thirty-seam incidence graph, a "
                "declared pointwise continuity update gives regional source "
                "minus outward flux, internal cancellation, and closed-graph "
                "conservation for zero-total source. Rational Gauss solutions "
                "exist exactly for neutral loads and form one translate of a "
                "nineteen-dimensional cycle kernel. For finite real linear "
                "state maps, all-state charge conservation is equivalent to "
                "the dual fixed-observable equation; an exact two-state "
                "counterexample shows that channel covariance alone does not "
                "imply conservation"
            ),
            "observed_counterpart": (
                "continuity, Gauss constraint, and protected-charge structure"
            ),
            "match": "exact finite precursor; physical Ward bridge open",
            "lean_declarations": {
                "RegionalContinuity": [
                    "regional_continuity",
                    "global_continuity",
                    "global_conservation_of_zero_total_source",
                    "FiniteContinuityWitness.regionalBalance",
                ],
                "DiscreteGauss": [
                    "gauss_solution_exists_iff_total_zero",
                    "rationalBoundarySection_is_gauss_solution",
                    "gauss_solution_iff_difference_is_cycle",
                    "gauss_cycle_space_finrank",
                ],
                "ProtectedCharge": [
                    "chargeExpectation_preserved_iff_dual_fixed",
                    "kernel_chargeExpectation_preserved_iff_pull_fixed",
                    "twoStateOddCharge_ne_zero",
                    "channel_covariance_does_not_imply_charge_conservation",
                ],
                "WardLimitManifest": [
                    "WardLimitManifest.wardIdentity",
                ],
            },
            "lean_receipts": _lean_receipt(
                "RegionalContinuity",
                "DiscreteGauss",
                "ProtectedCharge",
                "WardLimitManifest",
                declarations={
                    "RegionalContinuity": (
                        "regional_continuity",
                        "global_continuity",
                        "global_conservation_of_zero_total_source",
                        "FiniteContinuityWitness.regionalBalance",
                    ),
                    "DiscreteGauss": (
                        "gauss_solution_exists_iff_total_zero",
                        "rationalBoundarySection_is_gauss_solution",
                        "gauss_solution_iff_difference_is_cycle",
                        "gauss_cycle_space_finrank",
                    ),
                    "ProtectedCharge": (
                        "chargeExpectation_preserved_iff_dual_fixed",
                        "kernel_chargeExpectation_preserved_iff_pull_fixed",
                        "twoStateOddCharge_ne_zero",
                        "channel_covariance_does_not_imply_charge_conservation",
                    ),
                    "WardLimitManifest": (
                        "WardLimitManifest.wardIdentity",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the pointwise update is a declared finite premise; the load, "
                "current, update order, and charge have no physical identity. "
                "The guarded WardLimitManifest derives zero limiting residual "
                "from exact finite-residual vanishing and convergence on a "
                "shrinking scale with separating tests. Issue #694 must define "
                "that residual from finite continuity, identify it with physical "
                "distributional divergence, and supply the source, transport, "
                "chart, and common-tower evidence before continuum Ward use"
            ),
            "paper_ref": "screen-microphysics paper, finite conservation bridge",
        },
        {
            "id": "finite_fixed_word_locality_and_marginal_invariance",
            "statement": (
                "For the concrete single-site OPH localRepair, every fixed "
                "exogenous word of n sequential moves has an n-fold "
                "closed-neighborhood dependency upper bound. Separately, on "
                "a supplied finite bipartite split, row-normalized real maps "
                "preserve the remote algebraic marginal and Kraus-complete "
                "local matrix families preserve the remote partial trace"
            ),
            "observed_counterpart": (
                "finite propagation cones and classical or quantum "
                "no-signalling"
            ),
            "match": (
                "exact fixed-word and algebraic helpers; physical causality "
                "attachment open"
            ),
            "lean_declarations": {
                "DependencyCone": [
                    "localRepair_agree",
                    "applyWord_agree_on",
                    "no_influence_outside_ball",
                ],
                "NoSignalling": [
                    "sndMarginal_pushJoint_liftFst",
                    "remote_mass_changes_without_row_normalization",
                    "ptraceFst_local_kraus",
                ],
            },
            "lean_receipts": _lean_receipt(
                "DependencyCone",
                "NoSignalling",
                declarations={
                    "DependencyCone": (
                        "localRepair_agree",
                        "applyWord_agree_on",
                        "no_influence_outside_ball",
                    ),
                    "NoSignalling": (
                        "sndMarginal_pushJoint_liftFst",
                        "remote_mass_changes_without_row_normalization",
                        "ptraceFst_local_kraus",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the repair word is fixed externally and shared by both "
                "inputs; adaptive or globally state-dependent scheduling is "
                "not covered. The cone is an upper bound, not a minimal cone "
                "or graph-radius speed law. The bipartite split is supplied, "
                "the classical theorem permits signed arrays, and no OPH "
                "region-factor, spacelike, clock, stochastic-state, CPTP, or "
                "laboratory attachment is proved. A declared rich-fibre adapter "
                "conditionally attains coverage, but its source supplies no "
                "operators and its corners are not tensor factors. The row's "
                "sole live gate, #692, owns source-attached operator generation "
                "and a region-product, TensorSplitReceipt, or local-channel adapter. "
                "Source channel/adaptive-scheduler semantics (#693), physical "
                "clocks (#703), physical spacetime attachment (#694), and continuum causal/time-slice "
                "structure (#700) are "
                "downstream promotions outside this claim's gate. This row "
                "emits no prediction-ladder entry"
            ),
            "paper_ref": "consensus-protocol paper, finite locality boundary",
        },
        {
            "id": "conditional_adaptive_scheduler_locality_helper",
            "statement": (
                "For concrete localRepair, a supplied adaptive scheduler and "
                "supplied ConsultsOnly consultation region have the exact "
                "n-step cone bound ball(S union R,n); a one-site change outside "
                "ball(S,n) union ball(R,n) cannot change the probe readout. "
                "A two-cell control proves the consultation term is "
                "indispensable. Declared refinement maps and one-step "
                "intertwining imply cone-image inclusion and run/readback "
                "naturality"
            ),
            "observed_counterpart": (
                "adaptive finite update locality under an explicitly bounded "
                "consultation region"
            ),
            "match": (
                "exact conditional helper; source scheduler and physical "
                "channel attachment open"
            ),
            "lean_declarations": {
                "AdaptiveScheduler": [
                    "adaptiveRun_agree_on",
                    "adaptive_no_influence",
                    "consultation_region_not_droppable",
                    "ball_image",
                    "run_natural",
                    "readback_cone_bound",
                ],
            },
            "lean_receipts": _lean_receipt(
                "AdaptiveScheduler",
                declarations={
                    "AdaptiveScheduler": (
                        "adaptiveRun_agree_on",
                        "adaptive_no_influence",
                        "consultation_region_not_droppable",
                        "ball_image",
                        "run_natural",
                        "readback_cone_bound",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "sigma, R, ConsultsOnly, and every ConeRefinement map/law are "
                "supplied. The helper proves neither their source production nor "
                "fairness, liveness, positivity, normalized state/channel, CPTP, "
                "distance, clock, spacelike, continuum, or laboratory semantics. "
                "Issue #693 remains the live E2 source scheduler/channel gate; "
                "#703 and #700 retain the clock and continuum-causality "
                "attachments. This row emits no prediction-ladder entry"
            ),
            "paper_ref": "E2 adaptive-scheduler helper, finite locality boundary",
        },
        {
            "id": "finite_history_variational_helpers_and_bridge_obstruction",
            "statement": (
                "A supplied positive normalized exponential tilt on a finite "
                "history type obeys the finite information-projection "
                "Pythagorean and minimizer identities; inverse-noise mass on "
                "all strict nonminimizers tends to zero. Separately, a real "
                "local-action minimum under every single-site variation gives "
                "the scalar discrete Euler--Lagrange equation and local "
                "Noether transport, with a nonzero free-path witness. No "
                "finite real-path family contains every such variation. For "
                "the committed binary chain, the bilinear corner extension "
                "has no velocity solver, while an infinite positive-curvature "
                "family agrees on every source history and gives regular "
                "strictly convex Lagrangians and Hamiltonians; its checked "
                "curvatures one and two are distinct on both faces"
            ),
            "observed_counterpart": (
                "Gibbs path selection, least-action limits, discrete "
                "Euler--Lagrange equations, and Noether currents"
            ),
            "match": (
                "exact conditional helpers plus finite/real and real-enrichment "
                "non-identifiability boundaries; physical composition is open"
            ),
            "lean_declarations": {
                "PathGibbs": [
                    "pathGibbs_pythagorean",
                    "pathGibbs_minimizer",
                    "modal_path_least_action",
                    "noiseFamily_above_gap_mass_tendsto_zero",
                    "noiseFamily_nonminimal_mass_tendsto_zero",
                ],
                "DiscreteEulerLagrange": [
                    "stationary_localAction_discreteEulerLagrange",
                ],
                "DiscreteNoether": [
                    "noether_conserved",
                    "noether_current_constant_on_finite_chain",
                    "free_translation_nonzero_noether_witness",
                ],
                "FiniteHistoryBridge": [
                    "exists_real_site_variation_outside",
                    "not_all_real_site_variations_mem",
                ],
                "RealizedHistoryLegendreNoGo": [
                    "chainLogLagrangian_no_velocity_solver",
                    "chainCurvedLagrangian_realized_indistinguishable",
                    "chainCurvedLagrangian_one_two_midpoint_gap",
                    "chainCurved_legendreTransform",
                    "realizedHistory_legendre_nonidentifiability_receipt",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PathGibbs",
                "DiscreteEulerLagrange",
                "DiscreteNoether",
                "FiniteHistoryBridge",
                "RealizedHistoryLegendreNoGo",
                declarations={
                    "PathGibbs": (
                        "pathGibbs_pythagorean",
                        "pathGibbs_minimizer",
                        "modal_path_least_action",
                        "noiseFamily_above_gap_mass_tendsto_zero",
                        "noiseFamily_nonminimal_mass_tendsto_zero",
                    ),
                    "DiscreteEulerLagrange": (
                        "stationary_localAction_discreteEulerLagrange",
                    ),
                    "DiscreteNoether": (
                        "noether_conserved",
                        "noether_current_constant_on_finite_chain",
                        "free_translation_nonzero_noether_witness",
                    ),
                    "FiniteHistoryBridge": (
                        "exists_real_site_variation_outside",
                        "not_all_real_site_variations_mem",
                    ),
                    "RealizedHistoryLegendreNoGo": (
                        "chainLogLagrangian_no_velocity_solver",
                        "chainCurvedLagrangian_realized_indistinguishable",
                        "chainCurvedLagrangian_one_two_midpoint_gap",
                        "chainCurved_legendreTransform",
                        "realizedHistory_legendre_nonidentifiability_receipt",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the finite source packet and its log-transition corner action "
                "are attained, while the reference remains declared. A typed "
                "finite-to-real transfer exists only under an undercut receipt. "
                "The complete binary history law still cannot select the "
                "off-alphabet curvature of a regular Legendre system: the "
                "displayed family is constructed, not source-produced. Issue "
                "#683 remains open for a source-selected reference and real "
                "enrichment, saddle histories, physical action and clock, "
                "amplitudes, fields, continuum, and observable currents. This "
                "row emits no prediction-ladder entry"
            ),
            "paper_ref": "observers paper, conditional history boundary",
        },
    ]
    rows.extend(
        [
            {
                "id": "maxwell_classical_massless_kernel",
                "statement": (
                    "On the declared unbroken Maxwell action and deconfined "
                    "phase branch, the quadratic operator has zero hard mass "
                    "parameter and two transverse classical modes with "
                    "characteristic surface k^2=0"
                ),
                "observed_counterpart": (
                    "massless classical electromagnetic propagation"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "the Maxwell action, positive kinetic coefficient, field "
                    "content, and phase are supplied branch data; no photon "
                    "Hilbert space, positive-residue pole, or universal "
                    "zero-mass particle theorem is emitted"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "yang_mills_classical_massless_kernel",
                "statement": (
                    "On the declared pure Yang-Mills quadratic branch before "
                    "nonperturbative confinement, every color generator has "
                    "two transverse perturbative modes and zero hard "
                    "quadratic mass parameter"
                ),
                "observed_counterpart": (
                    "perturbative color-gauge kernel before confinement"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "this is not a free asymptotic-gluon claim and supplies "
                    "neither a continuum Yang-Mills gap nor a hadron mass"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "einstein_classical_massless_kernel",
                "statement": (
                    "On the declared pure Einstein-Hilbert linearization about "
                    "a suitable Ricci-flat background, the transverse-"
                    "traceless quadratic operator has zero hard mass parameter "
                    "and two classical modes with null characteristic"
                ),
                "observed_counterpart": (
                    "two massless classical gravitational-wave polarizations"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "the action and background are supplied branch data; no "
                    "graviton Hilbert space, quantum pole, or exclusion of "
                    "additional massive modes is emitted"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "simple_gut_xy_channel_absent",
                "statement": (
                    "The declared charged-double-triplet current fixture has "
                    "a direct-sum algebra with adjoint "
                    f"branch dimensions {color_adjoint_dimension}, "
                    f"{weak_adjoint_dimension}, and {abelian_dimension}. Its "
                    "adjoint therefore contains no mixed "
                    "(3,2,-5/6) (+) (bar3,2,+5/6) X/Y generator, so the "
                    "ordinary minimal simple-GUT X/Y exchange channel is absent"
                ),
                "observed_counterpart": (
                    "the Standard Model product adjoint contains no connected "
                    "simple-GUT X/Y generator"
                ),
                "match": "conditional algebraic channel exclusion",
                "artifact_ref": _rel("port_current"),
                "operator_census_ref": (
                    "code/a5_closure/receipts/"
                    "baryon_dimension_six_census.receipt.json"
                ),
                "derivation_kind": "direct_executable_algebraic_corollary",
                "adjoint_branching": adjoint_branching,
                "hypothesis_boundary": (
                    "the executable corollary applies to the declared "
                    "direct-sum matrix-current fixture. Its physical current "
                    "source gate is false, so the result is not a physical "
                    "current or proton-stability claim. General proton "
                    "stability does not follow. Conditional on the "
                    "declared one-generation matter table and baryon and "
                    "lepton labels, an exact dimension-six census admits "
                    "QQQL, QQUE, DUQL, and DUUE; the representatives remain "
                    "nonzero after the exterior-algebra relations. No "
                    "coefficient, physical decay amplitude, QCD matrix "
                    "element, or lifetime is supplied"
                ),
                "paper_ref": "Observers paper, gauge-channel boundary",
            },
        ]
    )
    return rows


def _alpha_rows(
    endpoint: dict[str, Any],
    bridge: dict[str, Any],
    alpha_hvp_verdict: dict[str, Any],
) -> list[dict[str, Any]]:
    ep = endpoint["endpoint"]
    co = endpoint["compare_only"]
    bridge_verdict = bridge["verdict"]
    if "reference_deficit_inside_certified_gap" not in bridge_verdict:
        raise SystemExit("anchor bridge artifact lacks the containment verdict field")
    cross_class = alpha_hvp_verdict["cross_class_agreement"]
    scope = alpha_hvp_verdict["scope"]
    if (
        alpha_hvp_verdict["verdict"]
        != "MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"
        or cross_class["recorded_accounting_replay_count"] != 1
        or cross_class["independently_evaluated_class_count"] != 0
        or cross_class["verdict"] != "NOT_EVALUABLE_NO_INDEPENDENT_CLASS"
        or scope["comparison_timing"] != "retrospective"
        or scope["prospective_freeze"] is not False
        or scope["physical_alpha_prediction_emitted"] is not False
    ):
        raise SystemExit("alpha/HVP verdict has left its retrospective audit boundary")
    containment = bridge_verdict["reference_deficit_inside_certified_gap"]
    return [
        {
            "id": "alpha_inv_thomson_endpoint",
            "value_central": float(ep["alpha_inv_central"]),
            "value_interval": [float(v) for v in ep["alpha_inv_interval"]],
            "measured": float(co["codata_alpha_inv"]),
            "measured_source": "CODATA 2022 via the endpoint artifact, compare-only",
            "payload_release": endpoint["inputs"]["payload_release"],
            "row_class": endpoint["row_class"],
            "tier": alpha_hvp_verdict["row_class"],
            "anchor_gap_interval": [float(v) for v in co["same_scheme_anchor_gap_interval_inv_alpha"]],
            "reference_deficit_inside_recorded_accounting_interval": containment,
            "audit_verdict": alpha_hvp_verdict["verdict"],
            "cross_class_agreement": cross_class,
            "reading": (
                "one retrospective KNT19 accounting row is arithmetically "
                "compatible with the recorded same-scheme interval. The "
                "multi-class HVP test is not evaluable because no independent "
                "frozen class is present. Containment does not identify the "
                "physical source of the gap, select one closure map, construct "
                "the same-quantity bridge, or close source-only transport"
            ),
            "artifact_refs": [
                _rel("endpoint"),
                _rel("anchor_bridge"),
                _rel("alpha_hvp_verdict"),
            ],
            "blocking_issues": [708, 696],
            "historical_blocking_issues": [425, 545],
        }
    ]


def _lepton_rows(
    surface: dict[str, Any],
    rectangle: dict[str, Any],
    coherent: dict[str, Any],
    koide: dict[str, Any],
) -> list[dict[str, Any]]:
    witness_point = rectangle["compare_only"].get("witness_point")
    if witness_point is None:
        raise SystemExit(
            "rectangle artifact lacks the witness_point block; rebuild the "
            "rectangle lane first"
        )
    width_floor = coherent.get("width_floor_audit")
    if width_floor is None:
        raise SystemExit(
            "coherent artifact lacks the width_floor_audit block; rebuild "
            "the coherent lane first"
        )
    witnesses = rectangle["compare_only"]["witness_masses_gev"]
    particles = [r["particle"] for r in rectangle["conditional_mass_rows"]]
    family = next(f for f in surface["families"] if f["family"] == "charged leptons")
    mcpr = next(r for r in family["rows"] if r["lane"].startswith("MCPR"))
    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "id": "charged_leptons_closure_target",
            "statement": (
                "The certified solve inverts exactly at the measured triple: "
                "one anchor-gap value closes the charged-lepton lane on the "
                "witness, that value lies inside the retrospective accounting "
                "interval, and its "
                "distance to the standard on-shell reference deficit is the "
                "scheme term of the resource-deferred historical anchor "
                "boundary #545. The live quantitative successor #696 does "
                "not discharge that missing source. The lepton "
                "scale is localized only under the recorded accounting "
                "packet. A "
                "source-emitted bridge value is a sharp falsification "
                "target: landing on the closure value satisfies the "
                "conditional lane on the witness, while landing outside the "
                "recorded interval refutes "
                "the decomposition."
            ),
            "witness_point": witness_point,
            "width_floor": width_floor["floor_attribution"],
            "tier": "T1_empirical_closure",
            "artifact_refs": [_rel("kappa_rectangle"), _rel("kappa_coherent")],
            "blocking_issues": [696, 697],
            "historical_blocking_issues": [425, 545],
        }
    )
    mcpr_masses = [float(m) / 1000.0 for m in mcpr["masses_MeV_display"]]
    rows.append(
        {
            "id": "charged_leptons_mcpr_conditional",
            "particles": particles,
            "values_gev": mcpr_masses,
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the kappa lane, compare-only",
            "relative_deltas": [v / w - 1.0 for v, w in zip(mcpr_masses, witnesses, strict=True)],
            "tier": mcpr["tier"],
            "row_class": mcpr["row_class"],
            "explanation": mcpr["explanation"],
            "artifact_ref": mcpr["artifact"],
            "blocking_objects": mcpr["blocking_objects"],
        }
    )
    for key, lane, artifact in (
        ("charged_leptons_kappa_rectangle", rectangle, "kappa_rectangle"),
        ("charged_leptons_kappa_coherent", coherent, "kappa_coherent"),
    ):
        mass_rows = lane["conditional_mass_rows"]
        containment = all(
            row["mass_interval"][0] < w < row["mass_interval"][1]
            for row, w in zip(mass_rows, witnesses, strict=True)
        )
        entry: dict[str, Any] = {
            "id": key,
            "particles": particles,
            "intervals_gev": [row["mass_interval"] for row in mass_rows],
            "centrals_gev": [row["mass_central"] for row in mass_rows],
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the lane, compare-only",
            "witness_inside_all_intervals": containment,
            "relative_half_widths": [
                (row["mass_interval"][1] - row["mass_interval"][0])
                / (2.0 * row["mass_central"])
                for row in mass_rows
            ],
            "logarithmic_half_width": (
                lane["kappa_interval"]["interval"][1]
                - lane["kappa_interval"]["interval"][0]
            )
            / 2.0,
            "one_sided_multiplicative_widths": {
                "lower": 1.0
                - math.exp(
                    -(
                        lane["kappa_interval"]["interval"][1]
                        - lane["kappa_interval"]["interval"][0]
                    )
                    / 2.0
                ),
                "upper": math.exp(
                    (
                        lane["kappa_interval"]["interval"][1]
                        - lane["kappa_interval"]["interval"][0]
                    )
                    / 2.0
                )
                - 1.0,
            },
            "tier": "T1_empirical_closure",
            "row_class": lane["row_class"],
            "epistemic_scope": lane["numerical_certificate"]["epistemic_scope"],
            "numerical_certificate": lane["numerical_certificate"],
            "artifact_ref": _rel(artifact),
            "blocking_issues": [696, 697],
            "historical_blocking_issues": [425, 545],
        }
        if key.endswith("coherent"):
            entry["width_reduction_factor"] = lane["kappa_interval"]["width_reduction_factor"]
            entry["premise"] = "payload-coherent anchor-gap premise, declared"
        rows.append(entry)

    tau_row = koide["conditional_tau"]
    balance = koide["balance_comparison"]
    rows.append(
        {
            "id": "charged_leptons_koide_conditional_tau",
            "premises": tau_row["premises"],
            "inputs": tau_row["inputs"],
            "tau_enclosure_mev_outward": tau_row["tau_enclosure_mev_outward"],
            "tau_central_mev": tau_row["tau_central_mev"],
            "measured_tau_mev": tau_row["measured_tau_mev"],
            "distance_sigma": tau_row["distance_sigma"],
            "spurious_root_excluded_by": tau_row["spurious_root_excluded_by"],
            "balance_target_inside_enclosure": balance["target_inside_enclosure"],
            "balance_distance_half_widths": balance["distance_in_enclosure_half_widths"],
            "forward_test": (
                "the enclosure is three orders of magnitude narrower than "
                "the measurement uncertainty; an improving tau-mass average "
                "outside the window refutes the balanced-circulant premise"
            ),
            "tier": "T2_conditional",
            "row_class": tau_row["row_class"],
            "artifact_ref": _rel("koide_balance"),
        }
    )
    return rows


def _ew_rows(conditional: dict[str, Any]) -> list[dict[str, Any]]:
    cc = conditional["comparison_compare_only"]
    rows = []
    for key in ("mH_gev", "mt_pole_gev", "MW_chart_gev", "MZ_chart_gev"):
        block = cc[key]
        row: dict[str, Any] = {
            "id": f"ew_{key}",
            "value_central": block["conditional_central"],
            "value_envelope": block["conditional_envelope"],
            "physical_comparison_status": block["physical_comparison_status"],
            "tier": "T2_conditional",
            "row_class": conditional["row_class"],
            "artifact_ref": _rel("conditional_ew"),
        }
        if block["physical_comparison_status"] == "COMPARE_ONLY":
            row.update(
                {
                    "measured": block["measured"],
                    "measured_sigma": block["measured_sigma"],
                    "measured_source": block["measured_source"],
                    "delta": block["delta"],
                    "delta_over_sigma": block["delta_over_sigma"],
                    "envelope_overlaps_one_sigma_band": block["envelope_overlaps_one_sigma_band"],
                }
            )
        else:
            row["reason"] = block["reason"]
        rows.append(row)
    return rows


def _quark_rows(
    obstruction: dict[str, Any],
    clebsch: dict[str, Any],
    selection: dict[str, Any],
) -> list[dict[str, Any]]:
    if obstruction["fork"] != "ii_fiber_survives":
        raise SystemExit(
            "fiber obstruction artifact is not on the survives fork; "
            "rebuild the quark section before aggregating"
        )
    compare = clebsch["compare_only"]
    predictions = clebsch["predictions"]
    return [
        {
            "id": "quark_absolute_masses_obstruction",
            "statement": (
                "No absolute quark mass is emitted: the two-modulus spread "
                "fiber survives every certified structure transport, so the "
                "six absolute masses are non-identifiable from the corpus, by "
                "theorem rather than by omission"
            ),
            "fork": obstruction["fork"],
            "fiber_cut_detected": obstruction["fiber_cut_detected"],
            "tier": obstruction["claim_tier"],
            "artifact_ref": _rel("fiber_obstruction"),
            "blocking_issues": [697],
            "historical_blocking_issues": obstruction["github_issues"],
        },
        {
            "id": "quark_down_type_clebsch_route_rejected",
            "values": predictions,
            "measured_references": compare["references"],
            "relative_deltas": {
                k: v for k, v in compare.items() if k.endswith("_relative")
            },
            "flag_2024_compare_only": compare["flag_2024"],
            "retrospective_flag_rejection": clebsch[
                "retrospective_flag_rejection"
            ],
            "permutation_scan": clebsch["permutation_scan"],
            "promotion_allowed": clebsch["promotion_allowed"],
            "status": clebsch["status"],
            "tier": "T2_conditional_rejected_candidate",
            "row_class": clebsch["row_class"],
            "premise": (
                "a cross-sector register relation, independent Yukawa "
                "coefficient identification, and a physical generation order; "
                "the pairing receipt supplies channel compatibility only"
            ),
            "selection_artifact_ref": _rel("clebsch_selection"),
            "selection_status": selection["status"],
            "artifact_ref": _rel("clebsch_lane"),
            "reading": (
                "The target-free F1/F2 scan fixes only the unordered multiset. "
                "All six assignments fail the retrospective conservative FLAG "
                "gate. The displayed GST value is sqrt(md/ms) under an assumed "
                "texture, not a derived CKM angle; a simultaneous diagonal "
                "mass ansatz would instead give the identity CKM matrix."
            ),
        },
    ]


def _hadron_rows(payload: dict[str, Any], standby: dict[str, Any]) -> list[dict[str, Any]]:
    integral = payload["integral"]
    norm = integral["normalization"]
    return [
        {
            "id": "hadronic_correction_engine",
            "delta_alpha_had_5_MZ": integral["value"],
            "uncertainty_total": integral["uncertainty"],
            "source_compilation": payload["source_compilation"]["id"],
            "pin_factor": norm["pin_factor"],
            "policy": (
                "The published-compilation payload is the correction engine of "
                "the fine-structure lane; source-only hadron rows stay "
                "suppressed. Historical issue #425 records the resource-deferred "
                "QCD backend; live bounded particle-output ownership is #697, "
                "and source-only QCD remains outside the available resources."
            ),
            "artifact_ref": _rel("hadron_payload"),
        },
        {
            "id": "qcd_solver_on_standby",
            "status": standby["status"],
            "invocation_gate": standby["policy"]["invocation_gate"],
            "artifact_ref": _rel("solver_standby"),
        },
    ]


def _principal_results(sections: dict[str, Any]) -> list[dict[str, Any]]:
    """Digest the five strongest structural and quantitative rows."""

    forced = {r["id"]: r for r in sections["forced_structure"]}
    leptons = {r["id"]: r for r in sections["charged_leptons"]}
    target = leptons["charged_leptons_closure_target"]
    alpha = sections["alpha"][0]
    wp = target["witness_point"]
    glo, ghi = alpha["anchor_gap_interval"]
    return [
        {
            "id": "intrinsic_rank_three_response_completion",
            "statement": (
                forced["intrinsic_rank_three_response_completion"]["statement"]
                + ". This is an exact intrinsic metric completion; physical "
                "position, scale, refinement, and gluing remain open."
            ),
        },
        {
            "id": "forced_gauge_structure",
            "statement": (
                "Complete compact port response from A1 and endogenous overlap "
                "transport from A2 force the abstract Lie type "
                "su(3)+su(2)+u(1) on the twelve-port carrier. Target-blind "
                "readback independently derives R=-J. Inside the declared "
                "exterior-response algebra, an exhaustive scan selects the "
                "charge-conjugate rank-15 chiral anomaly-free pair and its "
                "one-generation hypercharge multiset; its common central "
                "kernel is Z6. Source reconstruction of the matrix current and "
                "matter action, physical global-form selection, laboratory "
                "attachment, and continuum quantum field theory remain separate."
            ),
        },
        {
            "id": "carrier_class_dispersion_surface",
            "statement": (
                forced["carrier_class_dispersion_band"]["statement"]
                + ". The class theorem is exact; physical field attachment, "
                "finite scale, coherent frame, readout, and comparison remain open."
            ),
        },
        {
            "id": "koide_conditional_tau_window",
            "statement": (
                "Under the balanced-circulant and mass-ordering premises the "
                "measured electron and muon masses fix the tau mass inside "
                f"[{leptons['charged_leptons_koide_conditional_tau']['tau_enclosure_mev_outward'][0]}, "
                f"{leptons['charged_leptons_koide_conditional_tau']['tau_enclosure_mev_outward'][1]}] MeV, "
                f"{leptons['charged_leptons_koide_conditional_tau']['distance_sigma']} sigma from "
                "measurement. The balance premise was abstracted from the "
                "measured lepton triple, so this is a target-informed "
                "conditional postdiction with a frozen rejection rule."
            ),
        },
        {
            "id": "lepton_closure_target",
            "statement": (
                "The anchor-gap value "
                f"{wp['required_anchor_gap_at_witness_inv_alpha']:.4f} closes "
                "the charged-lepton lane exactly on the measured triple, "
                f"inside the retrospective accounting interval "
                f"[{glo:.4f}, {ghi:.4f}]; the "
                f"distance {wp['scheme_term_difference_inv_alpha']:+.4f} to "
                "the standard on-shell reference deficit "
                f"{wp['reference_deficit_inv_alpha']:.4f} is the scheme "
                "term of the resource-deferred historical anchor boundary "
                "#545; live quantitative successor #696 does not discharge "
                "that source requirement. The lepton "
                "scale is localized only under that recorded accounting "
                "packet. A "
                "source-emitted bridge value is a falsification target: the "
                "closure value would satisfy the conditional lane, while a "
                "value outside the interval refutes the declared decomposition."
            ),
        },
    ]


def build(
    out_path: Path = DEFAULT_OUT,
    md_path: Path | None = DEFAULT_MD,
    *,
    write: bool = True,
) -> dict[str, Any]:
    surface = _load("mass_surface")
    conditional = _load("conditional_ew")
    endpoint = _load("endpoint")
    bridge = _load("anchor_bridge")
    rectangle = _load("kappa_rectangle")
    coherent = _load("kappa_coherent")
    koide = _load("koide_balance")
    clebsch = _load("clebsch_lane")
    selection = _load("clebsch_selection")
    obstruction = _load("fiber_obstruction")
    matter = _load("matter_receipt")
    matter_menu = _load("matter_menu")
    port_current = _load("port_current")
    axis_center_descent = _load("axis_center_descent")
    carrier_modes = _load("carrier_modes")
    quantum_carrier_status = _load("quantum_carrier_status")
    carrier_class = _load("carrier_class_dispersion")
    carrier_frequency = _load("carrier_frequency_speed")
    gauge_kinetic = _load("gauge_kinetic_invariant_forms")
    oriented_face = _load("oriented_face_bracket_selector")
    invariant_metric = _load("invariant_metric_phase")
    alpha_hvp_verdict = _load("alpha_hvp_verdict")
    payload = _load("hadron_payload")
    standby = _load("solver_standby")

    sections = {
        "forced_structure": _forced_structure(
            matter,
            matter_menu,
            port_current,
            axis_center_descent,
            carrier_modes,
            carrier_class,
            carrier_frequency,
            gauge_kinetic,
            oriented_face,
            invariant_metric,
        ),
        "quantum_carrier_status": _quantum_carrier_status_row(
            quantum_carrier_status
        ),
        "alpha": _alpha_rows(endpoint, bridge, alpha_hvp_verdict),
        "charged_leptons": _lepton_rows(surface, rectangle, coherent, koide),
        "electroweak": _ew_rows(conditional),
        "quarks": _quark_rows(obstruction, clebsch, selection),
        "hadrons": _hadron_rows(payload, standby),
        "neutrinos": [
            {
                "id": "neutrino_dimensionless_pointer",
                "statement": (
                    "dimensionless PMNS and mass-splitting-ratio "
                    "comparisons live on the results status surface; the "
                    "absolute attachment stays compare-only"
                ),
                "artifact_ref": "code/particles/RESULTS_STATUS.md",
            }
        ],
    }
    result = {
        "artifact": "oph_postdiction_ledger",
        "generator": "code/particles/scripts/build_postdiction_ledger.py",
        "schema_version": 2,
        "row_class": "compare_only_postdiction_ledger",
        "guards": {
            "compare_only": True,
            "public_promotion_allowed": False,
            "changes_any_solve_path": False,
            "new_axiom_introduced": False,
            "hand_typed_measured_values": False,
        },
        "aggregation_policy": (
            "numeric values and measured references are read live from cited "
            "parents; structural rows are derived from validated Lean "
            "declarations, structured parents, or both, and identify direct "
            "algebraic corollaries explicitly; a missing or inconsistent "
            "receipt aborts the build"
        ),
        "principal_results": _principal_results(sections),
        "sections": sections,
    }
    if write:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if md_path is not None:
            md_path.write_text(_render_md(result), encoding="utf-8")
    return result


def _fmt(x: float, digits: int = 6) -> str:
    return f"{x:.{digits}g}"


def _render_md(ledger: dict[str, Any]) -> str:
    s = ledger["sections"]
    lines: list[str] = []
    add = lines.append
    add("# Postdiction Ledger")
    add("")
    add(
        "Generated deterministically by "
        "`code/particles/scripts/build_postdiction_ledger.py`; the JSON artifact "
        "is `code/particles/runs/status/postdiction_ledger.json`."
    )
    add("")
    add("Numeric values and measured references on this page are read live from "
        "the cited parent artifacts. Structural rows are derived from validated "
        "Lean declarations, structured parents, or both, and direct algebraic "
        "corollaries are identified. "
        "The ledger promotes nothing and changes no solve path. Interval rows "
        "report containment of the compare-only witness; conditional rows carry "
        "their declared premises; chart coordinates keep their NOT_EVALUABLE "
        "physical-comparison status.")
    add("")
    add("## Principal results")
    add("")
    for entry in ledger["principal_results"]:
        add(f"- {entry['statement']}")
    add("")
    add("## Forced structure")
    add("")
    add("These finite structural results precede or constrain numeric lanes. "
        "They include the icosahedral gauge packet and generic observer-law "
        "boundaries. Each row is checked in Lean, by a structured executable "
        "artifact, or by both, and records its own classical inputs and open "
        "physical premises.")
    add("")
    add("| Result | Observed counterpart | Match | Receipts |")
    add("| --- | --- | --- | --- |")
    for row in s["forced_structure"]:
        receipts = list(row.get("lean_receipts", []))
        if row.get("artifact_ref"):
            receipts.append(row["artifact_ref"])
        receipts.extend(row.get("artifact_refs", []))
        receipt_txt = ", ".join(f"`{r}`" for r in receipts)
        add(f"| {row['statement']} | {row['observed_counterpart']} | "
            f"`{row['match']}` | {receipt_txt} |")
    add("")
    add("Lean declaration bindings:")
    add("")
    for row in s["forced_structure"]:
        bindings = row.get("lean_declarations")
        if not bindings:
            continue
        rendered = "; ".join(
            f"`{module}`: " + ", ".join(f"`{name}`" for name in names)
            for module, names in bindings.items()
        )
        add(f"- `{row['id']}`: {rendered}")
    add("")
    add("Hypothesis boundaries:")
    add("")
    for row in s["forced_structure"]:
        add(f"- `{row['id']}`: {row['hypothesis_boundary']}")
    add("")
    add("## Quantum carrier gate")
    add("")
    carrier_status = s["quantum_carrier_status"]
    add(
        "The exact conditional four-dimensional propagating-mode vector is "
        "`(2, 16, 2)`: two Maxwell modes from one U(1) generator, sixteen "
        "perturbative color modes from eight SU(3) adjoint generators, and two "
        "Einstein transverse-traceless modes from one metric tensor field. "
        "The entries are neither particle counts nor one uniform gauge-algebra "
        "dimension vector."
    )
    add("")
    add("| Carrier | Quantum verdict | Blocking frontier |")
    add("| --- | --- | --- |")
    for row in carrier_status["rows"]:
        blockers = ", ".join(f"`{item}`" for item in row["blocking_frontier"])
        add(
            f"| `{row['carrier_id']}` | `{row['verdict']}` | {blockers} |"
        )
    add("")
    add(
        "The target-named status packet consumes no laboratory comparison value, "
        "permits no particle promotion, and is ineligible as a blind prediction. "
        "Its receipt is "
        f"`{carrier_status['artifact_ref']}`."
    )
    add("")
    add("## Fine-structure lane")
    add("")
    for row in s["alpha"]:
        lo, hi = row["value_interval"]
        glo, ghi = row["anchor_gap_interval"]
        add(f"- `alpha_em^-1` Thomson endpoint: `{_fmt(row['value_central'], 10)}` "
            f"in `[{_fmt(lo, 10)}, {_fmt(hi, 10)}]` against CODATA "
            f"`{_fmt(row['measured'], 10)}` (compare-only). Payload release "
            f"`{row['payload_release']}`.")
        inside = (
            "inside"
            if row["reference_deficit_inside_recorded_accounting_interval"]
            else "outside"
        )
        add(f"- Recorded retrospective same-scheme accounting interval "
            f"`[{_fmt(glo, 4)}, {_fmt(ghi, 4)}]` inverse-alpha units; the "
            f"standard reference deficit sits {inside} that interval.")
        add(f"- Independent-class verdict: `{row['audit_verdict']}`; evaluated "
            f"independent classes: "
            f"`{row['cross_class_agreement']['independently_evaluated_class_count']}`.")
        add(f"- Reading: {row['reading']}")
        add(
            "- Live blocking issues: "
            + ", ".join(f"#{i}" for i in row["blocking_issues"])
        )
        add(
            "- Historical resource-deferred boundaries: "
            + ", ".join(f"#{i}" for i in row["historical_blocking_issues"])
        )
    add("")
    add("## Charged leptons")
    add("")
    for row in s["charged_leptons"]:
        if row["id"].endswith("closure_target"):
            wp = row["witness_point"]
            add(f"- Closure target ({row['tier']}): the anchor-gap value "
                f"`{wp['required_anchor_gap_at_witness_inv_alpha']:.4f}` closes the "
                "lane exactly on the measured triple (inversion machine-checked); "
                f"the distance `{wp['scheme_term_difference_inv_alpha']:+.4f}` to the "
                f"on-shell reference deficit `{wp['reference_deficit_inv_alpha']:.4f}` "
                "is the live scheme term of the bridge. The certified width floor "
                "is the scheme-band ambiguity; no budget is shrunk without the "
                "source bridge.")
            continue
        if row["id"].endswith("koide_conditional_tau"):
            lo, hi = row["tau_enclosure_mev_outward"]
            measured, sigma = row["measured_tau_mev"]
            add(f"- Koide conditional tau ({row['tier']}): under the "
                "balanced-circulant and mass-ordering premises the measured "
                "electron and muon masses fix the tau mass inside "
                f"`[{lo}, {hi}]` MeV, `{row['distance_sigma']}` sigma from "
                f"the measured `{measured} +- {sigma}` MeV; the premise "
                "ancestry is declared and improving tau-mass averages test "
                "the premise directly.")
            continue
        if row["id"].endswith("mcpr_conditional"):
            deltas = ", ".join(
                f"{p} `{_fmt(d * 1e6, 3)} ppm`"
                for p, d in zip(row["particles"], row["relative_deltas"], strict=True)
            )
            add(f"- MCPR conditional triple ({row['tier']}): {deltas} against the "
                "PDG witness triple; the eight-register architecture is a "
                "declared model input.")
        else:
            kind = "coherent closure" if row["id"].endswith("coherent") else "rectangle"
            contained = "inside" if row["witness_inside_all_intervals"] else "OUTSIDE"
            one_sided = row["one_sided_multiplicative_widths"]
            add(f"- Kappa interval, {kind} ({row['tier']}): outward-rounded "
                "target-anchored diagnostic intervals with logarithmic half-width "
                f"`{_fmt(row['logarithmic_half_width'] * 100, 4)}%` and one-sided "
                f"multiplicative widths `-{_fmt(one_sided['lower'] * 100, 3)}%` / "
                f"`+{_fmt(one_sided['upper'] * 100, 3)}%`; the witness triple lies "
                f"{contained} every interval.")
            if "width_reduction_factor" in row:
                add(f"  - Width reduction over the rectangle: "
                    f"`{_fmt(row['width_reduction_factor'], 3)}x`; premise: "
                    f"{row['premise']}.")
    add("")
    add("## Electroweak sector")
    add("")
    add("| Quantity | Conditional central | Envelope | Measured | Delta/sigma | Status |")
    add("| --- | ---: | --- | --- | ---: | --- |")
    for row in s["electroweak"]:
        env = f"[{_fmt(row['value_envelope'][0], 8)}, {_fmt(row['value_envelope'][1], 8)}]"
        if row["physical_comparison_status"] == "COMPARE_ONLY":
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                f"`{row['measured']} +- {row['measured_sigma']}` ({row['measured_source']}) | "
                f"`{_fmt(row['delta_over_sigma'], 3)}` | compare-only |")
        else:
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                "chart coordinate | n/a | NOT_EVALUABLE |")
    add("")
    add("W/Z rows are running/tree chart coordinates. The strict one-loop "
        "consumer has a separate external fixture: interval receipts exclude "
        "scalar zeros in the declared principal-sheet boxes and isolate, for "
        "each of W and Z, one simple scalar zero with derivative and scalar-residue "
        "balls in its declared lower-half pole box on a channel-specific algebraic "
        "chart. They identify neither chart with the physical resonance sheet and "
        "prove no unique continuation, sign bridge, full-matrix Laurent "
        "residue, physical-current amplitude, or independent numerical replay. "
        "The fixture is not composed with the OPH chart, so no physical W/Z pole "
        "or mass comparison is defined. The Higgs and top rows are conditional "
        "on the declared selection axioms.")
    add("")
    add("## Quarks")
    add("")
    for row in s["quarks"]:
        if row["id"].endswith("obstruction"):
            add(f"- Absolute masses ({row['tier']}): {row['statement']} "
                f"(issues {', '.join(f'#{i}' for i in row['blocking_issues'])}).")
        else:
            vals = row["values"]
            flag = row["flag_2024_compare_only"]
            flag_refs = ", ".join(
                f"Nf={entry['nf']}: {_fmt(entry['reference_ms_over_md'], 4)}"
                for entry in flag
            )
            add(
                f"- Down-type register-Clebsch route, rejected "
                f"({row['tier']}): `ms/md = {_fmt(vals['ms_over_md'], 4)}` "
                f"against FLAG 2024 ({flag_refs}); all six generation "
                f"assignments are rejected by the retrospective conservative "
                f"gate. The diagnostic `sqrt(md/ms) = "
                f"{_fmt(vals['cabibbo_gst_sqrt_md_over_ms'], 4)}` is not a "
                f"derived Cabibbo angle. Premise: {row['premise']}. "
                f"{row['reading']}"
            )
    add("")
    add("## Hadrons")
    add("")
    for row in s["hadrons"]:
        if row["id"] == "hadronic_correction_engine":
            add(f"- Correction engine payload: `Delta alpha_had^(5)(M_Z^2) = "
                f"{row['delta_alpha_had_5_MZ']} +- {row['uncertainty_total']}` "
                f"from `{row['source_compilation']}` "
                f"(pin factor `{_fmt(row['pin_factor'], 7)}`). {row['policy']}")
        else:
            add(f"- QCD solver: `{row['status']}`; invocation is gated on the "
                "source-side parameter emissions recorded in the standby receipt.")
    add("")
    add("## Neutrinos")
    add("")
    for row in s["neutrinos"]:
        add(f"- {row['statement']} (`{row['artifact_ref']}`).")
    add("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed JSON or Markdown differs from a live rebuild",
    )
    args = parser.parse_args()
    result = build(args.out, args.md, write=not args.check)
    if args.check:
        expected_json = json.dumps(result, indent=2, sort_keys=True) + "\n"
        expected_md = _render_md(result)
        if not args.out.is_file() or args.out.read_text(encoding="utf-8") != expected_json:
            raise SystemExit(f"postdiction ledger JSON drift: {args.out}")
        if not args.md.is_file() or args.md.read_text(encoding="utf-8") != expected_md:
            raise SystemExit(f"postdiction ledger Markdown drift: {args.md}")
        print("postdiction ledger parity OK")
        return
    for name, rows in result["sections"].items():
        print(f"{name}: {len(rows)} rows")
    print(f"wrote {args.out}")
    print(f"wrote {args.md}")


if __name__ == "__main__":
    main()
