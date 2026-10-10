#!/usr/bin/env python3
"""Aggregate the certified postdiction rows into one ledger.

The ledger is a deterministic aggregator.  Numeric values and measured
references are read mechanically from their parent artifacts.  Structural rows are
derived from analytic paper proofs and validated finite evidence, with any direct algebraic
corollary identified as such.  A missing or inconsistent parent is a hard
failure, not a silently absent row.

Section one records the forced-structure layer: machine-checked finite
theorems, analytic paper proofs and executable certificates that precede or constrain numeric lanes,
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
from fractions import Fraction
import importlib.util
import hashlib
import json
import math
import re
import sys
import types
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
    "p_interval_contraction": RUNTIME
    / "p_interval_contraction_certificate_2026-07-14.json",
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
    "lambda_transmutation": RUNS / "qcd" / "lambda_qcd_source_transmutation.json",
    "nucleon_external_ratio": RUNS / "hadron" / "nucleon_mass_external_qcd_ratio.json",
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
    "SeamMaxwellContinuum": LEAN_SCREEN / "SeamMaxwellContinuum.lean",
    "SerialMaxwellReadout": LEAN_SCREEN / "SerialMaxwellReadout.lean",
    "ConeCochainBridge": LEAN_SCREEN / "ConeCochainBridge.lean",
    "WhitneyTimeBridge": LEAN_SCREEN / "WhitneyTimeBridge.lean",
    "WhitneyMaxwellDynamics": LEAN_SCREEN / "WhitneyMaxwellDynamics.lean",
    "WhitneyQuantumBridge": LEAN_SCREEN / "WhitneyQuantumBridge.lean",
    "WhitneyChargedMatter": LEAN_SCREEN / "WhitneyChargedMatter.lean",
    "WhitneySpatialConsistency": LEAN_SCREEN / "WhitneySpatialConsistency.lean",
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
    "WeylYukawaConventions": LEAN_SCREEN / "WeylYukawaConventions.lean",
    "TimeOrderLedger": REPO / "Lean" / "Time" / "TimeOrderLedger.lean",
    "ObserverHistory": REPO / "Lean" / "Time" / "ObserverHistory.lean",
    "ClockReadout": REPO / "Lean" / "Time" / "ClockReadout.lean",
    "WorldlineRealization": REPO / "Lean" / "Time" / "WorldlineRealization.lean",
    "ProperTimeCalibration": REPO / "Lean" / "Time" / "ProperTimeCalibration.lean",
    "ClockComparison": REPO / "Lean" / "Time" / "ClockComparison.lean",
    "SemanticEventProvenance": REPO / "Lean"
    / "ObserverPatchHolography" / "Provenance"
    / "SemanticEventProvenance.lean",
    "CausalInterval": REPO / "Lean" / "ObserverPatchHolography"
    / "Provenance" / "CausalInterval.lean",
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
    "SourceDerivedSpacetimeCarrier": REPO / "Lean" / "Geometry"
    / "SourceDerivedSpacetimeCarrier.lean",
    "SourceOrderFrameCompatibilityPacket": REPO / "Lean" / "Geometry"
    / "SourceOrderFrameCompatibilityPacket.lean",
    "RefiningLatticeCausalCone": REPO / "Lean" / "Geometry"
    / "RefiningLatticeCausalCone.lean",
    "SourceNetCausalCone": REPO / "Lean" / "Geometry" / "SourceNetCausalCone.lean",
    "SourceCountClock": REPO / "Lean" / "Time" / "SourceCountClock.lean",
    "SourceRecordProtection": REPO / "Lean" / "Geometry" / "SourceRecordProtection.lean",
    "SourcePopulationQuadrature": REPO / "Lean" / "Geometry" / "SourcePopulationQuadrature.lean",
    "GoldenSourceAssignment": REPO / "Lean" / "Geometry" / "GoldenSourceAssignment.lean",
    "GoldenSourceCountLimit": REPO / "Lean" / "Geometry" / "GoldenSourceCountLimit.lean",
    "SourceFeedbackTransport": REPO / "Lean" / "Geometry" / "SourceFeedbackTransport.lean",
    "SourceSeamPathTomography": REPO / "Lean" / "Geometry" / "SourceSeamPathTomography.lean",
    "ScalarRegionalTimeSlice": REPO / "Lean" / "QFT" / "ScalarRegionalTimeSlice.lean",
    "MetricKernelEnergy": REPO / "Lean" / "Geometry" / "MetricKernelEnergy.lean",
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
    "SourceOperatorGeneration": REPO / "Lean" / "QFT"
    / "SourceOperatorGeneration.lean",
    "JointSlotFactorisation": REPO / "Lean" / "QFT"
    / "JointSlotFactorisation.lean",
    "CPRestrictionNet": REPO / "Lean" / "QFT" / "CPRestrictionNet.lean",
    "TwoSlotCPNetWitness": REPO / "Lean" / "QFT"
    / "TwoSlotCPNetWitness.lean",
    "TowerAnchoredDiamond": REPO / "Lean" / "QFT"
    / "TowerAnchoredDiamond.lean",
    "SourceCorrelationCapstone": REPO / "Lean" / "QFT"
    / "SourceCorrelationCapstone.lean",
    "StructuralNetAdequacySurface": REPO / "Lean" / "QFT"
    / "StructuralNetAdequacySurface.lean",
    "TripleCarrierJoin": REPO / "Lean" / "QFT" / "TripleCarrierJoin.lean",
    "TripleCarrierOperatorJoin": REPO / "Lean" / "QFT" / "TripleCarrierOperatorJoin.lean",
    "LocalFaceMaxwellAction": LEAN_SCREEN / "LocalFaceMaxwellAction.lean",
    "CofinalSpectralTailFamily": REPO / "Lean" / "Thermodynamics"
    / "CofinalSpectralTailFamily.lean",
    "ModalMaxwellFactorizationBoundary": LEAN_SCREEN
    / "ModalMaxwellFactorizationBoundary.lean",
    "FiniteUnitaryScatteringNoGo": REPO / "Lean" / "QFT"
    / "FiniteUnitaryScatteringNoGo.lean",
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
    "SourceReferenceSelection": REPO / "Lean" / "InformationProjection"
    / "SourceReferenceSelection.lean",
    "SourceHistoryPacket": REPO / "Lean" / "InformationProjection"
    / "SourceHistoryPacket.lean",
    "LogTransitionAction": REPO / "Lean" / "InformationProjection"
    / "LogTransitionAction.lean",
    "StationarySaddleCoverage": REPO / "Lean" / "Variational"
    / "StationarySaddleCoverage.lean",
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
    "LayeredDiscreteGauss": LEAN_SCREEN / "LayeredDiscreteGauss.lean",
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
                rf"(?m)^\s*(?:noncomputable\s+)?(?:theorem|lemma|def)\s+"
                rf"{re.escape(declaration)}(?:\s|:)"
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
        "graviton": (
            "NOT_EVALUABLE_NO_PHYSICAL_SOURCE_CAUSAL_CONTINUUM_OR_QUANTUM_CARRIER"
        ),
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


def _refining_causal_control() -> dict[str, Any]:
    """Keep the supplied move-law control distinct from native source order."""
    source = "paper/tex_fragments/REFINING_CAUSAL_CONE.tex"
    label = "prop:refining-causal-cone"
    path = REPO / source
    if not path.is_file() or "\\label{" + label + "}" not in path.read_text(encoding="utf-8"):
        raise SystemExit("supplied causal-cone analytic theorem missing")
    verifier_path = CODE / "causal_refinement" / "verify_refining_cone.py"
    spec = importlib.util.spec_from_file_location("refining_causal_ledger_verifier", verifier_path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent supplied causal-cone verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    summary = verifier.verify(verifier.load())
    expected = {
        "verdict": "PASS_SUPPLIED_LAW_CAUSAL_REFINEMENT_CONTROL",
        "path_witnesses": 16, "largest_radius": 32, "largest_inner_speed": "29/32",
        "history_events": 81, "authenticated_edges": 794, "exact_history_width": 27,
        "native_physical_spacetime_selected": False,
        "continuum_claim": "analytic inner/outer bound; supplied grid, clock and density",
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("supplied causal-cone verifier scope/count mismatch")
    declarations = (
        "causal_cone_sandwich", "scaled_inner_cone_reachable",
        "scaled_reachable_outer_cone", "fixed_menu_missing_timelike_endpoint",
    )
    return {
        "source": source, "label": label,
        "lean_declarations": {"RefiningLatticeCausalCone": list(declarations)},
        "lean_receipts": _lean_receipt(
            "RefiningLatticeCausalCone",
            declarations={"RefiningLatticeCausalCone": declarations},
        ),
        "independent_verifier_result": summary,
        "analytic_controlled_causal_limit": True,
        "analytic_interval_volume_limit": True,
        "continuum_limit_formalized_in_lean": False,
        "native_OPH_law_selected": False,
        "physical_clock_calibrated": False,
        "observed_postdiction": False,
    }


def _source_net_causal_control(receipt_path: Path | None = None) -> dict[str, Any]:
    """Separate the finite source execution from its analytic declared-law limit."""
    source = "paper/tex_fragments/SOURCE_NET_CAUSAL_LIMIT.tex"
    labels = (
        "prop:source-net-causal-path", "prop:source-net-alexandrov-volume",
        "cor:source-net-general-diamond", "prop:golden-source-count-limit",
        "cor:source-record-ordering-fraction",
    )
    path = REPO / source
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    if any("\\label{" + label + "}" not in text for label in labels):
        raise SystemExit("source-net analytic theorem missing")
    declarations = (
        "covering_constructs_path", "reachable_outer", "same_layer_precedes_iff",
        "finite_antichain_card_le", "cone_speed_identity",
    )
    lean_refs = _lean_receipt(
        "SourceNetCausalCone", declarations={"SourceNetCausalCone": declarations},
    )
    directory = CODE / "causal_refinement"
    spec = importlib.util.spec_from_file_location(
        "source_net_causal_ledger_verifier", directory / "verify_source_net_causet.py",
    )
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent source-net verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt_path = receipt_path or directory / "source_net_causet_receipt.json"
    raw = receipt_path.read_bytes()
    try:
        summary = verifier.verify(verifier.load(receipt_path))
    except ValueError as exc:
        raise SystemExit(f"source-net independent replay failed: {exc}") from exc
    if receipt_path.read_bytes() != raw:
        raise SystemExit("source-net receipt changed during replay")
    expected = {
        "accepted": True,
        "scope": "CONDITIONAL_SOURCE_CODED_CAUSETS__SUPPLIED_POPULATION_READ_LAW_AND_CLOCK",
        "native_physical_spacetime_selected": False,
        "levels": [
            {"q": q, "source_records": q**3, "width": q**3, "height": height,
             "reads_per_full_trace": reads, "positive_inner_cone": q == 13,
             "inner_speed_lower": "75989/250000" if q == 13 else "0"}
            for q, height, reads in ((3, 3, 622), (5, 4, 8355),
                                     (8, 4, 99192), (13, 5, 1176764))
        ],
        "packet_sha256": hashlib.sha256(raw).hexdigest(),
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("source-net verifier scope/count mismatch")
    return {
        "source": source, "labels": list(labels),
        "lean_declarations": {"SourceNetCausalCone": list(declarations)},
        "lean_receipts": lean_refs,
        "receipt_pin": {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()},
        "independent_verifier_result": summary,
        "analytic_controlled_causal_limit": True,
        "analytic_declared_family_raw_count_limit": True,
        "analytic_contained_interval_ordering_fraction": "1/10",
        "volume_convention": "dt d^3x",
        "primitive_source_words_executed_as_repairs": False,
        "continuum_or_pair_count_limit_formalized_in_lean": False,
        "finite_runs_certify_asymptotic_limit": False,
        "native_OPH_population_or_law_selected": False,
        "physical_clock_calibrated": False,
        "field_action_or_quantum_continuum_identified": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "The source-axis cube, golden product population, every-neighbor read law, "
            "distinguished-event count and model tick are supplied. The analytic "
            "limit concerns bounded boundary-null regions and contained timelike "
            "diamonds with converging tips. Finite execution replays source records, "
            "local writer/value custody and interventions; it is not an asymptotic "
            "experiment or a physical dimension fit. Spatial populations are nested; "
            "changing radii and ticks do not give induced-prefix histories."
        ),
    }


def _fresh_structural_verifier(path: Path, dependency: str | None = None,
                              dependency_paths: dict[str, Path] | None = None):
    """Compile current verifier/helper bytes; restore existing import bindings."""
    dependencies = dict(dependency_paths or {})
    if dependency:
        dependencies[dependency] = path.parent / (dependency + ".py")
    names = list(dependencies)
    name = "postdiction_" + path.parent.name + "_" + path.stem
    names.append(name)
    saved = {key: sys.modules.get(key) for key in names}
    original_path = sys.path[:]
    try:
        for key in names:
            source = dependencies.get(key, path)
            module = types.ModuleType(key)
            module.__file__ = str(source)
            sys.modules[key] = module
            exec(compile(source.read_bytes(), str(source), "exec"), module.__dict__)
        return module
    finally:
        sys.path[:] = original_path
        for key, previous in saved.items():
            if previous is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = previous


def _structural_packet(package: str, verifier_name: str, receipt_name: str,
                       *, entry: str = "verify", dependency: str | None = None,
                       dependency_paths: dict[str, Path] | None = None,
                       receipt_path: Path | None = None):
    directory = CODE / package
    verifier = _fresh_structural_verifier(directory / verifier_name, dependency, dependency_paths)
    path = receipt_path or directory / receipt_name
    raw = path.read_bytes()
    try:
        packet = verifier.load(path)
        summary = getattr(verifier, entry)(packet)
    except (ValueError, TypeError, KeyError) as exc:
        raise SystemExit(f"{package} structural certificate rejected: {exc}") from exc
    if path.read_bytes() != raw:
        raise SystemExit(package + " receipt changed during verification")
    return packet, summary, {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def _protected_population_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_population", "verify_population.py", "runtime/population_receipt.json",
        entry="check", receipt_path=receipt_path)
    if (summary.get("accepted") is not True or summary.get("source_population_produced") is not False
            or summary.get("persistent_memory_supplied") is not True
            or any(summary.get(k) != value for k, value in packet["summary"].items())):
        raise SystemExit("protected population verification projection mismatch")
    declarations = ("linearRecord_repair_defect", "protects_iff_edge_constant",
                    "connected_protected_record_factors_through_total",
                    "active_distinguishing_seam_changes_record")
    return {
        "receipt": "code/source_population/runtime/population_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "mathematical_replay": True, "observed_postdiction": False,
        "lean_receipts": _lean_receipt("SourceRecordProtection",
            declarations={"SourceRecordProtection": declarations}),
        "lean_declarations": {"SourceRecordProtection": list(declarations)},
        "scope_boundary": (
            "The initial golden population and immutable retained address register are supplied. "
            "Exact pair-mean operations move live port readbacks; the protected address survives "
            "only the declared write footprint. This does not produce a persistent spatial "
            "population, select a routing law or exclude other observer dynamics."),
    }


def _source_population_quadrature_control() -> dict[str, Any]:
    declarations = ("cell_local_first_moment", "partition_uniform_quadrature",
                    "partition_quadrature_tendsto", "tent_field_lipschitz",
                    "constant_detects_mass_error", "misplaced_sample_counterexample")
    source = "paper/tex_fragments/SOURCE_POPULATION_QUADRATURE.tex"
    labels = ["prop:source-population-cell-transport", "prop:source-population-local-kernel"]
    golden = ("residue_bijective", "fibonacci_error", "golden_floor", "golden_orbit_eq",
              "golden_grid_error", "orbit_injective", "tensor_assignment", "site_card",
              "axisCell_disjoint", "axisCell_cover", "cell_disjoint", "cell_cover",
              "cell_mass", "golden_quadrature", "golden_quadrature_tendsto",
              "cube_ae_closed", "golden_quadrature_closed", "sample_outside_assigned_cell",
              "one_cell_width_counterexample")
    counts = ("partition_count_tendsto", "golden_count_tendsto", "timeCell_cover",
              "timeCell_disjoint", "eventCell_assignment", "event_weighted_count_tendsto",
              "sum_eventCell_mass", "event_uniform_weight_error", "event_count_tendsto",
              "sourceDelta_pos", "sourceDelta_tendsto", "source_event_count_tendsto",
              "goldenPosition_injective", "eventPosition_injective", "goldenCount_univ",
              "eventCount_univ", "aligned_endpoint_layer_control")
    text = (REPO/source).read_text(encoding="utf-8")
    if any("\\label{" + label + "}" not in text for label in labels):
        raise SystemExit("source population quadrature theorem missing")
    return {
        "source": source, "labels": labels,
        "lean_receipts": _lean_receipt("SourcePopulationQuadrature", "GoldenSourceAssignment",
            "GoldenSourceCountLimit", declarations={"SourcePopulationQuadrature": declarations,
                "GoldenSourceAssignment": golden, "GoldenSourceCountLimit": counts}),
        "lean_declarations": {"SourcePopulationQuadrature": list(declarations),
                              "GoldenSourceAssignment": list(golden),
                              "GoldenSourceCountLimit": list(counts)},
        "formal_actual_cell_integral_and_quadrature_limit": True,
        "formal_golden_partition_geometry": True,
        "formal_causal_pair_count_limit": False,
        "golden_assignment_bound": "2*sqrt(3)*L/q",
        "golden_cell_mass": "L^3/q^3",
        "golden_lipschitz_integral_error": "L^3*K*(2*sqrt(3)*L/q)",
        "fixed_globally_lipschitz_detector_convergence": True,
        "fixed_measurable_null_frontier_indicator_convergence": True,
        "formal_actual_source_time_grid_count_limit": True,
        "source_time_step": "Delta_q=L/(c*sqrt(q))",
        "source_event_layers": "j=0,...,floor(T/Delta_q), including aligned endpoint",
        "normalized_source_event_count": "Delta_q*(L^3/q^3)*cardinality",
        "last_temporal_cell_error_upper": "L^3*Delta_q",
        "volume_convention": "dt d^3x restricted to [0,T) x [0,L)^3",
        "arbitrary_varying_detector_convergence": False,
        "moving_causal_interval_limit_formalized": False,
        "actual_causal_cone_null_frontier_identified_in_lean": False,
        "equal_cell_weights_are_lumped_action_weights": False,
        "tent_integrand_lipschitz_coefficient": "K",
        "analytic_radial_quadrature_bound": "60*K*H*(epsilon+H)^3/epsilon^5",
        "observed_postdiction": False,
        "scope_boundary": (
            "The general quadrature theorem takes measured cells, masses and assignments "
            "as hypotheses. The new golden specialization constructs the actual fractional-part "
            "population, permuted half-open partition, exact Lebesgue cell masses and "
            "2*sqrt(3)*L/q assignment bound from Fibonacci arithmetic. A sample need not lie "
            "in its assigned cell. The Bochner integral estimate and limit fix one globally "
            "Lipschitz detector and its constant. Actual golden spatial and time-product counts "
            "now converge on any fixed measurable set with null frontier by dominated convergence, "
            "without an assumed count or weak-measure limit. Inclusive full-layer weights differ "
            "from the clipped final time cell by at most L^3*Delta_q. Actual causal-cone null "
            "frontiers, moving generated intervals, pair-count limits and arbitrary varying "
            "detectors remain separate. Dual-cell action and equal-volume "
            "weights are distinct finite measures; no physical population law follows."),
    }


def _source_feedback_transport_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_feedback_transport", "verify_transport.py", "transport_receipt.json",
        dependency_paths={"verify_routing": CODE / "source_routing" / "verify_routing.py"},
        receipt_path=receipt_path)
    expected = {"episodes": 16, "events": 4928, "hops": 740}
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("feedback transport complete replay projection mismatch")
    declarations = ("receive_error_identity", "exact_receive", "unit_signal_gain",
                    "route_error_identity", "route_error_bound", "bounded_hop_error",
                    "route_uniform_bound", "sharp_uniform_error", "repeated_exact_hops",
                    "semantic_projection_sound", "semantic_projection_complete",
                    "exact_induced_logical_order")
    source = "paper/tex_fragments/SOURCE_FEEDBACK_TRANSPORT.tex"
    labels = ["thm:source-feedback-transport", "thm:source-feedback-error"]
    text = (REPO/source).read_text(encoding="utf-8")
    if any("\\label{" + label + "}" not in text for label in labels):
        raise SystemExit("feedback transport analytic theorem missing")
    return {
        "receipt": "code/source_feedback_transport/transport_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": {"source": source, "labels": labels},
        "lean_receipts": _lean_receipt("SourceFeedbackTransport",
            declarations={"SourceFeedbackTransport": declarations}),
        "lean_declarations": {"SourceFeedbackTransport": list(declarations)},
        "mathematical_replay": True,
        "semantic_read_from_order_checked_for_all_logical_pairs": True,
        "all_primitive_events_and_immutable_versions_retained": True,
        "reusable_classical_transport": True,
        "copy_reset_factor_two_readout_and_schedule_supplied": True,
        "source_selected_feedback_or_quantum_operations": False,
        "full_q13_q21_routed_histories_executed": False,
        "physical_causal_order_capacity_or_clock_identified": False,
        "observed_postdiction": False,
        "noise_contract": packet["noise_contract"],
        "baseline_costs": [ep["costs"] for ep in packet["episodes"] if ep["variant"] == "baseline"],
        "compiler_census_not_executed": packet["compiler_targets_not_executed"],
        "scope_boundary": (
            "Immutable local record export, protected zero reset and factor-two local readout "
            "are supplied classical operations. Fixed non-interleaved six-event hops restore "
            "scratch ports and preserve old versions under the declared interventions. Exact "
            "writer ancestry gives the induced semantic read-from order, not physical resource, "
            "scheduler or custody-chain order. The sharp E0+d*(E+Z+2*M+2*R+D) error bound "
            "assumes separately bounded additive errors; the exact receipt measures no hardware "
            "noise. The q13/q21 census is an unicast compiler cost calculation, not execution "
            "of those full routed histories or a universal multicast cost bound. No native "
            "feedback selection, physical address channel, quantum copying or clock follows."),
    }


def _source_scalar_sequential_instrument_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_instruments", "verify_sequential_instrument.py",
        "sequential_instrument_receipt.json", receipt_path=receipt_path)
    expected = {"verified": True, **packet["summary"], "mathematical_replay": True,
                "parent_events_replayed": 5888, "physical_outcomes": False,
                "branch_conditioned_clock": False}
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("sequential instrument requires full fresh parent replay")
    if (summary["readout_times"] != 21 or summary["prior_pair_brackets"] != 210
            or len(summary["uniformly_resolved_steps"]) != 15):
        raise SystemExit("sequential instrument full readout census changed")
    source = "paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex"
    text = (REPO/source).read_text(encoding="utf-8")
    labels = ["thm:source-scalar-sequential-instrument", "prop:source-scalar-sequential-comparison"]
    if any("\\label{" + label + "}" not in text for label in labels):
        raise SystemExit("sequential instrument analytic theorem missing")
    return {
        "receipt": "code/source_scalar_instruments/sequential_instrument_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": {"source": source, "labels": labels},
        "mathematical_replay": True,
        "same_action_vacuum_preparation_and_recovered_clock_intervals": True,
        "reported_probabilities_are_unconditional_sequential_marginals": True,
        "reference_has_same_repeated_instrument": True,
        "backaction_and_all_prior_readout_pairs_included": True,
        "regional_ghz_gadget_on_scalar_action_graph": True,
        "regional_gadget": packet["regional_gadget"],
        "continuous_comparison_rows": packet["rows"],
        "energy_injected_per_readout_Qphi": packet["energy_injected_per_readout_Qphi"],
        "total_readout_injected_work_Qphi": packet["total_readout_injected_work_Qphi"],
        "instrument_and_preparation_laws_supplied": True,
        "sampled_quantum_outcome_histories": False,
        "branch_conditioned_clock_reconstructed": False,
        "pointer_action_graph_identified_with_W12_seams": False,
        "finite_duration_or_noisy_quantum_gates_certified": False,
        "instrument_theorems_formalized_in_lean": False,
        "physical_clock_or_observed_outcomes": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "A supplied centered binary quantum instrument defines a joint classical record "
            "algebra. All 21 unconditional marginals include backaction from earlier readouts "
            "and are compared with continuous evolution under the same repeated instrument. "
            "The exact commutator-cosine factors and inherited clock windows preserve 15 "
            "resolved slots. No sampled outcomes, independent-trial assumption or branch-clock "
            "reconstruction is made. A 31-edge GHZ pointer construction uses the scalar action "
            "graph; it is not W12 routing. Controlled kicks, pointer preparation/reset, "
            "Born readout, vacuum preparation and clock calibration are supplied. Ideal "
            "instantaneous gate counts do not enclose duration or noise, select quantum "
            "operations from OPH, or establish physical outcomes or an interacting continuum."),
    }


def _source_scalar_finite_instrument_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_finite_instrument", "verify_finite_instrument.py",
        "finite_instrument_receipt.json", receipt_path=receipt_path)
    expected = {"verified": True, **packet["summary"], "mathematical_replay": True,
                "parent_events_replayed": 5888, "physical_outcomes": False,
                "native_quantum_controls": False, "forced_oscillator_identities_verified": True}
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("finite instrument requires forced-oscillator proof and full fresh parent replay")
    if (summary["readouts"] != 21 or summary["prior_pairs"] != 210
            or summary["retained_operations"] != 2709
            or len(summary["resolved_steps"]) != 15
            or summary["all_parent_resolved_steps_survive"] is not True):
        raise SystemExit("finite instrument complete response and control census changed")
    source = "paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex"
    labels = ["prop:source-scalar-finite-pulse", "prop:source-scalar-finite-noise"]
    text = (REPO/source).read_text(encoding="utf-8")
    if any("\\label{"+label+"}" not in text for label in labels):
        raise SystemExit("finite instrument analytic theorem missing")
    return {
        "receipt": "code/source_scalar_finite_instrument/finite_instrument_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": {"source": source, "labels": labels},
        "mathematical_replay": True,
        "continuous_field_evolution_during_supplied_pointer_pulse": True,
        "finite_duration_and_declared_channel_noise_comparison": True,
        "baseline_and_intervention_errors_both_retained": True,
        "parameters": packet["parameters"], "comparison_rows": packet["rows"],
        "resource_contract": packet["resource_contract"],
        "physical_clock_or_observed_outcomes": False,
        "native_quantum_controls": False,
        "noisy_energy_bound": False,
        "branch_conditioned_clock": False,
        "original_equal_time_detector_region_preserved": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "A supplied finite rectangular field-pointer Hamiltonian has an exact forced-oscillator "
            "reduction on the ideal GHZ code, with a sinc-averaged smear. Complete repeated backaction, "
            "control windows and paired baseline/intervention channel errors preserve fifteen responses. "
            "Initial quantum states, control laws, entangling links, Born readout, precision and "
            "conditional reference-center times remain supplied. The scalar phase is code-restricted; "
            "generic noise need not preserve energy or first moments. Effective smear can spread beyond "
            "the local coupling region. No sampled outcomes, native quantum implementation, physical "
            "clock, noisy branch-clock reconstruction or physical time-slice attachment is established."),
    }


SOURCE_OPERATOR_JOIN_DECLARATIONS = (
    "operatorJoin88_injective", "operatorJoin247_injective", "operatorJoin86_injective",
    "operatorJoin_hinge", "operatorJoin_intersection", "operatorJoin_generate_top",
    "operatorJoin_spectator_locality", "operatorJoin_hinge_noncommuting",
    "operatorJoin88_retraction", "operatorJoin247_retraction",
    "operatorJoin88_counted_readout", "operatorJoin247_counted_readout",
    "operatorJoin88_evolve", "operatorJoin247_evolve",
    "operatorJoinSource88_generates", "operatorJoinSource247_generates",
    "operatorJoinAdmitted_generates", "operatorJoinAdmitted_evolved",
)
# This exact source was independently kernel checked.  A changed proof needs
# review before this pin changes; source-name checks alone do not verify Lean.
SOURCE_OPERATOR_JOIN_LEAN_SHA256 = "6325cd23dff82fe5ab8e37aa0d17df828a643f5c9aa1939fd8d3834bee4457fd"


def _source_operator_join_control(receipt_path: Path | None = None) -> dict[str, Any]:
    """Keep finite Python replay separate from the universal Lean theorems."""
    lean_path = LEAN_RECEIPTS["TripleCarrierOperatorJoin"]
    if not lean_path.is_file():
        raise SystemExit("operator join Lean source missing")
    lean_bytes = lean_path.read_bytes()
    if hashlib.sha256(lean_bytes).hexdigest() != SOURCE_OPERATOR_JOIN_LEAN_SHA256:
        raise SystemExit("operator join audited Lean source changed")
    lean_refs = _lean_receipt("TripleCarrierOperatorJoin", declarations={
        "TripleCarrierOperatorJoin": SOURCE_OPERATOR_JOIN_DECLARATIONS})
    packet, summary, pin = _structural_packet(
        "source_operator_join", "verify.py", "operator_join_packet.json",
        receipt_path=receipt_path)
    expected = {
        "verdict": "verified_source_operator_join_packet",
        "source_rows": 32, "actual_transitions": 31,
        "transition_carrier_points": 73346, "chronological_intervals": 496,
        "sparse_matrix_unit_witnesses": 18,
        "lean_kernel_replay_performed": False,
        "universal_operator_theorems_verified_by_python": False,
        "physical_adequacy": False, "observed_quantum_outcomes": False,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("operator join replay census or Python/Lean proof boundary changed")
    if lean_path.read_bytes() != lean_bytes:
        raise SystemExit("operator join Lean source changed during replay")
    counts = {tuple(row[:3]): row[3] for row in packet["state"]["triple_counts"]}
    permutations = packet["transitions"][0]["permutations"]
    transported = {tuple(permutations[i][q[i]] for i in range(3)): count
                   for q, count in counts.items()}
    changed = sorted(q for q in set(counts) | set(transported)
                     if counts.get(q, 0) != transported.get(q, 0))
    if not changed:
        raise SystemExit("operator join empirical-state nonstationarity control missing")
    witness = changed[0]
    return {
        "receipt": "code/source_operator_join/operator_join_packet.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "source_parent_pins": packet["parents"],
        "lean_receipts": lean_refs,
        "lean_declarations": {"TripleCarrierOperatorJoin": list(SOURCE_OPERATOR_JOIN_DECLARATIONS)},
        "lean_source_evidence": {
            "source": lean_refs[0], "bytes": len(lean_bytes),
            "sha256": hashlib.sha256(lean_bytes).hexdigest(),
            "kernel_replay_performed_by_ledger": False,
            "universal_proof_requires_separate_Lean_check": True,
        },
        "carrier": packet["carrier"],
        "finite_source_and_permutation_replay": True,
        "python_sparse_witnesses_prove_universal_generation": False,
        "observable_inclusions_are_unital_star_homomorphisms": True,
        "full_hinge_and_joint_generation_are_Lean_theorems": True,
        "state_marginals_use_unnormalized_partial_trace": True,
        "observable_retractions_use_normalized_partial_trace": True,
        "partial_trace_is_multiplicative": False,
        "counted_state_identified_with_uniform_tower_state": False,
        "counted_state_invariant_under_transport": False,
        "nonstationary_state_witness": {
            "transition": 0, "triple_label": list(witness),
            "denominator": packet["state"]["denominator"],
            "before_count": counts.get(witness, 0),
            "after_count": transported.get(witness, 0),
        },
        "coherent_operator_lifts_and_path_permutations_supplied": True,
        "source_selected_quantum_controls": False,
        "tower_or_refinement_morphism": False,
        "physical_time_slice_or_clock": False,
        "same_carrier_as_scalar64_instrument": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "The admitted complex operator lifts of counted transitions and field projectors "
            "generate the two pair algebras. Their identity-spectator inclusions have a full M13 "
            "intersection and generate M2366; these universal statements are proved in the pinned "
            "Lean module, not by the eighteen Python matrix-unit controls. Python independently "
            "replays all source rows, nonwrapping permutations and finite interval compositions. "
            "The shared hinge is noncommutative; only the exclusive spectator factors commute. "
            "The empirical state restricts to the counted pair states, differs from the supplied "
            "uniform tower states, and is not invariant under the row-dependent transports. "
            "Tensor assembly, coherent operator admission and path-read permutation gates remain "
            "supplied postprocessing. No tower/refinement morphism, physical region or clock, "
            "quantum outcome, physical time-slice or common scalar-instrument realization follows."),
    }


def _triple_observer_carrier_row() -> dict[str, Any]:
    declarations = {
        "TripleCarrierJoin": (
            "triplePath_projects_88", "triplePath_projects_247",
            "tripleCorrelationState_marginal_88", "tripleCorrelationState_marginal_247",
            "tripleCheckpoint_recovers_anchored_88", "tripleCheckpoint_recovers_anchored_247",
            "tripleMarginals_not_jointly_injective", "tripleCarrier_join_receipt",
        ),
        "TripleCarrierOperatorJoin": SOURCE_OPERATOR_JOIN_DECLARATIONS,
    }
    return {
        "id": "triple_observer_carrier_coupling",
        "statement": (
            "One explicit (86,88,247) carrier recovers both committed 32-row paths and counted "
            "pair states, both checkpoint partitions, and both marginal evolutions over exactly "
            "the 31 adjacent source transitions. The pair operator inclusions have full M13 "
            "intersection, generate M2366, and intertwine those transitions; admitted operator "
            "lifts of source generators retain full generation under the declared transport"),
        "observed_counterpart": "common coupling of the two committed correlation packets",
        "match": "exact conditional state/path/checkpoint coupling and finite operator join",
        "full_operator_join_control": _source_operator_join_control(),
        "lean_declarations": {key: list(value) for key, value in declarations.items()},
        "lean_receipts": _lean_receipt(*declarations, declarations=declarations),
        "artifact_refs": [
            "code/source_operator_join/operator_join_packet.json",
            "code/source_operator_join/verify.py",
        ],
        "hypothesis_boundary": (
            "The full step alignment, tensor slots and coherent operator lifts are declared, "
            "and pair marginals do not select the triple coupling uniquely. Observable inclusions "
            "are unital star homomorphisms; the opposite-direction partial traces are not "
            "multiplicative. The empirical state is not identified with uniform tower states "
            "or assumed transport invariant. Universal operator results belong to the Lean "
            "proof; finite Python controls are not its substitute. No source-selected quantum "
            "gate, regional-net/tower/refinement morphism, physical region, Cauchy/time-slice "
            "interpretation, physical clock or scalar64 instrument identification follows"),
        "paper_ref": "Consensus paper, triple-carrier coupling; observers paper, finite operator join",
    }


def _source_scalar_regional_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_regional", "verify_regional_time_slice.py", "regional_time_slice_receipt.json",
        dependency_paths={"verify_source_scalar_execution":
                          CODE / "source_scalar_execution/verify_source_scalar_execution.py"},
        receipt_path=receipt_path)
    expected = {"verified": True, "region_count": 8, "missing_row_counterexamples": 62,
                "regions": [{"region": r["name"], "sites": len(r["region_sites"]),
                    "coordinate_collar_sites": len(r["coordinate_collar_sites"]),
                    "minimal_linear_collar": r["minimal_linear_collar_dimension"]}
                    for r in packet["regions"]]}
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("regional scalar independent rank projection mismatch")
    declarations = ("reconstruct_velocity", "hidden_exterior_perturbation",
                    "determines_velocity_iff_kernel", "force_readout_suffices",
                    "no_collar_iff_zero_coupling")
    source = "paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex"
    label = "thm:source-scalar-regional-collar"
    if "\\label{" + label + "}" not in (REPO/source).read_text(encoding="utf-8"):
        raise SystemExit("regional scalar analytic theorem missing")
    return {
        "receipt": "code/source_scalar_regional/regional_time_slice_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": {"source": source, "label": label},
        "lean_receipts": _lean_receipt("ScalarRegionalTimeSlice",
            declarations={"ScalarRegionalTimeSlice": declarations}),
        "lean_declarations": {"ScalarRegionalTimeSlice": list(declarations)},
        "mathematical_replay": True, "exact_minimal_exterior_readout_ranks": True,
        "kernel_reconstruction_criterion_formalized": True,
        "weyl_algebra_and_rank_consequences_formalized": False,
        "all_real_weyl_parameters_are_supplied": True,
        "classical_record_trace_provides_weyl_access": False,
        "physical_regional_time_slice": False, "observed_quantum_outcomes": False,
        "source_action_population_or_physical_clock_selected": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "For the complete supplied finite scalar action, original regional canonical "
            "Weyl access follows from two field layers and exterior linear readouts C exactly "
            "when ker C is contained in ker A_R,E. The sharp minimum is rank A_R,E. "
            "All eight regions and 62 missing-readout controls are independently reconstructed. "
            "The two-layer fields have nonzero cross-time canonical brackets and define a "
            "common-unitary sheared net; they need not generate the original site algebra "
            "without exterior access. Rank/Weyl consequences are analytic; Lean proves the "
            "linear reconstruction criterion. Arbitrary real Weyl parameters, regions, action, "
            "quantization and physical attachment remain supplied. No classical record is "
            "promoted to an operational quantum outcome or physical Cauchy interface."),
    }


def _source_seam_routing_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_routing", "verify_routing.py", "runtime/path_tomography_receipt.json",
        receipt_path=receipt_path)
    if type(summary) is not list or len(summary) != 12:
        raise SystemExit("source routing requires all twelve independently replayed histories")
    expected = []
    for ep in packet["episodes"]:
        cost = ep["cost"]
        expected.append({"path_edges": cost["path_edges"], "variant": ep["variant"],
            "seam_means": cost["seam_mean_operations"], "events": cost["all_events"],
            "source_noise_gain": ep["conditioning"]["source_error_per_uniform_sample_error"],
            "provenance_edges": sum(len({r["writer"] for r in e["reads"]}) for e in ep["events"])})
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("source routing independent execution projection mismatch")
    if sum(r["events"] for r in summary) != 1449:
        raise SystemExit("source routing complete event census changed")
    for row in summary:
        previous, gain = 1, 1
        for _ in range(row["path_edges"]):
            previous, gain = gain, 2 * gain + previous
        if Fraction(row["source_noise_gain"]) != gain:
            raise SystemExit("source routing exact conditioning recurrence mismatch")
    declarations = ("invert_endpoint", "endpoint_difference", "invert_difference", "decode_encode",
                    "encode_observations_injective", "encode_lengths", "twice_meanCost", "canonical_hop_decode")
    source = "paper/tex_fragments/SOURCE_SEAM_PATH_TOMOGRAPHY.tex"
    label = "thm:source-seam-path-tomography"
    conditioning_label = "thm:source-seam-exactconditioning"
    gate_error_label = "cor:source-seam-gateerror"
    proof_text = (REPO/source).read_text(encoding="utf-8")
    if any("\\label{" + item + "}" not in proof_text
           for item in (label, conditioning_label, gate_error_label)):
        raise SystemExit("source routing analytic theorem missing")
    return {
        "receipt": "code/source_routing/runtime/path_tomography_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": {"source": source, "label": label},
        "conditioning_analytic_proof": {"source": source, "label": conditioning_label},
        "uniform_sample_error_gain": {
            "norm": "infinity", "observation_matrix_norm": 1,
            "recurrence": "g_0=1, g_1=3, g_d=2*g_(d-1)+g_(d-2)",
            "closed_form": "((1+sqrt(2))^(d+1)+(1-sqrt(2))^(d+1))/2",
            "initial_vector_error_bound": "g_d*epsilon",
            "sharp_for_independent_bounded_sample_errors": True,
            "all_depths_theorem_formalized": False,
        },
        "supplied_gate_sample_and_decoder_error_bound": {
            "analytic_proof": {"source": source, "label": gate_error_label},
            "mean_operations": "T=d*(d+1)/2",
            "initial_vector_error_bound": "g_d*(T*eta+sigma)+delta",
            "eta": "additive infinity-norm error per complete pair-mean map",
            "sigma": "sample and protected-memory retention error",
            "delta": "separately bounded decoding arithmetic error",
            "initial_state": "actual starting port vector",
            "sharp_gate_error_bound": False, "all_error_budgets_are_supplied": True,
            "attained_physical_precision": False, "formalized_in_lean": False,
        },
        "lean_receipts": _lean_receipt("SourceSeamPathTomography",
            declarations={"SourceSeamPathTomography": declarations}),
        "lean_declarations": {"SourceSeamPathTomography": list(declarations)},
        "mathematical_replay": True, "canonical_mean_only_destructive_readback": True,
        "destination_decoder_reads_remote_baselines": False,
        "all_transport_and_calibration_events_retained": True,
        "protected_destination_memory_supplied": True,
        "receipt_noise_bound_includes_dynamical_or_arithmetic_error": False,
        "full_metric_neighbor_order_refinement": False,
        "source_population_produced": False,
        "physical_capacity_or_clock_identified": False, "observed_postdiction": False,
        "baseline_costs": [ep["cost"] for ep in packet["episodes"] if ep["variant"] == "baseline"],
        "scope_boundary": (
            "Recursive suffix calibration uses only canonical pair means and destination-local "
            "protected samples. It reconstructs original path values and the changed suffix, "
            "with d(d+1)/2 means, d+1 samples and all 1449 events retained across twelve runs. "
            "The finite W12 support, isolated schedule, preparations and protected memory are "
            "declared. Local builder hashes identify the captured source bytes; the simulator "
            "base checkout commit does not claim those files were committed there. Register "
            "payload bits exclude metadata, temporary arithmetic and representation overhead. "
            "The analytic all-depth infinity-norm condition number is g_d, with g_0=1, "
            "g_1=3 and g_d=2*g_(d-1)+g_(d-2); finite inverse matrices independently verify "
            "the recorded cases. Its sharp sample-error gain assumes exact dynamics and "
            "arithmetic. A separate nonsharp analytic bound g_d*(T*eta+sigma)+delta "
            "allows supplied per-gate, sample/retention and decoding errors, relative to the "
            "actual initial ports. These all-depth conditioning/error statements are not "
            "formalized in Lean and assert no attained hardware precision. No repeated "
            "immutable-version read service, full metric-neighbour refinement, physical precision, "
            "population production or count-clock attachment follows."),
    }


def _common_source_scalar_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_packet", "verify_source_common_scalar.py", "source_common_scalar_receipt.json",
        dependency="source_scalar_intervals", receipt_path=receipt_path)
    last = packet["levels"][-1]
    expected = {
        "verified": True, "levels": len(packet["levels"]), "resolved_q": last["q"],
        "full_tensor_oscillators": last["dynamic_oscillators"],
        "largest_tridiagonal": last["tensor_algorithm_largest_matrix"],
        "uniform_window": packet["continuum"]["window"],
        "error_upper": last["graph_vs_compact_continuum_error_upper"][1],
        "resolved_signal_lower": last["resolved_graph_signal_lower"],
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("common scalar verification projection mismatch")
    if not 0 < Fraction(summary["error_upper"]) < Fraction(summary["resolved_signal_lower"]):
        raise SystemExit("common scalar response does not dominate its certified error")
    source = "paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex"
    labels = ["prop:common-scalar-full-leakage", "prop:common-scalar-continuum-limit",
              "thm:source-common-scalar-response"]
    text = (REPO/source).read_text(encoding="utf-8")
    if any("\\label{" + label + "}" not in text for label in labels):
        raise SystemExit("common scalar analytic theorem missing")
    return {
        "receipt": "code/source_scalar_packet/source_common_scalar_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "mathematical_replay": True, "observed_postdiction": False,
        "analytic_proof": {"source": source, "labels": labels},
        "analytic_free_scalar_detector_limit": True, "full_continuum_theorem_formalized_in_lean": False,
        "source_action_population_or_physical_clock_selected": False,
        "field_history_joined_to_causal_count_clock": False,
        "interacting_quantum_continuum": False,
        "scope_boundary": (
            "Protected golden source addresses, tensor mass-lumped action, Dirichlet cube, "
            "mass one, quantization and model time are declared. Full lattice leakage and "
            "sampled compact tails are bounded. The analytic limit concerns fixed-time "
            "smeared classical and coherent Weyl readouts for bounded Riemann-integrable "
            "smearings. The resolved q233 field edges fit the causal radius; that accuracy "
            "packet does not execute their field reads or attach a count clock. A separate "
            "q5 execution uses different preparation data. This action differs from the radial kernel."),
    }


def _source_scalar_execution_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_execution", "verify_source_scalar_execution.py",
        "source_scalar_execution_receipt.json", receipt_path=receipt_path)
    traces = packet["traces"]
    expected = {
        "verdict": "PASS", "mutable_sites": packet["model"]["mutable_sites"],
        "executed_times": packet["model"]["steps"], "traces_replayed": len(traces),
        "events_replayed": sum(len(t["events"]) for t in traces.values()),
        "dynamic_field_reads_replayed": sum(t["dynamic_field_reads"] for t in traces.values()),
        "full_mass_norm_error_final": packet["diagnostics"]["comparisons"][-1]["full_mass_norm_discretization_error"],
        "same_finite_operator_all_modes": True, "old_count_clock_applied": False,
        "quantum_covariance_or_probability_queried": False,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("scalar execution verification projection mismatch")
    return {
        "receipt": "code/source_scalar_execution/source_scalar_execution_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "mathematical_replay": True, "all_field_arithmetic_exact": True,
        "analytic_proof": "paper/tex_fragments/SOURCE_SCALAR_EXECUTION.tex",
        "observed_postdiction": False, "spatial_continuum_error_certified": False,
        "field_history_joined_to_causal_count_clock": False,
        "source_action_population_or_physical_clock_selected": False,
        "quantum_covariance_or_probability_queried": False,
        "scope_boundary": (
            "Protected q5 golden addresses, tensor action, Dirichlet boundary, polynomial "
            "preparation and model step are supplied. All substeps, exact values, immutable "
            "versions and actual writers are replayed under two schedules and two preparations. "
            "The error compares every mode of the same finite continuous-time action at the "
            "21 executed times. Its state-independent canonical commutator does not replace "
            "the original vacuum covariance. The q233 continuum accuracy and earlier "
            "point-event count-clock theorem are not transferred to this history."),
    }


def _source_scalar_quantum_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_quantum", "verify_quantum_probability.py",
        "quantum_probability_receipt.json", receipt_path=receipt_path)
    if (summary.get("verdict") != "PASS"
            or summary.get("parent_full_mathematical_replay") is not True
            or summary.get("parent_events_replayed") != 5888):
        raise SystemExit("scalar quantum projection requires full authenticated parent replay")
    expected = {
        "verdict": "PASS", "times": len(packet["rows"]),
        "resolved_steps": packet["summary"]["resolved_steps"],
        "uniform_added_variance_upper": packet["covariance"]["increase_upper"],
        "step16": packet["rows"][15],
        "quantum_probability_scope": "original finite-action vacuum and coherent momentum preparation; bounded sine effect",
        "spatial_continuum_error_transfer": False, "observed_outcomes": False,
        "parent_full_mathematical_replay": True, "parent_events_replayed": 5888,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("scalar quantum verification projection mismatch")
    row = summary["step16"]
    if not 0 < Fraction(row["probability_error_upper"]) < Fraction(row["split_response_abs_lower"]):
        raise SystemExit("scalar quantum response does not dominate its certified error")
    return {
        "receipt": "code/source_scalar_quantum/quantum_probability_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "mathematical_replay": True, "quantum_covariance_and_bounded_effect_probability": True,
        "analytic_proof": "paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex",
        "original_finite_action_vacuum_retained": True,
        "observed_postdiction": False, "observed_quantum_outcomes": False,
        "spatial_continuum_error_certified": False,
        "intermediate_time_error_enclosure": False,
        "field_history_joined_to_causal_count_clock": False,
        "source_action_population_or_physical_clock_selected": False,
        "scope_boundary": (
            "This separate extension uses all 64 q5 oscillators, the original finite-action "
            "vacuum and its coherent momentum displacement. The same symmetric quantum "
            "gates propagate both means and covariance; the squeezed baseline has sine-effect "
            "probability one half. Exact rational bounds compare the bounded sine-effect "
            "probability with continuous evolution of the same finite action at 21 stored "
            "times. The authenticated classical registers carry means, not an attested "
            "quantum state or outcome stream. The immutable earlier commutator receipt "
            "retains its own negative probability flag. No q233 continuum accuracy, "
            "intermediate-time bound, causal count clock or physical units are transferred."),
    }


def _source_scalar_clock_quantum_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_clock_quantum", "verify_clock_quantum.py",
        "clock_quantum_receipt.json", receipt_path=receipt_path)
    if (summary.get("verdict") != "PASS"
            or summary.get("full_clock_replay") is not True
            or summary.get("full_original_quantum_replay") is not True
            or type(summary.get("parent_events_replayed")) is not int
            or summary["parent_events_replayed"] != 5888):
        raise SystemExit("scalar clock/quantum projection requires both full parent proofs")
    expected = {
        "verdict": "PASS", **packet["summary"], "step16": packet["rows"][15],
        "clock_only_uncertainty": True, "physical_clock_identified": False,
        "parent_events_replayed": 5888, "full_clock_replay": True,
        "full_original_quantum_replay": True,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("scalar clock/quantum verification projection mismatch")
    row = summary["step16"]
    if not 0 < Fraction(row["total_probability_error_upper"]) < Fraction(row["split_response_abs_lower"]):
        raise SystemExit("scalar clock/quantum response must dominate the combined error")
    return {
        "receipt": "code/source_scalar_clock_quantum/clock_quantum_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "clock_receipt": packet["dependencies"]["reconstructed_clock"],
        "analytic_proofs": ["paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex",
                            "paper/tex_fragments/SOURCE_SCALAR_CLOCK_QUANTUM.tex"],
        "lean_receipts": ["Lean/Screen/ActionTimeGram.lean", "Lean/Screen/SourceActionTime.lean"],
        "mathematical_replay": True, "original_finite_action_vacuum_retained": True,
        "reference_probability_bound_for_every_time_in_each_interval": True,
        "clock_only_uncertainty": True,
        "hidden_common_stationary_history_required": True,
        "observed_postdiction": False, "observed_quantum_outcomes": False,
        "record_noise_certifies_quantum_state_error": False,
        "arbitrary_noisy_history_existence_proved": False,
        "continuous_split_trajectory_enclosed": False,
        "spatial_continuum_error_certified": False,
        "field_history_joined_to_causal_count_clock": False,
        "source_action_population_or_physical_clock_selected": False,
        "scope_boundary": (
            "Consumed-writer graph reconstruction supplies complete q5 field configurations. "
            "A declared trapezoidal action identifies one positive duration per whole "
            "configuration interval, with exact Gram checks and conditional stability. "
            "A common hidden stationary history and supplied M-norm record-error budget "
            "give elapsed-time intervals, not a proof that arbitrary noisy data admit such "
            "a history. The continuous original-vacuum detector is Lipschitz in time. "
            "Its comparison with the executed split detector holds at every time in each "
            "interval, with the Hamiltonian, vacuum, preparation and detector held fixed. "
            "All 21 rows are retained. No preparation/operator error, continuous split "
            "trajectory, regional time-slice, q233 spatial accuracy, count-volume clock "
            "or physical time is transported. Historical parent receipts retain their "
            "original scopes and flags."),
    }


def _source_scalar_time_refinement_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "source_scalar_time_refinement", "verify_source_scalar_time_refinement.py",
        "source_scalar_time_refinement_receipt.json", receipt_path=receipt_path)
    if (summary.get("verdict") != "PASS"
            or summary.get("full_spatial_parent_replayed") is not True
            or summary.get("observer_event_log_executed") is not False):
        raise SystemExit("scalar time refinement requires a full spatial proof, not an execution claim")
    last = packet["levels"][-1]
    expected = {"verdict": "PASS", "full_spatial_parent_replayed": True,
                "levels": len(packet["levels"]), "q": last["q"],
                "time_step": last["time_step"], "observer_event_log_executed": False}
    for key in ("temporal_probability_error_upper", "same_grid_total_error_upper",
                "split_response_lower", "nearby_continuum_error_upper"):
        expected[key] = last[key]
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("scalar time refinement verification projection mismatch")
    if not 0 < Fraction(last["nearby_continuum_error_upper"]) < Fraction(last["split_response_lower"]):
        raise SystemExit("scalar time refinement response must dominate all-window error")
    return {
        "receipt": "code/source_scalar_time_refinement/source_scalar_time_refinement_receipt.json",
        "receipt_pin": pin, "independent_verifier_result": summary,
        "analytic_proof": "paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex",
        "same_compact_spatial_parent": packet["parent"],
        "mathematical_replay": True, "joint_scalar_space_time_detector_limit": True,
        "all_modes_and_original_vacuum_covariance_retained": True,
        "asymptotic_prearrival_detector_convergence": True,
        "finite_q_prearrival_zero_response_certified": False,
        "nearest_in_window_update_readout": True,
        "complete_observer_event_log_executed": False,
        "q5_clock_error_or_preparation_transferred": False,
        "source_action_population_or_physical_clock_selected": False,
        "raw_ancestry_or_count_volume_identified": False,
        "full_fock_state_norm_convergence": False,
        "interacting_quantum_continuum": False, "observed_postdiction": False,
        "scope_boundary": (
            "The same protected golden population, tensor action, original vacuum, "
            "coherent momentum preparation and compact detector are supplied. Full-spectrum "
            "split errors and the existing spatial detector theorem give a joint "
            "observable limit for declared stable time refinements. The q233 certificate "
            "bounds mathematical local-update trajectories; no full event history is "
            "executed. Nearest in-window sampling extends the reference comparison to "
            "the whole target window. The pre-arrival statement is asymptotic and uses "
            "finite propagation of the limiting Dirichlet Klein-Gordon equation; it "
            "is neither finite-q zero response nor an identification of microscopic "
            "read ancestry. The separate q5 clock, field preparation and noise budget "
            "are not transported."),
    }


def _pauli_structure_rows() -> list[dict[str, Any]]:
    """Keep conditional exclusion and its countercontrols out of empirical counts."""
    specifications = (
        ("spin_exchange", "spin_exchange_selection_boundary",
         "SPIN_EXCHANGE_SELECTION_BOUNDARY.tex",
         "A supplied spin lift and conserved single-occupation current histories do not "
         "select exchange statistics. A separate same-channel two-particle comparison "
         "has exact coincidence probabilities 1, 49/625 and 337/625 for the declared "
         "fermionic, bosonic and distinguishable preparations.",
         "conditional statistics-identification boundary and exact finite comparison",
         "The same full internal state, preparation, ideal detector and error budgets "
         "are supplied. There are no native or laboratory outcomes; operator Gauss "
         "of the one-particle parent does not transfer to this new preparation."),
        ("pauli_stability", "conditional_pauli_scalar_q_selection",
         "PAULI_STABILITY_SELECTION.tex",
         "A common real scalar-q creation law and exact quadratic occupation/phase "
         "law force q in {1,-1}. Finite full representation dimension, or a stable "
         "fixed negative-energy ladder, excludes q=1 and yields Pauli exclusion "
         "and full canonical anticommutation by polarization.",
         "analytic conditional exclusion with independently checked finite controls",
         "The quantization class, unit phase law and finite full representation or "
         "stable quantum ladder are declared. The source response Hilbert space is "
         "not identified with the full matter representation. The theorem is spin "
         "independent; truncated bosons and other classes remain outside its premise."),
        ("pauli_source_selection", "conditional_pauli_operation_selection",
         "PAULI_SOURCE_SELECTION_BOUNDARY.tex",
         "On a finite positive Hilbert space, an unscaled contraction creation field "
         "with its exact quadratic unit-phase law squares to zero. Linearity gives "
         "creation anticommutation without a scalar-q premise. An exhaustive adjoint "
         "add/remove instrument additionally yields full CAR. A covariant symmetric "
         "Fock space capped at total occupation two preserves the supplied one-particle "
         "source dynamics and exact phase law but permits double occupation. Making "
         "its fields into Kraus contractions rescales the quadratic phase law.",
         "analytic conditional operation theorem and exact source-bound countercontrols",
         "Under the exact phase law contraction is equivalent to exclusion. "
         "The field/Kraus identification and exact unit normalization are supplied, "
         "not source-derived. General normalized instruments can have an idle outcome "
         "and do not force these hypotheses. The finite countercontrol is not a full "
         "A1-A3 model or impossibility theorem. Physical Spin, native operation and "
         "matter attachment remain unproved."),
    )
    rows = []
    for package, identifier, fragment, statement, match, boundary in specifications:
        packet, summary, pin = _structural_packet(package, "verify.py", "receipt.json")
        if summary.get("verified") is not True:
            raise SystemExit(package + " independent verification failed")
        scope = packet["scope"]
        if scope.get("physical_spin_statistics_theorem") is not False:
            raise SystemExit(package + " must retain its physical spin-statistics boundary")
        if package == "pauli_stability" and summary.get("source_premises_derived") is not False:
            raise SystemExit("Pauli selector cannot promote declared source premises")
        if package == "pauli_source_selection" and summary.get("source_premises_derived") is not False:
            raise SystemExit("operation selector cannot promote its field/Kraus identification")
        rows.append({
            "id": identifier, "statement": statement,
            "observed_counterpart": "Pauli exclusion (target only; no physical comparison)",
            "match": match, "observed_postdiction": False,
            "source_selected_statistics": False, "physical_spin_statistics_theorem": False,
            "artifact_refs": ["paper/tex_fragments/" + fragment,
                              "code/" + package + "/receipt.json",
                              "code/" + package + "/verify.py"],
            "receipt_pin": pin, "independent_verifier_result": summary,
            "hypothesis_boundary": boundary,
            "proof_kind": "analytic_operator_proof_with_finite_executable_controls",
            "paper_ref": "Standard Model gauge paper, exchange and exclusion sections",
        })
    return rows


def _fermion_source_current_control() -> dict[str, Any]:
    controls = {}
    for key, verifier, receipt in (
        ("mean_field", "verify_current.py", "current_receipt.json"),
        ("quantum_link", "verify_quantum_link.py", "quantum_link_receipt.json"),
        ("coupled", "verify_coupled.py", "coupled_receipt.json"),
    ):
        dependencies = ({"verify_current": CODE / "sm_fermion_current/verify_current.py"}
                        if key in {"quantum_link", "coupled"} else None)
        _, summary, pin = _structural_packet(
            "sm_fermion_current", verifier, receipt, dependency_paths=dependencies)
        if summary.get("verified") is not True or summary.get("physical_current_attached") is not False:
            raise SystemExit("fermion current must retain exact replay and physical boundary")
        if summary.get("operator_Gauss") is not (key == "quantum_link"):
            raise SystemExit("mean Gauss and dressed quantum Gauss must remain distinct")
        replay_key = "independent_quantum_state_replay" if key == "quantum_link" else "exact_rational_replay"
        if summary.get(replay_key) is not True:
            raise SystemExit("fermion current requires full independent mathematical replay")
        if key == "quantum_link" and summary.get("exact_continuous_quantum_edge_solution") is not True:
            raise SystemExit("quantum link requires its independently checked continuous solution")
        if key == "quantum_link" and summary.get("gauge_constraint_group") != "U1_hypercharge":
            raise SystemExit("quantum link certifies only the abelian Gauss constraint")
        if key == "coupled":
            for flag in ("same_declared_hybrid_action", "electric_to_link_to_matter_feedback",
                         "nonabelian_expected_first_variations_zero", "Higgs_Yukawa_expected_first_variations_zero"):
                if summary.get(flag) is not True:
                    raise SystemExit("coupled current requires replayed action and omitted mean forces")
            if (summary.get("electric_kinetic") != "fixed_lambda_convex_nonquadratic"
                    or summary.get("quadratic_Maxwell_kinetic") is not False
                    or summary.get("global_Z6_quotient_selected") is not False):
                raise SystemExit("coupled current must retain its regulated cover-action boundary")
        controls[key] = {
            "receipt": "code/sm_fermion_current/" + receipt,
            "receipt_pin": pin,
            "independent_verifier_result": summary,
        }
    return {
        "controls": controls, "mathematical_replay": True,
        "analytic_proof": "paper/tex_fragments/FERMION_SOURCE_CURRENT.tex",
        "coupled_action_proof": "paper/tex_fragments/FERMION_COUPLED_ACTION.tex",
        "observed_postdiction": False,
        "scope_boundary": (
            "The prepared q5 carrier supports exact finite CAR hypercharge transport. "
            "Classical electric feedback preserves mean Gauss in the delocalized "
            "Slater preparation, whose charge variance excludes the operator constraint. "
            "A separately supplied localized preparation and integer quantum link "
            "give an invariant flux-dressed sector, operator Gauss and an exact "
            "continuous one-link current history. Species weights, quantum state, "
            "action and time are declared. A separate regulated hybrid action adds "
            "all electric and Wilson factors and actual electric-to-matter feedback. "
            "Its invariant expectation sector checks all omitted Higgs/Yukawa and "
            "nonabelian first variations, with one occupied family and two vacua. "
            "The fixed nonquadratic electric law and direct-product cover do not "
            "supply the quadratic Maxwell action, global quotient, full quantum "
            "Gauss, excited nonabelian/Higgs execution, chiral continuum regulation, "
            "source selection or physical current attachment. None is an observed "
            "physical postdiction."),
    }


def _cartan_scalar_execution_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "sm_abelian_reduction", "verify_abelian.py", "abelian_receipt.json",
        receipt_path=receipt_path)
    if (summary.get("verified") is not True or summary.get("mathematical_replay") is not True
            or summary.get("continuous_history_error_enclosure") is not False
            or summary.get("producer_reverse_diagnostic_independently_replayed") is not False
            or summary.get("events") != sum(len(r["events"]) for r in packet["runs"])):
        raise SystemExit("Cartan execution must distinguish mathematical replay from custody and diagnostics")
    return {
        "receipt": "code/sm_abelian_reduction/abelian_receipt.json", "receipt_pin": pin,
        "independent_verifier_result": summary, "mathematical_replay": True,
        "analytic_proof": "paper/tex_fragments/CARTAN_SCALAR_REDUCTION.tex",
        "exact_symbolic_local_classical_reduction": True,
        "ideal_split_gauss_preservation": True, "trajectory_values_numerical": True,
        "continuous_history_error_enclosure": False, "observed_postdiction": False,
        "source_action_population_or_physical_clock_selected": False,
        "global_compact_subgroup_or_quantum_truncation": False,
        "scope_boundary": (
            "The kinetic-metric-aligned Cartan and upper-Higgs restriction is a local "
            "classical reduction of the declared action, with a principal-log plaquette "
            "chart in the discrete extension. Three q5 cohorts retain all local "
            "read/write operations. Exact symbolic reduction and ideal current/Gauss "
            "identities are separate from the 60-digit numerical replay of binary64 "
            "trajectories. The replay supplies no continuous-time error enclosure, "
            "q233 detector accuracy, global compact subgroup or quantum truncation."),
    }


def _cartan_scalar_continuous_readout_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, summary, pin = _structural_packet(
        "sm_abelian_readout", "verify_readout.py", "readout_receipt.json",
        dependency="rational", receipt_path=receipt_path)
    if (summary.get("verified") is not True or summary.get("mathematical_replay") is not True
            or type(summary.get("parent_events")) is not int or summary["parent_events"] != 5678):
        raise SystemExit("continuous Cartan readout requires full fresh parent replay")
    maximum = max(Fraction(r["error_upper"]) for r in packet["checkpoints"])
    expected = {"verified": True, "mathematical_replay": True,
        "continuous_classical_detector_bound": True,
        "current_lower": packet["current"]["lower_on_window"],
        "neighbor_response_lower_at_end": packet["neighbor_response"]["magnitude_lower_at_end"],
        "completed_checkpoints": 9, "max_checkpoint_intensity_error": str(maximum),
        "parent_events": 5678, "physical_premise_discharged": False, "observed_postdiction": False}
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("continuous Cartan independent readout projection mismatch")
    if not (0 < maximum < Fraction(6, 10**6)
            and Fraction(packet["current"]["lower_on_window"]) > 0
            and Fraction(packet["electric_feedback"]["lower_at_end"]) > 0):
        raise SystemExit("continuous Cartan readout error or signal budget failed")
    source = "paper/tex_fragments/CARTAN_SCALAR_READOUT.tex"
    label = "thm:cartan-scalar-readout"
    if "\\label{" + label + "}" not in (REPO/source).read_text(encoding="utf-8"):
        raise SystemExit("continuous Cartan analytic theorem missing")
    return {
        "receipt": "code/sm_abelian_readout/readout_receipt.json", "receipt_pin": pin,
        "independent_verifier_result": summary,
        "immutable_execution_parent_sha256": packet["parent_sha256"],
        "analytic_proof": {"source": source, "label": label},
        "mathematical_replay": True, "continuous_classical_current_and_readout_bound": True,
        "decoded_completed_checkpoint_errors": True,
        "window": packet["window"], "resolved_window": packet["resolved_window"],
        "current": packet["current"], "electric_feedback": packet["electric_feedback"],
        "neighbor_response": packet["neighbor_response"],
        "first_kick_is_exact_time_sample": False,
        "whole_state_numerical_trajectory_enclosed": False,
        "quantum_dynamics_or_born_outcomes": False, "full_fermion_generation_executed": False,
        "source_action_population_or_physical_clock_selected": False, "observed_postdiction": False,
        "scope_boundary": (
            "Exact rational first-exit and acceleration-difference estimates bound the nonlinear "
            "classical current, electric feedback and spatially distinct scalar-intensity response "
            "for the same supplied 64-site Cartan action through time 1/100. Every initial/completed "
            "readout is decoded from actual versioned writes and compared to the exact IVP; all "
            "nine intensity errors are below 6e-6. The ideal gauge copy is handled by exact gauge "
            "covariance, and binary64 initial values are distinguished from the exact preparation. "
            "The intermediate half-kick is not a full time sample. The older execution receipt "
            "retains its numerical-only scope. No full-state enclosure, physical clock, quantum "
            "outcomes, full fermion evolution or laboratory postdiction is supplied."),
    }


def _local_sm_action_control(receipt_path: Path | None = None) -> dict[str, Any]:
    packet, custody, pin = _structural_packet(
        "sm_local_action", "verify_local_action.py", "local_action_receipt.json",
        entry="validate_custody", receipt_path=receipt_path)
    if custody.get("accepted_custody") is not True or custody.get("mathematical_replay") is not False:
        raise SystemExit("local SM custody must not masquerade as mathematical replay")
    return {
        "receipt": "code/sm_local_action/local_action_receipt.json", "receipt_pin": pin,
        "custody_verifier_result": custody, "mathematical_replay": False,
        "full_mathematical_verifier": "code/sm_local_action/verify_local_action.py:verify",
        "observed_postdiction": False,
        "packet_inventory": {
            "action_sectors": len(packet["action_coefficients"]),
            "grassmann_generators": len(packet["generator_labels"]),
            "local_ward_directions": len(packet["ward_directions"]),
            "symmetric_second_parameter_jets": len(packet["symmetric_second_parameter_jets"]),
            "matter_currents": len(packet["matter_current_coefficients"]),
            "negative_controls": len(packet["negative_controls"]),
        },
        "scope_boundary": (
            "This ledger checks custody, exact encoding and declared hypotheses only. "
            "The separate full verifier reconstructs all action/current/Ward coefficients. "
            "The local flat classical action supplies connections, Spin convention, Higgs, "
            "three families and free parameters; it does not emit a source action, join the "
            "spatial population or construct a renormalized quantum theory."),
    }


def _source_count_clock_control(receipt_path: Path | None = None) -> dict[str, Any]:
    """Replay the finite decoder; retain the analytic/physical scope boundary."""
    source = "paper/tex_fragments/SOURCE_COUNT_CLOCK.tex"
    labels = ("thm:source-count-clock", "prop:source-count-clock-error",
              "prop:source-count-clock-uniqueness")
    prose = (REPO / source).read_text(encoding="utf-8")
    if any("\\label{" + label + "}" not in prose for label in labels):
        raise SystemExit("count-clock analytic theorem missing")
    directory = CODE / "causal_refinement"
    verifier_path = directory / "verify_source_count_clock.py"
    spec = importlib.util.spec_from_file_location("count_clock_ledger_verifier", verifier_path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent count-clock verifier")
    verifier = importlib.util.module_from_spec(spec)
    exec(compile(verifier_path.read_bytes(), str(verifier_path), "exec"), verifier.__dict__)
    receipt_path = receipt_path or directory / "source_count_clock_receipt.json"
    raw = receipt_path.read_bytes()
    try:
        summary = verifier.verify(verifier.load(receipt_path))
    except ValueError as exc:
        raise SystemExit(f"count-clock independent replay failed: {exc}") from exc
    if receipt_path.read_bytes() != raw:
        raise SystemExit("count-clock receipt changed during replay")
    expected = {
        "schema": "source-count-clock-v1", "source_q": 5,
        "authenticated_events": 500, "interval_counts": [2, 41, 80],
        "exact_clock_bounds": 3, "conditional_error_controls": 3,
        "proper_time_ratio_limit_is_analytic": True,
        "finite_clock_accuracy_certified": False, "physical_clock_identified": False,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("count-clock verifier scope/count mismatch")
    declarations = ("fourth_power_order", "fourth_root_unique", "common_weight_cancels",
                    "volume_ratio_enclosure", "clock_enclosure")
    return {
        "source": source, "labels": list(labels),
        "lean_receipts": _lean_receipt("SourceCountClock",
            declarations={"SourceCountClock": declarations}),
        "lean_declarations": {"SourceCountClock": list(declarations)},
        "receipt_pin": {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()},
        "independent_verifier_result": summary,
        "count_decoder_uses_timestamps_or_density": False,
        "analytic_same_family_proper_time_ratio_limit": True,
        "analytic_volume_only_inertial_clock_uniqueness": True,
        "analytic_resolved_timelike_worldline_accumulation": True,
        "count_volume_or_curve_limit_formalized_in_lean": False,
        "observed_postdiction": False,
        "scope_boundary": (
            "The geometric clock is retrospective and requires access to complete "
            "authenticated interval ancestry and a reference. The common order/count "
            "limit is the separately declared source-net law. Finite error controls "
            "check supplied arithmetic budgets, not actual q=5 timing accuracy. "
            "Native law selection, a running matter pointer and physical calibration "
            "are not derived. Lean checks finite scale and interval algebra only."
        ),
    }


def _protected_memory_control(receipt_path: Path | None = None) -> dict[str, Any]:
    """Replay the supplied stochastic law before projecting its exact evidence."""
    source = "paper/tex_fragments/PROTECTED_RECORD_MEMORY.tex"
    label = "prop:protected-record-memory"
    path = REPO / source
    if not path.is_file() or "\\label{" + label + "}" not in path.read_text(encoding="utf-8"):
        raise SystemExit("protected-memory analytic theorem missing")
    directory = CODE / "thermodynamics" / "protected_memory"
    verifier_path = directory / "verify_protected_memory.py"
    spec = importlib.util.spec_from_file_location("protected_memory_ledger_verifier", verifier_path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent protected-memory verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt_path = receipt_path or directory / "runtime" / "protected_memory_receipt.json"
    raw = receipt_path.read_bytes()
    try:
        summary = verifier.verify(verifier.load(receipt_path))
    except ValueError as exc:
        raise SystemExit(f"protected-memory independent replay failed: {exc}") from exc
    if receipt_path.read_bytes() != raw:
        raise SystemExit("protected-memory receipt changed during replay")
    expected = {
        "accepted": True,
        "scope": "EXACT_NEW_LOCAL_STOCHASTIC_LAW__PROTECTED_RECORD_AND_DECAYING_MEMORY",
        "states": 8, "protected_fibres": 2, "current_count": 2,
        "correlation_lags": 65, "distribution_samples": 8,
        "two_step_minorization": "7/96", "gap_lower": "7/192",
        "green_kubo_matrix": [["11417/9375", "1832/3125"], ["1832/3125", "3816/3125"]],
        "native_source_attachment": False,
        "historical_obstruction_overturned": False,
        "physical_clock_or_conductivity": False,
        "new_Lean_formalization": False,
    }
    if json.dumps(summary, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise SystemExit("protected-memory verifier scope/count mismatch")
    return {
        "source": source, "label": label,
        "receipt_pin": {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()},
        "independent_verifier_result": summary,
        "new_supplied_local_stochastic_law": True,
        "fibre_centered_currents_required": True,
        "observed_postdiction": False,
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
                "and same-current holonomy are not constructed"
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
                "same-source loop-to-kernel theorem are not constructed"
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
            "declared_local_action_control": _local_sm_action_control(),
            "cartan_scalar_execution_control": _cartan_scalar_execution_control(),
            "cartan_scalar_continuous_readout_control": _cartan_scalar_continuous_readout_control(),
            "fermion_source_current_control": _fermion_source_current_control(),
            "statement": (
                f"Inside the declared {component_count}-component "
                "exterior-response algebra, an exhaustive scan of all "
                f"{subsets_enumerated} subsets selects exactly one "
                "unordered charge-conjugate pair of nonempty chiral "
                "anomaly-free rank-"
                f"{matter['realized_package']['dimension']} projectors. "
                "Primitive determinant balance "
                "fixes the block charges up to conjugation, and the selected "
                f"representative has multiset {{{field_summary}}}. A separate declared "
                "flat classical action assembles gauge, Grassmann-Weyl, Higgs and symbolic "
                "Yukawa terms with local Ward and variational-current controls. This ledger "
                "checks that packet's custody; its separate full verifier checks the coefficients. "
                "A kinetic-aligned local Cartan/Higgs reduction additionally supplies a q5 "
                "charged execution with variational current and electric feedback, whose "
                "ideal algebra and numerical operation replay are independently checked. "
                "A separate exact consumer bounds that action's continuous classical current "
                "and spatially distinct intensity response, with every completed-checkpoint "
                "intensity error below 6e-6; this is no observed physical postdiction"
            ),
            "observed_counterpart": (
                "Standard Model one-generation hypercharge assignment"
            ),
            "realized_spectrum": spectrum,
            "match": "exact" if spectrum == sm_spectrum else "MISMATCH",
            "artifact_refs": [_rel("matter_receipt"), _rel("matter_menu"),
                "paper/tex_fragments/LOCAL_SM_JET_ACTION.tex", "Lean/Screen/LocalGaugeJetAction.lean",
                "code/sm_local_action/local_action_receipt.json", "code/sm_local_action/verify_local_action.py",
                "paper/tex_fragments/CARTAN_SCALAR_REDUCTION.tex", "Lean/Screen/CartanScalarReduction.lean",
                "code/sm_abelian_reduction/abelian_receipt.json", "code/sm_abelian_reduction/verify_abelian.py",
                "paper/tex_fragments/CARTAN_SCALAR_READOUT.tex",
                "code/sm_abelian_readout/readout_receipt.json", "code/sm_abelian_readout/verify_readout.py"],
            "subset_count": subsets_enumerated,
            "survivor_count": classification["survivor_count"],
            "survivor_dimension": matter["realized_package"]["dimension"],
            "derived_block_charges": matter_menu["declared_algebra"][
                "derived_block_charges"
            ],
            "lean_declarations": {
                "ExteriorSelection": list(exterior_declarations),
                "WeylYukawaConventions": ["conjugate_charge_dictionary", "weyl_up_neutral",
                    "weyl_down_neutral", "weyl_lepton_neutral", "dirac_up_neutral",
                    "dirac_down_neutral", "dirac_lepton_neutral", "mixed_convention_not_neutral"],
            },
            "lean_receipts": _lean_receipt(
                "ExteriorSelection", "WeylYukawaConventions",
                declarations={"ExteriorSelection": exterior_declarations,
                    "WeylYukawaConventions": ("conjugate_charge_dictionary", "weyl_up_neutral",
                        "weyl_down_neutral", "weyl_lepton_neutral", "dirac_up_neutral",
                        "dirac_down_neutral", "dirac_lepton_neutral", "mixed_convention_not_neutral")},
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
                "physical-instrument attachments not constructed"
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
                "and the publicization information chain are not constructed. "
                "Scientific owner #730 records those quantum-composition "
                "obligations; #739 records any residual premise-discharge "
                "obligation. "
                "No source rule selects that partition, a state, an observable, "
                "a detector algebra, or a public instrument, and no "
                "edge-to-partition identification is constructed"
            ),
            "paper_ref": "finite-event-algebra paper",
        },
        {
            "id": "modal_maxwell_factorization_boundary",
            "statement": (
                "On the complexified nonzero-momentum transverse fibre, "
                "the Fourier-modal block D=i[omega_a(k)/|k|](k cross -) "
                "has zero modal divergence and squares to the committed "
                "FZ-12 scalar spatial action. The opposite-sign pair "
                "G(E,B)=(D B,-D E) therefore squares to the existing "
                "second-order oscillator on both amplitudes; the same-sign "
                "mutation has the wrong second-order sign"
            ),
            "observed_counterpart": (
                "Fourier-modal algebra of the vacuum Maxwell curl pair; "
                "not a physical electromagnetic observation"
            ),
            "match": (
                "exact bounded modal factorization; local source-produced "
                "physical Maxwell theory not constructed"
            ),
            "lean_declarations": {
                "ModalMaxwellFactorizationBoundary": [
                    "modalCurlScale_sq_mul_dot_self",
                    "dot_modalCurl_zero",
                    "modalCurl_sq_on_transverse",
                    "complexMomentumDot_fourierCurl_zero",
                    "complexPhotonSpatialAction_complexifies",
                    "fourierCurl_sq_on_transverse",
                    "maxwellShapedModalGenerator_sq_wave",
                    "maxwellShapedModalGenerator_transverse",
                    "sameSignCurlMutation_sq_positive",
                    "sameSignCurlMutation_fails_wave",
                ]
            },
            "lean_receipts": _lean_receipt(
                "ModalMaxwellFactorizationBoundary",
                declarations={
                    "ModalMaxwellFactorizationBoundary": (
                        "modalCurlScale_sq_mul_dot_self",
                        "dot_modalCurl_zero",
                        "modalCurl_sq_on_transverse",
                        "complexMomentumDot_fourierCurl_zero",
                        "complexPhotonSpatialAction_complexifies",
                        "fourierCurl_sq_on_transverse",
                        "maxwellShapedModalGenerator_sq_wave",
                        "maxwellShapedModalGenerator_transverse",
                        "sameSignCurlMutation_sq_positive",
                        "sameSignCurlMutation_fails_wave",
                    )
                },
            ),
            "artifact_refs": [
                "code/electromagnetism/modal_maxwell_factorization.py",
                "code/electromagnetism/test_modal_maxwell_factorization.py",
                "paper/screen_microphysics_and_observer_synchronization.tex",
            ],
            "hypothesis_boundary": (
                "this is a pointwise Fourier-modal/pseudodifferential "
                "factorization of a committed mathematical "
                "oscillator. Its momentum-dependent multiplier supplies no "
                "local position-space operator, real-field assembly or "
                "opposite-momentum reality pairing, electric/magnetic "
                "identification, source-produced dynamics, U(1) potential "
                "or gauge quotient, Maxwell action, bridge from the finite "
                "Gauss receipts to this modal divergence, conserved physical "
                "current or source coupling, Lorentz covariance, continuum "
                "control, or laboratory readout. It removes no registered "
                "premise: PR-20, PR-21, and PR-22 remain declared inputs, "
                "while PR-53 and PR-54 name missing attachments. Scientific "
                "owner #733 records those obligations. This row "
                "emits no frozen prediction"
            ),
            "paper_ref": (
                "screen-microphysics paper, modal Maxwell factorization boundary"
            ),
        },
        {
            "id": "finite_unitary_scattering_limit_no_go",
            "statement": (
                "In every Hausdorff topological group, ordinary convergence "
                "of the natural powers g^n forces g to be the identity. Thus "
                "a supplied nonidentity finite-dimensional exact unitary U "
                "at fixed cutoff has no ordinary large-time limit in either "
                "the unitary subgroup or the ambient matrix topology, and "
                "hence has no full weak-operator limit at that same fixed "
                "finite dimension. The "
                "exact scope control (U^n)^(-1) U^n = 1 converges for every U"
            ),
            "observed_counterpart": (
                "fixed-cutoff scattering-construction boundary; not a "
                "physical observation or prediction"
            ),
            "match": (
                "exact direct-power obstruction; physical scattering and "
                "asymptotic comparison construction not constructed"
            ),
            "lean_declarations": {
                "FiniteUnitaryScatteringNoGo": [
                    "tendsto_powers_forces_identity",
                    "nontrivial_powers_have_no_limit",
                    "finite_unitary_powers_have_no_limit",
                    "finite_unitary_ambient_powers_have_no_limit",
                    "identical_relative_evolution_is_constant",
                    "identical_relative_evolution_tendsto",
                ]
            },
            "lean_receipts": _lean_receipt(
                "FiniteUnitaryScatteringNoGo",
                declarations={
                    "FiniteUnitaryScatteringNoGo": (
                        "tendsto_powers_forces_identity",
                        "nontrivial_powers_have_no_limit",
                        "finite_unitary_powers_have_no_limit",
                        "finite_unitary_ambient_powers_have_no_limit",
                        "identical_relative_evolution_is_constant",
                        "identical_relative_evolution_tendsto",
                    )
                },
            ),
            "artifact_refs": [
                "code/qft/finite_unitary_scattering_no_go.py",
                "code/qft/test_finite_unitary_scattering_no_go.py",
                "paper/tex_fragments/QFT_STRUCTURAL_INHERITANCE_STATUS.tex",
            ],
            "hypothesis_boundary": (
                "this is only a no-go for direct ordinary convergence of one "
                "fixed-cutoff power sequence U^n. Full weak-operator "
                "convergence at that same finite dimension is excluded too, "
                "because the standard finite-dimensional operator topologies "
                "coincide. It leaves relative or comparison dynamics, "
                "selected projected scalar or observable limits, infinite-"
                "dimensional weak limits, subsequential or Cesaro limits, "
                "continuum or infinite-volume limits, open-"
                "system evolution, and finite-time operational protocols "
                "unresolved by this result. It constructs no wave operator, S-matrix, cross "
                "section, pole, optical theorem, renormalization flow, or "
                "continuum QFT. Scientific owner #743 records those interacting-QFT, "
                "renormalization, and scattering obligations. This row emits "
                "no frozen prediction"
            ),
            "paper_ref": (
                "consensus paper, structural-QFT inheritance boundary"
            ),
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
                "matter attachment not constructed"
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
                "refinement, and gluing not constructed"
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
            "id": "layered_discrete_gauss_boundary",
            "statement": (
                "On a supplied finite directed graph with an integer layering, "
                "a finite depth, steady sourcing, shell equidistribution, and "
                "shell cardinality c n^2, exact divergence bookkeeping gives "
                "outward flux per shell edge Q/(c n^2) and makes flux times "
                "n^2 shell-independent. A three-vertex chain inhabits the "
                "premises. Two exact global controls show why the depth bound "
                "is essential: a globally steady drained source on a closed "
                "finite carrier has zero total load, and an n^2 shell law at "
                "every positive layer forces c=0"
            ),
            "observed_counterpart": (
                "a depth-bounded discrete inverse-square precursor with no "
                "physical radius, gravitational-field, or mass identification"
            ),
            "match": (
                "exact conditional finite theorem under PR-29, PR-30, and "
                "PR-31; physical Newtonian attachment not constructed"
            ),
            "lean_declarations": {
                "LayeredDiscreteGauss": [
                    "sum_vertexOutflow_eq_regionFlux",
                    "steady_total_source_eq_zero",
                    "unbounded_shellCard_forces_zero",
                    "regionFlux_eq_charge",
                    "shell_total_eq_charge",
                    "perEdgeFlux_eq",
                    "scale_free",
                    "chainWitness_charge",
                    "regionFlux_twelvePort",
                    "steady_witness_regionFlux_eq_source",
                ]
            },
            "lean_receipts": _lean_receipt(
                "LayeredDiscreteGauss",
                declarations={
                    "LayeredDiscreteGauss": (
                        "sum_vertexOutflow_eq_regionFlux",
                        "steady_total_source_eq_zero",
                        "unbounded_shellCard_forces_zero",
                        "regionFlux_eq_charge",
                        "shell_total_eq_charge",
                        "perEdgeFlux_eq",
                        "scale_free",
                        "chainWitness_charge",
                        "regionFlux_twelvePort",
                        "steady_witness_regionFlux_eq_source",
                    )
                },
            ),
            "hypothesis_boundary": (
                "PR-29 supplies the c n^2 shell count, PR-30 supplies equal "
                "outward flux within each shell, and PR-31 supplies steady "
                "sourcing through the declared depth. The exponent and "
                "isotropy are therefore conditional inputs. The theorem "
                "constructs no physical radius, field, mass density, "
                "continuum limit, or Einstein-branch join"
            ),
            "paper_ref": (
                "spacetime-recovery paper, layered discrete Gauss precursor"
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
                "readout, and comparison not constructed"
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
                "comparison not constructed"
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
            "match": "exact formal boundary; physical time realization not constructed",
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
                "physical clock not constructed"
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
            "match": "exact structural interface; source realization not constructed",
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
                "exact bounded conditional endpoint; source and limit not constructed"
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
                "four-dimensional ambient Lorentz module and observer-frame algebra"
            ),
            "match": (
                "exact ambient module algebra and coordinate bridge; source-causal "
                "event attachment and physical soldering remain outside C1"
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
                "rest spaces are algebraic. C1 alone supplies no source-event "
                "attachment or soldering. The separate source-derived causal-order "
                "packet and bounded C2 contract supply, respectively, an "
                "informational poset and algebraic overlap covariance from declared "
                "inputs. Faithful physical placement, calibrated density, "
                "manifoldlike refinement, rods, clocks, smooth spacetime, "
                "observable, decision rule, and prediction are not constructed"
            ),
            "paper_ref": "observers paper, canonical Lorentz module",
        },
        {
            "id": "source_derived_finite_one_three_causal_carrier",
            "supplied_law_refinement_control": _refining_causal_control(),
            "declared_source_net_control": _source_net_causal_control(),
            "declared_count_clock_control": _source_count_clock_control(),
            "protected_population_control": _protected_population_control(),
            "source_population_quadrature": _source_population_quadrature_control(),
            "common_free_scalar_control": _common_source_scalar_control(),
            "authenticated_scalar_execution_control": _source_scalar_execution_control(),
            "finite_scalar_quantum_probability_control": _source_scalar_quantum_control(),
            "reconstructed_clock_quantum_control": _source_scalar_clock_quantum_control(),
            "joint_scalar_time_refinement_control": _source_scalar_time_refinement_control(),
            "regional_scalar_time_slice_control": _source_scalar_regional_control(),
            "destination_local_seam_readback_control": _source_seam_routing_control(),
            "reusable_local_feedback_transport_control": _source_feedback_transport_control(),
            "sequential_scalar_quantum_instrument_control": _source_scalar_sequential_instrument_control(),
            "finite_duration_scalar_quantum_instrument_control": _source_scalar_finite_instrument_control(),
            "operational_cone_selection": {
                "analytic_proof": "paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex",
                "hypotheses": "nonzero closed convex pointed cone; finite irreducible carrier rotations; one continuous operational boost direction; time orientation",
                "conclusion": "unique future Lorentz cone",
                "native_operational_boost_covariance_derived": False,
                "cone_classification_formalized_in_lean": False,
                "observed_postdiction": False,
            },
            "metric_scalar_continuum": {
                "classification": "conditional_analytic_scalar_limit_not_observed_spacetime",
                "spatial_selection": "maximal separated finite conservative source-word menu in a metric ball",
                "action": "positive Voronoi-mass radial-kernel Klein-Gordon action",
                "operator_error": "(15/112)*epsilon^2*M4+960*h/epsilon^2*M1",
                "limit": "epsilon->0 and h/epsilon^2->0; smooth interior references",
                "physical_source_or_clock_selected": False,
                "quantum_continuum_claim": False,
                "analytic_proof": "paper/tex_fragments/SOURCE_METRIC_SCALAR_CONTINUUM.tex",
                "analytic_same_action_boost_detector_comparison": True,
                "boost_comparison_boundary": "spacetime-smearing and two transformed preparations supplied through a smooth continuum reference with common window/buffer; no finite-history boost map or raw-support covariance",
                "finite_algebra_receipt": _lean_receipt(
                    "MetricKernelEnergy",
                    declarations={"MetricKernelEnergy": (
                        "dirichlet_nonnegative", "dirichlet_constant",
                        "symmetric_second_moment", "weighted_degree_bound",
                        "massive_energy_nonnegative",
                    )},
                ),
            },
            "artifact_refs": [
                "Lean/Geometry/SourceRecordProtection.lean",
                "Lean/Geometry/SourcePopulationQuadrature.lean",
                "Lean/Geometry/GoldenSourceAssignment.lean",
                "Lean/Geometry/GoldenSourceCountLimit.lean",
                "Lean/Geometry/SourceFeedbackTransport.lean",
                "paper/tex_fragments/SOURCE_FEEDBACK_TRANSPORT.tex",
                "code/source_feedback_transport/transport_receipt.json",
                "code/source_feedback_transport/verify_transport.py",
                "paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex",
                "code/source_scalar_instruments/sequential_instrument_receipt.json",
                "code/source_scalar_instruments/verify_sequential_instrument.py",
                "paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex",
                "code/source_scalar_finite_instrument/finite_instrument_receipt.json",
                "code/source_scalar_finite_instrument/verify_finite_instrument.py",
                "Lean/Geometry/SourceSeamPathTomography.lean",
                "paper/tex_fragments/SOURCE_SEAM_PATH_TOMOGRAPHY.tex",
                "code/source_routing/runtime/path_tomography_receipt.json",
                "code/source_routing/verify_routing.py",
                "Lean/QFT/ScalarRegionalTimeSlice.lean",
                "paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex",
                "code/source_scalar_regional/regional_time_slice_receipt.json",
                "code/source_scalar_regional/verify_regional_time_slice.py",
                "paper/tex_fragments/SOURCE_RECORD_PROTECTION.tex",
                "paper/tex_fragments/SOURCE_POPULATION_QUADRATURE.tex",
                "code/source_population/runtime/population_receipt.json",
                "code/source_population/verify_population.py",
                "paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex",
                "code/source_scalar_packet/source_common_scalar_receipt.json",
                "code/source_scalar_packet/verify_source_common_scalar.py",
                "paper/tex_fragments/SOURCE_SCALAR_EXECUTION.tex",
                "code/source_scalar_execution/source_scalar_execution_receipt.json",
                "code/source_scalar_execution/verify_source_scalar_execution.py",
                "paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex",
                "code/source_scalar_quantum/quantum_probability_receipt.json",
                "code/source_scalar_quantum/verify_quantum_probability.py",
                "paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex",
                "Lean/Screen/ActionTimeGram.lean",
                "Lean/Screen/SourceActionTime.lean",
                "code/source_scalar_clock/source_scalar_clock_receipt.json",
                "code/source_scalar_clock/verify_source_scalar_clock.py",
                "paper/tex_fragments/SOURCE_SCALAR_CLOCK_QUANTUM.tex",
                "code/source_scalar_clock_quantum/clock_quantum_receipt.json",
                "code/source_scalar_clock_quantum/verify_clock_quantum.py",
                "paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex",
                "code/source_scalar_time_refinement/source_scalar_time_refinement_receipt.json",
                "code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py",
                "Lean/Geometry/MetricKernelEnergy.lean",
                "paper/tex_fragments/SOURCE_METRIC_SCALAR_CONTINUUM.tex",
                "Lean/Geometry/RefiningLatticeCausalCone.lean",
                "paper/tex_fragments/REFINING_CAUSAL_CONE.tex",
                "code/causal_refinement/refining_cone.py",
                "code/causal_refinement/verify_refining_cone.py",
                "code/causal_refinement/test_refining_cone.py",
                "code/causal_refinement/refining_cone_receipt.json",
                "Lean/Geometry/SourceNetCausalCone.lean",
                "paper/tex_fragments/SOURCE_NET_CAUSAL_LIMIT.tex",
                "code/causal_refinement/source_net_causet.py",
                "code/causal_refinement/verify_source_net_causet.py",
                "code/causal_refinement/test_source_net_causet.py",
                "code/causal_refinement/source_net_causet_receipt.json",
                "paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex",
                "paper/tex_fragments/SOURCE_COUNT_CLOCK.tex",
                "Lean/Time/SourceCountClock.lean",
                "code/causal_refinement/source_count_clock.py",
                "code/causal_refinement/verify_source_count_clock.py",
                "code/causal_refinement/test_source_count_clock.py",
                "code/causal_refinement/source_count_clock_receipt.json",
            ],
            "statement": (
                "Authenticated read-after-write provenance generates a finite "
                "locally finite partial order, and canonical source height obeys "
                "the exact longest-parent-path recursion. Independently of that "
                "order, an independent real axis and the proved rank-three source "
                "Gram quotient define an exact four-dimensional carrier with "
                "signature (+---); every source-unit spatial direction "
                "gives a future-null vector. A supplied positive time scale and "
                "spatial readback use source height only in the event placement. "
                "With the parent-edge "
                "speed bound, generated precedence maps into the carrier future "
                "cone. A separately certified faithful placement upgrades this to "
                "two-way order--cone equivalence and constructs the finite "
                "source-order/frame packet. The explicit four-event Boolean "
                "diamond is a non-chain control whose generated order agrees "
                "exactly with its carrier cone order. A separate supplied growing-menu "
                "lattice has exact finite inner/outer cone bounds and an analytic "
                "controlled cone and interval-volume limit; its local read/write "
                "control has 81 events, 794 authenticated edges and width 27. "
                "A separate conditional construction selects finite populations "
                "from the dense source metric and gives a positive scalar action "
                "with an analytic continuum field and detector error bound; "
                "physical selection of the coarsening and action is not derived. "
                "A declared local-read law on golden-orbit conservative source "
                "records has Lean-constructed finite cone bounds and analytic "
                "raw-count four-volume and contained-interval ordering-fraction "
                "limits, the latter equal to 1/10. Independent exact finite "
                "executions use 27, 125, 512 and 2197 sites with equal widths; "
                "the sharper q=13 fill gives inner speed at least 75989/250000"
                ". On this declared family, fourth roots of authenticated interval-count "
                "ratios recover geometric proper-time ratios, uniquely up to units among "
                "positive volume-only inertially additive readouts. Exact q=5 replay "
                "checks the decoder. A separate analytic finite-rotation/one-boost "
                "criterion selects the Lorentz cone under operational covariance. "
                "The same declared scalar action has controlled spacetime-smeared "
                "boost comparisons for two reference-supplied preparations. A separate "
                "protected-address experiment exactly distinguishes immutable records from "
                "moving live readbacks under pair means. Formal integrated cell estimates "
                "derive quadrature convergence from vanishing assignment distances. Fibonacci "
                "arithmetic constructs the actual golden permutation partition, its exact "
                "cell masses and fixed-globally-Lipschitz integral convergence in Lean. Actual "
                "spatial and declared-time-grid counts also converge on fixed measurable sets "
                "with null frontier; the inclusive final-layer error is at most L^3*Delta_q. On the "
                "same prepared golden coordinates, a different local massive Dirichlet "
                "scalar action has an analytic free-field detector limit and a finite full-Fock "
                "compact detector response whose exact certified error is smaller than its signal. "
                "For that action, stable local time refinement converges to the same detector "
                "limit; nearest-update readout preserves a resolved finite comparison, and "
                "the preparation-induced response vanishes asymptotically before continuum "
                "arrival. In a separate 64-site execution, graph-recovered configurations "
                "determine positive durations for a supplied action family; conditional "
                "timing intervals retain 15 of 21 resolved fixed-state quantum predictions. "
                "A sharp kernel criterion also identifies the minimum exterior linear readouts "
                "needed to reconstruct the original regional canonical algebra from two field "
                "layers; eight exact regional controls retain every exterior coupling. Separate "
                "canonical-mean path calibrations reconstruct remote initial data from destination "
                "samples alone, with all 1449 events and their precision costs retained. "
                "Their all-depth infinity-norm condition number obeys the exact recurrence "
                "g_0=1, g_1=3, g_d=2*g_(d-1)+g_(d-2); a separate sufficient bound propagates "
                "supplied gate, retained-sample and decoding errors. "
                "With supplied local copy/reset and factor-two readback, reusable version "
                "transport instead has unit per-hop signal gain and a sharp linear additive "
                "error bound; 4928 events retain exact semantic order and old-version reads. "
                "A separate centered quantum instrument includes all earlier-readout backaction "
                "in 21 unconditional marginals and retains 15 resolved comparisons against "
                "continuous evolution under that same repeated instrument. "
                "An exact finite rectangular pointer pulse keeps field evolution active; "
                "bounded pulse duration and declared channel noise retain all fifteen responses "
                "when both baseline and intervention errors and all control windows are counted. "
                "These results do not identify a physical clock or select native field dynamics"
            ),
            "observed_counterpart": (
                "finite causal-set-like order with an effective 1+3 Lorentz "
                "carrier and null-cone directions"
            ),
            "match": (
                "exact source-derived finite order and 1+3 carrier theorem "
                "stack; physical faithful placement and continuum manifold "
                "not constructed"
            ),
            "lean_declarations": {
                "CausalInterval": [
                    "finiteCausalSetAxioms",
                ],
                "SemanticEventProvenance": [
                    "sourceHeight_eq",
                ],
                "SourceDerivedSpacetimeCarrier": [
                    "sourceSpacetimeCarrier_finrank",
                    "sourceCarrier_one_three_signature",
                    "sourceUnitNullVector_futureCausal",
                    "generatedBefore_sourceCausalLE",
                    "generatedBeforeEq_iff_sourceCausalLE",
                    "exactDiamondConeOrder",
                ],
                "SourceOrderFrameCompatibilityPacket": [
                    "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement",
                    "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame",
                    "sourceOrderFramePacket_consequences",
                ],
            },
            "lean_receipts": _lean_receipt(
                "CausalInterval",
                "SemanticEventProvenance",
                "SourceDerivedSpacetimeCarrier",
                "SourceOrderFrameCompatibilityPacket",
                declarations={
                    "CausalInterval": (
                        "finiteCausalSetAxioms",
                    ),
                    "SemanticEventProvenance": (
                        "sourceHeight_eq",
                    ),
                    "SourceDerivedSpacetimeCarrier": (
                        "sourceSpacetimeCarrier_finrank",
                        "sourceCarrier_one_three_signature",
                        "sourceUnitNullVector_futureCausal",
                        "generatedBefore_sourceCausalLE",
                        "generatedBeforeEq_iff_sourceCausalLE",
                        "exactDiamondConeOrder",
                    ),
                    "SourceOrderFrameCompatibilityPacket": (
                        "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement",
                        "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame",
                        "sourceOrderFramePacket_consequences",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "source height is a canonical ordinal, not a physical clock. "
                "The generic cone-inclusion theorem supplies neither its "
                "positive time scale nor a spatial readback valued in the "
                "rank-three quotient or an edge-speed bound; two-way "
                "order--cone faithfulness is the "
                "additional converse field of a faithful placement. Event-frame "
                "transports are supplied, and the standard-frame constructor is "
                "a trivial gauge choice rather than a derivation of physical "
                "observer frames. The Boolean diamond is an exact mathematical "
                "non-chain control, not empirical evidence of manifoldlikeness. "
                "No physical event/link identification, Poisson sprinkling, "
                "calibrated density or count--volume law, independent dimension "
                "estimator, topology, manifoldlike refinement, uniqueness "
                "theorem, continuum convergence, smooth metric, curvature, "
                "Einstein equation, observable, decision rule, or prediction is "
                "constructed by this finite source theorem. The separate growing-menu "
                "control supplies the frame, local move law, mesh, model layer time "
                "and event-cell measure. Its cone and interval-volume limits are "
                "analytic conditional results, not a source-selected physical "
                "manifold, Poisson sprinkling or calibrated density. The finite "
                "Lean bounds and replay do not formalize the Riemann-volume argument. "
                "The separate metric scalar limit supplies its coarsening, clipped "
                "Voronoi masses, radial-kernel force law, inertia and canonical "
                "quantization. Its analytic estimate assumes smooth references "
                "with an interaction-radius boundary buffer and h/epsilon^2 "
                "tending to zero. Five Lean lemmas check finite energy algebra, "
                "not the analytic PDE limit or physical source selection. "
                "The golden population, cube, all-neighbor read law, model tick "
                "and distinguished-event counting are declared. Its raw-count "
                "limit is not a native physical calibration, Poisson sprinkling "
                "or a finite-run dimension fit. Accepted primitive repairs, "
                "field-action agreement and laboratory clock attachment remain "
                "separate; neither the analytic volume nor pair-count proof is "
                "formalized by the finite Lean path module. The later GoldenSourceAssignment "
                "module formalizes the spatial golden partition and fixed globally Lipschitz "
                "Bochner quadrature. GoldenSourceCountLimit derives actual fixed continuity-set "
                "counts for its declared time grid. GoldenSourceCausalLimit and GoldenSourcePairLimit "
                "prove moving generated-interval and strict-pair convergence; FlatDiamondNormalization "
                "proves the one-tenth limit under window containment. GoldenSourceVolumeLimit "
                "proves the vertical weighted-volume error and limit. The separate moving-tip "
                "finite volume-error estimate remains analytic. "
                "Regional Weyl reconstruction requires supplied access to all real Weyl "
                "parameters and proves no operational quantum outcome or physical time slice. "
                "Destination-local path tomography is destructive and assumes protected memory, "
                "an isolated declared schedule and exact arithmetic; its sample-error norm does "
                "not include imperfect dynamics or decoding. A separate analytic bound "
                "g_d*(T*eta+sigma)+delta includes supplied gate, sample/retention and decoding "
                "error budgets relative to actual starting ports, without claiming sharpness "
                "for gate error or attained physical precision. Register-payload counts "
                "exclude metadata and scratch memory. Reusable transport adds supplied local "
                "classical export/reset/readout and fixed non-interleaved hops; semantic "
                "read-from projection is not physical resource order. The repeated quantum "
                "instrument supplies controlled gates, vacuum preparation, Born readout and "
                "inherited clock intervals. Its marginals are unconditional, not independent "
                "samples or branch-conditioned clocks, and its ideal scalar-action-graph "
                "pointer gadget supplies no W12 quantum routing or gate-duration/noise bound"
            ),
            "paper_ref": (
                "flagship and spacetime papers, source-derived causal order and "
                "finite 1+3 carrier"
            ),
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
                "exact bounded algebraic contract; source and physical receipts not constructed"
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
                "identify that atlas family. The source-derived packet separately "
                "supplies finite informational events and order plus an ambient "
                "1+3 target under explicit conditions. It does not supply physical "
                "event and link attachment, two-way order-cone faithfulness, "
                "source-selected atlas coordinates, separation, open charts, "
                "calibrated count-volume density, manifoldlike refinement, topology, "
                "an operational clock, or smooth curvature. No physical spacetime, "
                "Einstein dynamics, observable, decision rule, or prediction follows"
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
                "two-by-two matrix corners at every window. A five-module "
                "continuation conditionally proves full observer-algebra "
                "generation from post-hoc counted transitions and field "
                "projectors, exact commuting factors on a declared Cartesian "
                "carrier, a conditional-expectation diamond with coverage and "
                "an explicit noninjective-correlation control, and transport to "
                "one constant A3 stage on the distinct 86/88 carrier. The 86/247 capstone kernel-checks fully "
                "disjoint source supports and an exact counted correlated state "
                "with the committed marginals on the instantiated diamond"
            ),
            "observed_counterpart": (
                "a causal local quantum-observable net with overlap descent"
            ),
            "match": (
                "substantial conditional finite interface, operator-generation, "
                "CP-diamond, and support-disjoint counted-correlation packet; "
                "source justification of the factor reading and a nonconstant "
                "realization are not constructed"
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
                "SourceOperatorGeneration": [
                    "sourceAlgebra_eq_top",
                    "obs86_sourceAlgebra_eq_top",
                    "obs88_sourceAlgebra_eq_top",
                    "obs247_sourceAlgebra_eq_top",
                    "obs384_sourceAlgebra_eq_top",
                    "obs86_projectors_only_ne_top",
                ],
                "JointSlotFactorisation": [
                    "obs86_lifted_eq_leftSlot",
                    "obs88_lifted_eq_rightSlot",
                    "slot_commute",
                    "obs86_checkpoint_pinch_no_signalling",
                    "rightSlotExpectation_fixes_right",
                    "rightSlotExpectation_posSemidef",
                ],
                "CPRestrictionNet": [
                    "restrictCP_comp",
                    "slot_ranges_generate_top",
                    "no_scalar_restriction_of_matrix_factor",
                ],
                "TwoSlotCPNetWitness": [
                    "slotExpectations_not_jointly_injective",
                    "checkpoint_pinch_fixes_right",
                    "twoSlot_left_no_scalar_hom",
                ],
                "TowerAnchoredDiamond": [
                    "anchoredTower_privateAlgebra",
                    "partition_members_in_left_region",
                    "anchoredCheckpointPartition_proper",
                    "anchoredNet_left_ne_top",
                ],
                "SourceCorrelationCapstone": [
                    "supports_disjoint",
                    "slotExpectations_erase_source_correlation",
                ],
                "StructuralNetAdequacySurface": [
                    "supportDisjoint_correlation_counted_state",
                    "supportDisjointNet_expect_recovers_marginals",
                    "supportDisjoint_composed_adequacy",
                    "supportDisjoint_swap_covariance",
                ],
            },
            "lean_receipts": _lean_receipt(
                "FiniteCausalObserverNet",
                "ObserverNetDescent",
                "RichFibreWitness",
                "RichFibreRegionalNet",
                "SourceOperatorGeneration",
                "JointSlotFactorisation",
                "CPRestrictionNet",
                "TwoSlotCPNetWitness",
                "TowerAnchoredDiamond",
                "SourceCorrelationCapstone",
                "StructuralNetAdequacySurface",
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
                    "SourceOperatorGeneration": (
                        "sourceAlgebra_eq_top",
                        "obs86_sourceAlgebra_eq_top",
                        "obs88_sourceAlgebra_eq_top",
                        "obs247_sourceAlgebra_eq_top",
                        "obs384_sourceAlgebra_eq_top",
                        "obs86_projectors_only_ne_top",
                    ),
                    "JointSlotFactorisation": (
                        "obs86_lifted_eq_leftSlot",
                        "obs88_lifted_eq_rightSlot",
                        "slot_commute",
                        "obs86_checkpoint_pinch_no_signalling",
                        "rightSlotExpectation_fixes_right",
                        "rightSlotExpectation_posSemidef",
                    ),
                    "CPRestrictionNet": (
                        "restrictCP_comp",
                        "slot_ranges_generate_top",
                        "no_scalar_restriction_of_matrix_factor",
                    ),
                    "TwoSlotCPNetWitness": (
                        "slotExpectations_not_jointly_injective",
                        "checkpoint_pinch_fixes_right",
                        "twoSlot_left_no_scalar_hom",
                    ),
                    "TowerAnchoredDiamond": (
                        "anchoredTower_privateAlgebra",
                        "partition_members_in_left_region",
                        "anchoredCheckpointPartition_proper",
                        "anchoredNet_left_ne_top",
                    ),
                    "SourceCorrelationCapstone": (
                        "supports_disjoint",
                        "slotExpectations_erase_source_correlation",
                    ),
                    "StructuralNetAdequacySurface": (
                        "supportDisjoint_correlation_counted_state",
                        "supportDisjointNet_expect_recovers_marginals",
                        "supportDisjoint_composed_adequacy",
                        "supportDisjoint_swap_covariance",
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
                "not tensor factors or a `TensorSplitReceipt`. The later "
                "operator-generation theorem transcribes source paths post hoc. "
                "The later 86/247 capstone grounds fully disjoint committed "
                "supports and an exact counted correlated state with the correct "
                "marginals. Its erasure theorem proves that coverage does not "
                "imply unique reconstruction. The two-observers-as-factors "
                "reading and region map remain declared. The constant-tower "
                "transport is on the separate 86/88 carrier. TripleCarrierJoin "
                "couples both pair paths, counted states, checkpoints, and "
                "adjacent-step marginals on one carrier, but supplies no "
                "regional-net or tower morphism; neither packet is a "
                "nonconstant source realization. Scientific owner #728 "
                "records the missing source-attached or otherwise justified "
                "regional construction and a nonconstant source composition; "
                "scientific owner #739 records a "
                "missing source realization or genuinely scoped no-go. "
                "No CP/CPTP channel, scheduler locality, spacetime causality, "
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
            "match": "substantial exact finite packet; global converse not constructed",
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
                "exact bounded algebraic packet; physical attachment not constructed"
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
                "the closed unit ball. The displayed finite unsharp battery also "
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
                "and normalized antipodal binary contexts do not derive "
                "affinity; the transverse cubic refutes the displayed finite "
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
                "typed convention. The full-effect theorem applies only "
                "after full coexistent-effect additivity is supplied. No physical "
                "Born derivation, observable, or prediction is emitted. Scientific "
                "owner #730 records the missing source-earned phase instrument, operational "
                "additivity, and public readback; scientific owner #739 records the "
                "remaining affinity-principle obligation"
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
                "committed source artifact, however, the transition action has a "
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
                "The exact obstruction is only a bounded negative result for "
                "two direct mechanisms on the certified "
                "artifact: a mixing-mode-retaining linear intertwiner into the "
                "idempotent heat bath and a deterministic empirical pushforward. "
                "It does not exclude stochastic, nonlinear, reverse-direction, "
                "dilated, or enriched-source constructions. Scientific owner "
                "#739 records the missing replacement common reference, collar, "
                "objective, genuinely varying refinement family, source-clock "
                "derivation, and repair-export decision; empirical energy-clock "
                "calibration is PR-15. Closed #732 records only the attained "
                "conditional composition milestone. "
                "The "
                "pinned 20-state collar table has an audit of all 15 syntactic "
                "coordinate maps obtained from nonempty subsets of the four "
                "committed fields, inducing only four distinct partitions: "
                "its repair-load count aggregation is an "
                "eight-state ergodic nonreversible H-theorem probe, but it "
                "fails the exact count-kernel strong-lumpability test, the "
                "declared record charge is constant, and the common reference "
                "is unidentified; the fine chain's only recurrent restriction "
                "is singleton freezeout. The audit does not exclude arbitrary "
                "partitions or nonlinear, statistical, stochastic, weakly "
                "lumpable, or history-dependent maps. The "
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
                "conservation. A separate exact eight-state supplied local "
                "resampling law preserves a nonconstant record while its "
                "fibre-centred currents have decaying memory, a gap of at least "
                "7/192, and a rational Green--Kubo matrix with certified tails"
            ),
            "observed_counterpart": (
                "Onsager symmetry, Green--Kubo response, Fick diffusion, and "
                "Fourier heat transport"
            ),
            "match": (
                "exact finite conditional transport structure; physical "
                "generator and coefficients not constructed"
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
            "artifact_refs": [
                "code/thermodynamics/protected_memory/protected_memory.py",
                "code/thermodynamics/protected_memory/verify_protected_memory.py",
                "code/thermodynamics/protected_memory/test_protected_memory.py",
                "code/thermodynamics/protected_memory/runtime/protected_memory_receipt.json",
                "paper/tex_fragments/PROTECTED_RECORD_MEMORY.tex",
            ],
            "protected_record_memory": _protected_memory_control(),
            "hypothesis_boundary": (
                "the reversible kernel, linear Poisson solver, graph, distance, "
                "clock increment, volumes, heat capacities, and conductances "
                "are declared finite inputs. No theorem identifies the "
                "Green--Kubo coefficient with graph conductance. The eight-state "
                "witness supplies a new faithful reference and lazy conditional "
                "resampling law; its equilibrium projection P differs from its "
                "transition T. The Poisson and tail bounds require centering "
                "within each protected fibre; a globally centred conserved "
                "record has persistent correlation. This is no native "
                "source attachment and does not overturn the historical "
                "idempotent-projector obstruction. Scientific "
                "owners #728, #729, #737, and #739 record the missing source "
                "evolution, physical equilibrium reference and conserved quantity, "
                "source-realized geometry, instrumentation, and source clock; "
                "empirical calibration is PR-15, and closed #732 is only a "
                "conditional composition milestone. The bounded algebraic C2 "
                "receipt is separate, "
                "and this row emits no prediction-ladder "
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
            "match": "exact finite precursor; physical Ward bridge not constructed",
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
                "shrinking scale with separating tests. Scientific owner #729 "
                "records the obligation to define "
                "that residual from finite continuity, identify it with physical "
                "distributional divergence, and supply the source, transport, "
                "chart, and common-tower evidence before continuum Ward use; "
                "scientific owner #737 records the missing instrument attachment, "
                "and #739 records deferred premise discharge"
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
                "attachment not constructed"
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
                "conditionally attains coverage. The later E1 packet transcribes "
                "operators post hoc and consumes this row's partial-trace helper "
                "on a declared Cartesian slot. The later 86/247 capstone grounds "
                "fully disjoint source supports and the exact counted correlated "
                "state with its marginals, while its erasure theorem proves that "
                "coverage cannot reconstruct that state. Scientific owner #728 "
                "records the missing source-attached or otherwise justified "
                "factor reading and nonconstant source composition; #739 owns the residual source "
                "realization or no-go. The missing source channel/adaptive-scheduler semantics "
                "are recorded by #728, physical clocks by #739, physical spacetime "
                "attachment by #729, instrumentation by #737, and continuum "
                "causal/time-slice structure by #730. Those scientific owners mark "
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
                "channel attachment not constructed"
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
                "Scientific owner #728 records the missing source scheduler/channel "
                "attachment; #739 and #730 record the clock and continuum-causality "
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
                "curvatures one and two are distinct on both faces. Supplied "
                "counting-reference and trivial-datum inputs conditionally "
                "realize the uniform transition kernel; scale and initial-law "
                "controls show this does not select the complete path reference. "
                "Two exact multiplier tilts have distinct mean actions, and "
                "the matching quadratic gives one unique positive parameter "
                "at the supplied empirical target. A "
                "concave three-record history is stationary but not minimal and "
                "is not a positive-Gibbs mode against one variation, giving only "
                "a scoped mode/minimizer control"
            ),
            "observed_counterpart": (
                "Gibbs path selection, least-action limits, discrete "
                "Euler--Lagrange equations, and Noether currents"
            ),
            "match": (
                "exact conditional helpers plus finite/real and real-enrichment "
                "non-identifiability boundaries; physical composition is not constructed"
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
                "SourceReferenceSelection": [
                    "heatBath_counting_trivial_eq_uniform",
                    "reference_realized_under_counting_trivial_inputs",
                    "nontrivial_datum_not_invariant",
                    "noncounting_reference_not_invariant",
                    "heatBath_scaledCounting_trivial_eq_uniform",
                    "uniform_transition_does_not_determine_path_reference",
                    "committed_tilts_have_distinct_mean_actions",
                ],
                "SourceHistoryPacket": [
                    "sourceMatchQuad_strictMonoOn_pos",
                    "sourcePositiveMeanMatch_unique",
                    "sourceMatchingPositiveParameter_existsUnique",
                ],
                "LogTransitionAction": [
                    "bare_log_action_multiplier_unique_of_nonconstant",
                ],
                "StationarySaddleCoverage": [
                    "stationaryMaximumHistory_stationary",
                    "stationaryMaximumHistory_not_minimal",
                    "gibbs_prefers_nonstationary",
                ],
            },
            "lean_receipts": _lean_receipt(
                "PathGibbs",
                "DiscreteEulerLagrange",
                "DiscreteNoether",
                "FiniteHistoryBridge",
                "RealizedHistoryLegendreNoGo",
                "SourceReferenceSelection",
                "SourceHistoryPacket",
                "LogTransitionAction",
                "StationarySaddleCoverage",
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
                    "SourceReferenceSelection": (
                        "heatBath_counting_trivial_eq_uniform",
                        "reference_realized_under_counting_trivial_inputs",
                        "nontrivial_datum_not_invariant",
                        "noncounting_reference_not_invariant",
                        "heatBath_scaledCounting_trivial_eq_uniform",
                        "uniform_transition_does_not_determine_path_reference",
                        "committed_tilts_have_distinct_mean_actions",
                    ),
                    "SourceHistoryPacket": (
                        "sourceMatchQuad_strictMonoOn_pos",
                        "sourcePositiveMeanMatch_unique",
                        "sourceMatchingPositiveParameter_existsUnique",
                    ),
                    "LogTransitionAction": (
                        "bare_log_action_multiplier_unique_of_nonconstant",
                    ),
                    "StationarySaddleCoverage": (
                        "stationaryMaximumHistory_stationary",
                        "stationaryMaximumHistory_not_minimal",
                        "gibbs_prefers_nonstationary",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the finite source packet and its log-transition corner action "
                "are attained, while the complete reference is not source-"
                "selected. Conditional repair of supplied trivial data under a "
                "supplied counting reference realizes only its uniform transition "
                "kernel; constant rescaling leaves that kernel unchanged, and "
                "distinct initial laws give distinct complete path references. The "
                "matching parameter exists uniquely at the supplied empirical "
                "target, but the constraint observable and level remain declared. A typed "
                "finite-to-real transfer exists only under an undercut receipt. "
                "The complete binary history law does not select the "
                "off-alphabet curvature of a regular Legendre system: the "
                "displayed family is constructed, not source-produced. The "
                "concave control concerns modes/minimizers only: constrained "
                "saddles, complex or signed stationary phase, and refinement "
                "routes are not excluded. Closed #731 records only the attained "
                "declared-bundle composition. Scientific owner #739 records source "
                "selection of the reference, real enrichment, stationary-phase "
                "mechanism, and source clock; #730 records the amplitude/interference "
                "interface, while physical fields and observable currents are not "
                "supplied here. This "
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
    rows.extend(
        [
            {
                "id": "local_face_maxwell_static_composition",
                "statement": (
                    "The committed oriented twenty-face incidence has rank 19 "
                    "and kernel exactly the port gradients. Its local seam "
                    "Hessian has five nonzero entries per row; gauge invariance "
                    "and stationary solvability are equivalent to conserved "
                    "seam current. A separate scalar action has the neutral "
                    "Green potential as its minimum modulo constants, and the "
                    "direct-product theorem retains port load and seam current "
                    "as distinct types"
                ),
                "observed_counterpart": (
                    "static Maxwell-shaped local curvature and Coulomb sectors"
                ),
                "match": "exact conditional finite structural composition",
                "lean_declarations": {
                    "LocalFaceMaxwellAction": [
                        "face_port_incidence_product_zero",
                        "ker_faceCurvature_eq_gradient",
                        "localKineticZ_five_per_row",
                        "localSourcedAction_gauge_invariant_iff",
                        "localStationary_solvable_iff",
                        "greenPotential_stationary",
                        "staticAction_global_minimum",
                    ]
                },
                "lean_receipts": _lean_receipt(
                    "LocalFaceMaxwellAction",
                    declarations={
                        "LocalFaceMaxwellAction": (
                            "face_port_incidence_product_zero",
                            "ker_faceCurvature_eq_gradient",
                            "localKineticZ_five_per_row",
                            "localSourcedAction_gauge_invariant_iff",
                            "localStationary_solvable_iff",
                            "greenPotential_stationary",
                            "staticAction_global_minimum",
                        )
                    },
                ),
                "artifact_ref": (
                    "code/electromagnetism/runtime/"
                    "local_face_maxwell_action_receipt.json"
                ),
                "hypothesis_boundary": (
                    "finite static carrier only: no temporal electromagnetic "
                    "evolution, charge-current continuity map, physical source, "
                    "spacetime/continuum attachment, or readout is supplied"
                ),
                "paper_ref": "Screen microphysics paper, local face action",
            },
            _triple_observer_carrier_row(),
            {
                "id": "cofinal_spectral_tail_four_law_composition",
                "statement": (
                    "One supplied ambient-cofinal spectral family with no "
                    "terminal regulator and strict finite-carrier growth carries "
                    "uniform and cofinal-tail concentration. A supplied "
                    "calibrated-stage energy identification yields exact "
                    "off-minimum-mass equality and the same envelope; the family "
                    "also has infinitely many stages and the directly "
                    "constructed finite four-law core"
                ),
                "observed_counterpart": (
                    "low-temperature control on a genuine regulator family"
                ),
                "match": "exact conditional thermodynamic structure",
                "lean_declarations": {
                    "CofinalSpectralTailFamily": [
                        "cofinal_tail_concentration",
                        "cofinalLadder_card_unbounded",
                        "constantCarrier_cofinalFamily_isEmpty",
                        "uniformGap_cofinalFamily_isEmpty",
                        "fourLaws_composed_cofinal",
                    ]
                },
                "lean_receipts": _lean_receipt(
                    "CofinalSpectralTailFamily",
                    declarations={
                        "CofinalSpectralTailFamily": (
                            "cofinal_tail_concentration",
                            "cofinalLadder_card_unbounded",
                            "constantCarrier_cofinalFamily_isEmpty",
                            "uniformGap_cofinalFamily_isEmpty",
                            "fourLaws_composed_cofinal",
                        )
                    },
                ),
                "hypothesis_boundary": (
                    "strict carrier growth is one sufficient refinement route; "
                    "the ambient regulator, energy, inverse temperature, repair, "
                    "source production, and physical continuum meaning remain "
                    "declared rather than derived"
                ),
                "paper_ref": "Observers paper, cofinal spectral-tail theorem",
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
            "scientific_owner_issues": [736],
        }
    ]


ALPHA_CLOSURE_MODES = (
    (
        "alpha_inv_closure_root",
        "thomson_structured_running",
        "root closure map, unified gauge width absent",
    ),
    (
        "alpha_inv_closure_gauge_width",
        "thomson_structured_running_plus_gauge_width",
        "closure map with the finite-screen unified gauge width on the inverse coupling",
    ),
)


def _alpha_closure_rows(
    endpoint: dict[str, Any],
    contraction: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare-only rows for the certified fixed points of the declared closure map.

    The certificate supplies interval-arithmetic existence and uniqueness for
    each declared map at its stated cutoffs. The rows carry the enclosure, the
    distance to the CODATA reference, and the owning issue. They consume no
    measured value on any solve path.
    """
    co = endpoint["compare_only"]
    measured = float(co["codata_alpha_inv"])
    if (
        contraction.get("exact_alpha_promoted") is not False
        or contraction.get("promotion_allowed") is not False
    ):
        raise SystemExit(
            "the P contraction certificate has left its compare-only boundary"
        )
    policy = contraction["consumer_policy"]
    if (
        policy["may_feed_compare_or_audit_surfaces"] is not True
        or policy["may_feed_live_particle_predictions"] is not False
        or policy["hidden_external_alpha_allowed"] is not False
        or policy["default_thomson_endpoint_allowed"] is not False
    ):
        raise SystemExit("the P contraction certificate consumer policy has changed")
    rows: list[dict[str, Any]] = []
    for row_id, mode, description in ALPHA_CLOSURE_MODES:
        block = contraction["modes"][mode]
        banach = block["banach"]
        if not (
            banach["existence"]
            and banach["contraction"]
            and banach["uniqueness_in_interval"]
            and banach["g_maps_interval_into_interior"]
        ):
            raise SystemExit(
                f"the {mode} fixed point lacks certified existence and uniqueness"
            )
        enclosure = block["certified_enclosure"]["alpha_inv"]
        point = block["fixed_point_point_estimate_display_only"]
        central = float(point["alpha_inv"])
        rows.append(
            {
                "id": row_id,
                "mode": mode,
                "closure_map_description": description,
                "map_definition": block["map_definition"],
                "value_central": central,
                "value_interval": [
                    float(enclosure["lo"]),
                    float(enclosure["hi"]),
                ],
                "value_interval_decimal": [
                    str(enclosure["lo"]),
                    str(enclosure["hi"]),
                ],
                "enclosure_width_decimal": str(enclosure["width"]),
                "p_central": float(point["P"]),
                "measured": measured,
                "measured_source": co["codata_source"],
                "deviation_inv_alpha": central - measured,
                "relative_deviation": (central - measured) / measured,
                "row_class": "certified_interval_fixed_point_of_declared_closure_map",
                "tier": "T2_conditional",
                "lipschitz_bound": float(banach["lipschitz_bound"]),
                "su2_cutoff": block["su2_cutoff"],
                "su3_cutoff": block["su3_cutoff"],
                "artifact_refs": [
                    _rel("p_interval_contraction"),
                    _rel("endpoint"),
                ],
                "scientific_owner_issues": [736],
            }
        )
    return rows


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
                "distance to the standard on-shell reference deficit is an "
                "unfixed scheme term. Scientific owner #736 records that "
                "missing source obligation. The lepton "
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
            "scientific_owner_issues": [736],
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
            "scientific_owner_issues": [736],
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
                "No absolute quark mass is emitted. The two-modulus spread "
                "fiber (R>0)^2, obtained by granting a candidate-only ordered "
                "shape law, survives the 2026-07 certified structure set "
                "(matter receipt #314, port receipt #566, twelve frozen "
                "selector candidates, with input hashes pinned at emission), so the "
                "six source-only absolute masses are non-identifiable from that "
                "set, by an explicit rescaling symmetry of the registered data. "
                "Live cuts: a Yukawa-typed source equation through the physical "
                "bindings or the family attachment, a new selector under the "
                "frozen discipline, or the conditional Higgs/top criticality "
                "coordinate"
            ),
            "fork": obstruction["fork"],
            "fiber_cut_detected": obstruction["fiber_cut_detected"],
            "tier": obstruction["claim_tier"],
            "artifact_ref": _rel("fiber_obstruction"),
            "scientific_owner_issues": [736],
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


def _hadron_rows(
    payload: dict[str, Any],
    standby: dict[str, Any],
    lambda_scale: dict[str, Any],
    nucleon: dict[str, Any],
) -> list[dict[str, Any]]:
    integral = payload["integral"]
    norm = integral["normalization"]
    if (
        lambda_scale.get("promotion_allowed") is not False
        or lambda_scale.get("checks_pass") is not True
        or nucleon.get("promotion_allowed") is not False
        or nucleon.get("checks_pass") is not True
    ):
        raise SystemExit(
            "a hadronic compare-only receipt has left its declared boundary"
        )
    lambda_central = float(lambda_scale["central"]["lambda3_gev"])
    lambda_interval = [float(v) for v in lambda_scale["lambda3_interval_gev"]]
    lambda_published = float(
        lambda_scale["machinery_validation"]["published_compare"]["lambda3_gev"]
    )
    nucleon_compare = nucleon["compare_only"]
    nucleon_central = float(nucleon["prediction"]["m_nucleon_gev_display"])
    nucleon_interval = [
        float(v) for v in nucleon["prediction"]["m_nucleon_interval_gev_display"]
    ]
    nucleon_measured = float(nucleon_compare["measured_m_proton_gev"])
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
                "suppressed. The resource-deferred QCD backend supplies no "
                "source-only result; "
                "scientific owner #736 records the bounded particle-output obligation, "
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
        {
            "id": "lambda_qcd_transmutation_scale",
            "value_central_gev": lambda_central,
            "value_interval_gev": lambda_interval,
            "published_central_gev": lambda_published,
            "published_source": (
                "published perturbative determinations of the three-flavor "
                "scale, compare-only"
            ),
            "relative_deviation": (lambda_central - lambda_published)
            / lambda_published,
            "row_class": lambda_scale["row_class"],
            "tier": "T2_conditional",
            "loop_order": lambda_scale["loop_order"],
            "declared_external_inputs": (
                "threshold locations are declared external quark scheme masses, "
                "and the interval is the swept threshold bracket"
            ),
            "claim_boundary": lambda_scale["claim_boundary"],
            "artifact_refs": [_rel("lambda_transmutation")],
            "scientific_owner_issues": [736],
        },
        {
            "id": "nucleon_mass_external_ratio",
            "value_central_gev": nucleon_central,
            "value_interval_gev": nucleon_interval,
            "measured_gev": nucleon_measured,
            "measured_source": "proton mass, compare-only",
            "relative_deviation": float(
                nucleon_compare["central_relative_difference"]
            ),
            "interval_contains_measured": bool(
                nucleon_compare["interval_contains_measured"]
            ),
            "row_class": nucleon["row_class"],
            "tier": "T2_conditional",
            "external_theory_factor": float(
                nucleon["external_theory_factor"]["R_nucleon_over_lambda3"]
            ),
            "external_theory_uncertainty": float(
                nucleon["external_theory_factor"]["uncertainty"]
            ),
            "claim_boundary": nucleon["claim_boundary"],
            "artifact_refs": [
                _rel("nucleon_external_ratio"),
                _rel("lambda_transmutation"),
            ],
            "scientific_owner_issues": [736],
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
                "position, scale, refinement, and gluing are not constructed."
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
                "finite scale, coherent frame, readout, and comparison are not constructed."
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
                "measurement. The balance premise is Koide's published "
                "relation, which the face circulant of the carrier holds as a "
                "finite structure, so this is a conditional postdiction with a "
                "frozen rejection rule and a declared premise ancestry."
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
                "unfixed scheme term; scientific owner #736 records that missing source "
                "requirement. The lepton "
                "scale is localized only under that recorded accounting "
                "packet. A "
                "source-emitted bridge value is a falsification target: the "
                "closure value would satisfy the conditional lane, while a "
                "value outside the interval refutes the declared decomposition."
            ),
        },
    ]


def _seam_maxwell_continuum_row() -> dict[str, Any]:
    """Derive a structural row only after the separate continuum replay."""
    path = CODE / "electromagnetism" / "verify_seam_maxwell_continuum.py"
    spec = importlib.util.spec_from_file_location("seam_continuum_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent continuum verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt = verifier.load_receipt(verifier.OUTPUT)
    count = verifier.verify(receipt)
    declarations = ("cosine_quartic_upper", "unit_seam_fourth_moment_eq",
                    "exact_symbol_quadratic_error", "exact_frequency_error",
                    "cosine_propagator_error", "sine_propagator_error")
    return {
        "id": "seam_maxwell_continuum",
        "statement": (
            "The complete seam symbol has global quadratic consistency error "
            "a^2|k|^4/20 and frequency error a^2|k|^3/20. On a supplied common "
            "Euclidean domain, its declared reversible curl pair assembles to "
            "real L2 fields converging strongly to local Maxwell evolution, "
            "uniformly on bounded times, with an O(a^2) H3-to-L2 rate and "
            "same-forcing Duhamel bound. Prescribed charge continuity propagates "
            "initial Gauss constraints. Lean checks scalar estimates; the L2 "
            "assembly and source/limit statements are proved in the paper"
        ),
        "observed_counterpart": "Classical continuum Maxwell equations; no measured data",
        "match": "analytic continuum field bridge on a supplied domain and dynamics",
        "lean_declarations": {"SeamMaxwellContinuum": list(declarations)},
        "lean_receipts": _lean_receipt("SeamMaxwellContinuum", declarations={
            "SeamMaxwellContinuum": declarations}),
        "artifact_refs": [
            "code/electromagnetism/runtime/seam_maxwell_continuum_receipt.json",
            "code/electromagnetism/seam_maxwell_continuum.py",
            "code/electromagnetism/verify_seam_maxwell_continuum.py",
            "paper/screen_microphysics_and_observer_synchronization.tex"],
        "universal_bounds": receipt["universal_bounds"],
        "independently_replayed_control_cases": count,
        "hypothesis_boundary": (
            "Supplied R3/Lebesgue field domain, complete seam symbol, reversible "
            "curl-pair evolution, common time and same data/forcing. No source-log "
            "refinement, operational clock, physical current, finite-scale local "
            "cone, decoded observer-field/action join or laboratory identification "
            "is constructed. The result discharges no physical premise and does "
            "not arm FZ-12. Owners #728 and #754 retain the physical maps"
        ),
        "paper_ref": "screen microphysics, Continuum Maxwell limit of the complete seam symbol",
    }


def _serial_maxwell_readout_row() -> dict[str, Any]:
    """Structural instrument evidence, with independent replay before inclusion."""
    path = CODE / "electromagnetism" / "verify_serial_maxwell_readout.py"
    spec = importlib.util.spec_from_file_location("serial_maxwell_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent serial Maxwell verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    result = verifier.verify(verifier.load())
    declarations = ("cycle_restores", "serial_restores", "response_after_prefix",
                    "decode_serial", "feedback_error", "decoded_fields",
                    "decoded_coupled_action", "decoded_joint_field_stationary")
    return {
        "id": "serial_maxwell_readout",
        "statement": (
            "Exact classical baseline/response feedback restores every pair and "
            "every finite overlapping probe word. Serial decoding recovers all "
            "42 typed potential coordinates and preserves electric/magnetic "
            "fields and the finite neutral-pair action. Decoded initial records "
            "and the charged path current advance a third rational slice"
        ),
        "observed_counterpart": "Finite Maxwell field and charged-path equations; no measured data",
        "match": "exact classical software instrument and decoded finite action join",
        "lean_declarations": {"SerialMaxwellReadout": list(declarations)},
        "lean_receipts": _lean_receipt("SerialMaxwellReadout", declarations={
            "SerialMaxwellReadout": declarations}),
        "artifact_refs": [
            "code/electromagnetism/runtime/serial_maxwell_readout_receipt.json",
            "code/electromagnetism/serial_maxwell_readout.py",
            "code/electromagnetism/verify_serial_maxwell_readout.py",
            "paper/observers_are_all_you_need.tex"],
        "events_per_gauge_execution": result["events"],
        "feedback_cycles_per_execution": result["cycles"],
        "independently_checked_field_derivatives": result["field_variations"],
        "exhaustive_path_replacements_per_execution": result["path_variations"],
        "hypothesis_boundary": (
            "Exact classical memory and writable ports, declared potential typing, "
            "initial slices, paths, finite action, h=1/2 and Lorentz unit tau=3. "
            "No source-selected physical clock, long noisy stability, quantum "
            "copying, spatial refinement into the continuum family or laboratory "
            "identification. No physical prediction is armed or premise discharged"
        ),
        "paper_ref": "observers synthesis, Serial feedback readout and finite Maxwell action",
    }


def _cone_whitney_bridge_row() -> dict[str, Any]:
    """Finite mathematical bridge; inclusion requires independent full replay."""
    path = CODE / "electromagnetism" / "verify_cone_whitney_bridge.py"
    spec = importlib.util.spec_from_file_location("cone_whitney_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent cone Whitney verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt = verifier.load()
    result = verifier.verify(receipt)
    declarations = {
        "ConeCochainBridge": (
            "extendOne_gradient", "extendOne_gauge_covariant",
            "coneCurl_extendOne", "coneCurl_gauge_invariant",
            "coneCurl_boundary_trace", "extended_curvature_closed",
            "extendTwo_curvature", "extendTwo_closed", "extendTwo_radial_unique",
            "closed_extension_iff_zero_flux", "constant_unit_flux_has_no_closed_extension",
            "zero_apex_constant_gauge_fails",
        ),
        "WhitneyTimeBridge": (
            "affine_quadratic_polynomial", "magnetic_window_action_defect",
            "magnetic_defect_divided_difference", "polarized_trapezoid_defect",
            "polarized_defect_divided_difference", "temporal_hat_mean",
            "temporal_hat_residual", "ampere_covector_decomposition",
            "scalar_nonzero_residual_control",
        ),
    }
    return {
        "id": "cone_whitney_bridge",
        "statement": (
            "Exact gauge-covariant cone cochain extension preserves decoded "
            "boundary fields, Faraday and zero magnetic divergence; closed "
            "two-form extension exists exactly at zero total boundary flux. "
            "Whitney interpolation supplies piecewise polynomial fields on a "
            "declared solid and time slabs. Exact temporal coefficient identities "
            "and numerical spatial action/residual controls distinguish the "
            "interpolated action from the counting action; all 68 finite-element "
            "free-coordinate derivatives are checked for each of two gauge histories"
        ),
        "observed_counterpart": (
            "Mathematical Maxwell cochain and variational identities; no observed data"
        ),
        "match": (
            "exact finite cochain and coefficient theorems; numerical geometric "
            "action/residual audit; no physical postdiction"
        ),
        "physical_comparison_status": "NOT_EVALUABLE",
        "observed_postdiction": False,
        "continuum_convergence_established": False,
        "lean_declarations": {key: list(value) for key, value in declarations.items()},
        "lean_receipts": _lean_receipt(*declarations, declarations=declarations),
        "artifact_refs": [
            "code/electromagnetism/runtime/cone_whitney_bridge_receipt.json",
            "code/electromagnetism/cone_whitney_bridge.py",
            "code/electromagnetism/verify_cone_whitney_bridge.py",
            "code/electromagnetism/test_cone_whitney_bridge.py",
            "paper/observers_are_all_you_need.tex",
            "paper/tex_fragments/CONE_WHITNEY_BRIDGE.tex",
        ],
        "certificate_scope": receipt["scope"],
        "numeric_policy": receipt["numeric_policy"],
        "tetrahedra": result["tetrahedra"],
        "independently_replayed_gauge_histories": result["gauge_histories"],
        "independently_checked_field_derivatives_per_history": result[
            "field_variations_per_history"],
        "full_volume_stationarity": receipt["controls"]["full_volume_stationarity"],
        "hypothesis_boundary": (
            "Supplied Euclidean embedding, apex, nonlocal radial extension, "
            "Whitney Hodge pairings, unit constitutive coefficients, uniform "
            "h=1/2 interpolation and the same finite source/clock contract. "
            "The Lorentz clock action is not a volume field term, and no physical "
            "clock or SI calibration is selected. Sources are finite-element "
            "covector impulses with the declared endpoint convention, not "
            "identified physical charge/current densities. Quadratic identities "
            "require symmetric pairings; the h-squared interior defect scaling "
            "requires uniform temporal regularity for a convergence reading, "
            "and the action has an explicit endpoint defect. Spatial Gram and "
            "action/residual values are float64 controls, not interval "
            "certificates. Whitney fields have tangential/normal conformity, "
            "not global vector continuity; full sourced weak residuals retain "
            "face jumps and temporal impulses and do not vanish here. Checking "
            "all finite-element coordinates does not establish arbitrary-test "
            "continuum stationarity. No refinement family, source-selected "
            "physical geometry, continuum convergence, measured postdiction, "
            "new prediction or premise discharge is supplied"
        ),
        "paper_ref": (
            "observers synthesis, Gauge-covariant cone reconstruction and metric/time action defects"
        ),
    }


def _whitney_dynamics_rows() -> list[dict[str, Any]]:
    """Same-geometry mathematical constructions, never observed postdictions."""
    path = CODE / "electromagnetism" / "verify_whitney_maxwell_dynamics.py"
    spec = importlib.util.spec_from_file_location("whitney_dynamics_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent Whitney dynamics verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt = verifier.load()
    verified = verifier.verify(receipt)
    specs = [
        ("whitney_maxwell_dynamics", "WhitneyMaxwellDynamics", (
            "two_slab_action_expansion", "scalar_potential_action_expansion",
            "temporal_gauge_ampere_iff_recurrence", "recurrence_next_unique",
            "gauss_residual_preserved", "pairEnergy_work", "modal_solution_bound",
            "criticalMode_unbounded_difference", "supercritical_unbounded_solution",
            "two_slab_endpoint_resonance"),
         "Full-metric prism action yields a unique next finite volume potential, preserves the initial Gauss residual under source continuity, and obeys an exact work identity. Positive-mode boundedness has the sharp window 0<h^2 lambda<12; stable initial-value dynamics can have a singular fixed-endpoint problem at h^2 lambda=3. Lean derives the exact Q(sqrt5) Whitney mass forms from the committed supplied coordinates and signed incidence, proves their assembled positivity and K<24M, and supplies the existing stability consumer. Certificate payloads are bound to a kernel witness, while independent numerical quadrature checks the signed assembled forms. A separately Gauss-projected numerical history has two 805-event gauge representatives, all 68 action derivatives verified, exact serial readback restoration and a separately identified 64-slab numerical continuation.",
         "exact finite action/stability theorems and algebraic matrix certificate; numerical trajectory and authenticated classical records; no observed postdiction",
         "Supplied Euclidean cone, constitutive interpretation, h=1/2, changed initial electric field, prescribed impulse source covectors, dense software evolution and exact writable classical ports. Exact mass positivity and the stability bound are derived for that geometry, not physical source selection. The local signed-coordinate certificate policy and global real face-order policy are separately proved; their bundle establishes no common local/global transport. Lean-free checks authenticate the committed witness payload and source closure plus numerical agreement; full kernel and consumer replay requires Lean. Float64 trajectory with 1e-9 comparison tolerances is not an interval trajectory certificate. Only three slices per gauge are instrumented. Finite-element stationarity is not arbitrary-test continuum stationarity; no physical clock, spatial refinement, laboratory evidence or new prediction.",
         "WHITNEY_MAXWELL_DYNAMICS.tex"),
        ("whitney_radiative_quantum", "WhitneyQuantumBridge", (
            "reconstruct_surjective", "same_action_normal_modes", "annihilation_creation",
            "canonical_position_momentum", "same_modes_quantized_energy",
            "quantumHamiltonian_monomial", "quantumHamiltonian_position", "quantumHamiltonian_momentum"),
         "The same geometric mass/stiffness action has a thirty-dimensional homogeneous transverse sector. Lean derives the concrete positive mass forms from the committed supplied cone geometry and constructs a complete mass-orthonormal positive spectral frame that pulls its source-free continuous-time action back to ordinary oscillators. The derived frame supplies the existing classical-action and polynomial-quantum consumers. Canonical bosonic quantization with supplied positive hbar gives exact polynomial CCR, occupation energies and Heisenberg equations. The paper constructs the factorial-weight Hilbert completion l2(N^30), self-adjoint diagonal Hamiltonian and strongly continuous unitary flow. A separate exact Legendre map identifies the stable prism-step modified oscillator Hamiltonian and a second-order classical temporal limit at fixed bounded spectrum.",
         "exact conditional algebraic quantum realization; analytic Hilbert/domain and temporal-limit proof in the paper; no observed postdiction",
         "Canonical quantization and hbar are imported. Lean constructs the complete positive normal frame from the certified mass forms of the supplied geometry. The local exact signed-coordinate certificate transport and global real face-order transport have separate naturality proofs; no common local/global policy is established. Hilbert completion, operator closures and temporal convergence are paper proofs. This is the zero-charge source-free radiative sector, not quantization of the prescribed-charge episode or interacting matter action. sqrt(lambda) and finite-step theta/h are distinct; no uniform operator-norm quantum error, spatial continuum limit, selected Born statistics or physical clock is supplied.",
         "WHITNEY_MAXWELL_DYNAMICS.tex"),
        ("whitney_charged_matter", "WhitneyChargedMatter", (
            "pathPhase_gauge", "interpolate_gauge", "interpolate_gauge_norm",
            "interpolate_vertex", "interpolate_face_trace"),
         "Real Whitney edge integrals dress nodal complex scalar fields by straight-segment U(1) phases, giving exact nodal interpolation, gauge covariance and matching face traces. Restriction of a supplied nonnegative-potential scalar-QED action to these fields and the same volume Maxwell variables gives a joint gauge-invariant action. The paper proves its Noether/Gauss identity, positive temporal-gauge velocity Hessian and global finite-dimensional classical evolution, and constructs nonzero locally charged but globally neutral exact initial data.",
         "exact finite interpolation gauge algebra; analytic coupled-action and global-existence proof; no observed postdiction",
         "Supplied scalar species, charge, mass, nonnegative quartic coefficient, Euclidean geometry, time and exact pulled-back scalar-QED Lagrangian. Real unwrapped edge integrals are needed. Differentiate the gauge-dependent scalar basis in all field variations; a naive projected scalar current omits terms. All thirteen scalar variations require total neutrality in the declared boundary convention. This algebra/global-existence result alone supplies no executed charged-matter episode, spatial error estimate or interacting quantum completion; those distinct constructions have separate rows. No source-selected matter content or Standard Model identification is supplied.",
         "WHITNEY_CHARGED_MATTER.tex"),
    ]
    rows = []
    for name, module, declarations, statement, match, boundary, fragment in specs:
        row = {"id": name, "statement": statement,
            "observed_counterpart": "Mathematical classical/quantum field structures; no observed data",
            "match": match, "physical_comparison_status": "NOT_EVALUABLE",
            "observed_postdiction": False, "continuum_convergence_established": False,
            "lean_declarations": {module: list(declarations)},
            "lean_receipts": _lean_receipt(module, declarations={module: declarations}),
            "artifact_refs": ["paper/observers_are_all_you_need.tex", "paper/tex_fragments/"+fragment],
            "hypothesis_boundary": boundary,
            "paper_ref": "observers synthesis, " + name.replace("_", " ")}
        if name in {"whitney_maxwell_dynamics", "whitney_radiative_quantum"}:
            row["artifact_refs"] += ["Lean/Screen/" + source + ".lean" for source in (
                "WhitneyFiniteCertificate", "WhitneyAlgebraicLDL",
                "WhitneySourceGeometry", "WhitneySourceAssembly",
                "WhitneyCertifiedConsumers", "WhitneySourceNaturality",
                "WhitneyGeneratedCertificate", "WhitneyOmittedCellCounterexample")]
        if name == "whitney_maxwell_dynamics":
            row.update({"certificate_scope": receipt["scope"], "numeric_policy": receipt["numeric_policy"],
                "independent_verifier_result": {key: value for key, value in verified.items()
                    if not key.endswith("_max_abs")},
                "finite_element_stationarity": "numerically verified for all 68 free coefficients per gauge history",
                "arbitrary_test_continuum_stationarity": False})
            row["artifact_refs"] += ["code/electromagnetism/"+x for x in (
                "runtime/whitney_maxwell_dynamics_receipt.json", "whitney_maxwell_dynamics.py",
                "verify_whitney_maxwell_dynamics.py", "test_whitney_maxwell_dynamics.py")]
        if name == "whitney_charged_matter":
            row["artifact_refs"].append("code/electromagnetism/test_whitney_charged_matter.py")
        rows.append(row)
    return rows


def _verify_whitney_coupled_parent(stem: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load and independently replay a current parent on every invocation."""
    if stem not in {"spatial_consistency", "charged_dynamics", "quantum_state"}:
        raise SystemExit("unknown coupled Whitney evidence parent")
    path = CODE / "electromagnetism" / f"verify_whitney_{stem}.py"
    spec = importlib.util.spec_from_file_location(f"whitney_{stem}_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"missing independent Whitney {stem} verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt = verifier.load()
    return receipt, verifier.verify(receipt)


def _whitney_coupled_rows() -> list[dict[str, Any]]:
    """Controlled action consistency, finite quantization and classical execution.

    The analytic theorems remain paper proofs. Numeric parents are replayed
    freshly, but only their stable categories, counts and exact rational
    strings enter this ledger.
    No floating-point replay residual becomes a platform-dependent projection.
    """
    spatial_receipt, spatial_verified = _verify_whitney_coupled_parent("spatial_consistency")
    charged_receipt, charged_verified = _verify_whitney_coupled_parent("charged_dynamics")
    quantum_receipt, quantum_verified = _verify_whitney_coupled_parent("quantum_state")
    for key in ("continuum_trajectory_claimed", "observer_history_claimed", "uniform_bound_certified_by_numerics"):
        if spatial_verified[key] is not False:
            raise SystemExit("coupled Whitney spatial evidence cannot promote " + key)
    if charged_verified["accepted"] is not True:
        raise SystemExit("coupled Whitney charged execution was not accepted")
    for key in ("observer_history", "quantum_state", "physical_continuum"):
        if charged_verified[key] is not False:
            raise SystemExit("coupled Whitney execution cannot promote " + key)
    for key in ("accepted", "initial_state_constructed", "initial_observables_provided"):
        if quantum_verified[key] is not True:
            raise SystemExit("coupled Whitney quantum initial evidence missing: " + key)
    for key in ("quantum_time_history", "observer_history", "physical_comparison",
                "analytic_operator_domain_proved_by_numeric_replay"):
        if quantum_verified[key] is not False:
            raise SystemExit("coupled Whitney quantum initial evidence cannot promote " + key)
    if (type(quantum_verified["real_configuration_dimension"]) is not int
            or quantum_verified["real_configuration_dimension"] != 56):
        raise SystemExit("coupled Whitney quantum configuration dimension must be exactly 56")
    if (type(quantum_verified["scope"]) is not str
            or quantum_verified["scope"] != quantum_receipt["scope"]):
        raise SystemExit("coupled Whitney quantum certificate scope mismatch")
    exact_quantum_keys = ("gaussian_sigma", "gaussian_norm_squared",
                          "matter_l2_coefficient", "matter_l4_coefficient")
    volume = quantum_verified["volume_in_Qsqrt5"]
    if type(volume) is not list or len(volume) != 2:
        raise SystemExit("coupled Whitney quantum exact quadratic-field volume missing")
    for value in [*(quantum_verified[key] for key in exact_quantum_keys), *volume]:
        try:
            canonical = type(value) is str and str(Fraction(value)) == value
        except (ValueError, ZeroDivisionError):
            canonical = False
        if not canonical:
            raise SystemExit("coupled Whitney quantum projection requires canonical exact rational strings")
    if (quantum_verified["gaussian_norm_squared"] != "1"
            or Fraction(quantum_verified["gaussian_sigma"]) <= 0):
        raise SystemExit("coupled Whitney quantum Gaussian must be normalized with positive width")
    spatial_keys = (
        "scope", "refinement_levels", "tetrahedra", "vertices",
        "finite_algebra_lean_declarations", "independent_derivative_stencil_points",
        "continuum_trajectory_claimed", "observer_history_claimed", "uniform_bound_certified_by_numerics",
    )
    charged_keys = (
        "accepted", "samples", "full_equations_per_sample", "gauss_equations_per_sample",
        "symmetry", "trajectory_error_status", "observer_history", "quantum_state", "physical_continuum",
    )
    quantum_keys = (
        "scope", "accepted", "initial_state_constructed", "initial_observables_provided",
        "quantum_time_history", "observer_history", "physical_comparison",
        "analytic_operator_domain_proved_by_numeric_replay",
        "real_configuration_dimension", *exact_quantum_keys, "volume_in_Qsqrt5",
    )
    specs = (
        ("whitney_spatial_consistency", "WHITNEY_SPATIAL_CONSISTENCY.tex", "thm:whitney-spatial-consistency",
         "On conforming shape-regular refinements of the same cone, the supplied dressed scalar/Whitney action and its full real first variation approximate the continuum scalar-electrodynamics action at O(delta), uniformly on bounded piecewise W^{1,infinity}_t W^{2,infinity}_x windows with compatible traces. The analytic proof controls the potential-dependent scalar interpolation and all dressing derivatives. Five finite Lean identities establish cancellation/error algebra; a manufactured-field refinement packet independently checks the implementation on three meshes.",
         "analytic smooth-window action/first-variation consistency; finite Lean algebra and numerical refinement checks; no observed postdiction",
         "Supplied geometry, smooth bounded fields, conforming shape-regular refinements, action time parameter and scalar-electrodynamics density. The uniform estimate is an analytic paper theorem, not a Lean convergence theorem or a bound certified by finite samples. Exact edge integration and the Whitney interior-path definition are required. No nonlinear trajectory convergence, arbitrary-H1 nodal estimate, observer source selection, calibrated physical clock or empirical comparison is established."),
        ("whitney_interacting_quantum", "WHITNEY_INTERACTING_QUANTUM.tex", "thm:whitney-interacting-quantum",
         "For the same interacting finite charged action with nonzero charge, the actual positive kinetic metric gives a smooth Schur metric on the global mean-zero-gauge Coulomb slice R^30 x C^13. Uniform fixed-mesh nodal coercivity gives two-sided polynomial metric bounds and proves that this configuration metric is complete. The declared Laplace--Beltrami operator with nonnegative potential is essentially self-adjoint on compactly supported smooth functions; its unique self-adjoint closure agrees with the Friedrichs realization for the chosen quantization measure and ordering. The residual U(1)-invariant subspace reduces the Hamiltonian and has a dense invariant operator core as well as a form core. An explicit metric-density-corrected Gaussian is exactly normalized, neutral and in the Hamiltonian domain. Independent exact initial-state moment identities provide mathematical observables; no quantum time history is computed.",
         "analytic gauge-reduction, metric-completeness and essential-self-adjointness proof; exact normalized initial state and initial observables; no observed postdiction",
         "Supplied fixed cone, scalar species/action, nonzero charge, nonnegative mass-squared/quartic coupling, hbar, quantization measure, Laplace--Beltrami ordering and Gaussian width. The Schur complement uses the full coupled kinetic metric, which includes every block beside the Maxwell one. Essential self-adjointness removes extension ambiguity for this declared operator, not the quantization choices themselves. The metric constants depend on the fixed mesh and couplings and are not uniform continuum-refinement estimates. The analytic proof is not formalized in Lean. The chosen Gaussian is not a derived vacuum or physical preparation. No unique quantization, equivalence to quantization before reduction, continuum interacting QFT, computed quantum-state history, physical Born statistics or laboratory identification is supplied."),
        ("whitney_charged_execution", "WHITNEY_CHARGED_EXECUTION.tex", "prop:whitney-charged-symmetric-lift",
         "A five-real-coordinate icosahedrally invariant sector of the same charged action evolves from the nonzero neutral initial data with all 68 full temporal-gauge Euler equations and all 13 Gauss equations checked at 81 samples. The analytic finite-group argument lifts the restricted equations to the full action. The independent verifier reconstructs unrestricted element Jacobians, re-integrates the path and RK4 controls, verifies degree-five quadrature, and rejects omitted-dressing and underintegration controls. Stored classical field readouts include electric cochains, matter fields and local charge.",
         "analytic full variational symmetry lift and numerical charged trajectory with independent replay; no observed postdiction",
         "Supplied dimensionless cone, scalar action, e=1/4, m^2=1/2, g=1/4, temporal gauge, initial data and action parameter on [0,2]. Exact-real symmetry and polynomial quadrature arguments are separate from floating-point samples. Sample residuals, independent integration and fixed-step comparisons do not enclose trajectory error rigorously. The trajectory occupies the fixed five-coordinate sector, with zero magnetic field but nonzero electric interaction and matter motion. No authenticated observer execution, actual quantum state, spatial continuum trajectory limit, physical clock or measured postdiction is supplied."),
    )
    rows = []
    for name, fragment, label, statement, match, boundary in specs:
        paper = "paper/tex_fragments/" + fragment
        path = REPO / paper
        if not path.is_file() or "\\label{" + label + "}" not in path.read_text(encoding="utf-8"):
            raise SystemExit("coupled Whitney analytic theorem missing: " + paper + "#" + label)
        row = {
            "id": name, "statement": statement,
            "observed_counterpart": "Mathematical classical/quantum field structures; no observed data",
            "match": match, "physical_comparison_status": "NOT_EVALUABLE",
            "observed_postdiction": False, "continuum_convergence_established": False,
            "spatial_action_consistency": name == "whitney_spatial_consistency",
            "analytic_paper_theorem": {"source": paper, "label": label},
            "lean_declarations": {}, "lean_receipts": [],
            "artifact_refs": ["paper/observers_are_all_you_need.tex", paper],
            "hypothesis_boundary": boundary,
            "paper_ref": "observers synthesis, " + name.replace("_", " "),
        }
        if name == "whitney_spatial_consistency":
            declarations = ("weighted_antisymmetric_cancellation", "weighted_path_phase_zero",
                            "polarized_cancellation", "centered_dressing_error", "centered_potential_variation")
            row.update({
                "lean_declarations": {"WhitneySpatialConsistency": list(declarations)},
                "lean_receipts": _lean_receipt("WhitneySpatialConsistency", declarations={"WhitneySpatialConsistency": declarations}),
                "lean_scope": "finite cancellation and error-decomposition algebra only; no formal spatial convergence theorem",
                "certificate_scope": spatial_receipt["scope"],
                "independent_verifier_result": {key: spatial_verified[key] for key in spatial_keys},
            })
            artifacts = ("runtime/whitney_spatial_consistency_receipt.json", "whitney_spatial_consistency.py",
                         "verify_whitney_spatial_consistency.py", "test_whitney_spatial_consistency.py")
        elif name == "whitney_interacting_quantum":
            supporting_labels = (
                "eq:whitney-interacting-global-metric-bound",
                "prop:whitney-interacting-gaussian-state",
                "eq:whitney-interacting-gaussian-state",
                "eq:whitney-interacting-gaussian-matter-moments",
                "eq:whitney-interacting-gaussian-magnetic-moment",
            )
            paper_text = path.read_text(encoding="utf-8")
            for supporting_label in supporting_labels:
                if "\\label{" + supporting_label + "}" not in paper_text:
                    raise SystemExit("coupled Whitney analytic theorem missing: " + paper + "#" + supporting_label)
            observable_source = "paper/tex_fragments/WHITNEY_COMMON_OBSERVABLES.tex"
            observable_labels = (
                "thm:whitney-common-observables",
                "eq:whitney-common-observables",
                "eq:whitney-common-observable-pvm",
                "eq:whitney-common-observable-law",
                "cor:whitney-continuum-detector-readouts",
            )
            observable_path = REPO / observable_source
            observable_text = observable_path.read_text(encoding="utf-8") if observable_path.is_file() else ""
            for observable_label in observable_labels:
                if "\\label{" + observable_label + "}" not in observable_text:
                    raise SystemExit("coupled Whitney analytic theorem missing: " + observable_source + "#" + observable_label)
            row["statement"] += (
                " The same reconstructed gauge-invariant scalar quadratic/quartic and magnetic smearings define self-adjoint multiplication operators on maximal neutral domains, a joint spectral PVM and a pushforward probability law for every normalized neutral state. Bounded Borel detector functions retain the identical classical configuration readout. Under the separately stated smooth real Neumann reference and Ritz initialization, scalar smearings and fixed Lipschitz detector responses converge classically at O(1/n)."
            )
            row["hypothesis_boundary"] += (
                " The shared-observable theorem uses real bounded scalar smearings and square-integrable magnetic smearings. Its continuum detector corollary retains the zero-current real sector, m^2>0, smooth Neumann reference, conforming refinement and Ritz data. These are analytic results, not numeric proofs of multiplication self-adjointness or quantum refinement. No ground-state assumption, quantum-classical state identification or physical detector calibration follows. General bounded Borel detectors have no asserted classical continuum convergence rate."
            )
            row["artifact_refs"].extend([observable_source, "paper/tex_fragments/WHITNEY_REAL_CONTINUUM.tex"])
            row.update({"hilbert_space_constructed": True, "computed_quantum_state_history": False,
                        "metric_completeness_established": True,
                        "essential_self_adjointness_established": True,
                        "unique_extension_given_ordering": True,
                        "initial_state_constructed": True,
                        "initial_observables_provided": True,
                        "shared_reconstructed_observable_algebra": {
                            "source": observable_source, "labels": list(observable_labels),
                            "self_adjoint_neutral_multipliers": True,
                            "joint_spectral_probability_law": True,
                            "conditional_real_sector_detector_convergence": True,
                            "ground_state_required": False,
                            "quantum_refinement_convergence": False,
                            "physical_detector_calibration": False,
                            "proved_by_numeric_parent_replay": False,
                        },
                        "analytic_supporting_results": {"source": paper, "labels": list(supporting_labels)},
                        "certificate_scope": quantum_receipt["scope"],
                        "independent_verifier_result": {key: quantum_verified[key] for key in quantum_keys},
                        "analytic_proof_formalized_in_lean": False})
            artifacts = ("whitney_interacting_quantum.py", "test_whitney_interacting_quantum.py",
                         "runtime/whitney_quantum_state_receipt.json", "whitney_quantum_state.py",
                         "verify_whitney_quantum_state.py", "test_whitney_quantum_state.py")
        else:
            row.update({"certificate_scope": charged_receipt["scope"],
                        "independent_verifier_result": {key: charged_verified[key] for key in charged_keys},
                        "rigorous_trajectory_error_enclosure": False,
                        "authenticated_observer_history": False, "computed_quantum_state_history": False})
            artifacts = ("runtime/whitney_charged_dynamics_receipt.json", "whitney_charged_dynamics.py",
                         "verify_whitney_charged_dynamics.py", "test_whitney_charged_dynamics.py")
        row["artifact_refs"].extend("code/electromagnetism/" + artifact for artifact in artifacts)
        rows.append(row)
    return rows


def _verify_whitney_completion_parent(stem: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Fresh replay of one new continuum, instrument, clock or trial parent."""
    if stem not in {"real_continuum", "charged_instrument", "ephemeris_clock", "quantum_history"}:
        raise SystemExit("unknown Whitney completion evidence parent")
    path = CODE / "electromagnetism" / f"verify_whitney_{stem}.py"
    spec = importlib.util.spec_from_file_location(f"whitney_{stem}_structural_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"missing independent Whitney {stem} verifier")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    receipt = verifier.load()
    return receipt, verifier.verify(receipt)


def _verify_whitney_checkpoint_parent(receipt_path: Path | None = None) -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    """Replay the separate checkpoint certificate from fresh verifier bytes."""
    directory = CODE / "electromagnetism"
    path = directory / "verify_whitney_charged_checkpoint.py"
    spec = importlib.util.spec_from_file_location("whitney_checkpoint_ledger_verifier", path)
    if spec is None or spec.loader is None:
        raise SystemExit("missing independent Whitney checkpoint verifier")
    verifier = importlib.util.module_from_spec(spec)
    exec(compile(path.read_bytes(), str(path), "exec"), verifier.__dict__)
    receipt_path = receipt_path or directory / "runtime/whitney_charged_checkpoint_receipt.json"
    raw = receipt_path.read_bytes()
    try:
        receipt = verifier.load(receipt_path)
        summary = verifier.verify(receipt)
    except ValueError as exc:
        raise SystemExit(f"Whitney checkpoint independent replay failed: {exc}") from exc
    if receipt_path.read_bytes() != raw:
        raise SystemExit("Whitney checkpoint receipt changed during replay")
    return raw, receipt, summary


def _whitney_checkpoint_control(receipt_path: Path | None = None) -> dict[str, Any]:
    """Keep exact checkpoint errors separate from the original instrument scope."""
    raw, receipt, verified = _verify_whitney_checkpoint_parent(receipt_path)
    required = {
        "accepted": True,
        "scope": "CERTIFIED_DECODED_CHARGED_CHECKPOINTS__SUPPLIED_ACTION_TIME",
        "events": 1782, "decoded_checkpoints": 81, "original_qv_dimension": 10,
        "exact_model_step": "1/40", "model_time_horizon": "2",
        "checkpoint_qv_error_certified": True, "parent_enclosure_freshly_verified": True,
        "observer_events_exactly_replayed": True,
        "continuous_observer_error_certified": False, "nonlinear_field_error_certified": False,
        "configuration_clock_error_certified": False, "physical_clock_calibrated": False,
        "external_signature_attestation": False, "quantum_history": False,
    }
    projection = {}
    for key, expected in required.items():
        actual = verified.get(key)
        if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
            raise SystemExit("Whitney checkpoint category/count mismatch: " + key)
        projection[key] = actual
    for key in ("historical_qv_error_upper", "maximum_decoded_reference_difference",
                "decoded_checkpoint_error_upper", "simple_checkpoint_error_upper"):
        value = verified.get(key)
        try:
            valid = (type(value) is str and str(Fraction(value)) == value
                     and Fraction(value) >= 0 and receipt.get("bounds", {}).get(key) == value)
        except (ValueError, ZeroDivisionError):
            valid = False
        if not valid:
            raise SystemExit("Whitney checkpoint exact bound attachment mismatch: " + key)
        projection[key] = value
    historical = Fraction(projection["historical_qv_error_upper"])
    difference = Fraction(projection["maximum_decoded_reference_difference"])
    transferred = Fraction(projection["decoded_checkpoint_error_upper"])
    simple = Fraction(projection["simple_checkpoint_error_upper"])
    if historical != Fraction(1, 10**10) or transferred != historical + difference or not (
            transferred <= simple == Fraction(10001, 10**14)):
        raise SystemExit("Whitney checkpoint triangle bound mismatch")
    return {
        "receipt": "code/electromagnetism/runtime/whitney_charged_checkpoint_receipt.json",
        "receipt_pin": {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()},
        "independent_verifier_result": projection,
        "observed_postdiction": False,
        "event_replay_scope": (
            "Exact writer/value, probe, retained response, feedback and decode replay. "
            "The frozen numerical-advance outputs are compared directly with certified "
            "historical samples; their floating solver is not recomputed by this consumer."
        ),
        "scope_boundary": (
            "The exact rational triangle inequality bounds ten original reduced q/v "
            "coordinates at the 81 decoded j/40 checkpoints only. Intermediate probe "
            "registers, continuous observer evolution, nonlinear fields, configuration-clock "
            "quadrature and physical calibration are not enclosed."
        ),
    }


def _whitney_completion_rows() -> list[dict[str, Any]]:
    """Project exact scopes and counts; never project recomputed float errors."""
    scopes = {
        "real_continuum": "UNIFORM_CONE_REFINEMENT__CONDITIONAL_REAL_QUARTIC_WAVE_CONTINUUM__FINITE_NUMERIC_CHECKS",
        "charged_instrument": "SELF_READING_COUPLED_CLASSICAL_EXECUTION__COMPUTATIONAL_PATCHES__MODEL_TIME_ONLY",
        "ephemeris_clock": "ACTION_DERIVED_MODEL_CLOCK__NUMERICAL_POLYLINE_READOUT",
        "quantum_history": "FULL_56D_NEUTRAL_POTENTIAL_PHASE_TRIAL__GLOBAL_NORM_BOUND__DECLARED_INPUTS",
    }
    required = {
        "real_continuum": {"accepted": True, "refinement_parameters": [1, 2, 4, 8],
            "tetrahedra": [20, 160, 1280, 10240], "uniform_shape_bound": "72",
            "conditional_real_sector_trajectory_bound": True, "full_charged_complex_trajectory_bound": False,
            "numerical_trajectory_error_certified": False, "physical_source_or_clock_selected": False,
            "formalized_in_lean": False},
        "charged_instrument": {"accepted": True, "patches": 5, "writable_real_registers": 10,
            "events": 1782, "decoded_samples": 81, "completed_repair_cycles": 405,
            "solver_advances": 80, "full_equations_per_sample": 68, "gauss_equations_per_sample": 13,
            "exact_record_restoration": True, "observer_software_history": True,
            "physical_clock_calibrated": False, "physical_observer_placement": False,
            "quantum_state_history": False, "spatial_trajectory_convergence": False,
            "rigorous_trajectory_enclosure": False, "empirical_prediction": False},
        "ephemeris_clock": {"accepted": True, "source_configurations": 81,
            "calibrated_physical_clock": False, "rigorous_numerical_enclosure": False,
            "quantum_clock_operator": False},
        "quantum_history": {"accepted": True, "real_configuration_dimension": 56,
            "state_samples": 5, "phase_configurations": 4, "global_time_coverage": True,
            "trial_history_computed": True, "exact_Hamiltonian_history_computed": False,
            "configuration_density_moves": False, "observer_history": False,
            "physical_state_preparation": False, "empirical_comparison": False,
            "ordinary_physics_benchmark": False, "analytic_proof_formalized_in_Lean": False,
            "numeric_quadrature_used_for_bound": False},
    }
    specs = (
        ("real_continuum", "WHITNEY_REAL_CONTINUUM.tex", "thm:whitney-real-continuum-trajectory",
         ["prop:whitney-uniform-cone-refinement", "lem:whitney-real-sector-invariant", "eq:whitney-real-trajectory-error"],
         "A conforming uniformly shape-regular refinement of the same solid supports the full charged action's invariant real scalar sector with A=phi=0. Given a C^2([0,T];H^2) real Neumann solution and mass-shifted Ritz initial position and velocity, the semidiscrete nonlinear wave trajectories converge at O(1/n) in H^1 position plus L^2 velocity, uniformly on the fixed time interval. The proof uses the actual invariant sector and retains its full-action equations and constraints.",
         "analytic conditional real-sector continuum trajectory bound; independent finite geometry/Ritz and wave checks; no observed postdiction",
         "Supplied cone and continuum time, scalar species/action, m^2>0, g>=0, real neutral sector, natural Neumann boundary conditions, and an existing C^2_t H^2_x continuum reference. The theorem requires Ritz-compatible initialization. The stored pulse trajectories instead use nodal initialization and verify numerical implementation, not the continuum error bound. No full charged-complex trajectory convergence, numerical interval enclosure, source-selected geometry/matter/clock or empirical comparison is established."),
        ("charged_instrument", "WHITNEY_CHARGED_INSTRUMENT.tex", "prop:whitney-charged-record-restoration", [],
         "Five computational observer-like patches execute 1782 events with local coordinate/velocity states, ring ports, destructive averaging probes, retained records and feedback. All 405 probe cycles restore their rational registers exactly. Eighty numerical action advances consume decoded states; 81 decoded frames reconstruct the same charged action's full fields, with all 68 configuration equations and 13 Gauss equations independently checked. A separate exact consumer freshly verifies the charged IVP enclosure and record replay: all ten decoded q/v coordinates at nominal j/40 checkpoints lie within 10001/10^14 of the exact trajectory. Historical samples remain comparison-only inputs, not instrument evolution inputs.",
         "exact software readback/restoration, replayed numerical fields and separately certified decoded checkpoint errors; no observed postdiction",
         "Supplied symmetry-coordinate patch placement, classical writable registers and records, numerical solver, cone, scalar action, initial data and action step. Hash-pinned replay proves internal software provenance without external attestation. The five patches are computational coordinates, not physical observer locations. Repair counts are operational events; their assignment to model time is supplied. The original numerical residuals are not rigorous trajectory enclosures. The separate checkpoint certificate transfers the historical sample bound by an exact rational triangle inequality; it does not enclose intermediate probe registers, continuous observer evolution, nonlinear field or configuration-clock readouts. No quantum history, continuum trajectory limit or laboratory clock calibration is attached."),
        ("ephemeris_clock", "WHITNEY_EPHEMERIS_CLOCK.tex", "prop:whitney-ephemeris-clock", ["eq:whitney-ephemeris-clock"],
         "The complete coupled kinetic metric, potential and supplied energy define the Jacobi-Maupertuis duration d_tau=sqrt(G[dq,dq]/(2(E-V))). On a regular nonturning stationary path of the fixed-energy Jacobi action, this timing recovers the natural action evolution and is invariant under positive reparameterization. A consumer integrates polygonal paths through 81 independently decoded configurations without using recorded velocities, timestamps or repair counts in its clock integral. A fixed smooth nonturning curve has an analytic O(delta^2) polygon-duration estimate.",
         "standard Jacobi timing attached to the same action and software records; independently replayed numerical readout; no observed postdiction",
         "Supplied full action, energy, temporal gauge, configuration-coordinate convention and ordered classical path. The analytic statement requires positive kinetic energy and E-V>0 and excludes turning points. The numerical polyline and quadrature comparisons do not certify segment-wide admissibility or integration/trajectory errors. Source authentication may inspect record metadata; the duration calculation consumes configurations only. This model-internal duration has no selected physical units, laboratory calibration, Lorentzian spacetime identification or quantum-clock operator."),
        ("quantum_history", "WHITNEY_QUANTUM_HISTORY.tex", "thm:whitney-quantum-trial-history", [],
         "The full 56-real-coordinate interacting Hilbert space admits the explicitly time-dependent neutral trial v(t)=exp(-itV/hbar)f_sigma. Exact global Gaussian moments and analytic coefficient bounds certify its Hilbert-norm distance from the exact Hamiltonian evolution for every time in a declared short interval. The packet supplies five times and four exact phase configurations while keeping the configuration probability density fixed; it does not restrict the quantum state to the classical five-coordinate trajectory.",
         "analytic trial-evolution error theorem and independently replayed exact global bound; no observed postdiction",
         "Supplied interacting action, Hilbert measure, operator ordering, hbar, Gaussian width and Hamiltonian time. This is a potential-phase trial with a conservative dimensionless horizon 2^-55 and norm-error target 1/10, not a computed exact Hamiltonian history or a practical ordinary-physics simulation. Global analytic envelopes and exact Gaussian moments avoid numerical quadrature or omitted Gaussian tails in the bound. No evolving configuration density, observer preparation, physical clock, continuum QFT or empirical comparison is established."),
    )
    rows = []
    for stem, fragment, label, supporting_labels, statement, match, boundary in specs:
        receipt, verified = _verify_whitney_completion_parent(stem)
        if (receipt.get("scope") != scopes[stem]
                or verified.get("scope", scopes[stem] if stem == "ephemeris_clock" else None) != scopes[stem]):
            raise SystemExit("Whitney completion certificate scope mismatch: "+stem)
        projection = {}
        for key, expected in required[stem].items():
            actual = verified.get(key)
            if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
                raise SystemExit("Whitney completion category/count mismatch: "+stem+"/"+key)
            projection[key] = actual
        projection["scope"] = scopes[stem]
        if stem == "quantum_history":
            for key in ("horizon", "target_norm_error", "squared_norm_error_upper"):
                value = verified.get(key)
                try:
                    valid = type(value) is str and str(Fraction(value)) == value and Fraction(value) > 0
                except (ValueError, ZeroDivisionError):
                    valid = False
                if not valid:
                    raise SystemExit("Whitney completion quantum bound needs a positive canonical rational: "+key)
                projection[key] = value
            if (Fraction(projection["horizon"]) != Fraction(1, 2**55)
                    or projection["target_norm_error"] != "1/10"
                    or Fraction(projection["squared_norm_error_upper"]) > Fraction(1, 100)):
                raise SystemExit("Whitney completion quantum certified interval/target mismatch")
            for key in ("horizon", "target_norm_error"):
                if receipt.get("error_certificate", {}).get(key) != projection[key]:
                    raise SystemExit("Whitney completion quantum certificate attachment mismatch")
            if receipt.get("error_certificate", {}).get("horizon_squared_error_upper") != projection["squared_norm_error_upper"]:
                raise SystemExit("Whitney completion quantum exact bound attachment mismatch")
        paper = "paper/tex_fragments/"+fragment
        text = (REPO/paper).read_text(encoding="utf-8")
        for named in [label, *supporting_labels]:
            if "\\label{"+named+"}" not in text:
                raise SystemExit("Whitney completion analytic theorem missing: "+paper+"#"+named)
        row = {"id": "whitney_"+stem, "statement": statement, "match": match,
            "observed_counterpart": "Mathematical field, software-instrument and quantum structures; no observed data",
            "observed_postdiction": False, "physical_comparison_status": "NOT_EVALUABLE",
            "source_selected_physical_continuum": False,
            "continuum_convergence_established": stem == "real_continuum",
            "continuum_convergence_scope": "conditional analytic invariant real-sector trajectory bound" if stem == "real_continuum" else "not established by this row",
            "certificate_scope": scopes[stem], "independent_verifier_result": projection,
            "analytic_paper_theorem": {"source": paper, "label": label},
            "analytic_supporting_results": {"source": paper, "labels": supporting_labels},
            "analytic_proof_formalized_in_lean": False, "lean_declarations": {}, "lean_receipts": [],
            "artifact_refs": ["paper/observers_are_all_you_need.tex", paper]+[
                "code/electromagnetism/"+name for name in (
                    f"runtime/whitney_{stem}_receipt.json", f"whitney_{stem}.py",
                    f"verify_whitney_{stem}.py", f"test_whitney_{stem}.py")],
            "hypothesis_boundary": boundary, "paper_ref": "observers synthesis, whitney "+stem.replace("_", " ")}
        if stem == "charged_instrument":
            row["decoded_checkpoint_control"] = _whitney_checkpoint_control()
            row["artifact_refs"].extend("code/electromagnetism/" + name for name in (
                "runtime/whitney_charged_checkpoint_receipt.json", "whitney_charged_checkpoint.py",
                "verify_whitney_charged_checkpoint.py", "test_whitney_charged_checkpoint.py"))
        rows.append(row)
    return rows


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
    contraction = _load("p_interval_contraction")
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
    lambda_scale = _load("lambda_transmutation")
    nucleon_ratio = _load("nucleon_external_ratio")

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
        "alpha": _alpha_rows(endpoint, bridge, alpha_hvp_verdict)
        + _alpha_closure_rows(endpoint, contraction),
        "charged_leptons": _lepton_rows(surface, rectangle, coherent, koide),
        "electroweak": _ew_rows(conditional),
        "quarks": _quark_rows(obstruction, clebsch, selection),
        "hadrons": _hadron_rows(
            payload, standby, lambda_scale, nucleon_ratio
        ),
        "neutrinos": [
            {
                "id": "neutrino_dimensionless_pointer",
                "statement": (
                    "dimensionless PMNS and mass-splitting-ratio "
                    "comparisons are recorded on the results surface; the "
                    "absolute attachment stays compare-only"
                ),
                "artifact_ref": "code/particles/RESULTS_STATUS.md",
            }
        ],
    }
    sections["forced_structure"].append(_seam_maxwell_continuum_row())
    sections["forced_structure"].append(_serial_maxwell_readout_row())
    sections["forced_structure"].append(_cone_whitney_bridge_row())
    sections["forced_structure"].extend(_whitney_dynamics_rows())
    sections["forced_structure"].extend(_whitney_coupled_rows())
    sections["forced_structure"].extend(_whitney_completion_rows())
    sections["forced_structure"].extend(_pauli_structure_rows())
    result = {
        "artifact": "oph_postdiction_ledger",
        "generator": "code/particles/scripts/build_postdiction_ledger.py",
        "schema_version": 3,
        "row_class": "compare_only_postdiction_ledger",
        "guards": {
            "compare_only": True,
            "public_promotion_allowed": False,
            "changes_any_solve_path": False,
            "new_axiom_introduced": False,
            "hand_typed_measured_values": False,
        },
        "aggregation_policy": (
            "numeric values and measured references are read mechanically from cited "
            "parents; structural rows distinguish analytic paper proofs, validated Lean "
            "declarations and structured parents, and identify direct "
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
    add("Numeric values and measured references on this page are read mechanically from "
        "the cited parent artifacts. Structural rows are derived from validated "
        "Lean declarations, analytic paper proofs, structured parents, or their combination, and direct algebraic "
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
    add("These structural results precede or constrain numeric lanes. "
        "They include the icosahedral gauge packet and generic observer-law "
        "boundaries. Each row distinguishes analytic paper proofs, finite Lean results "
        "and structured executable checks, and records its own classical inputs and missing "
        "physical attachments.")
    add("")
    add("| Result | Observed counterpart | Match | Receipts |")
    add("| --- | --- | --- | --- |")
    for row in s["forced_structure"]:
        receipts = list(row.get("lean_receipts", []))
        if row.get("artifact_ref"):
            receipts.append(row["artifact_ref"])
        receipts.extend(row.get("artifact_refs", []))
        receipt_txt = ", ".join(f"`{r}`" for r in receipts)
        add(f"| {_cell(row['statement'])} | {_cell(row['observed_counterpart'])} | "
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
        if row["id"] != "alpha_inv_thomson_endpoint":
            continue
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
            "- Scientific owner: "
            + ", ".join(f"#{i}" for i in row["scientific_owner_issues"])
        )
    closure_rows = [
        row for row in s["alpha"] if row["id"] != "alpha_inv_thomson_endpoint"
    ]
    if closure_rows:
        add("")
        add(
            "The closure map of the pixel lane has a certified fixed point in "
            "each declared mode, and the lane reads as a chain: the fixed point "
            "of the root map, the fixed point of the same map with the unified "
            "gauge width, and the term between that second fixed point and the "
            "reference value. Each row is compare-only: the CODATA reference "
            "sits outside every solve path, and the certificate permits no "
            "promotion."
        )
        add("")
        add(
            "| Closure map | Fixed point `alpha_em^-1` | Enclosure width | "
            "`P` | Distance to CODATA | Relative |"
        )
        add("| --- | ---: | ---: | ---: | ---: | ---: |")
        for row in closure_rows:
            width = float(row["enclosure_width_decimal"])
            add(
                f"| {row['closure_map_description']} | "
                f"`{_fmt(row['value_central'], 12)}` | "
                f"`{_fmt(width, 2)}` | "
                f"`{_fmt(row['p_central'], 12)}` | "
                f"`{row['deviation_inv_alpha']:+.6f}` | "
                f"`{row['relative_deviation']:+.2e}` |"
            )
        add("")
        for row in closure_rows:
            definition = row["map_definition"].rstrip(". ")
            add(
                f"- `{row['id']}`: {definition}. Existence and "
                "uniqueness of the fixed point are certified by interval "
                f"arithmetic with Lipschitz bound "
                f"`{_fmt(row['lipschitz_bound'], 4)}` at cutoffs "
                f"`{row['su2_cutoff']}` and `{row['su3_cutoff']}`, with the "
                "tails bounded."
            )
        gauge = closure_rows[-1]
        add(
            "- The distance from the gauge-width fixed point to the CODATA "
            f"reference, `{abs(gauge['deviation_inv_alpha']):.6f}` inverse-alpha "
            "units, is the open term of this lane. It carries the hadronic "
            "content that the Thomson-endpoint row above accounts for "
            "retrospectively, and closing it from the source side is work in "
            "progress under the scientific owner "
            + ", ".join(f"#{i}" for i in gauge["scientific_owner_issues"])
            + "."
        )
        add(
            "- Neither closure row is a frozen prediction, and neither is "
            "eligible as a blind prediction: both are retrospective "
            "comparisons of a certified fixed point against a reference value."
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
                "is the unfixed scheme term of the bridge. The certified width floor "
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
                f"(scientific owner {', '.join(f'#{i}' for i in row['scientific_owner_issues'])}).")
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
        elif row["id"] == "qcd_solver_on_standby":
            add(f"- QCD solver: `{row['status']}`; invocation is gated on the "
                "source-side parameter emissions recorded in the standby receipt.")
        elif row["id"] == "lambda_qcd_transmutation_scale":
            lo, hi = row["value_interval_gev"]
            add(f"- Transmutation scale ({row['tier']}): `Lambda_QCD^(3) = "
                f"{_fmt(row['value_central_gev'], 6)}` GeV in "
                f"`[{_fmt(lo, 6)}, {_fmt(hi, 6)}]` against the published "
                f"central `{_fmt(row['published_central_gev'], 4)}` GeV "
                f"(compare-only), `{row['relative_deviation']:+.2%}` relative. "
                f"{row['declared_external_inputs'].capitalize()}. "
                f"{row['claim_boundary']}")
        elif row["id"] == "nucleon_mass_external_ratio":
            lo, hi = row["value_interval_gev"]
            inside = "inside" if row["interval_contains_measured"] else "outside"
            add(f"- Nucleon mass ({row['tier']}): `{_fmt(row['value_central_gev'], 6)}` "
                f"GeV in `[{_fmt(lo, 6)}, {_fmt(hi, 6)}]` against the measured "
                f"proton mass `{_fmt(row['measured_gev'], 6)}` GeV "
                f"(compare-only), `{row['relative_deviation']:+.2%}` relative, "
                f"with the measured value {inside} the interval. The declared "
                f"external lattice-theory ratio is "
                f"`{row['external_theory_factor']}` with uncertainty "
                f"`{row['external_theory_uncertainty']}`. {row['claim_boundary']}")
    add("")
    add("## Neutrinos")
    add("")
    for row in s["neutrinos"]:
        add(f"- {row['statement']} (`{row['artifact_ref']}`).")
    add("")
    return "\n".join(lines)



def _cell(text: object) -> str:
    """Escape pipes so free text cannot break a Markdown table row."""
    return str(text).replace("|", "\\|")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed JSON or Markdown differs from a deterministic rebuild",
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
