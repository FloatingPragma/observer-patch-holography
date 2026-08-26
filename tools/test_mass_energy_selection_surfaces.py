"""Cross-surface gates for the V3.31 core-physics round: hopping-charge
minimal coupling, covariance and proper-time slope selection, the Gibbs
reference identification, modular-energy additivity, and carrier mode
oscillators."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

MONO_CLAIM = "OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING"
COV_CLAIM = "OPH-GEOMETRY-COMPOSITE-MOMENTUM-COVARIANCE"
SLOPE_CLAIM = "OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION"
OSC_CLAIM = "OPH-EM-CARRIER-MODE-OSCILLATORS"
GIBBS_CLAIM = "OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION"
BIND_CLAIM = "OPH-THERMO-MODULAR-ENERGY-ADDITIVITY"

MONO_LEAN = "Lean/Geometry/PortChargeMinimalCoupling.lean"
COV_LEAN = "Lean/Geometry/CompositeMomentumCovariance.lean"
SLOPE_LEAN = "Lean/Geometry/ProperTimeInternalAction.lean"
OSC_LEAN = "Lean/Screen/CarrierModeOscillators.lean"
GIBBS_LEAN = "Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean"
BIND_LEAN = "Lean/Thermodynamics/ModularEnergyAdditivity.lean"


def _registry() -> list[dict]:
    return yaml.safe_load((ROOT / "claims/claim_registry.yaml").read_text(
        encoding="utf-8"))["claims"]


def _claim(claim_id: str) -> dict:
    matches = [row for row in _registry() if row["claim_id"] == claim_id]
    assert len(matches) == 1, (claim_id, len(matches))
    return matches[0]


def _collapsed(relative_path: str) -> str:
    return " ".join((ROOT / relative_path).read_text(encoding="utf-8").split())


def test_claims_registered_with_gates() -> None:
    gates = {MONO_CLAIM: [740], COV_CLAIM: [736, 739], SLOPE_CLAIM: [736, 739],
             OSC_CLAIM: [736, 733], GIBBS_CLAIM: [729, 736],
             BIND_CLAIM: [736, 739]}
    for claim_id, expected in gates.items():
        assert _claim(claim_id)["gates"] == expected, claim_id


def test_every_claim_names_rows_and_discharges_none() -> None:
    for claim_id in (MONO_CLAIM, COV_CLAIM, SLOPE_CLAIM):
        claim = _claim(claim_id)
        assert claim["premise_dependencies"]["consumed"] == [], claim_id
        assert claim["premise_dependencies"]["open"] == ["PR-15", "PR-52", "PR-54"], claim_id
        assert "discharges none" in claim["statement"], claim_id


def test_monopole_claim_keeps_declared_scope() -> None:
    statement = _claim(MONO_CLAIM)["statement"]
    for token in ("total q at every step", "-(q/h) on the hopped seam",
                  "equivalent to the continuity equation on the interior steps",
                  "no energy of the charge and no equation of motion",
                  "kappa = -(q h)/(12 - 4 phi)",
                  "at the endpoints only",
                  "not derived from the Lorentz-module worldline"):
        assert token in statement, token


def test_covariance_claim_is_kinematic() -> None:
    statement = _claim(COV_CLAIM)["statement"]
    for token in ("exactly when lambda = 1 or E = 0", "5/4 and 3/4",
                  "(1 - lambda) E (m + lambda E) = 0",
                  "not about dynamics",
                  "does not derive the mass-energy identity"):
        assert token in statement, token


def test_slope_claim_names_both_declared_principles() -> None:
    statement = _claim(SLOPE_CLAIM)["statement"]
    for token in ("exactly when b = 0", "proper-time principle, declared",
                  "uniqueness of the declared member",
                  "inertial coefficient of the refinement-invariant form is exactly m + E",
                  "declared enrichment"):
        assert token in statement, token


def test_oscillator_claim_keeps_inferences_outside_theorems() -> None:
    statement = _claim(OSC_CLAIM)["statement"]
    for token in ("returns after six steps", "minimality of these periods is not stated",
                  "an inference outside the theorems", "h^2 < 2 / phi^2",
                  "7155 p / 61511", "disjoint scope", "discharges none"):
        assert token in statement, token
    assert _claim(OSC_CLAIM)["premise_dependencies"]["open"] == ["PR-15", "PR-52", "PR-53"]


def test_thermodynamic_claims_keep_declared_clauses() -> None:
    gibbs = _claim(GIBBS_CLAIM)["statement"]
    for token in ("Nothing derives that the record reference is Gibbs",
                  "nothing selects beta", "free exactly up to an additive constant",
                  "discharges neither"):
        assert token in gibbs, token
    assert _claim(GIBBS_CLAIM)["premise_dependencies"]["open"] == ["PR-15", "PR-52"]
    bind = _claim(BIND_CLAIM)["statement"]
    for token in ("correlated or not", "its sign is not fixed",
                  "-1/4 + log((3 + e)/4)", "a definitional composition",
                  "discharges none"):
        assert token in bind, token
    assert _claim(BIND_CLAIM)["premise_dependencies"]["open"] == ["PR-15", "PR-54"]


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        MONO_LEAN: ("hopping_continuity", "hoppingLoad_total",
                    "monopoleCoupledAction_eq_augmented",
                    "sourcePairing_gauge_invariant_iff",
                    "hopping_work_energy_forward", "bridge_iff_charge_fixed",
                    "bridge_off_endpoint_exhibit",
                    "charge_unconstrained_by_field_sector"),
        COV_LEAN: ("frameCovariant_slopeMomentum_iff", "slopeMomentum_mapFrame",
                   "lorentzQ_slopeMomentum", "slopeMomentum_shell_all_iff",
                   "standardBoost_moves_restVector",
                   "compositeFourMomentum_eq_slopeMomentum_one",
                   "covariance_selects_slope_one"),
        SLOPE_LEAN: ("properLength_refine", "clockAction_refine",
                     "refinementInvariant_iff_b_zero",
                     "restPhase_window_properLength",
                     "properTimeAgreement_iff_slope_one",
                     "lengthAction_hasDerivAt", "lengthStationary_iff_eom",
                     "length_inertial_coefficient", "slope_selection"),
        OSC_LEAN: ("cosHistory_ampere", "twoMode_eigen", "threeMode_eigen",
                   "threeMode_period_four", "window_iff_goldenRatio",
                   "projector_traces", "modeAngle_ratio_bounds",
                   "modeAngle_div_tendsto", "period_return_ratio_tick_free"),
        GIBBS_LEAN: ("modularEnergy_gibbs", "cap_firstLaw_gibbs", "slope_unique",
                     "gibbs_energy_unique_up_to_const", "screen_beta_unique",
                     "boundary_not_affine", "labInverseTemperature_not_forced"),
        BIND_LEAN: ("modularEnergy_productRef", "defect_gibbs",
                    "bindingConstant_nonneg", "defect_gibbs_le_const",
                    "two_point_defect_pos", "composite_inertial_ledger",
                    "faithful_is_gibbs"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path
        if relative_path in (COV_LEAN, SLOPE_LEAN):
            assert "dischargedRows_empty" in text, relative_path


def test_owner_paper_carries_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("PortChargeMinimalCoupling", "CompositeMomentumCovariance",
                  "ProperTimeInternalAction", "two declared principles",
                  "CarrierModeOscillators", "disjoint scope",
                  "ModularEnergyAdditivity", "sign of the full defect is not fixed"):
        assert token in observers, token
    gravity = _collapsed(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex")
    for token in ("GibbsReferenceEnergyIdentification", "Nothing derives that the record reference is Gibbs"):
        assert token in gravity, token


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert MONO_LEAN in rows["OL-N1"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"
    assert SLOPE_LEAN in rows["OL-H8"]["evidence"]
    assert OSC_LEAN in rows["OL-H8"]["evidence"]
    assert rows["OL-H8"]["status"] == "partial"
    assert GIBBS_LEAN in rows["OL-B1"]["evidence"]
    assert rows["OL-B1"]["status"] == "partial"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert MONO_LEAN in prows["PR-54"]["evidence"]
    assert COV_LEAN in prows["PR-54"]["evidence"]
    assert SLOPE_LEAN in prows["PR-15"]["evidence"]
    assert OSC_LEAN in prows["PR-15"]["evidence"]
    assert OSC_LEAN in prows["PR-53"]["evidence"]
    assert GIBBS_LEAN in prows["PR-15"]["evidence"]
    assert BIND_LEAN in prows["PR-54"]["evidence"]
    for pid in ("PR-15", "PR-52", "PR-54"):
        assert prows[pid]["disposition"] in ("remove", "axiomatize", "import"), pid
