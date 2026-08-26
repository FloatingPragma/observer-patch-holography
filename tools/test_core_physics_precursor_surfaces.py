"""Cross-surface gates for the recurrence-clock, charge-fixed interaction,
internal-energy inertia, and local energy balance rungs."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

CLK_CLAIM = "OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR"
CHG_CLAIM = "OPH-GEOMETRY-CHARGE-FIXED-INTERACTION"
INR_CLAIM = "OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR"
BAL_CLAIM = "OPH-EM-LOCAL-ENERGY-BALANCE"

CLK_LEAN = "Lean/QFT/SourceRecurrenceClock.lean"
CHG_LEAN = "Lean/Geometry/ChargeFixedInteraction.lean"
INR_LEAN = "Lean/Geometry/InternalEnergyInertia.lean"
BAL_LEAN = "Lean/Screen/LocalEnergyBalance.lean"


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
    gates = {CLK_CLAIM: [736], CHG_CLAIM: [740], INR_CLAIM: [736, 739],
             BAL_CLAIM: [729, 733]}
    for claim_id, expected in gates.items():
        assert _claim(claim_id)["gates"] == expected, claim_id


def test_every_claim_names_rows_and_discharges_none() -> None:
    for claim_id in (CLK_CLAIM, CHG_CLAIM, INR_CLAIM, BAL_CLAIM):
        claim = _claim(claim_id)
        assert claim["premise_dependencies"]["consumed"] == [], claim_id
        assert "open" in claim["status"], claim_id


def test_clock_claim_is_a_precursor_with_exact_rationals() -> None:
    statement = _claim(CLK_CLAIM)["statement"]
    for token in ("no cycle", "61511/7155", "61511/54356", "54356/7155",
                  "carries no tick", "standard Markov-chain period is one",
                  "no deterministic return interval exists",
                  "No physical duration is identified",
                  "forbids nothing here"):
        assert token in statement, token
    assert _claim(CLK_CLAIM)["premise_dependencies"]["open"] == ["PR-15"]


def test_charge_claim_states_polarization_and_conditionality() -> None:
    statement = _claim(CHG_CLAIM)["statement"]
    for token in ("kappa = -(q h)/(12 - 4 phi)", "polarization-type",
                  "no monopole identification",
                  "resting worldline induces no load",
                  "nonzero step, identification, and worldline shape are declared",
                  "2m times the forward second difference"):
        assert token in statement, token


def test_inertia_claim_is_conditional_on_a_declared_shape() -> None:
    statement = _claim(INR_CLAIM)["statement"]
    for token in ("nothing here derives the identity",
                  "inertial coefficient is exactly m + E",
                  "slope family m + lambda E",
                  "declared enrichments",
                  "no derived mass"):
        assert token in statement, token


def test_balance_claim_keeps_one_component_scope() -> None:
    statement = _claim(BAL_CLAIM)["statement"]
    for token in ("one energy component", "equal quarters",
                  "recovers the committed global energy balance exactly",
                  "no symmetric tensor, no covariance, no continuum limit",
                  "named open join"):
        assert token in statement, token


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        CLK_LEAN: ("acceptedStep_no_cycle", "periodic_iterate_is_fixed",
                   "kac_identity", "meanReturn_ratio_tick_free",
                   "no_single_period", "returnDuration_not_forced"),
        CHG_LEAN: ("coupledAction_reduction", "identification_iff_kappa",
                   "seam_crossing_charge_fixed", "inducedLoad_rest",
                   "coupled_worldline_equation",
                   "normalization_free_without_identification"),
        INR_LEAN: ("compositeStationary_iff_eom",
                   "composite_inertial_coefficient",
                   "shape_selection_discriminated",
                   "legendre_nonidentifiability_cited",
                   "dischargedRows_empty"),
        BAL_LEAN: ("faceEnergyDensity_sum", "seamFluxRate_conserved",
                   "local_energy_balance", "global_balance_from_local",
                   "fourDivergence_source_free"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path


def test_owner_papers_carry_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("SourceRecurrenceClock", "ChargeFixedInteraction",
                  "InternalEnergyInertia", "polarization"):
        assert token in observers, token
    gravity = _collapsed(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex")
    for token in ("LocalEnergyBalance", "one energy component"):
        assert token in gravity, token


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert CLK_LEAN in rows["OL-H8"]["evidence"]
    assert rows["OL-H8"]["status"] == "partial"
    assert CHG_LEAN in rows["OL-N1"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"
    assert BAL_LEAN in rows["OL-B1"]["evidence"]
    assert rows["OL-B1"]["status"] == "partial"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert CLK_LEAN in prows["PR-15"]["evidence"]
    assert CHG_LEAN in prows["PR-54"]["evidence"]
    assert BAL_LEAN in prows["PR-52"]["evidence"]
    for pid in ("PR-15", "PR-52", "PR-54"):
        assert prows[pid]["disposition"] in ("remove", "axiomatize", "import"), pid
