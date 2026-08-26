"""Cross-surface gates for the V3.34 fourth wave: the proper-length-clocked
chain, the golden-sector irreducibility, the timelike-class force law, and
the carrier evolution flow."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

CHAIN_CLAIM = "OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN"
IRR_CLAIM = "OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY"
BLOCK_CLAIM = "OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW"
FLOW_CLAIM = "OPH-EM-CARRIER-EVOLUTION-FLOW"

CHAIN_LEAN = "Lean/Geometry/ProperLengthClockedChain.lean"
IRR_LEAN = "Lean/Screen/GoldenSectorIrreducibility.lean"
BLOCK_LEAN = "Lean/Geometry/TimelikeClassForceLaw.lean"
FLOW_LEAN = "Lean/Screen/CarrierEvolutionFlow.lean"


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
    gates = {CHAIN_CLAIM: [736, 739], IRR_CLAIM: [733, 728], BLOCK_CLAIM: [740],
             FLOW_CLAIM: [730, 733]}
    for claim_id, expected in gates.items():
        claim = _claim(claim_id)
        assert claim["gates"] == expected, claim_id
        assert claim["premise_dependencies"]["consumed"] == [], claim_id


def test_claims_keep_declared_scope() -> None:
    chain = _claim(CHAIN_CLAIM)["statement"]
    for token in ("increments in {0, 1}", "exactly d and below one",
                  "declared enrichment of the join", "nothing selects u from the source",
                  "discharges none"):
        assert token in chain, token
    irr = _claim(IRR_CLAIM)["statement"]
    for token in ("dimension nine", "irreducible real representations",
                  "This module proves real irreducibility only",
                  "The later dedicated GoldenSectorComplexIrreducibility module transports the span certificates",
                  "proves complex irreducibility, inequivalence, and real-type commutants",
                  "external inference", "discharges neither"):
        assert token in irr, token
    block = _claim(BLOCK_CLAIM)["statement"]
    for token in ("shifted later ports", "zero clock difference",
                  "Along the respective scaled Ampere evolutions",
                  "h times the field-energy transfer", "the unit drops out",
                  "discharges none"):
        assert token in block, token
    flow = _claim(FLOW_CLAIM)["statement"]
    for token in ("determinant one and trace 2 - h^2 lam", "Phi_h = T",
                  "exactly the gradient space",
                  "This module does not itself package the energy as a positive Hilbert form",
                  "the later CurlSectorEigenbasis module supplies the nineteen nonzero curl modes",
                  "Reading t as physical time is declared",
                  "no identification with that flow is claimed", "discharges none"):
        assert token in flow, token


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        CHAIN_LEAN: ("clockedIndex_rest", "clockedIndex_unit_increments",
                     "clockedIndex_uniform_rate", "clockedReturns_eq_returnCount",
                     "clockedIndexPath_refine", "no_stepDuration_matches_uniform",
                     "unit_not_selected", "clocked_vs_index_scope"),
        IRR_LEAN: ("goldenPlus_irreducible", "goldenMinus_irreducible",
                   "finrank_WPlus", "finrank_SPlus", "Phi_S_surjective",
                   "W_invariant", "conj_certificates"),
        BLOCK_LEAN: ("transportedAction_shifted_variation", "shifted_variation_rejoin",
                     "delay_action_difference", "position_stationary_iff",
                     "static_field_drift", "blockMomentum_delay",
                     "exchange_action_eq_h_mul_transfer"),
        FLOW_LEAN: ("ampere_iff_stepMatrix", "stepMatrix_det", "stepMatrix_eq_conj",
                    "rotFlow_add", "rotFlow_nat_step", "modeForm_rotFlow",
                    "shearFlow_step", "assembled_flow", "carrier_flow_five_modes",
                    "carrier_flow_stone_reading"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path


def test_owner_paper_carries_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("ProperLengthClockedChain", "GoldenSectorIrreducibility",
                  "TimelikeClassForceLaw", "CarrierEvolutionFlow"):
        assert token in observers, token


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert CHAIN_LEAN in rows["OL-H8"]["evidence"]
    assert IRR_LEAN in rows["OL-F2"]["evidence"]
    assert BLOCK_LEAN in rows["OL-N1"]["evidence"]
    assert FLOW_LEAN in rows["OL-C2"]["evidence"]
    assert rows["OL-C2"]["status"] == "partial"
    assert rows["OL-N1"]["status"] == "owed"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert CHAIN_LEAN in prows["PR-15"]["evidence"]
    assert FLOW_LEAN in prows["PR-15"]["evidence"]
    assert IRR_LEAN in prows["PR-53"]["evidence"]
    assert BLOCK_LEAN in prows["PR-54"]["evidence"]
    assert prows["PR-15"]["disposition"] == "import"
