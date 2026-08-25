"""Cross-surface gates for the defect-mechanism, premise-derivation,
cap-generator-split, source-binding, and Maxwell-clock-join rungs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

MECH_CLAIM = "OPH-DM-SPARSE-DEFECT-MECHANISM"
DERIV_CLAIM = "OPH-DM-COLLAR-PREMISE-DERIVATION"
SPLIT_CLAIM = "OPH-DM-FINITE-CAP-GENERATOR-SPLIT"
BIND_CLAIM = "OPH-QUANTUM-SOURCE-BOUND-INSTRUMENT-INTERFACE"
JOIN_CLAIM = "OPH-GEOMETRY-COMMON-WORLD-MAXWELL-CLOCK-JOIN"
ALL_CLAIMS = (MECH_CLAIM, DERIV_CLAIM, SPLIT_CLAIM, BIND_CLAIM, JOIN_CLAIM)

MECH_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/SparseRecordDefectWitness.lean"
DERIV_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean"
SPLIT_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean"
BIND_LEAN = "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean"
JOIN_LEAN = "Lean/Geometry/CommonWorldMaxwellClockJoin.lean"

DARK_PAPER = "cosmology/oph_dark_matter_paper.tex"
OBSERVERS_PAPER = "paper/observers_are_all_you_need.tex"


def _registry() -> list[dict]:
    return yaml.safe_load((ROOT / "claims/claim_registry.yaml").read_text(
        encoding="utf-8"))["claims"]


def _claim(claim_id: str) -> dict:
    matches = [row for row in _registry() if row["claim_id"] == claim_id]
    assert len(matches) == 1, (claim_id, len(matches))
    return matches[0]


def _csv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _collapsed(relative_path: str) -> str:
    return " ".join((ROOT / relative_path).read_text(encoding="utf-8").split())


def test_all_five_claims_registered_with_gates_and_class() -> None:
    gates = {
        MECH_CLAIM: [742, 751],
        DERIV_CLAIM: [742, 751],
        SPLIT_CLAIM: [742, 751],
        BIND_CLAIM: [730],
        JOIN_CLAIM: [740],
    }
    for claim_id, expected in gates.items():
        claim = _claim(claim_id)
        assert claim["gates"] == expected, (claim_id, claim["gates"])
        assert claim["claim_class"] == "conditional_implication", claim_id


def test_mechanism_claim_carries_exact_defect_law() -> None:
    statement = _claim(MECH_CLAIM)["statement"]
    for token in ("(1 - p)/2", "15/(16 pi^2)", "p > 1 - 2 eps",
                  "modeling premise", "p <= 1"):
        assert token in statement, token
    assert MECH_LEAN in _claim(MECH_CLAIM)["evidence"]


def test_derivation_claim_carries_both_derivations() -> None:
    statement = _claim(DERIV_CLAIM)["statement"]
    for token in ("pairwise independent", "not claimed necessary",
                  "arithmetic progression", "r/s - 1 <= count <= r/s + 1",
                  "at most one per-cut amplitude", "Irregular or ergodic"):
        assert token in statement, token
    assert DERIV_LEAN in _claim(DERIV_CLAIM)["evidence"]


def test_split_claim_carries_uniqueness_and_quotient_invariance() -> None:
    statement = _claim(SPLIT_CLAIM)["statement"]
    for token in ("K = 2 pi B + K_anom + c 1", "bulk-weighted mean",
                  "non-unique", "exp(-t) w", "Gibbs", "carried entirely",
                  "nonzero boost"):
        assert token in statement, token
    assert SPLIT_LEAN in _claim(SPLIT_CLAIM)["evidence"]


def test_binding_claim_carries_split_and_continuum() -> None:
    statement = _claim(BIND_CLAIM)["statement"]
    for token in ("placeholder schema", "expected-frequency table",
                  "-(3/256) t", "propositional extensionality",
                  "0 <= t <= 1", "only at t = 0",
                  "No run, public outcome, authenticated export"):
        assert token in statement, token
    deps = _claim(BIND_CLAIM)["premise_dependencies"]
    assert deps["consumed"] == ["PR-02", "PR-04"]
    assert deps["open"] == ["PR-64", "PR-65"]
    assert deps["boundary"] == ["PR-03"]


def test_join_claim_carries_formal_product_and_open_physical_join() -> None:
    statement = _claim(JOIN_CLAIM)["statement"]
    for token in ("formal same-index product", "inherited componentwise",
                  "No theorem identifies a port", "not a physical Maxwell-clock join",
                  "physical carrier map", "interacting common action"):
        assert token in statement, token
    deps = _claim(JOIN_CLAIM)["premise_dependencies"]
    assert deps["classification"] == "explicit_edges"
    assert "PR-53" in deps["open"] and "PR-54" in deps["open"]


def test_dependency_graph_carries_nodes_and_edges() -> None:
    graph = json.loads((ROOT / "claims/dependency_graph.json").read_text(
        encoding="utf-8"))
    for claim_id in ALL_CLAIMS:
        assert claim_id in graph["nodes"], claim_id
    edges = {(edge["from"], edge["to"]) for edge in graph["edges"]}
    expected = {
        ("OPH-DM-CONT", MECH_CLAIM),
        ("OPH-DM-COLLAR-SATURATION", MECH_CLAIM),
        ("OPH-DM-PER-CUT-COMPOSITION", DERIV_CLAIM),
        ("OPH-DM-DEEP-PROFILE-CHARACTERIZATION", DERIV_CLAIM),
        ("OPH-DM-CONT", SPLIT_CLAIM),
        ("OPH-QUANTUM-OPERATIONAL-ADDITIVITY-BOUNDARY", BIND_CLAIM),
        ("OPH-GEOMETRY-COMMON-WORLD-INSTRUMENT-JOIN", JOIN_CLAIM),
        ("OPH-EM-CERTIFIED-SCALED-STEP-INSTRUMENT", JOIN_CLAIM),
    }
    missing = expected - edges
    assert not missing, missing


def test_matrix_rows_exist_for_every_new_claim() -> None:
    fals = {row["claim_id"] for row in _csv("claims/falsification_matrix.csv")}
    nov = {row["claim_id"] for row in _csv("claims/novelty_matrix.csv")}
    for claim_id in ALL_CLAIMS:
        assert claim_id in fals, claim_id
        assert claim_id in nov, claim_id


def test_lean_modules_carry_their_headline_declarations() -> None:
    expectations = {
        MECH_LEAN: ("recordDefect_isLeast_randomized",
                    "unrecordedCut_error_ge_half", "recordDefect_lt_iff",
                    "sparseRecordRemainder_stress_eq", "RecordedCollar"),
        DERIV_LEAN: ("compoundVariance_disjUnion", "cut_eq",
                     "count_linear_deviation", "composed_deep_law",
                     "correlated_pair_breaks_variance_additivity",
                     "density_necessary"),
        SPLIT_LEAN: ("modularGenerator_split", "split_unique",
                     "split_nonunique_without_normalization",
                     "anomalousPart_rescale",
                     "anomalousPart_supported_on_collar",
                     "anomalousPart_collar_boost_free", "boostedState"),
        BIND_LEAN: ("SourceBoundDeterminedData",
                    "SourceBoundInstrumentBinding",
                    "committed_corpus_does_not_determine_binding",
                    "affineRunValuation_gap",
                    "affineRunValuation_isEffectValuation_iff"),
        JOIN_LEAN: ("MaxwellClockJoinedArchitecture", "join_injective",
                    "joined_conservation", "joinSuppliesNoCarrierMap",
                    "fourFifthsJoined"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)


def test_owner_papers_carry_the_new_results() -> None:
    dark = _collapsed(DARK_PAPER)
    for token in ("erasure-record model", "(1-p)/2",
                  "arithmetic progression with spacing",
                  "machine-checked finite counterpart",
                  "bulk-weighted mean"):
        assert token in dark, token
    observers = _collapsed(OBSERVERS_PAPER)
    for token in ("-(3/256)", "placeholder binding schema",
                  "CommonWorldMaxwellClockJoin", "formal same-index product",
                  "not an exhaustive theorem"):
        assert token in observers, token


def test_ledger_rows_cite_the_new_modules_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    for path in (MECH_LEAN, DERIV_LEAN, SPLIT_LEAN):
        assert path in rows["OL-I3"]["evidence"], path
    assert rows["OL-I3"]["status"] == "owed"
    assert BIND_LEAN in rows["OL-C5"]["evidence"]
    assert rows["OL-C5"]["status"] == "partial"
    assert JOIN_LEAN in rows["OL-N1"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"


def test_premise_register_rows_cite_the_new_modules_as_open() -> None:
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in register["rows"]}
    for premise_id in ("PR-65", "PR-03", "PR-64"):
        assert BIND_LEAN in rows[premise_id]["evidence"], premise_id
    for premise_id in ("PR-15", "PR-53", "PR-54"):
        assert JOIN_LEAN in rows[premise_id]["evidence"], premise_id
