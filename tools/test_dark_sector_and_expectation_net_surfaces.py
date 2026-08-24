"""Cross-surface gates for the collar-saturation, deep-profile, per-cut,
expectation-net, and certified-step rungs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

COLLAR_CLAIM = "OPH-DM-COLLAR-SATURATION"
PROFILE_CLAIM = "OPH-DM-DEEP-PROFILE-CHARACTERIZATION"
PERCUT_CLAIM = "OPH-DM-PER-CUT-COMPOSITION"
NET_CLAIM = "OPH-QFT-SOURCE-HISTORY-STATE-PRESERVING-RECORD-NET"
STEP_CLAIM = "OPH-EM-CERTIFIED-SCALED-STEP-INSTRUMENT"
ALL_CLAIMS = (COLLAR_CLAIM, PROFILE_CLAIM, PERCUT_CLAIM, NET_CLAIM, STEP_CLAIM)

COLLAR_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/CollarScaleSaturation.lean"
PROFILE_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean"
PERCUT_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean"
NET_LEAN = "Lean/QFT/SourceHistoryExpectationNet.lean"
STEP_LEAN = "Lean/Screen/CertifiedScaledStepInstrument.lean"

DARK_PAPER = "cosmology/oph_dark_matter_paper.tex"
CONSENSUS_PAPER = "paper/reality_as_consensus_protocol.tex"
OBSERVERS_PAPER = "paper/observers_are_all_you_need.tex"
SCREEN_PAPER = "paper/screen_microphysics_and_observer_synchronization.tex"
FLAGSHIP_PAPER = "flagship/from_observer_consensus_to_standard_physics.tex"


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


def test_all_five_claims_are_registered_with_their_gates() -> None:
    gates = {
        COLLAR_CLAIM: [742, 751],
        PROFILE_CLAIM: [742, 751],
        PERCUT_CLAIM: [742, 751],
        NET_CLAIM: [730, 743],
        STEP_CLAIM: [737, 739],
    }
    for claim_id, expected in gates.items():
        claim = _claim(claim_id)
        assert claim["gates"] == expected, (claim_id, claim["gates"])
        assert claim["claim_class"] == "conditional_implication", claim_id


def test_collar_saturation_claim_carries_rate_and_delimitation() -> None:
    statement = _claim(COLLAR_CLAIM)["statement"]
    for token in ("kappa * B * exp", "rho^(1/3)", "15/(4 pi^2)",
                  "identically zero defect", "stays a premise"):
        assert token in statement, token
    assert COLLAR_LEAN in _claim(COLLAR_CLAIM)["evidence"]


def test_profile_claim_carries_uniqueness_and_load_bearing_receipts() -> None:
    statement = _claim(PROFILE_CLAIM)["statement"]
    for token in ("unique constant a0 > 0", "r * sqrt(Mb * a0 / G)",
                  "Cauchy lemma", "matching radius",
                  "additive law", "r^2 law",
                  "value of a0 is not determined"):
        assert token in statement, token
    assert PROFILE_LEAN in _claim(PROFILE_CLAIM)["evidence"]


def test_percut_claim_carries_dictionary_and_non_identifiability() -> None:
    statement = _claim(PERCUT_CLAIM)["statement"]
    for token in ("a0 = G n^2 c", "quadrature closure",
                  "theorem rather than a premise", "n^2 c",
                  "no numeric value of n, c, or a0"):
        assert token in statement, token
    assert PERCUT_LEAN in _claim(PERCUT_CLAIM)["evidence"]


def test_net_claim_carries_exact_witnesses_and_obstruction() -> None:
    statement = _claim(NET_CLAIM)["statement"]
    for token in ("383/415", "2/89", "94/1754", "103/1754", "197/1754",
                  "no complex-linear map", "diagonal embedding"):
        assert token in statement, token
    deps = _claim(NET_CLAIM)["premise_dependencies"]
    assert deps["classification"] == "explicit_edges"
    assert deps["open"] == ["PR-15", "PR-52", "PR-54", "PR-58"]


def test_step_claim_carries_certificate_and_no_run_boundary() -> None:
    statement = _claim(STEP_CLAIM)["statement"]
    for token in ("h = 4/5", "84/25 < 4", "sqrt 5 < 9/4", "25/2",
                  "declared selection", "No run is authorized"):
        assert token in statement, token
    deps = _claim(STEP_CLAIM)["premise_dependencies"]
    assert deps["consumed"] == ["PR-66"]
    assert deps["open"] == ["PR-15", "PR-53", "PR-54"]


def test_dependency_graph_carries_the_new_nodes_and_edges() -> None:
    graph = json.loads((ROOT / "claims/dependency_graph.json").read_text(
        encoding="utf-8"))
    for claim_id in ALL_CLAIMS:
        assert claim_id in graph["nodes"], claim_id
    edges = {(edge["from"], edge["to"]) for edge in graph["edges"]}
    expected = {
        ("OPH-GR-D2-COLLAR-RECOVERY", COLLAR_CLAIM),
        ("OPH-DM-CONT", COLLAR_CLAIM),
        ("OPH-DM-CONT", PROFILE_CLAIM),
        (PROFILE_CLAIM, PERCUT_CLAIM),
        ("OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS", NET_CLAIM),
        ("OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT", NET_CLAIM),
        ("OPH-EM-SCALED-MAXWELL-STABILITY-ACTION", STEP_CLAIM),
    }
    missing = expected - edges
    assert not missing, missing


def test_matrix_rows_exist_for_every_new_claim() -> None:
    fals = {row["claim_id"] for row in _csv("claims/falsification_matrix.csv")}
    nov = {row["claim_id"] for row in _csv("claims/novelty_matrix.csv")}
    for claim_id in ALL_CLAIMS:
        assert claim_id in fals, claim_id
        assert claim_id in nov, claim_id


def test_lean_modules_exist_and_carry_their_headline_declarations() -> None:
    expectations = {
        COLLAR_LEAN: ("eta_and_stress_tendsto_zero_of_density",
                      "envelopeBound_lt_iff",
                      "recoveryEnvelope_admits_zero_defect",
                      "sparseCutRemainder_stress_pos"),
        PROFILE_LEAN: ("deepProfile_characterization",
                       "linear_on_pos_of_additive_monotone",
                       "anomalous_dominates_iff",
                       "additiveLaw_not_quadratureClosed",
                       "sqrtLaw_not_scaleCovariant"),
        PERCUT_LEAN: ("varianceSlope", "perCut_a0_dictionary",
                      "perCut_nonidentifiability"),
        NET_LEAN: ("regionalExpectation_tower",
                   "recordMean_regionalExpectation",
                   "no_statePreserving_bimodule_projection_onto_interval01",
                   "383 / 415"),
        STEP_LEAN: ("stepCertificate_strict", "fourFifthsInstrument",
                    "step_selection_not_forced", "InstrumentRunHandoff"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "sorry" not in text.replace("assert_no_sorry", ""), relative_path


def test_owner_papers_carry_the_new_results() -> None:
    dark = _collapsed(DARK_PAPER)
    for token in ("Compact-source saturation", "Profile characterization",
                  "Quadrature composition", "Deep-regime scale covariance",
                  "a_0=G\\,n^2c", "milgrom2009", "verlinde2016"):
        assert token in dark, token
    consensus = _collapsed(CONSENSUS_PAPER)
    for token in ("SourceHistoryExpectationNet", "383/415", "2/89"):
        assert token in consensus, token
    observers = _collapsed(OBSERVERS_PAPER)
    for token in ("383/415", "2/89"):
        assert token in observers, token
    screen = _collapsed(SCREEN_PAPER)
    for token in ("CertifiedScaledStepInstrument", "84/25"):
        assert token in screen, token
    flagship = _collapsed(FLAGSHIP_PAPER)
    for token in ("The dark sector and the deep galaxy law", "a_0=Gn^2c"):
        assert token in flagship, token


def test_ledger_rows_cite_the_new_modules_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert NET_LEAN in rows["OL-C6"]["evidence"]
    assert rows["OL-C6"]["status"] == "partial"
    for path in (COLLAR_LEAN, PROFILE_LEAN, PERCUT_LEAN):
        assert path in rows["OL-I3"]["evidence"], path
    assert rows["OL-I3"]["status"] == "owed"
    assert STEP_LEAN in rows["OL-F2"]["evidence"]
    assert rows["OL-F2"]["status"] == "partial"


def test_premise_register_rows_cite_the_new_modules_as_open() -> None:
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in register["rows"]}
    for premise_id in ("PR-15", "PR-52", "PR-54", "PR-58"):
        assert NET_LEAN in rows[premise_id]["evidence"], premise_id
    for premise_id in ("PR-15", "PR-53", "PR-54", "PR-66"):
        assert STEP_LEAN in rows[premise_id]["evidence"], premise_id
