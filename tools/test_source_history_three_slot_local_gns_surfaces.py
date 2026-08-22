"""Cross-surface gates for the source-history three-slot local GNS rung."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM_ID = "OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS"
LEAN_PATH = "Lean/QFT/SourceHistoryThreeSlotLocalGNS.lean"
ASSUMPTION = "tensor_local_history_coordinate_reading_adapter"
SOURCE_GNS_PARENT = "OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT"
ISING_PARENT = "OPH-QFT-FINITE-TWO-SITE-ISING-HILBERT-DYNAMICS"
CONSENSUS_PAPER = "paper/reality_as_consensus_protocol.tex"
OBSERVERS_PAPER = "paper/observers_are_all_you_need.tex"


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _row(rows: list[dict], key: str, value: str) -> dict:
    matches = [row for row in rows if row[key] == value]
    assert len(matches) == 1, (key, value, len(matches))
    return matches[0]


def _csv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _collapsed(relative_path: str) -> str:
    return " ".join((ROOT / relative_path).read_text(encoding="utf-8").split())


def test_claim_records_exact_scope_and_open_physical_exits() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    assert claim["owner_paper"] == CONSENSUS_PAPER
    assert claim["assumptions"] == [ASSUMPTION]
    assert claim["evidence"] == [
        LEAN_PATH,
        "Lean/QFT/SourceHistoryGNSDynamics.lean",
        "Lean/QFT/FiniteTwoSiteIsingField.lean",
        CONSENSUS_PAPER,
        OBSERVERS_PAPER,
    ]
    assert claim["gates"] == [730, 743]
    assert claim["premise_dependencies"] == {
        "classification": "explicit_edges",
        "consumed": [],
        "open": ["PR-15", "PR-52", "PR-54", "PR-58"],
        "boundary": [],
    }
    statement = claim["statement"]
    for text in [
        "A0, A1, and A2",
        "H = H01 + H12",
        "47/877",
        "103/1754",
        "normalized positive functional",
        "compatible restriction of one global state",
        "not a construction of conditional expectations",
        "issues 730 and 743 remain open",
    ]:
        assert text in statement


def test_claim_has_exactly_the_two_scientific_parents() -> None:
    graph = _json("claims/dependency_graph.json")
    incoming = {
        edge["from"] for edge in graph["edges"] if edge["to"] == CLAIM_ID
    }
    assert incoming == {SOURCE_GNS_PARENT, ISING_PARENT}
    assert graph["nodes"].count(CLAIM_ID) == 1


def test_assumption_marks_tensor_local_reading_as_adapter() -> None:
    dictionary = (ROOT / "claims/assumption_dictionary.md").read_text(
        encoding="utf-8"
    )
    rows = [line for line in dictionary.splitlines() if f"`{ASSUMPTION}`" in line]
    assert len(rows) == 1
    row = rows[0]
    assert "mathematical adapter" in row
    assert "does not source-select the full matrix-factor interpretation" in row
    assert "state-preserving regional conditional expectations" in row


def test_novelty_and_falsification_rows_are_unique_and_bounded() -> None:
    novelty = _row(_csv("claims/novelty_matrix.csv"), "claim_id", CLAIM_ID)
    assert novelty["novelty_type"] == (
        "source_history_three_slot_tensor_local_bond_decomposition_and_"
        "selected_GNS_transport"
    )
    assert "normalized positive compatible state restrictions" in novelty[
        "oph_specific_delta"
    ]

    falsification = _row(
        _csv("claims/falsification_matrix.csv"), "claim_id", CLAIM_ID
    )
    physical = falsification["physical_identification_falsifier"]
    assert "compatible restrictions of one global state" in physical
    assert "PR-15, PR-52, PR-54, or PR-58" in physical
    assert "no phenomenological prediction" in falsification[
        "phenomenological_falsifier"
    ]


def test_lean_packet_exposes_locality_propagation_and_gns_transport() -> None:
    source = (ROOT / LEAN_PATH).read_text(encoding="utf-8")
    assert "native_decide" not in source
    for declaration in [
        "def historyEquiv",
        "def algebra0",
        "def algebra1",
        "def algebra2",
        "def algebra01",
        "def algebra12",
        "def algebra012",
        "theorem algebra0_algebra2_locality",
        "theorem sourceHistoryHamiltonian_eq_bond01_add_bond12",
        "theorem sourceHistoryDensity_bond01_expectation",
        "theorem sourceHistoryDensity_bond12_expectation",
        "theorem restrictedSourceState_nonneg",
        "theorem restrictedSourceState_compatible",
        "theorem restrictedSourceState_receipt",
        "theorem sourceHistoryGenerator_maps_algebra0_to_algebra01",
        "theorem sourceHistoryGenerator_maps_algebra2_to_algebra12",
        "theorem represented_endpoint_locality",
        "theorem sourceHistoryThreeSlotLocalGNSAttachment",
    ]:
        assert declaration in source


def test_claim_does_not_promote_compatible_restriction_to_expectation() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    combined = " ".join(
        [claim["statement"], claim["oph_specific_delta"], claim["status"]]
    )
    assert "state-preserving regional expectation" in combined
    assert "unconstructed" in combined or "not a construction" in combined
    assert claim["claim_class"] == "conditional_implication"


def test_qft_umbrella_imports_packet_once_and_states_its_boundary() -> None:
    umbrella = (ROOT / "Lean/QFT.lean").read_text(encoding="utf-8")
    assert umbrella.splitlines().count(
        "import QFT.SourceHistoryThreeSlotLocalGNS"
    ) == 1
    collapsed_umbrella = " ".join(umbrella.split())
    for text in [
        "three binary tensor slots",
        "two positive adjacent domain-wall bonds",
        "source functional restricts compatibly",
        "not Lorentzian regions or calibrated times",
    ]:
        assert text in collapsed_umbrella


def test_all_three_papers_state_the_local_rung_and_physical_boundary() -> None:
    consensus = _collapsed(CONSENSUS_PAPER)
    for text in [
        "The same carrier has an exact tensor-local history decomposition.",
        r"H_{\rm hist}=H_{01}+H_{12}",
        "first-commutator finite-support statements",
        "does not assert a source-state-preserving conditional expectation",
        "not a physical QFT",
    ]:
        assert text in consensus

    observers = _collapsed(OBSERVERS_PAPER)
    for text in [
        "The eight histories also reindex exactly as three binary tensor slots.",
        r"H_{\rm hist}=H_{01}+H_{12}",
        "source functional restricts compatibly along all these algebra inclusions",
        "not calibrated physical times or Lorentzian regions",
    ]:
        assert text in observers


def test_observation_and_premise_surfaces_keep_the_physical_exits_open() -> None:
    observation = _row(
        _json("tracking/observation_ledger.json")["rows"], "id", "OL-C6"
    )
    assert observation["status"] == "partial"
    assert observation["open_premises"] == ["PR-15", "PR-52", "PR-58"]
    assert LEAN_PATH in observation["evidence"]
    assert "SourceHistoryThreeSlotLocalGNS identifies" in observation["notes"]
    assert "normalized compatible functionals" in observation["notes"]

    premises = _json("tracking/premise_register.json")["rows"]
    for premise_id in ["PR-15", "PR-52", "PR-54", "PR-58"]:
        premise = _row(premises, "id", premise_id)
        assert LEAN_PATH in premise["evidence"]
        assert premise["evidence_roles"][LEAN_PATH] == "statement"
    assert _row(premises, "id", "PR-15")["disposition"] == "import"
    for premise_id in ["PR-52", "PR-54", "PR-58"]:
        assert _row(premises, "id", premise_id)["disposition"] == "remove"
