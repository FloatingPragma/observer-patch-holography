"""Cross-surface scope gates for the source-history finite GNS attachment."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM_ID = "OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT"
LEAN_PATH = "Lean/QFT/SourceHistoryGNSDynamics.lean"
HISTORY_CLAIM_ID = "OPH-FINITE-HISTORY-VARIATIONAL-HELPERS"
GNS_CLAIM_ID = "OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT"


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _row(rows: list[dict], key: str, value: str) -> dict:
    matches = [row for row in rows if row[key] == value]
    assert len(matches) == 1, (key, value, len(matches))
    return matches[0]


def test_claim_records_the_source_counted_attachment_and_open_exits() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["assumptions"] == ["source_counted_history_quantum_extension"]
    assert claim["evidence"] == [
        LEAN_PATH,
        "Lean/InformationProjection/SourceHistoryPacket.lean",
        "Lean/QFT/ColimitSelectedStateGNS.lean",
        "paper/reality_as_consensus_protocol.tex",
        "paper/observers_are_all_you_need.tex",
        "extra/compact_proof_of_oph.tex",
    ]
    assert claim["gates"] == [730, 739, 743]
    assert claim["premise_dependencies"] == {
        "classification": "explicit_edges",
        "consumed": [],
        "open": ["PR-15", "PR-52", "PR-54", "PR-58"],
        "boundary": [],
    }
    statement = claim["statement"]
    for text in [
        "window count divided by 1754",
        "exact energy expectation is 197/1754",
        "does not annihilate the cyclic unit",
        "off-diagonal matrix unit E_01",
        "full matrix algebra",
        "No regional net",
        "issues 730, 739, and 743 are open",
    ]:
        assert text in statement


def test_dependency_graph_has_only_the_two_scientific_parents() -> None:
    graph = _json("claims/dependency_graph.json")
    incoming = {
        edge["from"] for edge in graph["edges"] if edge["to"] == CLAIM_ID
    }
    assert incoming == {HISTORY_CLAIM_ID, GNS_CLAIM_ID}


def test_lean_packet_keeps_source_and_quantum_layers_distinct() -> None:
    source = (ROOT / LEAN_PATH).read_text(encoding="utf-8")
    for declaration in [
        "theorem sourceHistoryDensity_isState",
        "theorem sourceHistoryDensity_eq_windowCount",
        "theorem sourceHistoryHamiltonian_eq_changes",
        "theorem sourceHistoryDensity_meanEnergy",
        "theorem sourceHistoryProbeGenerator_ne_zero",
        "theorem stageRepresentation_injective",
        "theorem cyclicUnit_hamiltonian_matrixCoefficient",
        "theorem representedHamiltonian_cyclicUnit_ne_zero",
        "theorem representedHeisenberg_intertwines",
        "theorem sourceHistoryFiniteGNSAttachment",
    ]:
        assert declaration in source
    assert "the source data only fixes its diagonal" in source
    assert "has no physical clock calibration" in source
    umbrella = (ROOT / "Lean/QFT.lean").read_text(encoding="utf-8")
    assert umbrella.splitlines().count("import QFT.SourceHistoryGNSDynamics") == 1
    assert "exact mean energy `197 / 1754`" in umbrella


def test_papers_state_the_result_and_the_non_vacuum_boundary() -> None:
    consensus = (ROOT / "paper/reality_as_consensus_protocol.tex").read_text(
        encoding="utf-8"
    )
    assert "A source-counted finite history Hamiltonian on the selected GNS" in consensus
    assert r"\operatorname{Tr}(\rho_{\rm hist}H_{\rm hist})" in consensus
    assert r"\pi_{\rm hist}(H_{\rm hist})\Omega_{\rm hist}\ne0" in consensus
    assert "The retained run was hash-pinned" in consensus
    assert "jointly preregistered" in consensus
    observers = (ROOT / "paper/observers_are_all_you_need.tex").read_text(
        encoding="utf-8"
    )
    assert "The retained length-three history packet also has an exact finite operator" in observers
    assert "The cyclic unit has nonzero" in observers
    compact = (ROOT / "extra/compact_proof_of_oph.tex").read_text(
        encoding="utf-8"
    )
    assert "source-counted operator attachment is exact" in compact
    assert "unit has nonzero energy and is not a vacuum" in compact


def test_observation_and_premise_rows_stay_open() -> None:
    observation = _row(
        _json("tracking/observation_ledger.json")["rows"], "id", "OL-C6"
    )
    assert observation["status"] == "partial"
    assert observation["open_premises"] == ["PR-15", "PR-52", "PR-58"]
    assert LEAN_PATH in observation["evidence"]
    assert "this carrier has no regional net" in observation["notes"]

    premises = _json("tracking/premise_register.json")["rows"]
    for premise_id in ["PR-15", "PR-52", "PR-54", "PR-58"]:
        row = _row(premises, "id", premise_id)
        assert LEAN_PATH in row["evidence"]
        assert row["evidence_roles"][LEAN_PATH] == "statement"
    assert _row(premises, "id", "PR-15")["disposition"] == "import"
    for premise_id in ["PR-52", "PR-54", "PR-58"]:
        assert _row(premises, "id", premise_id)["disposition"] == "remove"


def test_no_observation_or_prediction_is_promoted() -> None:
    frozen = _json("claims/frozen_prediction_register.json")
    for row in frozen["rows"]:
        encoded = json.dumps(row, sort_keys=True)
        assert CLAIM_ID not in encoded
        assert LEAN_PATH not in encoded
