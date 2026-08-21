"""Cross-surface scope gates for the finite two-site Ising benchmark (#743)."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM_ID = "OPH-QFT-FINITE-TWO-SITE-ISING-HILBERT-DYNAMICS"
LEAN_PATH = "Lean/QFT/FiniteTwoSiteIsingField.lean"


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _row(rows: list[dict], key: str, value: str) -> dict:
    matches = [row for row in rows if row[key] == value]
    assert len(matches) == 1, (key, value, len(matches))
    return matches[0]


def test_claim_is_exact_finite_and_keeps_physical_exits_open() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["evidence"] == [
        LEAN_PATH,
        "paper/reality_as_consensus_protocol.tex",
    ]
    assert claim["gates"] == [743]
    assert claim["premise_dependencies"] == {
        "classification": "explicit_edges",
        "consumed": [],
        "open": ["PR-15", "PR-52", "PR-54", "PR-58"],
        "boundary": [],
    }
    statement = claim["statement"]
    assert "not any sum of independent one-site Hamiltonians" in statement
    assert "H_del = (1 / 2) I" in statement
    assert "not a Hamiltonian attachment" in statement
    assert "ground space is degenerate" in statement
    assert "issue 743 remains open" in statement


def test_observation_row_stays_partial_and_records_the_non_edge() -> None:
    row = _row(_json("tracking/observation_ledger.json")["rows"], "id", "OL-C6")
    assert row["status"] == "partial"
    assert row["lane_issue"] == 730
    assert row["open_premises"] == ["PR-15", "PR-52", "PR-58"]
    assert LEAN_PATH in row["evidence"]
    notes = row["notes"]
    assert "non-reducibility to any sum of one-site Hamiltonians" in notes
    assert "not composed with the selected-state tower GNS representation" in notes
    assert "PR-54" in notes


def test_paper_states_exact_model_and_continuum_boundary() -> None:
    paper = (ROOT / "paper/reality_as_consensus_protocol.tex").read_text(
        encoding="utf-8"
    )
    assert r"\paragraph{An exact finite two-site Hamiltonian benchmark.}" in paper
    assert r"H_{\mathrm I}=\frac12\bigl(1-Z\otimes Z\bigr)" in paper
    assert r"H_L\otimes1+1\otimes H_R" in paper
    assert r"H_{\mathrm{del}}=\tfrac12 1" in paper
    assert "not an interacting\nquantum field theory" in paper
    assert "not yet a Hamiltonian attachment" in paper
    assert "lane~\\#743 remains open" in paper


def test_lean_surface_and_umbrella_keep_the_controls_visible() -> None:
    source = (ROOT / LEAN_PATH).read_text(encoding="utf-8")
    for declaration in [
        "theorem isingHamiltonian_not_sum_one_site",
        "theorem isingGroundSpace_degenerate",
        "theorem isingGroundDensity_stationary",
        "theorem left_right_equal_time_locality",
        "theorem uncoupledHeisenberg_eq",
        "theorem isingLeftPauliX_dynamics_nontrivial",
        "theorem isingInteraction_load_bearing",
    ]:
        assert declaration in source
    assert "The ground space is\ndegenerate" in source
    umbrella = (ROOT / "Lean/QFT.lean").read_text(encoding="utf-8")
    assert umbrella.splitlines().count("import QFT.FiniteTwoSiteIsingField") == 1
    assert "not attached to the tower\nGNS representation" in umbrella


def test_dependency_graph_has_exact_inputs_and_no_false_gns_edge() -> None:
    edges = _json("claims/dependency_graph.json")["edges"]
    pairs = {(edge["from"], edge["to"]) for edge in edges}
    assert ("OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE", CLAIM_ID) in pairs
    assert ("OPH-EVENTALGEBRA-PRODUCT-SPLIT-SLOT-LOCAL-CHSH", CLAIM_ID) in pairs
    assert ("OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT", CLAIM_ID) not in pairs


def test_no_frozen_prediction_is_emitted() -> None:
    frozen = _json("claims/frozen_prediction_register.json")
    for row in frozen["rows"]:
        encoded = json.dumps(row, sort_keys=True)
        assert CLAIM_ID not in encoded
        assert LEAN_PATH not in encoded
