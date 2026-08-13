"""Cross-surface scope gates for the fixed-cutoff unitary no-go (#743)."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM_ID = "OPH-FINITE-UNITARY-SCATTERING-LIMIT-NO-GO"
POSTDICTION_ID = "finite_unitary_scattering_limit_no_go"
LEAN_PATH = "Lean/QFT/FiniteUnitaryScatteringNoGo.lean"

CLAIM_EVIDENCE = [
    "paper/reality_as_consensus_protocol.tex",
    "paper/tex_fragments/QFT_STRUCTURAL_INHERITANCE_STATUS.tex",
    LEAN_PATH,
    "code/qft/finite_unitary_scattering_no_go.py",
    "code/qft/test_finite_unitary_scattering_no_go.py",
]

OBSERVATION_EVIDENCE = [
    LEAN_PATH,
    "code/qft/finite_unitary_scattering_no_go.py",
    "code/qft/test_finite_unitary_scattering_no_go.py",
    "paper/tex_fragments/QFT_STRUCTURAL_INHERITANCE_STATUS.tex",
]

POSTDICTION_ARTIFACTS = [
    "code/qft/finite_unitary_scattering_no_go.py",
    "code/qft/test_finite_unitary_scattering_no_go.py",
    "paper/tex_fragments/QFT_STRUCTURAL_INHERITANCE_STATUS.tex",
]

LEAN_DECLARATIONS = [
    "tendsto_powers_forces_identity",
    "nontrivial_powers_have_no_limit",
    "finite_unitary_powers_have_no_limit",
    "finite_unitary_ambient_powers_have_no_limit",
    "identical_relative_evolution_is_constant",
    "identical_relative_evolution_tendsto",
]


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _row(rows: list[dict], key: str, value: str) -> dict:
    matches = [row for row in rows if row[key] == value]
    assert len(matches) == 1, (key, value, len(matches))
    return matches[0]


def test_claim_is_exact_but_records_missing_scattering_construction() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["evidence"] == CLAIM_EVIDENCE
    assert claim["gates"] == [743]
    assert "equivalence_of_standard_operator_topologies_in_finite_dimension" in (
        claim["imported_results"]
    )
    assert claim["status"] == (
        "exact_direct_power_limit_obstruction_attained__physical_scattering_open"
    )
    statement = claim["statement"]
    assert "full weak-operator convergence" in statement
    assert "finite-dimensional operator topologies coincide" in statement
    assert "selected projected scalar or observable limits" in statement
    assert "infinite-dimensional weak limits" in statement
    assert "No wave operator, S-matrix" in statement


def test_observation_row_remains_owed_with_exact_evidence_and_premise() -> None:
    row = _row(_json("tracking/observation_ledger.json")["rows"], "id", "OL-K2")
    assert row["status"] == "owed"
    assert row["lane_issue"] == 743
    assert row["premises"] == []
    assert row["open_premises"] == ["PR-54"]
    assert row["evidence"] == OBSERVATION_EVIDENCE
    notes = row["notes"]
    assert "full weak-operator convergence" in notes
    assert "selected projected scalar or observable limit" in notes
    assert "infinite-dimensional weak limit" in notes
    assert "a weaker or projected limit" not in notes


def test_detailed_paper_states_fixed_cutoff_scope_without_issue_status() -> None:
    fragment = (ROOT / "paper/tex_fragments/QFT_STRUCTURAL_INHERITANCE_STATUS.tex").read_text(
        encoding="utf-8"
    )
    prose = " ".join(fragment.split())
    assert "Any fixed-cutoff direct-power construction" in prose
    assert "Any asymptotic construction must" not in prose
    assert "full weak-operator convergence does not evade" in prose
    assert "selected projected scalar or observable readouts" in prose
    assert "infinite-dimensional weak limit" in prose
    assert "interacting-QFT/RG/scattering construction remains unresolved" in prose
    assert r"it does not close issue~\#743" not in prose


def test_no_frozen_prediction_row_is_emitted() -> None:
    frozen = _json("claims/frozen_prediction_register.json")
    for row in frozen["rows"]:
        encoded = json.dumps(row, sort_keys=True)
        assert CLAIM_ID not in encoded
        assert POSTDICTION_ID not in encoded
        assert LEAN_PATH not in encoded


def test_postdiction_row_matches_theorem_and_boundary() -> None:
    rows = _json("code/particles/runs/status/postdiction_ledger.json")["sections"][
        "forced_structure"
    ]
    row = _row(rows, "id", POSTDICTION_ID)
    assert row["observed_counterpart"] == (
        "fixed-cutoff scattering-construction boundary; not a "
        "physical observation or prediction"
    )
    assert row["match"] == (
        "exact direct-power obstruction; physical scattering and "
        "asymptotic comparison construction not constructed"
    )
    assert row["lean_declarations"]["FiniteUnitaryScatteringNoGo"] == (
        LEAN_DECLARATIONS
    )
    assert row["artifact_refs"] == POSTDICTION_ARTIFACTS
    assert set(POSTDICTION_ARTIFACTS).issubset(CLAIM_EVIDENCE)
    boundary = row["hypothesis_boundary"]
    assert "Full weak-operator convergence" in boundary
    assert "finite-dimensional operator topologies coincide" in boundary
    assert "selected projected scalar or observable limits" in boundary
    assert "infinite-dimensional weak limits" in boundary
    assert "Scientific owner #743 records" in boundary
    assert "no frozen prediction" in boundary


def test_qft_umbrella_imports_the_kernel_checked_module() -> None:
    umbrella = (ROOT / "Lean/QFT.lean").read_text(encoding="utf-8")
    import_line = "import QFT.FiniteUnitaryScatteringNoGo"
    assert umbrella.splitlines().count(import_line) == 1
