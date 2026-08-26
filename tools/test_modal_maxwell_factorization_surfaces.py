"""Cross-surface scope gates for the modal Maxwell factorization (#733)."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM_ID = "OPH-MODAL-MAXWELL-FACTORIZATION-BOUNDARY"
POSTDICTION_ID = "modal_maxwell_factorization_boundary"
LEAN_PATH = "Lean/Screen/ModalMaxwellFactorizationBoundary.lean"

CLAIM_EVIDENCE = [
    "paper/screen_microphysics_and_observer_synchronization.tex",
    "Lean/Screen/SeamCurrentFreePhotonLift.lean",
    LEAN_PATH,
    "code/electromagnetism/modal_maxwell_factorization.py",
    "code/electromagnetism/test_modal_maxwell_factorization.py",
]

OBSERVATION_EVIDENCE = [
    "Lean/Screen/LightSignalAdequacySurface.lean",
    LEAN_PATH,
    "code/electromagnetism/modal_maxwell_factorization.py",
    "code/electromagnetism/test_modal_maxwell_factorization.py",
    "paper/screen_microphysics_and_observer_synchronization.tex",
    "Lean/Screen/LightSignalMaxwellComposition.lean",
    "Lean/Screen/SeamU1HolonomyClassification.lean",
    "Lean/Screen/PositionSpaceMaxwellAction.lean",
    "Lean/Screen/LocalFaceMaxwellAction.lean",
    "code/electromagnetism/runtime/local_face_maxwell_action_receipt.json",
    "Lean/Screen/TemporalMaxwellEvolution.lean",
    "Lean/Screen/ScaledMaxwellStability.lean",
    "Lean/Screen/SeamChargeContinuity.lean",
]

POSTDICTION_ARTIFACTS = [
    "code/electromagnetism/modal_maxwell_factorization.py",
    "code/electromagnetism/test_modal_maxwell_factorization.py",
    "paper/screen_microphysics_and_observer_synchronization.tex",
]

LEAN_DECLARATIONS = [
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

EXCLUDED_PHYSICAL_CONTENT = (
    "local position-space",
    "real field",
    "U(1)",
    "Maxwell action",
    "Gauss",
    "physical current",
    "Lorentz covariance",
    "continuum",
    "laboratory readout",
)


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _row(rows: list[dict], key: str, value: str) -> dict:
    matches = [row for row in rows if row[key] == value]
    assert len(matches) == 1, (key, value, len(matches))
    return matches[0]


def _csv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_claim_is_exactly_bounded_and_records_missing_attachments() -> None:
    claim = _row(_json("claims/claim_registry.yaml")["claims"], "claim_id", CLAIM_ID)
    assert claim["owner_paper"] == (
        "paper/screen_microphysics_and_observer_synchronization.tex"
    )
    assert claim["evidence"] == CLAIM_EVIDENCE
    assert claim["gates"] == [733]
    assert claim["status"] == (
        "exact_modal_factorization_attained__local_source_produced_physical_Maxwell_open"
    )
    assert claim["claim_class"] == "conditional_implication"
    statement = claim["statement"]
    assert "modal/pseudodifferential factorization" in statement
    assert "removes no registered premise" in statement
    assert "emits no frozen prediction" in statement
    assert "attachments remain unresolved" in statement
    assert "does not close issue" not in statement
    assert "same-sign mutation has the oscillator sign on an amplitude with nonzero spatial action" in claim["falsifier"]
    for excluded in EXCLUDED_PHYSICAL_CONTENT:
        assert excluded in statement


def test_observation_row_advances_only_the_modal_partial_rung() -> None:
    row = _row(_json("tracking/observation_ledger.json")["rows"], "id", "OL-F1")
    assert row["status"] == "partial"
    assert row["lane_issue"] == 733
    assert row["premises"] == ["PR-20", "PR-21", "PR-22", "PR-66"]
    assert row["open_premises"] == ["PR-15", "PR-53", "PR-54"]
    assert row["evidence"] == OBSERVATION_EVIDENCE
    notes = row["notes"]
    assert "modal package supplies an exact pseudodifferential first-order factorization" in notes
    assert "LocalFaceMaxwellAction supplies the local finite operator H=C^T C" in notes
    assert "TemporalMaxwellEvolution supplies the unit-step identities conditional on the declared update" in notes
    assert "ScaledMaxwellStability relocates that declaration" in notes
    assert "Euler-Lagrange equations of one declared discrete action" in notes
    assert "SeamChargeContinuity proves the exact declared-source consistency criterion" in notes
    assert "neutral-partner or neutralising-background requirement" in notes
    assert "declared unit step is unstable" in notes
    assert "PR-15, PR-53, and PR-54 stay open" in notes
    assert "does not close issue" not in notes
    for required_boundary in (
        "nonnegative under h^2 Lambda <= 4",
        "sharp constant Lambda = 3 + sqrt(5)",
        "step index is not physical time",
        "sources are declared",
    ):
        assert required_boundary in notes


def test_detailed_paper_states_exact_result_and_full_boundary() -> None:
    paper = (
        ROOT / "paper/screen_microphysics_and_observer_synchronization.tex"
    ).read_text(encoding="utf-8")
    prose = " ".join(paper.split())
    assert r"\label{prop:modal-maxwell-factorization-boundary}" in paper
    assert r"\mathcal G_{a,k}^{\,2}(E,B)" in paper
    assert "same-sign mutation" in prose
    assert "generally pseudodifferential" in prose
    assert "The result depends on every stated premise" in prose
    assert "implies no empirical prediction" in prose
    assert (
        "no physical local-field, gauge-action, source-current, continuum, or readout "
        "attachment" in prose
    )
    assert r"issue~\#733 remains open" not in prose
    for excluded in (
        "local position-space curl operator",
        "assembled real field",
        r"\(\mathrm U(1)\) potential or gauge quotient",
        "Maxwell action",
        "seam-incidence Gauss receipts",
        "conserved physical current",
        "Lorentz covariance",
        "continuum control",
        "laboratory readout",
    ):
        assert excluded in prose


def test_postdiction_row_matches_theorem_and_boundary() -> None:
    rows = _json("code/particles/runs/status/postdiction_ledger.json")["sections"][
        "forced_structure"
    ]
    row = _row(rows, "id", POSTDICTION_ID)
    assert row["match"] == (
        "exact bounded modal factorization; local source-produced "
        "physical Maxwell theory not constructed"
    )
    assert row["lean_declarations"]["ModalMaxwellFactorizationBoundary"] == (
        LEAN_DECLARATIONS
    )
    assert row["artifact_refs"] == POSTDICTION_ARTIFACTS
    assert set(POSTDICTION_ARTIFACTS).issubset(CLAIM_EVIDENCE)
    boundary = row["hypothesis_boundary"]
    assert "PR-20, PR-21, and PR-22 remain declared inputs" in boundary
    assert "PR-53 and PR-54 name missing attachments" in boundary
    assert "Scientific owner #733 records" in boundary
    assert "emits no frozen prediction" in boundary


def test_claim_matrices_and_dependency_edge_are_exact() -> None:
    assert len([r for r in _csv("claims/novelty_matrix.csv") if r["claim_id"] == CLAIM_ID]) == 1
    falsification_rows = [
        r for r in _csv("claims/falsification_matrix.csv") if r["claim_id"] == CLAIM_ID
    ]
    assert len(falsification_rows) == 1
    assert (
        "same-sign pair has the oscillator sign on an amplitude with nonzero spatial action"
        in falsification_rows[0]["mathematical_falsifier"]
    )
    graph = _json("claims/dependency_graph.json")
    assert graph["nodes"].count(CLAIM_ID) == 1
    matches = [
        edge
        for edge in graph["edges"]
        if edge["from"] == "OPH-A5-SEAM-CURRENT-EDGE30"
        and edge["to"] == CLAIM_ID
    ]
    assert len(matches) == 1
    assert "supplies no local real field" in matches[0]["role"]


def test_no_frozen_prediction_or_marketing_surface_is_emitted() -> None:
    frozen = _json("claims/frozen_prediction_register.json")
    for row in frozen["rows"]:
        encoded = json.dumps(row, sort_keys=True)
        assert CLAIM_ID not in encoded
        assert POSTDICTION_ID not in encoded
        assert LEAN_PATH not in encoded

    forbidden_marker = "modal-maxwell-factorization-boundary"
    for relative_path in (
        "README.md",
        "README_FR.md",
        "flagship/from_observer_consensus_to_standard_physics.tex",
    ):
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        assert CLAIM_ID not in text
        assert forbidden_marker not in text.lower()


def test_screen_umbrella_imports_the_kernel_checked_module() -> None:
    umbrella = (ROOT / "Lean/Screen/OPHScreen.lean").read_text(encoding="utf-8")
    import_line = "import ModalMaxwellFactorizationBoundary"
    assert umbrella.splitlines().count(import_line) == 1
    lakefile = (ROOT / "Lean/lakefile.lean").read_text(encoding="utf-8")
    assert lakefile.count("`ModalMaxwellFactorizationBoundary") == 1
