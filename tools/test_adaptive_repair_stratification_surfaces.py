"""Cross-surface regression gates for the Issue #750 adaptive repair result."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

STRATIFICATION = "OPH-CONS-ADAPTIVE-REPAIR-STRATIFICATION"
PUBLIC_ENDPOINT = "OPH-CONS-ADAPTIVE-PUBLIC-ENDPOINT"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_claim_contracts_keep_assumptions_and_evidence_explicit() -> None:
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }

    stratification = claims[STRATIFICATION]
    assert stratification["owner_paper"] == "extra/observable_normal_forms.tex"
    assert stratification["assumptions"] == [
        "pathwise_weak_fair_adaptive_scheduler",
        "work_conserving_adaptive_scheduler",
    ]
    assert stratification["premise_dependencies"]["consumed"] == ["PR-77", "PR-78"]
    assert stratification["evidence"] == [
        "extra/observable_normal_forms.tex",
        "Lean/ObserverPatchHolography/Execution/AdaptiveRunStratification.lean",
        "Lean/ObserverPatchHolography/Execution/AdaptiveRunCounterexamples.lean",
    ]
    assert "selects no scheduler" in stratification["statement"]

    endpoint = claims[PUBLIC_ENDPOINT]
    assert endpoint["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert endpoint["assumptions"] == [
        "pathwise_weak_fair_adaptive_scheduler",
        "repair_completeness",
        "local_diamond",
    ]
    assert endpoint["premise_dependencies"]["consumed"] == ["PR-01", "PR-77"]
    assert endpoint["evidence"] == [
        "extra/observable_normal_forms.tex",
        "paper/reality_as_consensus_protocol.tex",
        "flagship/from_observer_consensus_to_standard_physics.tex",
        "Lean/ObserverPatchHolography/Execution/AdaptiveRunStratification.lean",
        "Lean/Tower/AdaptiveFixedPointEndpoint.lean",
        "Lean/docs/A4_PUBLIC_WORLD_ENDPOINT.md",
    ]
    assert "does not derive fairness" in endpoint["statement"]


def test_dependency_edges_keep_provider_and_consumer_load_bearing() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]) for row in graph["edges"]}
    assert {
        ("OPH-CONS-D1", STRATIFICATION),
        ("OPH-ADAPTIVE-SCHEDULER-LOCALITY-HELPER", STRATIFICATION),
        (STRATIFICATION, PUBLIC_ENDPOINT),
        ("OPH-PUBLIC-WORLD-ENDPOINT-FINITE", PUBLIC_ENDPOINT),
    } <= edges


def test_premises_keep_scheduler_selection_open() -> None:
    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    fairness = rows["PR-77"]
    work_conservation = rows["PR-78"]

    assert fairness["consuming_lanes"] == [750]
    assert fairness["disposition"] == "axiomatize"
    assert "not derived or selected by canonical repair" in fairness["statement"]
    assert work_conservation["consuming_lanes"] == [750]
    assert work_conservation["disposition"] == "axiomatize"
    assert "not required for premise-free eventual constancy" in work_conservation["statement"]


def test_public_imports_and_manuscripts_expose_the_same_boundary() -> None:
    assert _text("Lean/ObserverPatchHolography.lean").count(
        "import ObserverPatchHolography.Execution.AdaptiveRunCounterexamples"
    ) == 1
    assert _text("Lean/Tower.lean").count(
        "import Tower.AdaptiveFixedPointEndpoint"
    ) == 1

    component = _text("extra/observable_normal_forms.tex")
    consensus = _text("paper/reality_as_consensus_protocol.tex")
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")

    assert "Adaptive canonical repair: stutter, fairness, and endpoint promotion" in component
    assert "canonical adaptive repair is eventually constant without a fairness" in component
    assert "Adaptive repair attempts: stabilization, fairness, and work conservation" in consensus
    assert "Canonical endpoint for weak-fair adaptive repair" in consensus
    assert "These results select no scheduler" in flagship
    assert "provide no stochastic hitting" in flagship
