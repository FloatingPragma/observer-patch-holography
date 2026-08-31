"""Cross-surface gates for the Issue #750 cumulative attempt-capacity result."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CAPACITY = "OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY"
ENDPOINT = "OPH-CONS-CAPACITY-BOUNDED-PUBLIC-ENDPOINT"
ADAPTIVE = "OPH-CONS-ADAPTIVE-REPAIR-STRATIFICATION"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_claims_keep_capacity_assumptions_and_nonclaims_explicit() -> None:
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }

    capacity = claims[CAPACITY]
    assert capacity["owner_paper"] == "extra/observable_normal_forms.tex"
    assert capacity["assumptions"] == [
        "bounded_waste_adaptive_scheduler",
        "sufficient_cumulative_attempt_budget",
    ]
    assert capacity["premise_dependencies"]["consumed"] == ["PR-79", "PR-80"]
    assert "including equality stutter" in capacity["statement"]
    assert "selects no scheduler" in capacity["statement"]

    endpoint = claims[ENDPOINT]
    assert endpoint["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert endpoint["premise_dependencies"]["consumed"] == [
        "PR-01",
        "PR-79",
        "PR-80",
    ]
    assert "no larger than B" in endpoint["statement"]
    assert "supplies no scheduler" in endpoint["statement"]


def test_dependency_edges_keep_capacity_provider_and_consumer_load_bearing() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]) for row in graph["edges"]}
    assert {
        (ADAPTIVE, CAPACITY),
        (CAPACITY, ENDPOINT),
        ("OPH-PUBLIC-WORLD-ENDPOINT-FINITE", ENDPOINT),
    } <= edges


def test_premises_keep_scheduler_and_budget_sources_open() -> None:
    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    bounded = rows["PR-79"]
    budget = rows["PR-80"]

    assert bounded["disposition"] == "axiomatize"
    assert bounded["consuming_lanes"] == [750]
    assert "does not derive q or select a scheduler" in bounded["statement"]
    assert budget["disposition"] == "axiomatize"
    assert budget["consuming_lanes"] == [750]
    assert "no physical time" in budget["statement"]


def test_lean_provider_consumer_and_umbrellas_are_wired() -> None:
    provider = _text(
        "Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacity.lean"
    )
    examples = _text(
        "Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacityExamples.lean"
    )
    consumer = _text("Lean/Tower/CumulativeCapacityEndpoint.lean")

    assert "def BoundedWaste" in provider
    assert "theorem cumulativeGenuineChangeCost_le_initialMismatch" in provider
    assert "theorem boundedWaste_eventually_normal" in provider
    assert "theorem delayed_normalizing_attempt_no_go" in examples
    assert "theorem delayThenProbe_pathwiseWeakFair" in examples
    assert "theorem firstBroken_attempt_threshold_iff" in examples
    assert "theorem sharpInstances_patchCardinality_differs" in examples
    assert "theorem boundedWaste_public_endpoint_exists_unique" in consumer

    assert _text("Lean/ObserverPatchHolography.lean").count(
        "import ObserverPatchHolography.Execution.CumulativeAttemptCapacityExamples"
    ) == 1
    assert _text("Lean/Tower.lean").count(
        "import Tower.CumulativeCapacityEndpoint"
    ) == 1


def test_manuscripts_and_proof_indices_expose_the_same_boundary() -> None:
    component = _text("extra/observable_normal_forms.tex")
    consensus = _text("paper/reality_as_consensus_protocol.tex")
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    flagship_words = " ".join(flagship.split())
    proof_index = _text("Lean/docs/PROOF_INDEX.md")
    standalone_index = _text("Lean/ObservableNormalForms/PROOF_INDEX.md")

    assert "Cumulative attempt capacity and sharpness" in component
    assert "uniformly" in component
    assert "Capacity-bounded canonical endpoint" in consensus
    assert "bounded-waste upper horizon" in flagship_words
    assert "work-conserving case exactly" in flagship_words
    assert "Cumulative attempt capacity charges equality stutters" in proof_index
    assert "No stochastic scheduler, rate, physical clock, resource" in proof_index
    assert "Cumulative attempt-capacity classification" not in standalone_index
