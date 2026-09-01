"""Cross-surface gates for the Issue #750 fixed-federation execution classification."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CLASSIFICATION = "OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION"
PROGRESS = "OPH-COMPUTATION-FIXED-FEDERATION-PROGRESS"
CAPACITY = "OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_claim_and_premises_keep_the_exact_execution_boundary() -> None:
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }
    claim = claims[CLASSIFICATION]
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["premise_dependencies"]["consumed"] == ["PR-84", "PR-85"]
    assert "reverse is vacuous" in claim["statement"]
    assert "order-sharp Theta(n^2)" in claim["statement"]
    assert "Weak fairness alone supplies no uniform finite attempt horizon" in claim["statement"]
    assert "not an exact constant match" in claim["statement"]

    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    bounded = rows["PR-85"]
    assert bounded["disposition"] == "axiomatize"
    assert bounded["consuming_lanes"] == [750]
    assert "not implied by pathwise weak fairness" in bounded["statement"]
    assert "fixedRoundRobinScheduler" in bounded["notes"]


def test_dependency_graph_keeps_code_reuse_out_of_scientific_parentage() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]) for row in graph["edges"]}
    assert (PROGRESS, CLASSIFICATION) in edges
    assert (CAPACITY, CLASSIFICATION) not in edges


def test_lean_provider_controls_and_actual_consumers_are_wired() -> None:
    execution = _text("Lean/Computation/FixedFederationExecution.lean")
    complexity = _text("Lean/Computation/FixedFederationComplexity.lean")
    examples = _text("Lean/Computation/FixedFederationExecutionExamples.lean")
    generic = _text(
        "Lean/ObserverPatchHolography/Execution/RankedAttemptCapacity.lean"
    )
    adaptive = _text(
        "Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacity.lean"
    )
    fixed_tower = _text("Lean/Tower/FixedFederationExecutionEndpoint.lean")
    adaptive_tower = _text("Lean/Tower/CumulativeCapacityEndpoint.lean")

    for name in (
        "stableConsensusTail_nodePathwiseWeakFair",
        "nodePathwiseWeakFair_iff_eventuallyStableConsensus",
        "nodeMemberRecurrent_iff_nodeSiteRecurrent",
        "exists_mathematical_fair_scheduler",
        "roundRobinScheduler_selects_within",
    ):
        assert f"theorem {name}" in execution

    for name in (
        "fixedProgram_atMostOneDownstreamConsumer",
        "canonicalAcceptedStep_linearDefectRank_lt",
        "fixedProgram_acceptedSteps_quadratic",
        "fixedRoundRobin_nodeBoundedWaste",
        "fixedProgram_nodeBoundedWaste_triangle",
    ):
        assert f"theorem {name}" in complexity

    for name in (
        "recurrence_strictly_stronger_than_tail_fairness",
        "weakFair_no_uniform_attempt_bound",
        "negComb_quadratic_lower",
        "fixedProgram_sharp_quadratic_certificates",
    ):
        assert f"theorem {name}" in examples

    provider_call = "OPH.RankedAttempt.boundedWaste_eventually_quiescent"
    assert f"theorem boundedWaste_eventually_quiescent" in generic
    assert "OPH.RankedAttempt.boundedWaste_exists_quiescent_by_rank" in adaptive
    assert provider_call in fixed_tower
    assert provider_call in adaptive_tower
    assert "theorem canonicalAcceptedSteps_within_quadratic" in fixed_tower
    assert "fixedProgram_acceptedSteps_quadratic phi hsteps" in fixed_tower
    assert "theorem roundRobin_reaches_bounded_unique_output" in fixed_tower


def test_manuscripts_and_indices_keep_nonclaims_and_order_sharpness() -> None:
    component = " ".join(_text("extra/observable_normal_forms.tex").split())
    consensus = " ".join(_text("paper/reality_as_consensus_protocol.tex").split())
    flagship = " ".join(
        _text("flagship/from_observer_consensus_to_standard_physics.tex").split()
    )
    guide = _text("Lean/docs/LIBRARY_GUIDE.md")
    index = _text("Lean/docs/PROOF_INDEX.md")

    assert "Fixed-federation execution classification" in component
    assert "vacuous after stable consensus" in consensus
    assert "order-sharp" in consensus
    assert "bounded-waste" in flagship
    assert "FixedFederationExecution.lean" in guide
    assert "FixedFederationExecutionEndpoint.lean" in index
    assert "finite `PublicWorld`" in index


def test_umbrellas_register_each_new_module_once() -> None:
    root = _text("Lean/ObserverPatchHolography.lean")
    tower = _text("Lean/Tower.lean")
    lake = _text("Lean/lakefile.lean")
    assert root.count("import Computation.FixedFederationExecutionExamples") == 1
    assert root.count(
        "import ObserverPatchHolography.Execution.RankedAttemptCapacity"
    ) == 1
    assert tower.count("import Tower.FixedFederationExecutionEndpoint") == 1
    assert lake.count("Computation.FixedFederationExecutionExamples") == 1
