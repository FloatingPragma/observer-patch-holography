"""Cross-surface gates for the Issue #750 fixed-federation result."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PROGRESS = "OPH-COMPUTATION-FIXED-FEDERATION-PROGRESS"
ENDPOINT = "OPH-COMPUTATION-FIXED-FEDERATION-ENDPOINT"
BOUNDARY = "OPH-COMPUTATION-REPAIR-PROGRAM-BOUNDARY"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_claims_and_premise_keep_the_exact_fairness_boundary() -> None:
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }
    progress = claims[PROGRESS]
    endpoint = claims[ENDPOINT]

    assert (
        progress["novelty_type"]
        == "domain_specific_application_delta_with_sharp_weak_step_no_go"
    )
    assert progress["premise_dependencies"]["consumed"] == ["PR-81"]
    assert endpoint["premise_dependencies"]["consumed"] == ["PR-81"]
    assert "historical weak RepairStep is unchanged" in progress["statement"]
    assert "no finite PublicWorld presentation" in endpoint["statement"]

    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    fairness = rows["PR-81"]
    assert fairness["disposition"] == "axiomatize"
    assert fairness["consuming_lanes"] == [750]
    assert "does not cause accepted-step termination" in fairness["statement"]


def test_dependency_edges_connect_boundary_provider_and_consumer() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]) for row in graph["edges"]}
    assert {(BOUNDARY, PROGRESS), (PROGRESS, ENDPOINT)} <= edges
    assert ("OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY", PROGRESS) not in edges


def test_provider_counterexamples_consumer_and_umbrellas_are_wired() -> None:
    provider = _text("Lean/Computation/FixedFederationProgress.lean")
    controls = _text("Lean/Computation/FixedFederationCounterexamples.lean")
    consumer = _text("Lean/Tower/FixedFederationEndpoint.lean")

    for name in (
        "canonicalAcceptedStep_defectRank_lt",
        "canonicalAcceptedStep_wellFounded",
        "attemptRun_eventually_constant",
        "nodePathwiseWeakFair_eventually_consensus",
        "fixed_federation_fair_universality",
        "fixedObserverEndpointUniqueOutput",
    ):
        assert f"theorem {name}" in provider

    for name in (
        "weakNodeAttempt_to_repairStep",
        "weakAttemptRun_isRepairStepRun",
        "weak_fair_stuttering_no_go",
        "singletonStart_not_consensus",
        "singletonCanonical_positive",
        "fanout_defect_count_increases",
        "fanoutRoot_globallyEnabledFair",
        "globallyEnabledFair_no_go",
        "fanoutRoot_not_pathwiseWeakFair",
    ):
        assert f"theorem {name}" in controls

    assert "theorem fairAttempt_reaches_unique_output" in consumer
    assert _text("Lean/ObserverPatchHolography.lean").count(
        "import Computation.FixedFederationCounterexamples"
    ) == 1
    assert _text("Lean/Tower.lean").count(
        "import Tower.FixedFederationEndpoint"
    ) == 1


def test_manuscripts_and_proof_indices_expose_the_same_nonclaims() -> None:
    component = _text("extra/observable_normal_forms.tex")
    consensus = _text("paper/reality_as_consensus_protocol.tex")
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    proof_index = _text("Lean/docs/PROOF_INDEX.md")
    standalone_index = _text("Lean/ObservableNormalForms/PROOF_INDEX.md")

    assert "Fixed computation federations" in component
    assert "fixed-federation continuation" in consensus.lower()
    assert "generated formula nodes form one federation independent of input" in flagship.lower()
    assert "historical patch-frame relation" in component.lower()
    assert "Fixed computation-federation progress" in proof_index
    assert "finite `PublicWorld`" in proof_index
    assert "Fixed computation-federation progress" not in standalone_index
    assert "Fixed computation output endpoint" not in standalone_index
