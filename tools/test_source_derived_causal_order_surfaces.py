"""Cross-surface gates for the Issue #763 source-derived causal-order result."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CLAIM = "OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_claim_keeps_assumptions_and_nonclaims_explicit() -> None:
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }

    claim = claims[CLAIM]
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["assumptions"] == [
        "certified_causal_read_support",
        "append_only_ancestry_rank",
    ]
    assert claim["premise_dependencies"]["consumed"] == ["PR-81", "PR-82"]
    assert "least strict transitive relation" in claim["statement"]
    assert "writer-blind" in claim["statement"]
    assert "informational" in claim["statement"]


def test_premises_keep_support_and_rank_sources_open() -> None:
    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    support = rows["PR-81"]
    rank = rows["PR-82"]

    assert support["disposition"] == "axiomatize"
    assert support["consuming_lanes"] == [763]
    assert "nothing derives the support from executor traces" in support["statement"]
    assert rank["disposition"] == "axiomatize"
    assert rank["consuming_lanes"] == [763]
    assert "append-only" in rank["statement"]


def test_dependency_edge_keeps_transactional_provider() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]) for row in graph["edges"]}
    assert ("OPH-CONS-D1", CLAIM) in edges


def test_lean_package_and_umbrellas_are_wired() -> None:
    core = _text(
        "Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean"
    )
    mismatch = _text(
        "Lean/ObserverPatchHolography/Provenance/MismatchProvenance.lean"
    )
    seam = _text("Lean/ObserverPatchHolography/Provenance/SeamDeltaAggregation.lean")
    interval = _text("Lean/ObserverPatchHolography/Provenance/CausalInterval.lean")
    bridge = _text("Lean/QFT/SourceDerivedEventPrecedence.lean")

    history = _text(
        "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean"
    )
    quotient = _text(
        "Lean/ObserverPatchHolography/Provenance/QuotientInvariance.lean"
    )
    refinement = _text(
        "Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean"
    )

    assert "def generatedRecordOrder" in core
    assert "theorem generatedBefore_le_of_transitive" in core
    assert "theorem eq_generatedBefore_of_exact" in core
    assert "theorem orientation_not_determined" in mismatch
    assert "theorem repairs_without_reading" in mismatch
    assert "theorem aggregate_unique_of_touches" in seam
    assert "theorem separate_writers" in seam
    assert "theorem responses_incomparable" in interval
    assert "theorem interval_injection_answer_eq_univ" in interval
    assert "theorem chainWitnessPrecedenceAdapter_exact" in bridge
    assert "theorem swap_execParents" in history
    assert "theorem eqvGen_independentSwap_invariant" in history
    assert "theorem execParents_rank_lt" in history
    assert "theorem unconsumed_mismatch_persists" in history
    assert "theorem dependent_swap_changes_edges" in history
    assert "theorem hiddenSpec_no_visible_edge" in quotient
    assert "theorem hiddenSpec_not_changes" in quotient
    assert "theorem generatedBefore_natural" in refinement
    assert "theorem reversed_map_not_natural" in refinement

    assert _text("Lean/ObserverPatchHolography.lean").count(
        "import ObserverPatchHolography.Provenance.CausalInterval"
    ) == 1
    assert _text("Lean/QFT.lean").count(
        "import QFT.SourceDerivedEventPrecedence"
    ) == 1


def test_manuscript_and_indices_expose_the_same_boundary() -> None:
    consensus = _text("paper/reality_as_consensus_protocol.tex")
    consensus_words = " ".join(consensus.split())
    proof_index = _text("Lean/docs/PROOF_INDEX.md")
    guide = _text("Lean/docs/LIBRARY_GUIDE.md")

    assert "Source-derived event precedence from mismatch provenance" in consensus
    assert "a verifier rather than a choice" in consensus_words
    assert "determines no causal arrow" in consensus_words
    assert "carry their own receipts and are not supplied here" in consensus_words
    assert "Source-derived event precedence" in proof_index
    assert "source-derived causal-order package" in guide


def test_composition_surfaces_state_the_generated_order() -> None:
    shared = _text("paper/tex_fragments/PAPER.tex")
    shared_words = " ".join(shared.split())
    flagship = " ".join(
        _text("flagship/from_observer_consensus_to_standard_physics.tex").split()
    )
    observers = " ".join(_text("paper/observers_are_all_you_need.tex").split())
    consensus = " ".join(
        _text("paper/reality_as_consensus_protocol.tex").split()
    )
    micro = " ".join(
        _text("paper/screen_microphysics_and_observer_synchronization.tex").split()
    )
    dark = " ".join(
        _text("cosmology/oph_dark_matter_paper.tex").split()
    )

    assert "finite halves of" in shared_words
    assert "That commit identification is a named premise" in shared_words
    assert "provenance-disciplined" in shared_words
    assert "an exact adapter is its verifier rather than a choice of it" in consensus
    assert "bombelli1987" in consensus
    assert "lamport1978" in consensus
    assert "reflexive transitive closure of certified read-from parenthood" in observers
    assert "generated transitive closure of certified read-from provenance" in flagship
    assert "causal set generated from record provenance rather than posited" in flagship
    assert "provenance-disciplined" in micro
    assert "long-lived by theorem" in dark


def test_registry_carries_the_commit_identification_premise() -> None:
    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    base = rows["PR-83"]
    assert base["disposition"] == "axiomatize"
    assert base["consuming_lanes"] == [763]
    assert "Nothing derives this identification" in base["statement"]
    claims = {
        row["claim_id"]: row for row in _json("claims/claim_registry.yaml")["claims"]
    }
    manifold = claims["OPH-GR-D4B-EVENT-MANIFOLD"]
    assert "PR-83" in manifold["premise_dependencies"]["consumed"]
    assert "commit_populated_event_base" in manifold["assumptions"]
