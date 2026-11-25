"""Cross-surface gates for the source-derived causal-order and 1+3 result."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
ORDER_CLAIM = "OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE"
SPACETIME_CLAIM = "OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM"
POSET_PACKAGE = "evidence/source_net_causal_poset"


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _claims() -> dict[str, dict]:
    return {
        row["claim_id"]: row
        for row in _json("claims/claim_registry.yaml")["claims"]
    }


def test_source_order_claim_keeps_its_finite_boundary_explicit() -> None:
    claim = _claims()[ORDER_CLAIM]
    assert claim["owner_paper"] == "paper/reality_as_consensus_protocol.tex"
    assert claim["assumptions"] == [
        "certified_causal_read_support",
        "fresh_duplicate_free_complete_executed_event_carrier",
    ]
    assert claim["premise_dependencies"]["consumed"] == ["PR-81"]
    assert claim["premise_dependencies"]["boundary"] == ["PR-82"]
    assert "least strict transitive relation" in claim["statement"]
    assert "writer-blind" in claim["statement"]
    assert "causal-set-like" in claim["statement"]
    assert "count-to-volume" in claim["statement"]
    assert "attained maximum authenticated-parent-chain length" in claim["statement"]
    assert "finite causet compiler" in claim["statement"]
    assert (
        "Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean"
        in claim["evidence"]
    )




def test_dependency_graph_uses_the_new_two_input_carrier() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]): row["role"] for row in graph["edges"]}
    assert ("OPH-CONS-D1", ORDER_CLAIM) in edges
    assert (ORDER_CLAIM, SPACETIME_CLAIM) in edges
    assert ("OPH-GR-PORT-RESPONSE-COMPLETION", SPACETIME_CLAIM) in edges
    assert "OPH-GR-SOURCE-CAUSAL-HISTORY-FAMILY" not in graph["nodes"]
    assert ("OPH-GR-D3-CAP-H3", SPACETIME_CLAIM) not in edges
    assert "OPH-GR-FCC-CAUSET-COMPATIBILITY-CONTROL" not in graph["nodes"]


def test_lean_package_and_umbrellas_are_wired() -> None:
    provenance = _text(
        "Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean"
    )
    compiler = _text(
        "Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean"
    )
    cover_compiler = _text(
        "Lean/ObserverPatchHolography/Provenance/FiniteCausetCoverCompiler.lean"
    )
    history = _text(
        "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean"
    )
    geometry = _text("Lean/Geometry/EventPopulationChartInterface.lean")
    carrier = _text("Lean/Geometry/SourceDerivedSpacetimeCarrier.lean")
    frame_packet = _text("Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean")
    einstein = _text("Lean/Geometry/SourceOrderEinsteinComposition.lean")

    assert "def AuthenticatedDirectSemanticParent" in provenance
    assert "def generatedRecordOrder" in provenance
    assert "theorem generatedBefore_le_of_transitive" in provenance
    assert "theorem eq_generatedBefore_of_exact" in provenance
    assert "inductive ParentChainTo" in provenance
    assert "theorem sourceHeight_chain_attained" in provenance
    assert "theorem sourceHeight_eq_max_parentChainLength" in provenance
    assert "theorem authenticatedParent_commit_iff" in compiler
    assert "theorem generatedBefore_iff_transGen" in compiler
    assert "def predecessorRank" in compiler
    assert "theorem finiteStrictTransitiveRelation_realized" in compiler
    assert "def CoverRelation" in cover_compiler
    assert "theorem rel_iff_transGen_coverRelation" in cover_compiler
    assert "theorem coverLog_parentEdge_iff" in cover_compiler
    assert "theorem coverLog_generatedBefore_iff" in cover_compiler
    assert "theorem finiteStrictTransitiveRelation_cover_realized" in cover_compiler
    assert "theorem semanticEventLogOfExecution_generatedBefore_iff" in history
    assert "theorem execParents_rank_lt" in history

    assert "structure SourceDerivedCausalChartInterface" in geometry
    source_interface = geometry.split(
        "structure SourceDerivedCausalChartInterface", 1
    )[1].split("namespace SourceDerivedCausalChartInterface", 1)[0]
    assert "semanticLog" in source_interface
    assert "generated_cone_base" in source_interface
    assert "populate" not in source_interface
    assert "prec :" not in source_interface
    assert "separation :" not in source_interface
    assert "theorem separation_base" in geometry
    assert "theorem generatedBefore_iff_causalLE_and_ne" in geometry
    assert "structure EventPopulationChartInterface" not in geometry
    assert "SourceDerivedOrderEventPopulationChartInterface" not in geometry
    assert "def ofLegacy" not in geometry

    for declaration in (
        "abbrev SourceSpacetimeCarrier := ℝ × FrameQuotient",
        "theorem sourceSpacetimeCarrier_finrank",
        "theorem sourceCarrier_one_three_signature",
        "theorem sourceUnitNullVector_futureCausal",
        "structure RankSpatialCausalPlacement",
        "theorem generatedBefore_sourceCausalLE",
        "theorem eventPoint_injective_of_sameHeightSpatialInjective",
        "namespace EnumeratedCausalPlacement",
        "theorem finiteLog_has_separated_forwardCausalPlacement",
        "structure FaithfulRankSpatialCausalPlacement",
        "noncomputable def ofIncomparableSpacelike",
        "theorem generatedBeforeEq_iff_sourceCausalLE",
        "theorem eventPoint_injective",
        "def toSourceDerivedCausalChartInterface",
        "def faithfulPlacement",
        "theorem exactDiamondConeOrder",
        "theorem responses_spacelike",
        "theorem parentEdge_null",
    ):
        assert declaration in carrier

    assert "sourceCharts : SourceDerivedCausalChartInterface" in frame_packet
    assert "def SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement" in frame_packet
    assert "structure SourceDirectionEinsteinShapePremises" in einstein
    shape = einstein.split("structure SourceDirectionEinsteinShapePremises", 1)[1]
    shape = shape.split("namespace SourceDirectionEinsteinShapePremises", 1)[0]
    for removed in (
        "SourceOrderFrameCompatibilityPacket",
        "entropyStress",
        "universalCoupling",
        "newton",
        "referenceLambda",
        "vacuumReference",
        "physicalScale",
    ):
        assert removed not in shape
    assert "theorem constantMetricAmbiguity" in einstein
    assert "structure SourceIndexedEinsteinPremises" in einstein

    root_umbrella = _text("Lean/ObserverPatchHolography.lean")
    assert root_umbrella.splitlines().count(
        "import ObserverPatchHolography.Provenance.FiniteCausetCoverCompiler"
    ) == 1

    umbrella = _text("Lean/Geometry.lean")
    assert umbrella.count("import Geometry.SourceDerivedSpacetimeCarrier") == 1
    assert umbrella.count("import Geometry.SourceOrderFrameCompatibilityPacket") == 1
    assert umbrella.count("import Geometry.SourceOrderEinsteinComposition") == 1


def test_registry_carries_the_finite_carrier_and_conditional_continuum() -> None:
    claims = _claims()
    evidence_paths = {
        path
        for claim in claims.values()
        for path in claim.get("evidence", [])
    }
    for retired in (
        "code/geometry/event_manifold_reconstruction.py",
        "code/geometry/realized_event_receipts.py",
        "code/geometry/runs/realized_event_receipt_report.json",
        "code/geometry/bulk_depth_receipts.py",
        "code/geometry/runs/bulk_depth_receipt_report.json",
        "code/geometry/realized_branch_receipts.py",
        "code/geometry/runs/realized_branch_receipt_report.json",
    ):
        assert retired not in evidence_paths
    carrier = claims[SPACETIME_CLAIM]
    assert carrier["assumptions"][0] == "source_causal_continuum_certificate"
    assert "commit_populated_event_base" not in carrier["assumptions"]
    assert "SourceDerivedSpacetimeCarrier.lean" in " ".join(carrier["evidence"])
    assert "dimension 4" in carrier["statement"]
    assert "Lorentz inertia (1,3)" in carrier["statement"]
    assert "rank-four chart" in carrier["statement"]
    assert "calibrated count-volume" in carrier["statement"]
    assert "attained maximum authenticated-parent-chain length" in carrier["statement"]
    assert "enumeration-dependent injective placement" in carrier["statement"]
    assert "cross-height incomparable events can become cone-related" in carrier["statement"]
    assert (
        "Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean"
        in carrier["evidence"]
    )

    assert "OPH-GR-SOURCE-CAUSAL-HISTORY-FAMILY" not in claims
    assert "OPH-GR-FCC-CAUSET-COMPATIBILITY-CONTROL" not in claims
    assert f"{POSET_PACKAGE}/source_net_causal_limit_receipt.json" in carrier["evidence"]
    assert f"{POSET_PACKAGE}/carrier_source_net_receipt.json" in carrier["evidence"]
    assert "Lean/Geometry/SourceNetLayeredOrder.lean" in carrier["evidence"]
    assert "Lean/Geometry/OrderingFractionFourDimensional.lean" in carrier["evidence"]

    charts = claims["OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE"]
    assert "no population map" in charts["statement"]
    assert "FaithfulRankSpatialCausalPlacement" in charts["statement"]
    assert "Separation is therefore a theorem rather than a field" in charts["statement"]
    assert "declared_exact_base_order_cone_clause" in charts["assumptions"]
    assert "declared_separation_and_cone_clauses" not in charts["assumptions"]
    assert "source_native_finite_causal_chart_interface" in charts["novelty_type"]

    shape = claims["OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION"]
    assert "SourceDirectionEinsteinShapePremises" in shape["statement"]
    assert "no source-order packet" in shape["statement"]
    assert "separate SourceIndexedEinsteinPremises adapter" in shape["statement"]




def test_causal_poset_package_carries_the_finite_readings() -> None:
    family = _json(f"{POSET_PACKAGE}/source_net_causal_limit_receipt.json")
    assert family["schema"] == "oph.exact.source-net-causal-limit.v1"
    assert [level["q"] for level in family["levels"]] == [5, 8, 13, 21, 34, 55]
    current = _json(f"{POSET_PACKAGE}/source_net_causal_limit_receipt_2026-09-24.json")
    assert [level["q"] for level in current["levels"]] == [5, 8, 13, 21, 34, 55, 89]
    assert current["rer_cross_check"]["all_agree"] is True
    for level in current["levels"]:
        for fam in level["families"]:
            assert fam["vertical_pair_counting"] == "exact_all_pairs"
    q89 = [fam for fam in current["levels"][-1]["families"] if fam["dimension"] == 3][0]
    assert 4.0 < q89["vertical_intervals"][-1]["myrheim_meyer_dimension"] < 4.2
    current_carrier = _json(f"{POSET_PACKAGE}/carrier_source_net_receipt_2026-09-24.json")
    assert [level["q"] for level in current_carrier["levels"]] == [5, 8, 13, 21, 34, 55, 89]
    for level in current_carrier["levels"]:
        assert level["provenance"]["derived_rank_equals_round"] is True
        assert level["intervention"]["equals_future_cone_all_rounds"] is True
    assert family["rer_cross_check"]["all_agree"] is True
    assert family["scope"]["native_repair_selected"] is False
    assert family["scope"]["physical_clock_or_spacetime_identified"] is False
    three = {
        level["q"]: [fam for fam in level["families"] if fam["dimension"] == 3][0]
        for level in family["levels"]
    }
    for q, low, high in ((13, 4.0, 4.3), (34, 4.0, 4.2), (55, 4.0, 4.2)):
        assert low < three[q]["vertical_intervals"][-1]["myrheim_meyer_dimension"] < high
    carrier = _json(f"{POSET_PACKAGE}/carrier_source_net_receipt.json")
    assert carrier["schema"] == "oph.exact.carrier-source-net.v1"
    assert carrier["readback_metric"]["scale_to_paper_position_s"] == "1"
    for level in carrier["levels"]:
        assert level["neighbours"]["equals_source_net_digest"] is True
        assert level["provenance"]["derived_rank_equals_round"] is True
        assert level["intervention"]["equals_future_cone_all_rounds"] is True
    readme = _text(f"{POSET_PACKAGE}/README.md")
    assert "layerPrec_isPartialOrder" in readme
    assert "orderingFraction_eq" in readme
    assert (ROOT / POSET_PACKAGE / "verify_causal_poset_archive.py").is_file()
    assert (ROOT / POSET_PACKAGE / "build_causal_poset.py").is_file()


def test_scientific_registry_cites_the_causal_poset_package() -> None:
    claims = _claims()
    carrier_evidence = claims[SPACETIME_CLAIM]["evidence"]
    assert f"{POSET_PACKAGE}/verify_causal_poset_archive.py" in carrier_evidence
    local_domain_evidence = claims["OPH-GR-FINITE-LOCAL-DOMAIN"]["evidence"]
    assert f"{POSET_PACKAGE}/source_net_causal_limit_receipt.json" in local_domain_evidence
    assert not any("causet_likeness" in item for item in local_domain_evidence)
    assert not (ROOT / "evidence/causet_likeness").exists()


@pytest.fixture
def causal_archive_checker():
    path = ROOT / POSET_PACKAGE / "verify_causal_poset_archive.py"
    spec = importlib.util.spec_from_file_location("causal_archive_checker", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_causal_archive_fresh_metadata_and_poset_replay(causal_archive_checker) -> None:
    assert causal_archive_checker.main() == 0


@pytest.mark.parametrize("mutation", (
    "family_digest", "carrier_digest", "producer_digest", "missing_producer",
    "extra_producer", "duplicate_research", "inventory_digest", "total_bytes",
    "duplicate_inventory", "floating_size",
))
def test_causal_archive_rejects_contradictory_metadata(
    mutation, causal_archive_checker, tmp_path, monkeypatch,
) -> None:
    manifest = _json(f"{POSET_PACKAGE}/archive_manifest.json")
    if mutation in ("family_digest", "carrier_digest"):
        key = "source_net_receipt_sha256" if mutation == "family_digest" else "carrier_receipt_sha256"
        manifest["result"][key] = "sha256:" + "0" * 64
    elif mutation == "producer_digest":
        manifest["source"]["producer_files"]["tests/test_exact_carrier_source_net.py"] = "0" * 64
    elif mutation == "missing_producer":
        del manifest["source"]["producer_files"]["tests/test_exact_carrier_source_net.py"]
    elif mutation == "extra_producer":
        manifest["source"]["producer_files"]["oph_exact/invented.py"] = "0" * 64
    elif mutation == "duplicate_research":
        paths = manifest["source"]["research_files_pinned_by_the_receipts"]
        paths.append(paths[0])
    elif mutation == "inventory_digest":
        manifest["curated_archive"]["inventory_sha256"] = "0" * 64
    elif mutation == "total_bytes":
        manifest["curated_archive"]["total_bytes"] += 1
    elif mutation == "duplicate_inventory":
        manifest["inventory"].append(manifest["inventory"][0].copy())
    else:
        manifest["inventory"][0]["bytes"] = float(manifest["inventory"][0]["bytes"])
    candidate = tmp_path / "manifest.json"
    candidate.write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr(causal_archive_checker, "MANIFEST", candidate)
    with pytest.raises(SystemExit, match="FAIL:"):
        causal_archive_checker.main()


@pytest.mark.parametrize("mutation", ("conflicting_shared_source", "family_attachment", "log_attachment"))
def test_causal_archive_rejects_conflicting_receipt_metadata(
    mutation, causal_archive_checker,
) -> None:
    manifest = _json(f"{POSET_PACKAGE}/archive_manifest.json")
    family = _json(f"{POSET_PACKAGE}/source_net_causal_limit_receipt.json")
    carrier = _json(f"{POSET_PACKAGE}/carrier_source_net_receipt.json")
    key = {
        "conflicting_shared_source": "oph_exact/source_net.py",
        "family_attachment": "data/exact/source_net_causal_limit_receipt.json",
        "log_attachment": "data/exact/carrier_source_net_logs/q5_event_log.json.gz",
    }[mutation]
    carrier["source_pins"][key] = "0" * 64
    generation = manifest["historical_generations"][0]
    with pytest.raises(SystemExit, match="FAIL:"):
        causal_archive_checker.check_metadata(
            generation["source"], generation["result"], family, carrier, generation["family"], generation["carrier"],
        )


def test_causal_archive_rejects_duplicate_json_metadata(causal_archive_checker, tmp_path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"commit":"first","commit":"second"}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate key"):
        causal_archive_checker.strict_json(path)


def test_publication_surfaces_reject_the_retired_gluing_story() -> None:
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    gravity = _text(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"
    )
    combined = "\n".join((flagship, gravity))
    assert "thm:source-causal-precursor" in flagship
    assert "thm:source-causal-continuum" in flagship
    assert "source-derived finite" in combined.lower()
    assert "source-causal continuum" in combined.lower()
    for stale in (
        "Sufficient gluing criterion",
        "Identifying commits with record germs is a named premise",
        "event_manifold_receipts_E1_E6",
        "commit_populated_event_base",
        "OPH-GR-FCC-CAUSET-COMPATIBILITY-CONTROL",
    ):
        assert stale not in combined

    publication_surfaces = tuple(
        ROOT.glob("flagship/*.tex")
    ) + tuple(ROOT.glob("paper/*.tex")) + tuple(ROOT.glob("cosmology/*.tex")) + (
        ROOT / "paper/archive/observer_patch_holography_as_string_vacuum_selector.tex",
    )
    all_publication_text = "\n".join(
        path.read_text(encoding="utf-8") for path in publication_surfaces
    )
    for stale in (
        "Sufficient gluing criterion",
        "populated record-germ",
        "event_manifold_receipts_E1_E6",
        "commit_populated_event_base",
        "SOURCE_NATIVE_EDGE_LOCAL_RANK3_NULL_CONE_DECORATION_RECEIPT",
        "OPH-GR-FCC-CAUSET-COMPATIBILITY-CONTROL",
        "parent-free",
    ):
        assert stale not in all_publication_text
    for tracking_phrase in (
        "github.com/FloatingPragma/observer-patch-holography/issues/",
        "issue #",
        "completion plan",
        "work in progress",
    ):
        assert tracking_phrase.lower() not in all_publication_text.lower()
