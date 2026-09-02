"""Cross-surface gates for the source-derived causal-order and 1+3 result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORDER_CLAIM = "OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE"
SPACETIME_CLAIM = "OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM"
HISTORY_CLAIM = "OPH-GR-SOURCE-CAUSAL-HISTORY-FAMILY"


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


def test_premise_register_retires_the_old_event_base_and_exposes_continuum() -> None:
    rows = {row["id"]: row for row in _json("tracking/premise_register.json")["rows"]}
    assert "PR-83" not in rows
    support = rows["PR-81"]
    rank = rows["PR-82"]
    continuum = rows["PR-16"]
    assert support["disposition"] == "axiomatize"
    assert "does not infer support completeness" in support["statement"]
    assert rank["disposition"] == "axiomatize"
    assert "does not import this rank" in rank["statement"]
    assert "not a repair schedule" in rank["statement"]
    assert continuum["name"] == (
        "source-causal continuum and stable-causality certificate"
    )
    for phrase in (
        "physical event and link identification",
        "two-way order-cone-faithful placements",
        "calibrated event count to spacetime volume",
        "manifoldlikeness",
        "stable-topology",
        "convergent curvature estimators",
    ):
        assert phrase in continuum["statement"]
    assert "code/geometry/event_manifold_reconstruction.py" not in continuum["evidence"]


def test_dependency_graph_uses_the_new_two_input_carrier_and_history_control() -> None:
    graph = _json("claims/dependency_graph.json")
    edges = {(row["from"], row["to"]): row["role"] for row in graph["edges"]}
    assert ("OPH-CONS-D1", ORDER_CLAIM) in edges
    assert (ORDER_CLAIM, SPACETIME_CLAIM) in edges
    assert ("OPH-GR-PORT-RESPONSE-COMPLETION", SPACETIME_CLAIM) in edges
    assert (ORDER_CLAIM, HISTORY_CLAIM) in edges
    assert (SPACETIME_CLAIM, HISTORY_CLAIM) in edges
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

    history = claims[HISTORY_CLAIM]
    assert "24 to 384 events" in history["statement"]
    assert "Width is exactly 2" in history["statement"]
    assert "0.992221496953873" in history["statement"]
    assert "no admissible global scale" in history["statement"]
    assert history["claim_class"] == "emitted_artifact"
    assert (
        "evidence/source_causal_history_family/source_causal_history_family_publication_projection.json"
        in history["evidence"]
    )
    assert "OPH-GR-FCC-CAUSET-COMPATIBILITY-CONTROL" not in claims

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


def test_observation_ledger_carries_exact_lean_gains_without_physical_promotion() -> None:
    rows = {
        row["id"]: row for row in _json("tracking/observation_ledger.json")["rows"]
    }
    compiler = (
        "Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean"
    )
    order = rows["OL-A3"]
    carrier = rows["OL-A4"]
    assert compiler in order["evidence"]
    assert compiler in carrier["evidence"]
    assert "attained maximum authenticated-parent-chain length" in order["notes"]
    assert "grammar expressivity" in order["notes"]
    assert "Exact base order-cone agreement now derives chart separation" in order["notes"]
    assert "separated one-way realization" in carrier["notes"]
    assert "not source-selected or order-reflecting" in carrier["notes"]
    assert "calibrated count-volume" in carrier["notes"]
    assert carrier["status"] == "partial"


def test_history_receipt_records_positive_custody_and_negative_placement() -> None:
    receipt = _json(
        "evidence/source_causal_history_family/source_causal_history_family_publication_projection.json"
    )
    assert receipt["schema"] == (
        "oph.source-causal-history-family-publication-projection.v1"
    )
    assert receipt["round_cutoffs"] == [4, 8, 16, 32, 64]
    assert receipt["INFORMATIONAL_INDEPENDENT_CUTOFF_GENERATION_RECEIPT"] is True
    assert receipt["all_cutoffs_independently_generated"] is True
    assert len(receipt["cutoff_run_evidence_sha256s"]) == 5
    scaling = receipt["scaling_diagnostic"]
    assert scaling["event_counts"] == [24, 48, 96, 192, 384]
    assert scaling["widths"] == [2, 2, 2, 2, 2]
    assert scaling["ordering_fractions"][-1] == 0.992221496953873
    assert all(
        level["generated_from_own_cutoff_capture"] is True
        and level["independent_cutoff_run_evidence_sha256"]
        == evidence_sha256
        for level, evidence_sha256 in zip(
            receipt["levels"],
            receipt["cutoff_run_evidence_sha256s"],
            strict=True,
        )
    )
    cone = receipt["prescribed_single_frame_source_port_placement"]
    assert cone["event_count"] == 24
    assert cone["comparable_pair_count"] == 244
    assert cone["incomparable_pair_count"] == 32
    assert cone["injective_four_coordinate_map"] is False
    assert cone["precedence_iff_future_causal"] is False
    assert cone["global_time_scale_interval_nonempty"] is False
    assert cone["inter_carrier_frame_gluing_source_derived"] is False
    assert cone["consumed_record_barycentre_rule_source_derived"] is False
    assert cone["other_source_selected_placements_excluded"] is False
    assert cone["physical_no_go_for_other_source_selected_placements"] is False
    assert "not a no-go" in cone["interpretation"]
    assert receipt["promotion_and_nonclaim_flags"]["physical_promotion_allowed"] is False
    full_path = ROOT / (
        "evidence/source_causal_history_family/"
        "source_causal_history_family_receipt.json"
    )
    full_raw = full_path.read_bytes()
    full_receipt = json.loads(full_raw.decode("ascii"))
    assert receipt["full_receipt_file_sha256"] == (
        "sha256:" + hashlib.sha256(full_raw).hexdigest()
    )
    assert receipt["full_receipt_report_sha256"] == full_receipt["report_sha256"]
    verifier = _text(
        "evidence/source_causal_history_family/verify_source_causal_history_family_projection.py"
    )
    assert "reconstructs every cutoff separately" in " ".join(verifier.split())
    assert "EXPECTED_ORDERING_FRACTIONS" in verifier


def test_every_simulation_claiming_surface_points_to_current_rer_evidence() -> None:
    source_history_surfaces = (
        "flagship/from_observer_consensus_to_standard_physics.tex",
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex",
        "paper/screen_microphysics_and_observer_synchronization.tex",
        "paper/observers_are_all_you_need.tex",
        "paper/reality_as_consensus_protocol.tex",
        "paper/tex_fragments/PAPER.tex",
        "paper/tex_fragments/EINSTEIN_DERIVATION_SCOPE.tex",
        "paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex",
        "extra/observer_patch_holography_as_string_vacuum_selector.tex",
        "book/08-chapter-07.md",
        "book/35-chapter-34.md",
        "book/45-appendix-c.md",
        "README.md",
        "README_FR.md",
    )
    causet_control_surfaces = source_history_surfaces + (
        "cosmology/oph_cosmology_finite_source_cmb_program.tex",
        "cosmology/oph_dark_matter_paper.tex",
        "cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.tex",
    )
    for relative in source_history_surfaces:
        assert "evidence/source_causal_history_family" in _text(relative), relative
    for relative in causet_control_surfaces:
        assert "evidence/causet_likeness" in _text(relative), relative

    claims = _claims()
    history_evidence = claims[HISTORY_CLAIM]["evidence"]
    assert (
        "evidence/source_causal_history_family/"
        "source_causal_history_family_receipt.json"
    ) in history_evidence
    assert (
        "evidence/source_causal_history_family/"
        "verify_source_causal_history_family_projection.py"
    ) in history_evidence
    local_domain_evidence = claims["OPH-GR-FINITE-LOCAL-DOMAIN"]["evidence"]
    assert "evidence/causet_likeness/causet_likeness_receipt.json" in local_domain_evidence
    assert (
        "evidence/causet_likeness/verify_causet_likeness_receipt.py"
        in local_domain_evidence
    )


def test_publication_surfaces_reject_the_retired_gluing_story() -> None:
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    shared = _text("paper/tex_fragments/PAPER.tex")
    gravity = _text(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"
    )
    combined = "\n".join((flagship, shared, gravity))
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
        ROOT / "extra/observer_patch_holography_as_string_vacuum_selector.tex",
        ROOT / "paper/tex_fragments/PAPER.tex",
        ROOT / "paper/tex_fragments/EINSTEIN_DERIVATION_SCOPE.tex",
        ROOT / "paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex",
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
