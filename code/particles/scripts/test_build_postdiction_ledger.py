"""Tests for the postdiction ledger aggregator."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import build_postdiction_ledger as ledger


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ledger")
    return ledger.build(tmp / "postdiction_ledger.json", tmp / "POSTDICTION_LEDGER.md")


@pytest.fixture(scope="module")
def whitney_rows():
    rows = ledger._whitney_dynamics_rows()
    assert len(rows) == 3
    return {row["id"]: row for row in rows}


def test_guards_are_compare_only(result):
    guards = result["guards"]
    assert guards["compare_only"] is True
    assert guards["public_promotion_allowed"] is False
    assert guards["changes_any_solve_path"] is False
    assert guards["hand_typed_measured_values"] is False


def test_all_sections_present(result):
    expected = {
        "forced_structure",
        "quantum_carrier_status",
        "alpha",
        "charged_leptons",
        "electroweak",
        "quarks",
        "hadrons",
        "neutrinos",
    }
    assert set(result["sections"]) == expected
    for rows in result["sections"].values():
        assert rows


def test_quantum_carrier_status_preserves_classical_count_without_promotion(result):
    status = result["sections"]["quantum_carrier_status"]
    assert status["classical_mode_vector_order"] == [
        "photon",
        "gluon",
        "graviton",
    ]
    assert status["classical_mode_vector"] == [2, 16, 2]
    assert status["artifact_ref"] == (
        "code/particles/runs/status/quantum_carrier_status.json"
    )
    rows = {row["carrier_id"]: row for row in status["rows"]}
    assert rows["photon"]["verdict"] == (
        "NOT_EVALUABLE_NO_SOURCE_SELECTED_MAXWELL_QUANTUM_SECTOR"
    )
    assert rows["gluon"]["verdict"] == "NOT_EVALUABLE_NO_QCD"
    assert rows["graviton"]["verdict"] == (
        "NOT_EVALUABLE_NO_PHYSICAL_SOURCE_CAUSAL_CONTINUUM_OR_QUANTUM_CARRIER"
    )
    assert rows["photon"]["blocking_frontier"] == [
        "source_selected_unbroken_u1_quantum_maxwell_sector",
        "finite_source_to_lorentzian_quantum_eft_construction",
    ]
    assert rows["gluon"]["blocking_frontier"] == [
        "finite_source_to_lorentzian_quantum_eft_construction",
        "source_derived_qcd_physical_spectral_sector",
    ]
    assert rows["graviton"]["blocking_frontier"] == [
        "source_selected_faithful_3p1_causal_manifold_and_smooth_einstein_limit",
        "finite_source_to_lorentzian_linearized_quantum_carrier",
    ]


def test_forced_structure_receipts_exist(result):
    for row in result["sections"]["forced_structure"]:
        for ref in row.get("lean_receipts", []):
            assert (ledger.REPO / ref).exists(), ref
        if "artifact_ref" in row:
            assert (ledger.REPO / row["artifact_ref"]).exists()
        for ref in row.get("artifact_refs", []):
            assert (ledger.REPO / ref).exists()


def test_finite_unitary_limit_row_records_missing_scattering_scope(result):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    row = rows["finite_unitary_scattering_limit_no_go"]
    assert row["match"] == (
        "exact direct-power obstruction; physical scattering and "
        "asymptotic comparison construction not constructed"
    )
    assert row["lean_declarations"]["FiniteUnitaryScatteringNoGo"] == [
        "tendsto_powers_forces_identity",
        "nontrivial_powers_have_no_limit",
        "finite_unitary_powers_have_no_limit",
        "finite_unitary_ambient_powers_have_no_limit",
        "identical_relative_evolution_is_constant",
        "identical_relative_evolution_tendsto",
    ]
    boundary = row["hypothesis_boundary"]
    for route in (
        "relative or comparison dynamics",
        "selected projected scalar or observable limits",
        "infinite-dimensional weak limits",
        "Cesaro limits",
        "continuum or infinite-volume limits",
        "open-system evolution",
        "finite-time operational protocols",
    ):
        assert route in boundary
    assert "only a no-go for direct ordinary convergence" in boundary
    assert "Full weak-operator convergence" in boundary
    assert "finite-dimensional operator topologies coincide" in boundary
    assert "constructs no wave operator, S-matrix" in boundary
    assert "Scientific owner #743 records" in boundary
    assert "no frozen prediction" in boundary


def test_modal_maxwell_row_records_missing_physical_em_attachments(result):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    row = rows["modal_maxwell_factorization_boundary"]
    assert row["match"] == (
        "exact bounded modal factorization; local source-produced "
        "physical Maxwell theory not constructed"
    )
    assert row["lean_declarations"]["ModalMaxwellFactorizationBoundary"] == [
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
    boundary = row["hypothesis_boundary"]
    for excluded in (
        "local position-space operator",
        "opposite-momentum reality pairing",
        "U(1) potential or gauge quotient",
        "Maxwell action",
        "finite Gauss receipts",
        "Lorentz covariance",
        "laboratory readout",
    ):
        assert excluded in boundary
    assert "PR-20, PR-21, and PR-22 remain declared inputs" in boundary
    assert "PR-53 and PR-54 name missing attachments" in boundary
    assert "Scientific owner #733" in boundary
    assert "emits no frozen prediction" in boundary


def test_adaptive_scheduler_helper_keeps_its_inputs_and_v3_gate_explicit(result):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    row = rows["conditional_adaptive_scheduler_locality_helper"]
    assert row["match"] == (
        "exact conditional helper; source scheduler and physical "
        "channel attachment not constructed"
    )
    assert row["lean_declarations"]["AdaptiveScheduler"] == [
        "adaptiveRun_agree_on",
        "adaptive_no_influence",
        "consultation_region_not_droppable",
        "ball_image",
        "run_natural",
        "readback_cone_bound",
    ]
    boundary = row["hypothesis_boundary"]
    assert "sigma, R, ConsultsOnly" in boundary
    assert "Scientific owner #728 records the missing source scheduler/channel" in boundary
    assert "no prediction-ladder entry" in boundary


def test_e1_conditional_packet_keeps_source_regional_correlation_gate(result):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    row = rows["finite_causal_observer_net_interface"]
    assert row["match"] == (
        "substantial conditional finite interface, operator-generation, "
        "CP-diamond, and support-disjoint counted-correlation packet; "
        "source justification of the factor reading and a nonconstant "
        "realization are not constructed"
    )
    assert row["lean_declarations"]["TwoSlotCPNetWitness"] == [
        "slotExpectations_not_jointly_injective",
        "checkpoint_pinch_fixes_right",
        "twoSlot_left_no_scalar_hom",
    ]
    assert "two-observers-as-factors reading and region map remain declared" in row[
        "hypothesis_boundary"
    ]
    assert "coverage does not imply unique reconstruction" in row[
        "hypothesis_boundary"
    ]
    assert "constant-tower transport is on the separate 86/88 carrier" in row[
        "hypothesis_boundary"
    ]
    assert "Scientific owner #728 records" in row["hypothesis_boundary"]


def test_b7_reference_and_stationary_controls_do_not_overclose(result):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    row = rows["finite_history_variational_helpers_and_bridge_obstruction"]
    assert row["lean_declarations"]["SourceReferenceSelection"] == [
        "heatBath_counting_trivial_eq_uniform",
        "reference_realized_under_counting_trivial_inputs",
        "nontrivial_datum_not_invariant",
        "noncounting_reference_not_invariant",
        "heatBath_scaledCounting_trivial_eq_uniform",
        "uniform_transition_does_not_determine_path_reference",
        "committed_tilts_have_distinct_mean_actions",
    ]
    assert row["lean_declarations"]["StationarySaddleCoverage"] == [
        "stationaryMaximumHistory_stationary",
        "stationaryMaximumHistory_not_minimal",
        "gibbs_prefers_nonstationary",
    ]
    assert row["lean_declarations"]["SourceHistoryPacket"] == [
        "sourceMatchQuad_strictMonoOn_pos",
        "sourcePositiveMeanMatch_unique",
        "sourceMatchingPositiveParameter_existsUnique",
    ]
    assert row["lean_declarations"]["LogTransitionAction"] == [
        "bare_log_action_multiplier_unique_of_nonconstant",
    ]
    boundary = row["hypothesis_boundary"]
    assert "supplied counting reference" in boundary
    assert "distinct initial laws" in boundary
    assert "exists uniquely" in boundary
    assert "modes/minimizers only" in boundary
    assert "saddles, complex or signed stationary phase" in boundary
    assert "Closed #731 records only the attained declared-bundle composition" in boundary
    assert "Scientific owner #739 records source selection" in boundary
    assert "superseded historical gate" not in boundary


def test_hypercharge_spectrum_matches_receipt(result):
    row = next(
        r
        for r in result["sections"]["forced_structure"]
        if r["id"] == "hypercharge_spectrum"
    )
    receipt = json.loads(ledger.PARENTS["matter_receipt"].read_text(encoding="utf-8"))
    menu = json.loads(ledger.PARENTS["matter_menu"].read_text(encoding="utf-8"))
    assert row["realized_spectrum"] == receipt["realized_package"]["charge_spectrum"]
    assert row["match"] == "exact"
    assert row["subset_count"] == menu["subset_classification"]["subsets_enumerated"]
    assert row["survivor_count"] == menu["subset_classification"]["survivor_count"]
    assert row["survivor_dimension"] == receipt["realized_package"]["dimension"]
    assert row["lean_declarations"]["ExteriorSelection"] == menu[
        "subset_classification"
    ]["lean_cross_reference"]["theorems"]
    assert any(
        path.endswith("/ExteriorSelection.lean")
        for path in row["lean_receipts"]
    )


def test_classical_carriers_and_xy_boundary_are_visible(result):
    rows = {r["id"]: r for r in result["sections"]["forced_structure"]}
    for row_id in (
        "maxwell_classical_massless_kernel",
        "yang_mills_classical_massless_kernel",
        "einstein_classical_massless_kernel",
    ):
        assert rows[row_id]["match"] == "conditional structural"
        assert rows[row_id]["artifact_ref"] == (
            "code/particles/runs/status/carrier_mode_acceptance.json"
        )
    xy = rows["simple_gut_xy_channel_absent"]
    assert xy["match"] == "conditional algebraic channel exclusion"
    assert "(3,2,-5/6) (+) (bar3,2,+5/6)" in xy["statement"]
    assert xy["derivation_kind"] == "direct_executable_algebraic_corollary"
    assert xy["artifact_ref"] == (
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
    )
    assert xy["operator_census_ref"] == (
        "code/a5_closure/receipts/baryon_dimension_six_census.receipt.json"
    )
    assert xy["adjoint_branching"]["mixed_xy_bifundamental_dimension"] == 0
    assert "lean_receipts" not in xy
    assert "general proton stability does not follow" in xy[
        "hypothesis_boundary"
    ].lower()
    assert "physical current source gate is false" in xy[
        "hypothesis_boundary"
    ].lower()
    assert xy["observed_counterpart"] == (
        "the Standard Model product adjoint contains no connected "
        "simple-GUT X/Y generator"
    )
    assert "no observed proton decay" not in xy["observed_counterpart"].lower()


def test_lie_type_and_conditional_z6_descent_are_current(result):
    rows = {r["id"]: r for r in result["sections"]["forced_structure"]}
    gauge = rows["gauge_lie_algebra"]
    assert gauge["artifact_ref"] == (
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
    )
    assert gauge["match"] == (
        "axiom-forced abstract Lie type; conditional matrix witness"
    )
    assert "Complete compact port response" in gauge["statement"]
    assert "ordered source histories" in gauge["hypothesis_boundary"]
    assert "A2HolonomyBridge" in gauge["lean_declarations"]
    assert "issues 567 and 599" not in gauge["hypothesis_boundary"]

    global_form = rows["global_form_z6"]
    assert global_form["match"] == (
        "exact conditional kernel and maximal faithful image"
    )
    assert (
        "code/a5_closure/receipts/axis_center_descent_reference.receipt.json"
        in global_form["artifact_refs"]
    )
    assert "Z6Descent" in global_form["lean_declarations"]
    assert "sixAxisToKernel_range" in global_form["lean_declarations"]["Z6Descent"]
    assert "character completeness" in global_form["hypothesis_boundary"]
    assert "laboratory attachment" in global_form["hypothesis_boundary"]


def test_rank_three_completion_and_carrier_class_are_visible(result):
    rows = {r["id"]: r for r in result["sections"]["forced_structure"]}

    completion = rows["intrinsic_rank_three_response_completion"]
    assert "rank-three Gram quotient" in completion["statement"]
    assert "same abstract three-dimensional Euclidean completion" in completion[
        "statement"
    ]
    assert completion["match"].startswith("exact intrinsic metric completion")
    assert "frameQuotient_finrank" in completion["lean_declarations"][
        "PrimitivePortFrameQuotient"
    ]
    assert "d6Position_denseRange" in completion["lean_declarations"][
        "SeamCurrentCarrierQuotient"
    ]

    carrier = rows["carrier_class_dispersion_band"]
    assert "B0/C4^2 at least 10/21" in carrier["statement"]
    assert "5 D6 B0 = 12 B6 D0" in carrier["statement"]
    assert "D6/D0 = (12/5)(B6/B0)" in carrier["statement"]
    assert "zero-anisotropy mixture" in carrier["statement"]
    assert "multi-radius members retain radial-moment dependence" in carrier[
        "statement"
    ]
    assert carrier["artifact_ref"] == (
        "code/a5_fingerprint/runtime/carrier_class_dispersion_receipt.json"
    )
    assert "cross_order_lock" in carrier["lean_declarations"][
        "A5CarrierClassBand"
    ]
    assert "cross_order_polynomial" in carrier["lean_declarations"][
        "A5CarrierClassBand"
    ]
    assert "multi_radius_negative_control" in carrier["lean_declarations"][
        "A5CarrierClassBand"
    ]


def test_layered_discrete_gauss_row_is_premise_complete_and_bounded(result):
    rows = {r["id"]: r for r in result["sections"]["forced_structure"]}
    row = rows["layered_discrete_gauss_boundary"]
    assert "Q/(c n^2)" in row["statement"]
    assert "three-vertex chain inhabits" in row["statement"]
    assert "forces c=0" in row["statement"]
    assert "PR-29" in row["match"] and "PR-31" in row["match"]
    assert "physical radius" in row["hypothesis_boundary"]
    assert "perEdgeFlux_eq" in row["lean_declarations"]["LayeredDiscreteGauss"]


def test_bounded_time_row_does_not_promote_order_to_time(result):
    row = next(
        item
        for item in result["sections"]["forced_structure"]
        if item["id"] == "bounded_observer_time_calibration"
    )
    assert "supplied strictly increasing natural-number rank" in row["statement"]
    assert "order-compatible scalar readout" in row["statement"]
    assert "affine consistency is equivalent to a cross-product equation" in row["statement"]
    assert "event and both readings differ from the anchors" in row["statement"]
    assert "ordinal observer time" not in row["observed_counterpart"]
    assert row["match"].startswith(
        "exact bounded conditional algebra with finite controls"
    )


def test_source_derived_finite_one_three_row_binds_exact_stack_without_continuum_promotion(
    result,
):
    rows = {r["id"]: r for r in result["sections"]["forced_structure"]}
    row = rows["source_derived_finite_one_three_causal_carrier"]
    assert row["match"] == (
        "exact source-derived finite order and 1+3 carrier theorem "
        "stack; physical faithful placement and continuum manifold "
        "not constructed"
    )
    assert row["lean_declarations"] == {
        "CausalInterval": ["finiteCausalSetAxioms"],
        "SemanticEventProvenance": ["sourceHeight_eq"],
        "SourceDerivedSpacetimeCarrier": [
            "sourceSpacetimeCarrier_finrank",
            "sourceCarrier_one_three_signature",
            "sourceUnitNullVector_futureCausal",
            "generatedBefore_sourceCausalLE",
            "generatedBeforeEq_iff_sourceCausalLE",
            "exactDiamondConeOrder",
        ],
        "SourceOrderFrameCompatibilityPacket": [
            "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement",
            "SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame",
            "sourceOrderFramePacket_consequences",
        ],
    }
    assert "exact longest-parent-path recursion" in row["statement"]
    assert "exact four-dimensional carrier with signature (+---)" in row[
        "statement"
    ]
    assert "an independent real axis and the proved rank-three source" in row[
        "statement"
    ]
    assert "Independently of that order" in row["statement"]
    assert "use source height only in the event placement" in row["statement"]
    assert "Adjoining that ordinal axis" not in row["statement"]
    assert "two-way order--cone equivalence" in row["statement"]
    assert "four-event Boolean diamond" in row["statement"]
    boundary = row["hypothesis_boundary"]
    for excluded in (
        "not a physical clock",
        "spatial readback valued in the rank-three quotient",
        "additional converse field of a faithful placement",
        "not empirical evidence of manifoldlikeness",
        "calibrated density or count--volume law",
        "independent dimension estimator",
        "continuum convergence",
        "smooth metric",
        "Einstein equation",
    ):
        assert excluded in boundary


def test_thermodynamic_receipt_owners_are_separate(result):
    row = next(
        item
        for item in result["sections"]["forced_structure"]
        if item["id"] == "thermodynamic_four_law_package"
    )
    boundary = row["hypothesis_boundary"]
    assert "is only a bounded negative result" in boundary
    assert "two direct mechanisms" in boundary
    assert "mixing-mode-retaining linear intertwiner" in boundary
    assert "deterministic empirical pushforward" in boundary
    assert "does not exclude stochastic, nonlinear, reverse-direction" in boundary
    assert "Scientific owner #739 records the missing replacement common reference" in boundary
    assert "empirical energy-clock calibration is PR-15" in boundary
    assert "Closed #732 records only the attained conditional composition milestone" in boundary
    assert "superseded historical owners" not in boundary
    assert "five receipts stay open under issue #688" not in boundary


def test_born_frame_row_keeps_post_hoc_orientation_out_of_source_selection(result):
    row = next(
        item
        for item in result["sections"]["forced_structure"]
        if item["id"] == "finite_born_frame_rank_gap"
    )
    assert "source-attached real S3 algebraic contexts" in row["statement"]
    assert "post-hoc raw-count product-gap diagnostic" in row["statement"]
    assert "emits no source selection or validation" in row["statement"]
    assert "bornWeight_re_matrixConj" in row["lean_declarations"]["ConjugationGauge"]
    assert "designatedCycle_normalized_products" in row["lean_declarations"][
        "RepairCurrentOrientation"
    ]
    assert "oriented_born_capstone" in row["lean_declarations"][
        "SourceOrientedCompletion"
    ]
    assert "statistic and designation rule were not preregistered" in row[
        "hypothesis_boundary"
    ]
    assert "phase pairing is an arbitrary typed convention" in row[
        "hypothesis_boundary"
    ]
    assert (
        "code/thermodynamics/repair_current_orientation/verify_repair_current_orientation.py"
        in row["artifact_refs"]
    )


def test_alpha_row_values_match_endpoint(result):
    row = result["sections"]["alpha"][0]
    endpoint = json.loads(ledger.PARENTS["endpoint"].read_text(encoding="utf-8"))
    assert row["value_central"] == pytest.approx(
        float(endpoint["endpoint"]["alpha_inv_central"])
    )
    assert row["measured"] == pytest.approx(
        float(endpoint["compare_only"]["codata_alpha_inv"])
    )
    verdict = json.loads(
        ledger.PARENTS["alpha_hvp_verdict"].read_text(encoding="utf-8")
    )
    assert row["reference_deficit_inside_recorded_accounting_interval"] is True
    assert row["audit_verdict"] == verdict["verdict"]
    assert row["cross_class_agreement"]["independently_evaluated_class_count"] == 0
    assert "does not identify the physical source" in row["reading"]
    assert row["scientific_owner_issues"] == [736]
    assert "historical_context_issues" not in row


def test_lepton_rows_match_parents_and_contain_witness(result):
    rows = {r["id"]: r for r in result["sections"]["charged_leptons"]}
    coherent = json.loads(ledger.PARENTS["kappa_coherent"].read_text(encoding="utf-8"))
    row = rows["charged_leptons_kappa_coherent"]
    assert row["witness_inside_all_intervals"] is True
    assert row["intervals_gev"] == [
        r["mass_interval"] for r in coherent["conditional_mass_rows"]
    ]
    assert row["width_reduction_factor"] == coherent["kappa_interval"][
        "width_reduction_factor"
    ]
    assert rows["charged_leptons_kappa_rectangle"]["witness_inside_all_intervals"] is True
    for row_id in (
        "charged_leptons_closure_target",
        "charged_leptons_kappa_rectangle",
        "charged_leptons_kappa_coherent",
    ):
        assert rows[row_id]["scientific_owner_issues"] == [736]
        assert "historical_context_issues" not in rows[row_id]


def test_ew_rows_preserve_comparison_status(result):
    rows = {r["id"]: r for r in result["sections"]["electroweak"]}
    assert rows["ew_mH_gev"]["physical_comparison_status"] == "COMPARE_ONLY"
    assert rows["ew_MW_chart_gev"]["physical_comparison_status"] == "NOT_EVALUABLE"
    assert "measured" not in rows["ew_MW_chart_gev"]
    parent = json.loads(ledger.PARENTS["conditional_ew"].read_text(encoding="utf-8"))
    assert rows["ew_mH_gev"]["delta_over_sigma"] == parent[
        "comparison_compare_only"
    ]["mH_gev"]["delta_over_sigma"]


def test_quark_section_is_obstruction_plus_conditional_texture(result):
    rows = {r["id"]: r for r in result["sections"]["quarks"]}
    obstruction = rows["quark_absolute_masses_obstruction"]
    assert obstruction["fork"] == "ii_fiber_survives"
    assert obstruction["fiber_cut_detected"] is False
    texture = rows["quark_down_type_clebsch_route_rejected"]
    assert texture["tier"] == "T2_conditional_rejected_candidate"
    assert texture["promotion_allowed"] is False
    assert texture["retrospective_flag_rejection"][
        "all_six_permutations_rejected"
    ] is True
    assert texture["permutation_scan"]["retrospective_metric"][
        "target_informed"
    ] is True
    assert "cabibbo_gst_sqrt_md_over_ms" in texture["values"]
    assert obstruction["scientific_owner_issues"] == [736]
    assert "historical_context_issues" not in obstruction


def test_hadron_row_carries_pinned_payload(result):
    rows = {r["id"]: r for r in result["sections"]["hadrons"]}
    engine = rows["hadronic_correction_engine"]
    payload = json.loads(ledger.PARENTS["hadron_payload"].read_text(encoding="utf-8"))
    assert engine["delta_alpha_had_5_MZ"] == payload["integral"]["value"]
    assert engine["uncertainty_total"] == payload["integral"]["uncertainty"]


def test_fail_closed_on_missing_parent(tmp_path, monkeypatch):
    monkeypatch.setitem(
        ledger.PARENTS, "endpoint", tmp_path / "absent_endpoint.json"
    )
    with pytest.raises(SystemExit, match="parent missing"):
        ledger.build(tmp_path / "out.json", None)


def test_markdown_rendered(tmp_path):
    md = tmp_path / "ledger.md"
    ledger.build(tmp_path / "out.json", md)
    text = md.read_text(encoding="utf-8")
    assert "# Postdiction Ledger" in text
    assert "## Forced structure" in text
    assert "## Quantum carrier gate" in text
    assert "NOT_EVALUABLE" in text
    assert "Recorded retrospective same-scheme accounting interval" in text
    assert "Certified same-scheme anchor gap" not in text
    assert "`code/particles/scripts/build_postdiction_ledger.py`" in text
    assert "`code/particles/runs/status/postdiction_ledger.json`" in text
    assert "Each row is checked in Lean, by a structured executable artifact" in text
    assert "Every step is machine checked in the Lean workspace" not in text
    assert "Scientific owner: #736" in text
    assert "Historical context:" not in text


def test_principal_results_prioritize_strong_structural_rows(result):
    principal = result["principal_results"]
    assert principal[0]["id"] == "intrinsic_rank_three_response_completion"
    assert {p["id"] for p in principal} == {
        "intrinsic_rank_three_response_completion",
        "forced_gauge_structure",
        "carrier_class_dispersion_surface",
        "koide_conditional_tau_window",
        "lepton_closure_target",
    }
    assert "target-informed conditional postdiction" in principal[3]["statement"]
    wp = next(
        r
        for r in result["sections"]["charged_leptons"]
        if r["id"] == "charged_leptons_closure_target"
    )["witness_point"]
    assert f"{wp['required_anchor_gap_at_witness_inv_alpha']:.4f}" in principal[4]["statement"]


def test_closure_target_row_reads_lane_artifacts(result):
    row = next(
        r
        for r in result["sections"]["charged_leptons"]
        if r["id"] == "charged_leptons_closure_target"
    )
    parent = json.loads(ledger.PARENTS["kappa_rectangle"].read_text(encoding="utf-8"))
    assert row["witness_point"] == parent["compare_only"]["witness_point"]
    assert "545" in row["width_floor"]


def test_build_is_deterministic_and_has_no_wall_clock_field(tmp_path):
    first = ledger.build(
        tmp_path / "ignored-a.json",
        tmp_path / "ignored-a.md",
        write=False,
    )
    second = ledger.build(
        tmp_path / "ignored-b.json",
        tmp_path / "ignored-b.md",
        write=False,
    )
    assert first == second
    assert ledger._render_md(first) == ledger._render_md(second)
    assert "generated_utc" not in first
    assert first["schema_version"] == 3


def test_fail_closed_on_missing_lean_declaration(tmp_path, monkeypatch):
    menu = json.loads(ledger.PARENTS["matter_menu"].read_text(encoding="utf-8"))
    menu["subset_classification"]["lean_cross_reference"]["theorems"].append(
        "fabricated_exterior_theorem"
    )
    path = tmp_path / "matter-menu.json"
    path.write_text(json.dumps(menu), encoding="utf-8")
    monkeypatch.setitem(ledger.PARENTS, "matter_menu", path)
    monkeypatch.setattr(ledger, "_rel", lambda key: f"test/{key}.json")
    with pytest.raises(SystemExit, match="Lean declaration missing"):
        ledger.build(tmp_path / "out.json", None)


def test_fail_closed_on_inconsistent_port_current(tmp_path, monkeypatch):
    receipt = json.loads(
        ledger.PARENTS["port_current"].read_text(encoding="utf-8")
    )
    receipt["closure"]["derived_block_dimensions"]["even_block_su3"] = 9
    path = tmp_path / "port-current.json"
    path.write_text(json.dumps(receipt), encoding="utf-8")
    monkeypatch.setitem(ledger.PARENTS, "port_current", path)
    monkeypatch.setattr(ledger, "_rel", lambda key: f"test/{key}.json")
    with pytest.raises(SystemExit, match="product algebra"):
        ledger.build(tmp_path / "out.json", None)


def test_whitney_rows_remain_mathematical_not_observed_or_continuum(whitney_rows):
    assert set(whitney_rows) == {
        "whitney_maxwell_dynamics",
        "whitney_radiative_quantum",
        "whitney_charged_matter",
    }
    for row in whitney_rows.values():
        assert row["physical_comparison_status"] == "NOT_EVALUABLE"
        assert row["observed_postdiction"] is False
        assert row["continuum_convergence_established"] is False
        assert "no observed data" in row["observed_counterpart"]
        assert "no observed postdiction" in row["match"]
        assert not {"measured", "delta_over_sigma"} & row.keys()


def test_whitney_rows_are_in_the_generated_forced_structure(result, whitney_rows):
    rows = {row["id"]: row for row in result["sections"]["forced_structure"]}
    for name, expected in whitney_rows.items():
        assert rows[name] == expected


def test_whitney_dynamics_keeps_exact_numeric_and_instrumented_scopes(whitney_rows):
    row = whitney_rows["whitney_maxwell_dynamics"]
    verified = row["independent_verifier_result"]
    assert {
        key: verified[key]
        for key in (
            "gauge_histories",
            "events_per_history",
            "instrumented_slices_per_history",
            "probe_cycles_per_history",
            "field_variations_per_history",
            "continuation_slabs",
            "exact_global_stiffness_bound",
        )
    } == {
        "gauge_histories": 2,
        "events_per_history": 805,
        "instrumented_slices_per_history": 3,
        "probe_cycles_per_history": 252,
        "field_variations_per_history": 68,
        "continuation_slabs": 64,
        "exact_global_stiffness_bound": 24,
    }
    # Residuals are checked by the verifier but omitted from the generated row
    # so platform-dependent rounding does not change the ledger's bytes.
    assert not {"gauss_max_abs", "ampere_max_abs"} & verified.keys()
    assert row["certificate_scope"] == (
        "NUMERIC_STATIONARY_FINITE_VOLUME_MAXWELL__"
        "EXACT_CLASSICAL_READBACK__EXACT_LOCAL_STABILITY_BOUND"
    )
    policy = row["numeric_policy"]
    assert policy["atol"] == policy["rtol"] == 1e-9
    assert "float64" in policy["evolution"]
    assert "no interval or exact-trajectory certificate" in policy["evolution"]
    assert "exact pair-average restoration" in policy["registers"]
    assert "exact Q(sqrt(5)) LDL certificate" in policy["stability"]
    assert row["arbitrary_test_continuum_stationarity"] is False
    assert row["finite_element_stationarity"] == (
        "numerically verified for all 68 free coefficients per gauge history"
    )
    assert "separately identified 64-slab numerical continuation" in row["statement"]
    boundary = row["hypothesis_boundary"]
    for qualification in (
        "changed initial electric field",
        "prescribed impulse source covectors",
        "Only three slices per gauge are instrumented",
        "not arbitrary-test continuum stationarity",
    ):
        assert qualification in boundary
    assert {
        "recurrence_next_unique",
        "gauss_residual_preserved",
        "pairEnergy_work",
        "modal_solution_bound",
        "two_slab_endpoint_resonance",
    } <= set(row["lean_declarations"]["WhitneyMaxwellDynamics"])
    assert {
        "code/electromagnetism/runtime/whitney_maxwell_dynamics_receipt.json",
        "code/electromagnetism/verify_whitney_maxwell_dynamics.py",
    } <= set(row["artifact_refs"])


def test_whitney_quantum_keeps_paper_hilbert_proof_and_free_sector(whitney_rows):
    row = whitney_rows["whitney_radiative_quantum"]
    assert set(row["lean_declarations"]) == {"WhitneyQuantumBridge"}
    assert set(row["lean_declarations"]["WhitneyQuantumBridge"]) == {
        "reconstruct_surjective",
        "same_action_normal_modes",
        "annihilation_creation",
        "canonical_position_momentum",
        "same_modes_quantized_energy",
        "quantumHamiltonian_monomial",
        "quantumHamiltonian_position",
        "quantumHamiltonian_momentum",
    }
    assert "exact polynomial CCR" in row["statement"]
    assert "paper constructs the factorial-weight Hilbert completion l2(N^30)" in row[
        "statement"
    ]
    assert "analytic Hilbert/domain and temporal-limit proof in the paper" in row["match"]
    boundary = row["hypothesis_boundary"]
    for qualification in (
        "Canonical quantization and hbar are imported",
        "Lean assumes a complete positive normal frame",
        "Hilbert completion, operator closures and temporal convergence are paper proofs",
        "zero-charge source-free radiative sector",
        "not quantization of the prescribed-charge episode or interacting matter action",
        "sqrt(lambda) and finite-step theta/h are distinct",
        "no uniform operator-norm quantum error",
    ):
        assert qualification in boundary
    assert "independent_verifier_result" not in row


def test_whitney_matter_keeps_gauge_algebra_distinct_from_analytic_evolution(whitney_rows):
    row = whitney_rows["whitney_charged_matter"]
    assert set(row["lean_declarations"]) == {"WhitneyChargedMatter"}
    assert set(row["lean_declarations"]["WhitneyChargedMatter"]) == {
        "pathPhase_gauge",
        "interpolate_gauge",
        "interpolate_gauge_norm",
        "interpolate_vertex",
        "interpolate_face_trace",
    }
    assert "exact finite interpolation gauge algebra" in row["match"]
    assert "analytic coupled-action and global-existence proof" in row["match"]
    assert "paper proves its Noether/Gauss identity" in row["statement"]
    assert "global finite-dimensional classical evolution" in row["statement"]
    boundary = row["hypothesis_boundary"]
    for qualification in (
        "nonnegative quartic coefficient",
        "Real unwrapped edge integrals are needed",
        "Differentiate the gauge-dependent scalar basis in all field variations",
        "All thirteen scalar variations require total neutrality",
        "No executed charged-matter episode",
        "source-selected matter content",
        "interacting quantum completion",
    ):
        assert qualification in boundary
    assert "independent_verifier_result" not in row
    assert "code/electromagnetism/test_whitney_charged_matter.py" in row["artifact_refs"]


@pytest.mark.parametrize("row_id", [
    "whitney_maxwell_dynamics",
    "whitney_radiative_quantum",
    "whitney_charged_matter",
])
def test_whitney_named_theorem_cannot_be_replaced_by_empty_provider(
    whitney_rows, row_id, tmp_path, monkeypatch
):
    row = whitney_rows[row_id]
    module, declarations = next(iter(row["lean_declarations"].items()))
    empty = tmp_path / f"{module}.lean"
    empty.write_text("-- No declared theorem is present.\n", encoding="utf-8")
    monkeypatch.setitem(ledger.LEAN_RECEIPTS, module, empty)
    with pytest.raises(SystemExit, match=f"Lean declaration missing: {module}"):
        ledger._lean_receipt(module, declarations={module: tuple(declarations)})


def test_whitney_builder_rejects_altered_dynamics_provider(monkeypatch):
    original = Path.read_bytes
    provider = ledger.REPO / "Lean/Screen/WhitneyMaxwellDynamics.lean"
    seen = []

    def altered(path):
        data = original(path)
        if path == provider:
            seen.append(path)
            return data + b"\n-- independent provider mutation\n"
        return data

    monkeypatch.setattr(Path, "read_bytes", altered)
    with pytest.raises(
        ValueError, match="provider pin Lean/Screen/WhitneyMaxwellDynamics.lean"
    ):
        ledger._whitney_dynamics_rows()
    assert seen == [provider]
