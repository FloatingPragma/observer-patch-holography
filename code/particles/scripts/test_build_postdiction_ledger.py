"""Tests for the postdiction ledger aggregator."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import pytest

import build_postdiction_ledger as ledger


def test_finite_spectrum_comparison_projection_and_rendering():
    rows = ledger._finite_spectrum_comparison_rows()
    parent = json.loads((ledger.REPO / "model/particle_spectrum_audit/results.json").read_text())
    assert rows["charged_leptons"]["ratios"] == parent["lepton_ratios"]
    assert rows["quarks"]["inverse_coordinates"] == parent["quark_exponent_diagnostics"]
    assert rows["neutrinos"]["comparison"] == parent["neutral"]
    assert rows["neutrinos"]["prospective_FZ03_verdict"] is None
    retained = json.loads(ledger.DEFAULT_OUT.read_text())
    for section, row in rows.items():
        retained["sections"][section] = [row] + [old for old in retained["sections"][section]
                                                if old["id"] != row["id"]]
    text = ledger._render_md(retained)
    assert "105.6509313" in text and "+0.07249 ppm" in text
    assert "u=6.63023" in text and "0.004988" in text
    assert "MCPR conditional triple" in text  # The separate branch is retained.
    assert "no prospective FZ-03 verdict" in rows["neutrinos"]["hypothesis_boundary"]
    changed = deepcopy(rows["charged_leptons"])
    changed["particles"][0]["value_GeV"] = 0.001
    assert "masses are 1," in "\n".join(ledger._render_finite_spectrum_comparison(changed))


@pytest.mark.parametrize("mutation", ["custody", "promotion"])
def test_finite_spectrum_comparison_rejects_invalid_parent(monkeypatch, mutation):
    original = ledger._public_receipt

    def changed(relative, schema):
        packet, checksum = original(relative, schema)
        if schema == "oph.particle_spectrum_audit.result.v1":
            packet = deepcopy(packet)
            if mutation == "custody":
                packet["source_sha256"]["derive.py"] = "0"*64
            else:
                packet["scope"]["physical_spectrum_identity_established"] = True
        return packet, checksum

    monkeypatch.setattr(ledger, "_public_receipt", changed)
    with pytest.raises(SystemExit, match="finite spectrum comparison"):
        ledger._finite_spectrum_comparison_rows()


def test_retained_public_record_projection_is_current():
    """Check live receipt projection without rebuilding unrelated proof rows."""
    retained = json.loads(ledger.DEFAULT_OUT.read_text(encoding="utf-8"))
    assert retained["retrospective_public_record_comparisons"] == ledger._retrospective_public_record_references()


def test_spectral_cross_scale_projection_keeps_actual_join_boundary():
    summary = ledger._spectral_cross_scale_reference()
    assert summary["coupling"]["within_one_reported_standard_uncertainty"] is True
    assert summary["tail_weight_interval"] == [3114, 3114]
    assert summary["joint_status"]["declared_chiral_prescriptions_passing_both_lepton_ratios"] == []
    assert summary["joint_status"]["joint_likelihood_or_significance"] is None


@pytest.mark.parametrize("mutation", ["promotion", "custody"])
def test_spectral_cross_scale_projection_rejects_false_join(monkeypatch, mutation):
    original = ledger._public_receipt
    def changed(relative, schema):
        packet, checksum = original(relative, schema)
        if schema == "oph.empirical_diagnostics.fixed_spectral_cross_scale.v1":
            packet = deepcopy(packet)
            if mutation == "promotion":
                packet["joint_status"]["all_public_datasets_reproduced"] = True
            else:
                pins = packet["custody"]["input_sha256"]
                pins[next(iter(pins))] = "0"*64
        return packet, checksum
    monkeypatch.setattr(ledger, "_public_receipt", changed)
    with pytest.raises(SystemExit, match="[Ss]pectral cross-scale"):
        ledger._spectral_cross_scale_reference()


def test_complete_occupation_projection_keeps_same_data_and_scope():
    summary = ledger._occupation_control_reference()
    assert summary['observations'] == 125
    assert summary['same_observations_as_anharmonic_transfer'] is True
    assert summary['new_independent_data'] is False
    assert summary['full_uncertainty_agreement_established'] is False
    assert summary['statistical_rejection_claimed'] is False
    assert summary['certified'] == sum(row['certified'] for row in summary['per_case'])
    assert Fraction(summary['minimum_common_remaining_allowance']) >= Fraction(summary['minimum_common_remaining_allowance_floor_5'])


@pytest.mark.parametrize('mutation', ['fit', 'agreement', 'significance', 'duplicate', 'summary', 'custody'])
def test_complete_occupation_projection_rejects_promotion_or_incomplete_coverage(monkeypatch, mutation):
    original = ledger._public_receipt
    def changed(relative, schema):
        packet, checksum = original(relative, schema)
        if schema == 'oph.occupation_control.transfer_results.v1':
            packet = deepcopy(packet)
            if mutation == 'fit':
                packet['prediction']['settings']['new_parameter_fit'] = True
            elif mutation == 'agreement':
                packet['comparison']['full_uncertainty_agreement_established'] = True
            elif mutation == 'significance':
                packet['comparison']['statistical_rejection_claimed'] = True
            elif mutation == 'duplicate':
                rows = packet['comparison']['runs'][0]['held_out_rows']
                rows[1] = deepcopy(rows[0])
            elif mutation == 'summary':
                packet['comparison']['per_case'][0]['minimum_common_remaining_allowance'] = '0'
            else:
                pins = packet['prediction']['source_bindings']
                pins[next(iter(pins))] = '0'*64
        return packet, checksum
    monkeypatch.setattr(ledger, '_public_receipt', changed)
    with pytest.raises(SystemExit, match='occupation comparison'):
        ledger._occupation_control_reference()


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ledger")
    result = ledger.build(tmp / "postdiction_ledger.json", tmp / "POSTDICTION_LEDGER.md")
    assert (tmp / "POSTDICTION_LEDGER.md").read_text(encoding="utf-8") == ledger._render_md(result)
    return result


@pytest.fixture(scope="module")
def whitney_rows():
    rows = ledger._whitney_dynamics_rows()
    assert len(rows) == 3
    return {row["id"]: row for row in rows}


@pytest.fixture(scope="module")
def coupled_whitney_rows(result):
    # New initial-state evidence strengthens its existing quantum row.
    assert len(result["sections"]["forced_structure"]) == 59
    names = {"whitney_spatial_consistency", "whitney_interacting_quantum", "whitney_charged_execution"}
    rows = {row["id"]: row for row in result["sections"]["forced_structure"] if row["id"] in names}
    assert set(rows) == names
    return rows


def test_guards_are_compare_only(result):
    guards = result["guards"]
    assert guards["compare_only"] is True
    assert guards["public_promotion_allowed"] is False
    assert guards["changes_any_solve_path"] is False
    assert guards["hand_typed_measured_values"] is False
    assert "separately supplied apparatus models and calibrations" in result["aggregation_policy"]
    assert "do not establish joint agreement for one fixed physical construction" in result["aggregation_policy"]


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


def test_retrospective_tilt_keeps_exact_observed_tail_and_no_promotion():
    comparison = ledger._retrospective_tilt_reference()
    assert comparison["pooled"]["positions"]["edge_center"][
        "posterior_mass_le_target_exact"] == "3114/508985"
    assert comparison["analysis"]["FZ15_prospective_verdict"] is None
    assert comparison["analysis"]["new_independent_data"] is False
    assert comparison["prediction"]["native_observer_model_derivation_established"] is False


@pytest.mark.parametrize("section,field,value", [
    ("analysis", "new_independent_data", True),
    ("analysis", "FZ15_prospective_verdict", "COMPATIBLE"),
    ("analysis", "tail_mass_is_p_value", True),
    ("prediction", "native_observer_model_derivation_established", True),
])
def test_retrospective_tilt_rejects_promotion(tmp_path, section, field, value):
    source = ledger.CODE / "cosmology/tilt_posterior/results.json"
    packet = json.loads(source.read_text())
    packet[section][field] = value
    path = tmp_path / "changed-tilt.json"
    path.write_text(json.dumps(packet))
    with pytest.raises(SystemExit, match="cannot promote"):
        ledger._retrospective_tilt_reference(path)


def test_retrospective_tilt_rejects_broken_exact_tail(tmp_path):
    source = ledger.CODE / "cosmology/tilt_posterior/results.json"
    packet = json.loads(source.read_text())
    packet["pooled"]["positions"]["edge_center"]["weight_strictly_below"] = "3115"
    path = tmp_path / "changed-tail.json"
    path.write_text(json.dumps(packet))
    with pytest.raises(SystemExit, match="tail accounting mismatch"):
        ledger._retrospective_tilt_reference(path)


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
    assert "Historical lane #743 (closed; live owner #1134) records" in boundary
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
    assert "Historical lane #733 (closed; live owner #1134)" in boundary
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
    assert "Historical lane #728 (closed; live owner #1134) records the missing source scheduler/channel" in boundary
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
    assert "Historical lane #728 (closed; live owner #1134) records" in row["hypothesis_boundary"]


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
    assert "Historical lane #739 (closed; live owner #1134) records source selection" in boundary
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
    control = row["supplied_law_refinement_control"]
    assert control["analytic_controlled_causal_limit"] is True
    assert control["analytic_interval_volume_limit"] is True
    assert all(control[key] is False for key in (
        "continuum_limit_formalized_in_lean", "native_OPH_law_selected",
        "physical_clock_calibrated", "observed_postdiction",
    ))
    assert control["independent_verifier_result"]["authenticated_edges"] == 794
    assert control["independent_verifier_result"]["exact_history_width"] == 27
    assert control["lean_declarations"] == {"RefiningLatticeCausalCone": [
        "causal_cone_sandwich", "scaled_inner_cone_reachable",
        "scaled_reachable_outer_cone", "fixed_menu_missing_timelike_endpoint",
    ]}
    assert control["source"] in row["artifact_refs"]
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


@pytest.mark.parametrize("field,value", [
    ("native_physical_spacetime_selected", True),
    ("history_events", True),
    ("history_events", 80),
    ("authenticated_edges", 793),
    ("largest_inner_speed", "1"),
    ("continuum_claim", "source-selected physical spacetime"),
])
def test_refining_causal_projection_rejects_changed_scope_and_counts(monkeypatch, field, value):
    original_spec = ledger.importlib.util.spec_from_file_location
    def spec(name, path, *args, **kwargs):
        item = original_spec(name, path, *args, **kwargs)
        if name == "refining_causal_ledger_verifier":
            original_exec = item.loader.exec_module
            def execute(module):
                original_exec(module)
                original_verify = module.verify
                def changed(packet):
                    result = original_verify(packet)
                    result[field] = value
                    return result
                module.verify = changed
            item.loader.exec_module = execute
        return item
    monkeypatch.setattr(ledger.importlib.util, "spec_from_file_location", spec)
    with pytest.raises(SystemExit, match="scope/count mismatch"):
        ledger._refining_causal_control()


def test_refining_causal_projection_requires_analytic_theorem(monkeypatch):
    original = Path.read_text
    target = ledger.REPO / "paper/tex_fragments/REFINING_CAUSAL_CONE.tex"
    def erased(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        return text.replace("\\label{prop:refining-causal-cone}", "") if path == target else text
    monkeypatch.setattr(Path, "read_text", erased)
    with pytest.raises(SystemExit, match="analytic theorem missing"):
        ledger._refining_causal_control()


@pytest.mark.parametrize("name", ["causal_cone_sandwich", "scaled_inner_cone_reachable",
                                 "scaled_reachable_outer_cone", "fixed_menu_missing_timelike_endpoint"])
def test_refining_causal_projection_requires_finite_lean_bounds(monkeypatch, tmp_path, name):
    original = ledger.LEAN_RECEIPTS["RefiningLatticeCausalCone"]
    text = original.read_text(encoding="utf-8")
    target = tmp_path / "RefiningLatticeCausalCone.lean"
    target.write_text(text.replace("theorem " + name + " ", "theorem erased_" + name + " "), encoding="utf-8")
    monkeypatch.setitem(ledger.LEAN_RECEIPTS, "RefiningLatticeCausalCone", target)
    with pytest.raises(SystemExit, match="Lean declaration missing"):
        ledger._refining_causal_control()


def test_source_net_projection_separates_finite_analytic_and_physical_results(result):
    rows = result["sections"]["forced_structure"]
    assert len(rows) == 59
    row = next(r for r in rows if r["id"] == "source_derived_finite_one_three_causal_carrier")
    control = row["declared_source_net_control"]
    assert control["analytic_controlled_causal_limit"] is True
    assert control["analytic_declared_family_raw_count_limit"] is True
    assert control["analytic_contained_interval_ordering_fraction"] == "1/10"
    assert control["volume_convention"] == "dt d^3x"
    for key in (
        "primitive_source_words_executed_as_repairs",
        "continuum_or_pair_count_limit_formalized_in_lean",
        "finite_runs_certify_asymptotic_limit", "native_OPH_population_or_law_selected",
        "physical_clock_calibrated", "field_action_or_quantum_continuum_identified",
        "observed_postdiction",
    ):
        assert control[key] is False
    summary = control["independent_verifier_result"]
    assert [r["source_records"] for r in summary["levels"]] == [27, 125, 512, 2197]
    assert [r["positive_inner_cone"] for r in summary["levels"]] == [False, False, False, True]
    assert summary["levels"][-1]["inner_speed_lower"] == "75989/250000"
    receipt = ledger.CODE / "causal_refinement/source_net_causet_receipt.json"
    raw = receipt.read_bytes()
    assert control["receipt_pin"] == {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
    assert control["source"] in row["artifact_refs"]
    assert "not constructed" in row["match"]


@pytest.mark.parametrize("field,value", [
    ("accepted", 1), ("native_physical_spacetime_selected", True),
    ("scope", "PHYSICAL_SPACETIME_DERIVED"), ("packet_sha256", "0"*64),
    ("q", True), ("width", 2197.0), ("height", 99),
    ("inner_speed_lower", "1"), ("positive_inner_cone", 1),
    ("reads_per_full_trace", 0),
])
def test_source_net_projection_rejects_verifier_summary_drift(result, monkeypatch, field, value):
    row = next(r for r in result["sections"]["forced_structure"]
               if r["id"] == "source_derived_finite_one_three_causal_carrier")
    summary = deepcopy(row["declared_source_net_control"]["independent_verifier_result"])
    target = summary if field in summary else summary["levels"][-1]
    target[field] = value
    original = ledger.importlib.util.spec_from_file_location
    def spec(name, path, *args, **kwargs):
        item = original(name, path, *args, **kwargs)
        if name == "source_net_causal_ledger_verifier":
            execute = item.loader.exec_module
            def changed(module):
                execute(module)
                module.verify = lambda packet: deepcopy(summary)
            item.loader.exec_module = changed
        return item
    monkeypatch.setattr(ledger.importlib.util, "spec_from_file_location", spec)
    with pytest.raises(SystemExit, match="scope/count mismatch"):
        ledger._source_net_causal_control()


@pytest.mark.parametrize("label", [
    "prop:source-net-causal-path", "prop:source-net-alexandrov-volume",
    "cor:source-net-general-diamond", "prop:golden-source-count-limit",
    "cor:source-record-ordering-fraction",
])
def test_source_net_projection_requires_each_analytic_statement(monkeypatch, label):
    original = Path.read_text
    target = ledger.REPO / "paper/tex_fragments/SOURCE_NET_CAUSAL_LIMIT.tex"
    def erased(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        return text.replace("\\label{" + label + "}", "") if path == target else text
    monkeypatch.setattr(Path, "read_text", erased)
    with pytest.raises(SystemExit, match="source-net analytic theorem missing"):
        ledger._source_net_causal_control()


@pytest.mark.parametrize("name", [
    "covering_constructs_path", "reachable_outer", "same_layer_precedes_iff",
    "finite_antichain_card_le", "cone_speed_identity",
])
def test_source_net_projection_requires_constructed_lean_paths(monkeypatch, tmp_path, name):
    source = ledger.LEAN_RECEIPTS["SourceNetCausalCone"]
    target = tmp_path / "SourceNetCausalCone.lean"
    text = source.read_text(encoding="utf-8")
    assert "theorem " + name in text
    target.write_text(text.replace("theorem " + name, "theorem erased_" + name), encoding="utf-8")
    monkeypatch.setitem(ledger.LEAN_RECEIPTS, "SourceNetCausalCone", target)
    with pytest.raises(SystemExit, match="Lean declaration missing"):
        ledger._source_net_causal_control()


@pytest.mark.parametrize("forgery", ["source_pin", "late_value", "count_includes_reads", "scope"])
def test_source_net_projection_replays_resealed_false_receipt(tmp_path, forgery):
    source = ledger.CODE / "causal_refinement/source_net_causet_receipt.json"
    packet = json.loads(source.read_text(encoding="utf-8"))
    first = packet["levels"][0]
    if forgery == "source_pin":
        packet["source_pins"]["code/causal_refinement/source_net_causet.py"] = "0"*64
    elif forgery == "late_value":
        first["forward_execution"]["layer_value_sums"][-1] += 1
    elif forgery == "count_includes_reads":
        first["center_intervals"][-1]["inclusive_event_count"] += first["forward_execution"]["authenticated_read_count"]
    else:
        packet["scope"]["finite_runs_demonstrate_fixed_time_asymptotic_limit"] = True
    target = tmp_path / "resealed.json"
    target.write_text(json.dumps(packet, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="source-net independent replay failed"):
        ledger._source_net_causal_control(target)


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
    assert "Historical lane #739 (closed; live owner #1134) records the missing replacement common reference" in boundary
    assert "empirical energy-clock calibration is PR-15" in boundary
    assert "Closed #732 records only the attained conditional composition milestone" in boundary
    assert "superseded historical owners" not in boundary
    assert "five receipts stay open under issue #688" not in boundary


def test_protected_memory_projection_replays_exact_packet_and_bytes(tmp_path):
    source = ledger.CODE / "thermodynamics/protected_memory/runtime/protected_memory_receipt.json"
    raw = source.read_bytes()
    copied = tmp_path / source.name
    copied.write_bytes(raw)
    control = ledger._protected_memory_control(copied)
    assert control == ledger._protected_memory_control()
    assert control["receipt_pin"] == {
        "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
    }
    report = control["independent_verifier_result"]
    assert report["accepted"] is True
    assert report["gap_lower"] == "7/192"
    assert report["green_kubo_matrix"] == [
        ["11417/9375", "1832/3125"], ["1832/3125", "3816/3125"],
    ]
    assert control["fibre_centered_currents_required"] is True
    assert control["new_supplied_local_stochastic_law"] is True
    assert control["observed_postdiction"] is False


def test_protected_memory_projection_rejects_rehashed_false_memory(tmp_path):
    source = ledger.CODE / "thermodynamics/protected_memory/runtime/protected_memory_receipt.json"
    packet = json.loads(source.read_text(encoding="utf-8"))
    gk = packet["green_kubo"]
    # Forge instantaneous mixing while retaining internally matching sum and
    # remainder summaries. All provider pins remain valid; the ledger hashes
    # the new file itself, so a stale receipt hash cannot be the rejection.
    gk["correlations"][1:] = [[["0", "0"], ["0", "0"]] for _ in range(64)]
    gk["matrix"] = deepcopy(gk["equal_time"])
    gk["partial_sum"] = deepcopy(gk["equal_time"])
    gk["exact_remainder"] = [["0", "0"], ["0", "0"]]
    forged = tmp_path / "rehashed_false_memory.json"
    forged.write_text(json.dumps(packet, sort_keys=True) + "\n", encoding="utf-8")
    assert hashlib.sha256(forged.read_bytes()).digest() != hashlib.sha256(source.read_bytes()).digest()
    with pytest.raises(SystemExit, match="protected-memory independent replay failed"):
        ledger._protected_memory_control(forged)


@pytest.mark.parametrize("field,value", [
    ("accepted", 1),
    ("states", True),
    ("gap_lower", "1"),
    ("native_source_attachment", True),
    ("historical_obstruction_overturned", True),
    ("physical_clock_or_conductivity", True),
    ("new_Lean_formalization", True),
])
def test_protected_memory_projection_rejects_provider_promotion(monkeypatch, field, value):
    original_spec = ledger.importlib.util.spec_from_file_location
    def spec(name, path, *args, **kwargs):
        item = original_spec(name, path, *args, **kwargs)
        if name == "protected_memory_ledger_verifier":
            original_exec = item.loader.exec_module
            def execute(module):
                original_exec(module)
                original_verify = module.verify
                def changed(packet):
                    result = original_verify(packet)
                    result[field] = value
                    return result
                module.verify = changed
            item.loader.exec_module = execute
        return item
    monkeypatch.setattr(ledger.importlib.util, "spec_from_file_location", spec)
    with pytest.raises(SystemExit, match="protected-memory verifier scope/count mismatch"):
        ledger._protected_memory_control()


def test_protected_memory_projection_requires_fresh_provider_bytes(monkeypatch):
    target = ledger.CODE / "thermodynamics/protected_memory/protected_memory.py"
    original = Path.read_bytes
    def altered(path):
        raw = original(path)
        return raw + b"\n# changed provider\n" if path == target else raw
    monkeypatch.setattr(Path, "read_bytes", altered)
    with pytest.raises(SystemExit, match="protected-memory independent replay failed: source pins"):
        ledger._protected_memory_control()


def test_protected_memory_row_keeps_classification_and_lean_declarations(result):
    row = next(item for item in result["sections"]["forced_structure"]
               if item["id"] == "finite_green_kubo_graph_transport")
    assert row["protected_record_memory"] == ledger._protected_memory_control()
    assert row["match"] == (
        "exact finite conditional transport structure; physical generator and coefficients not constructed"
    )
    assert set(row["lean_declarations"]) == {"GreenKubo", "GraphDiffusion"}
    assert len(row["lean_declarations"]["GreenKubo"]) == 8
    assert len(row["lean_declarations"]["GraphDiffusion"]) == 11
    assert "paper/tex_fragments/PROTECTED_RECORD_MEMORY.tex" in row["artifact_refs"]
    assert "centering within each protected fibre" in row["hypothesis_boundary"]
    assert "does not overturn the historical idempotent-projector obstruction" in row["hypothesis_boundary"]
    assert "this row emits no prediction-ladder entry" in row["hypothesis_boundary"]


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
    assert row["historical_owner_issues"] == [736]
    assert row["live_owner_issues"] == [1134]
    assert "historical_context_issues" not in row


ALPHA_CLOSURE_ROWS = (
    ("alpha_inv_closure_root", "thomson_structured_running"),
    (
        "alpha_inv_closure_gauge_width",
        "thomson_structured_running_plus_gauge_width",
    ),
)


def test_alpha_closure_rows_match_the_contraction_certificate(result):
    rows = {row["id"]: row for row in result["sections"]["alpha"]}
    assert set(rows) == {
        "alpha_inv_thomson_endpoint",
        "alpha_inv_closure_root",
        "alpha_inv_closure_gauge_width",
    }
    certificate = json.loads(
        ledger.PARENTS["p_interval_contraction"].read_text(encoding="utf-8")
    )
    endpoint = json.loads(ledger.PARENTS["endpoint"].read_text(encoding="utf-8"))
    measured = float(endpoint["compare_only"]["codata_alpha_inv"])
    for row_id, mode in ALPHA_CLOSURE_ROWS:
        row = rows[row_id]
        block = certificate["modes"][mode]
        point = block["fixed_point_point_estimate_display_only"]
        enclosure = block["certified_enclosure"]["alpha_inv"]
        assert row["mode"] == mode
        assert row["value_central"] == pytest.approx(float(point["alpha_inv"]))
        assert row["p_central"] == pytest.approx(float(point["P"]))
        assert row["value_interval_decimal"] == [
            str(enclosure["lo"]),
            str(enclosure["hi"]),
        ]
        assert row["enclosure_width_decimal"] == str(enclosure["width"])
        assert row["measured"] == pytest.approx(measured)
        assert row["deviation_inv_alpha"] == pytest.approx(
            row["value_central"] - measured
        )
        assert row["relative_deviation"] == pytest.approx(
            (row["value_central"] - measured) / measured
        )
        assert row["tier"] == "T2_conditional"
        assert row["historical_owner_issues"] == [736]
        assert row["live_owner_issues"] == [1134]
        assert ledger._rel("p_interval_contraction") in row["artifact_refs"]


def test_alpha_closure_chain_orders_the_two_maps(result):
    rows = {row["id"]: row for row in result["sections"]["alpha"]}
    root = rows["alpha_inv_closure_root"]
    gauge = rows["alpha_inv_closure_gauge_width"]
    assert abs(root["relative_deviation"]) > abs(gauge["relative_deviation"])
    assert abs(gauge["relative_deviation"]) < 3e-6
    text = ledger._render_md(result)
    assert "`alpha_inv_closure_root`" in text
    assert "`alpha_inv_closure_gauge_width`" in text
    assert f"{abs(gauge['deviation_inv_alpha']):.6f}" in text
    assert "is open work owned by [#1134](https://github.com/FloatingPragma/observer-patch-holography/issues/1134) (historical lane #736, closed)" in text


def test_alpha_closure_rows_fail_closed_when_the_certificate_promotes():
    endpoint = json.loads(ledger.PARENTS["endpoint"].read_text(encoding="utf-8"))
    certificate = json.loads(
        ledger.PARENTS["p_interval_contraction"].read_text(encoding="utf-8")
    )
    for mutation in (
        {"promotion_allowed": True},
        {"exact_alpha_promoted": True},
    ):
        payload = deepcopy(certificate)
        payload.update(mutation)
        with pytest.raises(SystemExit):
            ledger._alpha_closure_rows(endpoint, payload)
    payload = deepcopy(certificate)
    payload["consumer_policy"]["may_feed_live_particle_predictions"] = True
    with pytest.raises(SystemExit):
        ledger._alpha_closure_rows(endpoint, payload)
    payload = deepcopy(certificate)
    payload["modes"]["thomson_structured_running"]["banach"][
        "uniqueness_in_interval"
    ] = False
    with pytest.raises(SystemExit):
        ledger._alpha_closure_rows(endpoint, payload)


def test_lepton_rows_match_parents_and_contain_witness(result):
    rows = {r["id"]: r for r in result["sections"]["charged_leptons"]}
    ratio_parent = json.loads(ledger.PARENTS["koide_identity"].read_text(encoding="utf-8"))
    mcpr = rows["charged_leptons_mcpr_conditional"]
    ratios = mcpr["dimensionless_ratio_comparison"]
    parent_ratios = ratio_parent["compare_only_PDG_2026_ratio_envelope"]
    assert ratios["m_mu_over_m_e"] == parent_ratios["independent_ratios"]["m_mu_over_m_e"]
    assert ratios["m_mu_over_m_e"]["inside_componentwise_one_sigma_interval"] is False
    assert ratios["m_tau_over_m_e_inside_componentwise_one_sigma_interval"] is True
    assert ratios["common_scale_under_identity_mass_scheme_necessary_condition_passed"] is False
    assert mcpr["tier"] == "T2"
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
        assert rows[row_id]["historical_owner_issues"] == [736]
        assert rows[row_id]["live_owner_issues"] == [1134]
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
    assert obstruction["historical_owner_issues"] == [736]
    assert obstruction["live_owner_issues"] == [1134]
    assert "historical_context_issues" not in obstruction
    rscc = rows["quark_conditional_rscc_comparison"]
    parent = json.loads(ledger.PARENTS["conditional_quark_comparison"].read_text(encoding="utf-8"))
    branch = parent["branches"]["full_rscc"]
    assert rscc["flag_derived_ms_over_md"][0]["computed"] == branch["flag_derived_ms_over_md"][0]["computed"]
    assert rscc["new_independent_evidence_claimed"] is False
    assert rscc["pdg_2026_individual_marginals"] == {"count": 6, "all_inside": True}
    pdg = rscc["pdg_2026_ms_over_mud"]
    assert pdg["inside_quoted_marginal_range"] is False
    assert float(rscc["pdg_2026_ratio_from_separately_displayed_mass_centers"]) < float(pdg["quoted_lower"])


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


def test_markdown_rendered(result):
    # The module fixture verifies the written file against this renderer;
    # prose assertions do not require another full numerical parent replay.
    text = ledger._render_md(result)
    assert "# Postdiction Ledger" in text
    assert "necessary single-common-scale condition" in text
    assert "Conditional RSCC ratios" in text
    assert "## Forced structure" in text
    assert "## Quantum carrier gate" in text
    assert "NOT_EVALUABLE" in text
    assert "Recorded retrospective same-scheme accounting interval" in text
    assert "Certified same-scheme anchor gap" not in text
    assert "`code/particles/scripts/build_postdiction_ledger.py`" in text
    assert "`code/particles/runs/status/postdiction_ledger.json`" in text
    assert "Each row distinguishes analytic paper proofs, finite Lean results" in text
    assert "Every step is machine checked in the Lean workspace" not in text
    assert "Live owner: [#1134](https://github.com/FloatingPragma/observer-patch-holography/issues/1134) (historical lane #736, closed)" in text
    assert "Historical context:" not in text
    assert "not identified as restrictions of one fixed physical construction" in text
    assert "no pooled all-scale agreement verdict" in text


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
    koide_statement = principal[3]["statement"]
    assert "conditional postdiction" in koide_statement
    assert "declared premise ancestry" in koide_statement
    assert "Koide's published relation" in koide_statement
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
        "Lean constructs the complete positive normal frame",
        "certified mass forms of the supplied geometry",
        "separate naturality proofs",
        "no common local/global policy is established",
        "Hilbert completion, operator closures and temporal convergence are paper proofs",
        "zero-charge source-free radiative sector",
        "not quantization of the prescribed-charge episode or interacting matter action",
        "sqrt(lambda) and finite-step theta/h are distinct",
        "no uniform operator-norm quantum error",
    ):
        assert qualification in boundary
    assert {
        "Lean/Screen/WhitneySourceGeometry.lean",
        "Lean/Screen/WhitneySourceAssembly.lean",
        "Lean/Screen/WhitneyCertifiedConsumers.lean",
        "Lean/Screen/WhitneySourceNaturality.lean",
    } <= set(row["artifact_refs"])
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
        "This algebra/global-existence result alone supplies no executed charged-matter episode",
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


def test_coupled_whitney_rows_separate_analytic_numeric_and_empirical_scopes(coupled_whitney_rows):
    rows = coupled_whitney_rows
    for row in rows.values():
        assert row["physical_comparison_status"] == "NOT_EVALUABLE"
        assert row["observed_postdiction"] is False
        assert row["continuum_convergence_established"] is False
        assert "no observed postdiction" in row["match"]
        assert "no observed data" in row["observed_counterpart"]
        assert not {"measured", "delta_over_sigma"} & row.keys()
        theorem = row["analytic_paper_theorem"]
        assert "\\label{"+theorem["label"]+"}" in (ledger.REPO/theorem["source"]).read_text(encoding="utf-8")
    spatial = rows["whitney_spatial_consistency"]
    assert spatial["spatial_action_consistency"] is True
    assert "O(delta)" in spatial["statement"]
    assert "no formal spatial convergence theorem" in spatial["lean_scope"]
    assert len(spatial["lean_declarations"]["WhitneySpatialConsistency"]) == 5
    assert spatial["independent_verifier_result"]["tetrahedra"] == [20, 160, 1280]
    assert spatial["independent_verifier_result"]["uniform_bound_certified_by_numerics"] is False
    quantum = rows["whitney_interacting_quantum"]
    assert quantum["hilbert_space_constructed"] is True
    assert quantum["computed_quantum_state_history"] is False
    assert quantum["analytic_proof_formalized_in_lean"] is False
    assert quantum["lean_receipts"] == [] and quantum["lean_declarations"] == {}
    for key in ("metric_completeness_established", "essential_self_adjointness_established",
                "unique_extension_given_ordering", "initial_state_constructed",
                "initial_observables_provided"):
        assert quantum[key] is True
    summary = quantum["independent_verifier_result"]
    assert summary["accepted"] is True
    assert summary["initial_state_constructed"] is True
    assert summary["initial_observables_provided"] is True
    assert all(summary[key] is False for key in ("quantum_time_history", "observer_history", "physical_comparison",
                                               "analytic_operator_domain_proved_by_numeric_replay"))
    assert summary["real_configuration_dimension"] == 56
    assert {key: summary[key] for key in ("gaussian_sigma", "gaussian_norm_squared",
                                        "matter_l2_coefficient", "matter_l4_coefficient")} == {
        "gaussian_sigma": "1/2", "gaussian_norm_squared": "1",
        "matter_l2_coefficient": "1/5", "matter_l4_coefficient": "3/35",
    }
    assert summary["volume_in_Qsqrt5"] == ["10", "10/3"]
    support = quantum["analytic_supporting_results"]
    assert support["labels"] == [
        "eq:whitney-interacting-global-metric-bound",
        "prop:whitney-interacting-gaussian-state",
        "eq:whitney-interacting-gaussian-state",
        "eq:whitney-interacting-gaussian-matter-moments",
        "eq:whitney-interacting-gaussian-magnetic-moment",
    ]
    paper = (ledger.REPO / support["source"]).read_text(encoding="utf-8")
    assert all("\\label{" + label + "}" in paper for label in support["labels"])
    assert "code/electromagnetism/runtime/whitney_quantum_state_receipt.json" in quantum["artifact_refs"]
    assert "code/electromagnetism/verify_whitney_quantum_state.py" in quantum["artifact_refs"]
    shared = quantum["shared_reconstructed_observable_algebra"]
    assert shared == {
        "source": "paper/tex_fragments/WHITNEY_COMMON_OBSERVABLES.tex",
        "labels": ["thm:whitney-common-observables", "eq:whitney-common-observables",
                   "eq:whitney-common-observable-pvm", "eq:whitney-common-observable-law",
                   "cor:whitney-continuum-detector-readouts"],
        "self_adjoint_neutral_multipliers": True,
        "joint_spectral_probability_law": True,
        "conditional_real_sector_detector_convergence": True,
        "ground_state_required": False,
        "quantum_refinement_convergence": False,
        "physical_detector_calibration": False,
        "proved_by_numeric_parent_replay": False,
    }
    assert shared["source"] in quantum["artifact_refs"]
    assert "paper/tex_fragments/WHITNEY_REAL_CONTINUUM.tex" in quantum["artifact_refs"]
    assert "shared_reconstructed_observable_algebra" not in summary
    assert "General bounded Borel detectors have no asserted classical continuum convergence rate" in quantum["hypothesis_boundary"]
    for text in ("actual positive kinetic metric", "Schur metric", "R^30 x C^13", "Friedrichs", "form core", "complete", "essentially self-adjoint", "operator core", "Hamiltonian domain"):
        assert text in quantum["statement"]
    for text in ("hbar", "quantization measure", "ordering", "not formalized in Lean", "unique quantization"):
        assert text in quantum["hypothesis_boundary"]
    charged = rows["whitney_charged_execution"]
    assert charged["rigorous_trajectory_error_enclosure"] is False
    assert charged["authenticated_observer_history"] is False
    assert charged["computed_quantum_state_history"] is False
    assert charged["lean_receipts"] == [] and charged["lean_declarations"] == {}
    summary = charged["independent_verifier_result"]
    assert summary["accepted"] is True and summary["samples"] == 81
    assert summary["full_equations_per_sample"] == 68 and summary["gauss_equations_per_sample"] == 13
    assert summary["symmetry"]["real_configuration_fixed_dimension"] == 5
    assert all(summary[key] is False for key in ("observer_history", "quantum_state", "physical_continuum"))


def test_coupled_projection_is_stable_when_replayed_float_diagnostics_change(coupled_whitney_rows, monkeypatch):
    rows = coupled_whitney_rows
    by_parent = {"spatial_consistency": rows["whitney_spatial_consistency"],
                 "charged_dynamics": rows["whitney_charged_execution"],
                 "quantum_state": rows["whitney_interacting_quantum"]}
    calls = []
    def replay(stem):
        calls.append(stem)
        row = by_parent[stem]
        summary = deepcopy(row["independent_verifier_result"])
        # These are deliberately new values, not patched expected row bytes.
        # A replay checks its own tolerances; the ledger must not serialize
        # its platform-dependent recomputed diagnostic floats.
        summary.update({"numeric_diagnostics": [{"action": 0.123456789}],
                        "full_euler_max_abs": 9.1e-14, "energy_drift": 1.1e-13,
                        "rk4_endpoint_errors": [6.1e-6, 4.1e-7, 2.1e-8]})
        return {"scope": row["certificate_scope"]}, summary
    monkeypatch.setattr(ledger, "_verify_whitney_coupled_parent", replay)
    actual = {row["id"]: row for row in ledger._whitney_coupled_rows()}
    assert calls == ["spatial_consistency", "charged_dynamics", "quantum_state"]
    assert actual == rows


@pytest.mark.parametrize("parent,key,value", [
    ("spatial_consistency", "uniform_bound_certified_by_numerics", True),
    ("spatial_consistency", "continuum_trajectory_claimed", True),
    ("spatial_consistency", "observer_history_claimed", True),
    ("charged_dynamics", "accepted", False),
    ("charged_dynamics", "quantum_state", True),
    ("charged_dynamics", "observer_history", True),
    ("charged_dynamics", "physical_continuum", True),
    ("quantum_state", "accepted", False),
    ("quantum_state", "accepted", 1),
    ("quantum_state", "initial_state_constructed", False),
    ("quantum_state", "initial_observables_provided", False),
    ("quantum_state", "quantum_time_history", True),
    ("quantum_state", "quantum_time_history", 0),
    ("quantum_state", "observer_history", True),
    ("quantum_state", "physical_comparison", True),
    ("quantum_state", "analytic_operator_domain_proved_by_numeric_replay", True),
    ("quantum_state", "real_configuration_dimension", 56.0),
    ("quantum_state", "real_configuration_dimension", 55),
    ("quantum_state", "scope", "OTHER_SCOPE"),
    ("quantum_state", "gaussian_sigma", "0"),
    ("quantum_state", "gaussian_sigma", "-1/2"),
    ("quantum_state", "gaussian_sigma", 0.5),
    ("quantum_state", "gaussian_sigma", "1/0"),
    ("quantum_state", "gaussian_norm_squared", "2"),
    ("quantum_state", "matter_l2_coefficient", "2/10"),
    ("quantum_state", "matter_l4_coefficient", "nan"),
    ("quantum_state", "volume_in_Qsqrt5", ["10", "20/6"]),
    ("quantum_state", "volume_in_Qsqrt5", ["10"]),
])
def test_coupled_builder_rejects_promotion_flags_even_from_provider(
    coupled_whitney_rows, monkeypatch, parent, key, value
):
    by_parent = {"spatial_consistency": coupled_whitney_rows["whitney_spatial_consistency"],
                 "charged_dynamics": coupled_whitney_rows["whitney_charged_execution"],
                 "quantum_state": coupled_whitney_rows["whitney_interacting_quantum"]}
    def replay(stem):
        row = by_parent[stem]
        summary = deepcopy(row["independent_verifier_result"])
        if stem == parent:
            summary[key] = value
        return {"scope": row["certificate_scope"]}, summary
    monkeypatch.setattr(ledger, "_verify_whitney_coupled_parent", replay)
    with pytest.raises(SystemExit, match="coupled Whitney"):
        ledger._whitney_coupled_rows()


@pytest.mark.parametrize("parent,provider", [
    ("spatial_consistency", "Lean/Screen/WhitneySpatialConsistency.lean"),
    ("charged_dynamics", "code/electromagnetism/verify_whitney_charged_dynamics.py"),
    ("quantum_state", "code/electromagnetism/verify_whitney_quantum_state.py"),
])
def test_coupled_fresh_parent_replay_rejects_changed_provider_bytes(monkeypatch, parent, provider):
    original = Path.read_bytes
    path = ledger.REPO/provider
    seen = []
    def altered(target):
        raw = original(target)
        if target == path:
            seen.append(target)
            return raw+b"\n# independent changed-provider mutation\n"
        return raw
    monkeypatch.setattr(Path, "read_bytes", altered)
    with pytest.raises(ValueError, match="source pin"):
        ledger._verify_whitney_coupled_parent(parent)
    assert seen == [path]


def test_coupled_named_spatial_algebra_cannot_be_replaced_by_empty_provider(
    coupled_whitney_rows, tmp_path, monkeypatch
):
    row = coupled_whitney_rows["whitney_spatial_consistency"]
    declarations = row["lean_declarations"]["WhitneySpatialConsistency"]
    empty = tmp_path/"WhitneySpatialConsistency.lean"
    empty.write_text("-- no finite algebra declarations\n", encoding="utf-8")
    monkeypatch.setitem(ledger.LEAN_RECEIPTS, "WhitneySpatialConsistency", empty)
    with pytest.raises(SystemExit, match="Lean declaration missing: WhitneySpatialConsistency"):
        ledger._lean_receipt("WhitneySpatialConsistency", declarations={"WhitneySpatialConsistency": tuple(declarations)})


@pytest.mark.parametrize("missing_label", [
    "thm:whitney-interacting-quantum",
    "eq:whitney-interacting-global-metric-bound",
    "prop:whitney-interacting-gaussian-state",
    "eq:whitney-interacting-gaussian-state",
    "eq:whitney-interacting-gaussian-matter-moments",
    "eq:whitney-interacting-gaussian-magnetic-moment",
    "thm:whitney-common-observables",
    "eq:whitney-common-observables",
    "eq:whitney-common-observable-pvm",
    "eq:whitney-common-observable-law",
    "cor:whitney-continuum-detector-readouts",
])
def test_coupled_missing_analytic_quantum_theorem_is_not_silently_promoted(
    coupled_whitney_rows, monkeypatch, missing_label
):
    by_parent = {"spatial_consistency": coupled_whitney_rows["whitney_spatial_consistency"],
                 "charged_dynamics": coupled_whitney_rows["whitney_charged_execution"],
                 "quantum_state": coupled_whitney_rows["whitney_interacting_quantum"]}
    def replay(stem):
        row = by_parent[stem]
        return {"scope": row["certificate_scope"]}, deepcopy(row["independent_verifier_result"])
    monkeypatch.setattr(ledger, "_verify_whitney_coupled_parent", replay)
    original = Path.read_text
    quantum_path = ledger.REPO/"paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex"
    observable_path = ledger.REPO/"paper/tex_fragments/WHITNEY_COMMON_OBSERVABLES.tex"
    def empty(path, *args, **kwargs):
        source = original(path, *args, **kwargs)
        if path in (quantum_path, observable_path):
            return source.replace("\\label{" + missing_label + "}", "")
        return source
    monkeypatch.setattr(Path, "read_text", empty)
    with pytest.raises(SystemExit, match="analytic theorem missing"):
        ledger._whitney_coupled_rows()


@pytest.fixture(scope="module")
def completion_whitney_rows(result):
    names = {"whitney_real_continuum", "whitney_charged_instrument", "whitney_ephemeris_clock", "whitney_quantum_history"}
    rows = {row["id"]: row for row in result["sections"]["forced_structure"] if row["id"] in names}
    assert set(rows) == names
    return rows


def completion_parent_fixture(rows, stem):
    row = rows["whitney_"+stem]
    summary = deepcopy(row["independent_verifier_result"])
    packet = {"scope": row["certificate_scope"]}
    if stem == "quantum_history":
        packet["error_certificate"] = {"horizon": summary["horizon"], "target_norm_error": summary["target_norm_error"],
                                       "horizon_squared_error_upper": summary["squared_norm_error_upper"]}
    return packet, summary


def test_new_bridges_preserve_precise_conditional_scopes(completion_whitney_rows):
    rows = completion_whitney_rows
    for row in rows.values():
        assert row["observed_postdiction"] is False
        assert row["physical_comparison_status"] == "NOT_EVALUABLE"
        assert row["source_selected_physical_continuum"] is False
        assert row["analytic_proof_formalized_in_lean"] is False
        assert row["lean_receipts"] == [] and row["lean_declarations"] == {}
        assert row["independent_verifier_result"]["accepted"] is True
        proof = row["analytic_paper_theorem"]
        text = (ledger.REPO/proof["source"]).read_text(encoding="utf-8")
        assert "\\label{"+proof["label"]+"}" in text
        for label in row["analytic_supporting_results"]["labels"]:
            assert "\\label{"+label+"}" in text
        assert all((ledger.REPO/path).is_file() for path in row["artifact_refs"])
    real = rows["whitney_real_continuum"]
    assert real["continuum_convergence_established"] is True
    assert "conditional analytic invariant real-sector" in real["continuum_convergence_scope"]
    assert real["independent_verifier_result"]["full_charged_complex_trajectory_bound"] is False
    assert "Ritz" in real["hypothesis_boundary"] and "nodal initialization" in real["hypothesis_boundary"]
    observer = rows["whitney_charged_instrument"]["independent_verifier_result"]
    assert observer["events"] == 1782 and observer["decoded_samples"] == 81
    assert observer["observer_software_history"] is True and observer["physical_observer_placement"] is False
    assert observer["rigorous_trajectory_enclosure"] is False
    checkpoint = rows["whitney_charged_instrument"]["decoded_checkpoint_control"]
    assert checkpoint["independent_verifier_result"]["checkpoint_qv_error_certified"] is True
    assert checkpoint["independent_verifier_result"]["decoded_checkpoints"] == observer["decoded_samples"]
    assert checkpoint["independent_verifier_result"]["events"] == observer["events"]
    assert checkpoint["independent_verifier_result"]["simple_checkpoint_error_upper"] == "10001/100000000000000"
    assert checkpoint["independent_verifier_result"]["continuous_observer_error_certified"] is False
    assert checkpoint["observed_postdiction"] is False
    assert "floating solver is not recomputed" in checkpoint["event_replay_scope"]
    clock = rows["whitney_ephemeris_clock"]["independent_verifier_result"]
    assert clock["source_configurations"] == 81 and clock["calibrated_physical_clock"] is False
    assert "duration" not in clock and "comparison_error" not in clock
    history = rows["whitney_quantum_history"]["independent_verifier_result"]
    assert history["horizon"] == "1/36028797018963968"
    assert history["target_norm_error"] == "1/10"
    assert history["global_time_coverage"] is True and history["trial_history_computed"] is True
    assert history["exact_Hamiltonian_history_computed"] is False
    assert history["ordinary_physics_benchmark"] is False
    assert history["configuration_density_moves"] is False


def test_new_bridge_projection_has_no_float_diagnostic_leak(completion_whitney_rows, monkeypatch):
    calls = []
    def replay(stem):
        calls.append(stem)
        packet, summary = completion_parent_fixture(completion_whitney_rows, stem)
        summary.update(numeric_diagnostics={"changed": .123456789}, duration=1.2345, comparison_error=.01)
        return packet, summary
    monkeypatch.setattr(ledger, "_verify_whitney_completion_parent", replay)
    monkeypatch.setattr(ledger, "_whitney_checkpoint_control", lambda: deepcopy(
        completion_whitney_rows["whitney_charged_instrument"]["decoded_checkpoint_control"]))
    actual = {row["id"]: row for row in ledger._whitney_completion_rows()}
    assert actual == completion_whitney_rows
    assert calls == ["real_continuum", "charged_instrument", "ephemeris_clock", "quantum_history"]
    def no_float(value):
        if isinstance(value, dict):
            return all(no_float(v) for v in value.values())
        if isinstance(value, list):
            return all(no_float(v) for v in value)
        return type(value) is not float
    assert all(no_float(row["independent_verifier_result"]) for row in actual.values())


@pytest.mark.parametrize("stem,key,value", [
    ("real_continuum", "accepted", 1),
    ("real_continuum", "refinement_parameters", [1, 2, 4, 8.0]),
    ("real_continuum", "uniform_shape_bound", "144/2"),
    ("real_continuum", "conditional_real_sector_trajectory_bound", False),
    ("real_continuum", "full_charged_complex_trajectory_bound", True),
    ("real_continuum", "numerical_trajectory_error_certified", True),
    ("real_continuum", "physical_source_or_clock_selected", True),
    ("charged_instrument", "events", 1782.0),
    ("charged_instrument", "decoded_samples", 80),
    ("charged_instrument", "observer_software_history", False),
    ("charged_instrument", "physical_clock_calibrated", True),
    ("charged_instrument", "physical_observer_placement", True),
    ("charged_instrument", "quantum_state_history", True),
    ("charged_instrument", "rigorous_trajectory_enclosure", True),
    ("ephemeris_clock", "accepted", 1),
    ("ephemeris_clock", "source_configurations", 81.0),
    ("ephemeris_clock", "calibrated_physical_clock", True),
    ("ephemeris_clock", "rigorous_numerical_enclosure", True),
    ("quantum_history", "exact_Hamiltonian_history_computed", True),
    ("quantum_history", "configuration_density_moves", True),
    ("quantum_history", "ordinary_physics_benchmark", True),
    ("quantum_history", "global_time_coverage", 1),
    ("quantum_history", "numeric_quadrature_used_for_bound", True),
    ("quantum_history", "horizon", "1/18014398509481984"),
    ("quantum_history", "horizon", "2/72057594037927936"),
    ("quantum_history", "target_norm_error", .1),
    ("quantum_history", "squared_norm_error_upper", "1/99"),
])
def test_new_bridge_provider_promotions_and_types_fail(completion_whitney_rows, monkeypatch, stem, key, value):
    def replay(name):
        packet, summary = completion_parent_fixture(completion_whitney_rows, name)
        if name == stem:
            summary[key] = value
        return packet, summary
    monkeypatch.setattr(ledger, "_verify_whitney_completion_parent", replay)
    monkeypatch.setattr(ledger, "_whitney_checkpoint_control", lambda: deepcopy(
        completion_whitney_rows["whitney_charged_instrument"]["decoded_checkpoint_control"]))
    with pytest.raises(SystemExit, match="Whitney completion"):
        ledger._whitney_completion_rows()


@pytest.mark.parametrize("stem", ["real_continuum", "charged_instrument", "ephemeris_clock", "quantum_history"])
def test_new_bridge_source_bytes_are_replayed_fresh(monkeypatch, stem):
    original = Path.read_bytes
    path = ledger.CODE/"electromagnetism"/f"verify_whitney_{stem}.py"
    seen = []
    def altered(target):
        data = original(target)
        if target == path:
            seen.append(target)
            return data+b"\n# altered independent provider\n"
        return data
    monkeypatch.setattr(Path, "read_bytes", altered)
    with pytest.raises(ValueError, match="pin"):
        ledger._verify_whitney_completion_parent(stem)
    assert seen


@pytest.mark.parametrize("stem", ["real_continuum", "charged_instrument", "ephemeris_clock", "quantum_history"])
def test_new_bridge_missing_theorem_label_fails(completion_whitney_rows, monkeypatch, stem):
    monkeypatch.setattr(ledger, "_verify_whitney_completion_parent",
                        lambda name: completion_parent_fixture(completion_whitney_rows, name))
    monkeypatch.setattr(ledger, "_whitney_checkpoint_control", lambda: deepcopy(
        completion_whitney_rows["whitney_charged_instrument"]["decoded_checkpoint_control"]))
    proof = completion_whitney_rows["whitney_"+stem]["analytic_paper_theorem"]
    original = Path.read_text
    def removed(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        if path == ledger.REPO/proof["source"]:
            return text.replace("\\label{"+proof["label"]+"}", "")
        return text
    monkeypatch.setattr(Path, "read_text", removed)
    with pytest.raises(SystemExit, match="analytic theorem missing"):
        ledger._whitney_completion_rows()


def checkpoint_parent_fixture(rows):
    path = ledger.CODE / "electromagnetism/runtime/whitney_charged_checkpoint_receipt.json"
    raw = path.read_bytes()
    summary = deepcopy(rows["whitney_charged_instrument"]["decoded_checkpoint_control"]["independent_verifier_result"])
    return raw, json.loads(raw), summary


def test_checkpoint_projection_is_exact_and_preserves_parent_specs(completion_whitney_rows, monkeypatch):
    raw, packet, summary = checkpoint_parent_fixture(completion_whitney_rows)
    summary["numeric_diagnostics"] = {"forbidden_in_projection": 0.1}
    monkeypatch.setattr(ledger, "_verify_whitney_checkpoint_parent", lambda path: (raw, packet, summary))
    control = ledger._whitney_checkpoint_control()
    assert control == completion_whitney_rows["whitney_charged_instrument"]["decoded_checkpoint_control"]
    assert control["receipt_pin"] == {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
    assert len(completion_whitney_rows) == 4
    assert "numeric_diagnostics" not in control["independent_verifier_result"]
    assert all(type(x) is not float for x in control["independent_verifier_result"].values())


@pytest.mark.parametrize("key,value", [
    ("events", 1782.0), ("decoded_checkpoints", 80), ("original_qv_dimension", 9),
    ("checkpoint_qv_error_certified", False), ("parent_enclosure_freshly_verified", False),
    ("observer_events_exactly_replayed", 1), ("continuous_observer_error_certified", True),
    ("nonlinear_field_error_certified", True), ("configuration_clock_error_certified", True),
    ("physical_clock_calibrated", True), ("external_signature_attestation", True),
    ("quantum_history", True), ("exact_model_step", "1/80"),
    ("historical_qv_error_upper", "1/100000000000"),
    ("decoded_checkpoint_error_upper", "1/10000000000"),
    ("maximum_decoded_reference_difference", 0.0),
    ("simple_checkpoint_error_upper", "100010/1000000000000000"),
])
def test_checkpoint_projection_rejects_error_flag_or_bound_promotions(completion_whitney_rows, monkeypatch, key, value):
    raw, packet, summary = checkpoint_parent_fixture(completion_whitney_rows)
    summary[key] = value
    if key in packet["bounds"]:
        packet["bounds"][key] = value
    monkeypatch.setattr(ledger, "_verify_whitney_checkpoint_parent", lambda path: (raw, packet, summary))
    with pytest.raises(SystemExit, match="Whitney checkpoint"):
        ledger._whitney_checkpoint_control()


@pytest.mark.parametrize("mutation", ["coherent_frame_substitution", "understated_error", "continuous_flag"])
def test_checkpoint_ledger_replays_actual_forged_receipt(tmp_path, mutation):
    source = ledger.CODE / "electromagnetism/runtime/whitney_charged_checkpoint_receipt.json"
    packet = json.loads(source.read_text(encoding="utf-8"))
    if mutation == "coherent_frame_substitution":
        historical = json.loads((ledger.CODE / "electromagnetism/runtime/whitney_charged_dynamics_receipt.json").read_bytes())
        row = packet["checkpoints"][79]
        row["decoded_qv"] = deepcopy(packet["checkpoints"][78]["decoded_qv"])
        old = historical["samples"][79]
        differences = [abs(ledger.Fraction(x) - ledger.Fraction(y)) for x, y in
                       zip(row["decoded_qv"], old["q_reduced"] + old["v_reduced"], strict=True)]
        row["absolute_reference_differences"] = list(map(str, differences))
        row["checkpoint_error_upper"] = str(ledger.Fraction(1, 10**10) + max(differences))
        maximum = max(ledger.Fraction(x) for r in packet["checkpoints"] for x in r["absolute_reference_differences"])
        packet["bounds"]["maximum_decoded_reference_difference"] = str(maximum)
        packet["bounds"]["decoded_checkpoint_error_upper"] = str(ledger.Fraction(1, 10**10) + maximum)
    elif mutation == "understated_error":
        packet["bounds"]["decoded_checkpoint_error_upper"] = "1/10000000000"
    else:
        packet["contract"]["continuous_observer_error_certified"] = True
    path = tmp_path / "checkpoint.json"
    path.write_text(json.dumps(packet), encoding="utf-8")
    with pytest.raises(SystemExit, match="Whitney checkpoint independent replay failed"):
        ledger._whitney_checkpoint_control(path)


def test_checkpoint_ledger_executes_fresh_source_bytes(monkeypatch):
    original = Path.read_bytes
    verifier_path = ledger.CODE / "electromagnetism/verify_whitney_charged_checkpoint.py"
    seen = []
    def altered(path):
        data = original(path)
        if path == verifier_path:
            seen.append(path)
            return data + b"\n# altered checkpoint provider\n"
        return data
    monkeypatch.setattr(Path, "read_bytes", altered)
    with pytest.raises(SystemExit, match="Whitney checkpoint independent replay failed"):
        ledger._whitney_checkpoint_control()
    assert len(seen) >= 2


def test_count_clock_is_structural_and_independently_replayed(result):
    row = next(r for r in result["sections"]["forced_structure"]
               if r["id"] == "source_derived_finite_one_three_causal_carrier")
    clock = row["declared_count_clock_control"]
    assert clock["independent_verifier_result"]["interval_counts"] == [2, 41, 80]
    assert clock["independent_verifier_result"]["authenticated_events"] == 500
    assert clock["count_decoder_uses_timestamps_or_density"] is False
    assert clock["count_volume_or_curve_limit_formalized_in_lean"] is False
    assert clock["observed_postdiction"] is False
    assert clock["independent_verifier_result"]["finite_clock_accuracy_certified"] is False
    assert clock["independent_verifier_result"]["physical_clock_identified"] is False
    assert row["operational_cone_selection"]["native_operational_boost_covariance_derived"] is False


@pytest.mark.parametrize("mutation", ["count", "physical_promotion"])
def test_count_clock_ledger_rejects_forged_evidence(tmp_path, mutation):
    source = ledger.CODE / "causal_refinement/source_count_clock_receipt.json"
    item = json.loads(source.read_text(encoding="utf-8"))
    if mutation == "count":
        item["interval_counts"][1] += 1
    else:
        item["scope"]["physical_clock_identified"] = True
    path = tmp_path / "clock.json"
    path.write_text(json.dumps(item), encoding="utf-8")
    with pytest.raises(SystemExit, match="count-clock independent replay failed"):
        ledger._source_count_clock_control(path)


@pytest.fixture(scope="module")
def common_scalar_control():
    return ledger._common_source_scalar_control()


def test_source_scalar_control_uses_full_replay_and_exact_packet_bounds(common_scalar_control):
    control = common_scalar_control
    packet = json.loads((ledger.REPO / control['receipt']).read_bytes())
    last = packet['levels'][-1]
    summary = control['independent_verifier_result']
    assert summary['resolved_q'] == last['q']
    assert summary['full_tensor_oscillators'] == last['dynamic_oscillators']
    assert summary['error_upper'] == last['graph_vs_compact_continuum_error_upper'][1]
    assert summary['resolved_signal_lower'] == last['resolved_graph_signal_lower']
    assert ledger.Fraction(summary['error_upper']) < ledger.Fraction(summary['resolved_signal_lower'])
    assert control['mathematical_replay'] is True
    assert control['analytic_free_scalar_detector_limit'] is True
    for key in ('observed_postdiction', 'source_action_population_or_physical_clock_selected',
                'field_history_joined_to_causal_count_clock', 'interacting_quantum_continuum',
                'full_continuum_theorem_formalized_in_lean'):
        assert control[key] is False


def test_source_scalar_control_rejects_coherent_error_and_signal_forgery(tmp_path):
    packet = json.loads((ledger.CODE / 'source_scalar_packet/source_common_scalar_receipt.json').read_bytes())
    row = packet['levels'][-1]
    row['graph_vs_compact_continuum_error_upper'] = ['0', '0']
    row['resolved_graph_signal_lower'] = '1/2'
    row['graph_compact_induced_response'] = ['-1/2', '-1/2']
    row['error_smaller_than_resolved_signal'] = True
    path = tmp_path / 'coherent-scalar-forgery.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        ledger._common_source_scalar_control(path)


def test_source_quadrature_control_preserves_exact_formal_scope():
    control = ledger._source_population_quadrature_control()
    assert control['formal_actual_cell_integral_and_quadrature_limit'] is True
    assert 'formal_golden_partition_geometry_or_causal_pair_count_limit' not in control
    assert control['formal_golden_partition_geometry'] is True
    assert control['formal_causal_pair_count_limit'] is False
    assert control['golden_cell_mass'] == 'L^3/q^3'
    assert control['golden_assignment_bound'] == '2*sqrt(3)*L/q'
    assert control['fixed_globally_lipschitz_detector_convergence'] is True
    assert 'arbitrary_varying_detector_or_indicator_convergence' not in control
    assert control['fixed_measurable_null_frontier_indicator_convergence'] is True
    assert control['formal_actual_source_time_grid_count_limit'] is True
    assert control['arbitrary_varying_detector_convergence'] is False
    assert control['moving_causal_interval_limit_formalized'] is False
    assert control['actual_causal_cone_null_frontier_identified_in_lean'] is False
    assert control['source_time_step'] == 'Delta_q=L/(c*sqrt(q))'
    assert control['last_temporal_cell_error_upper'] == 'L^3*Delta_q'
    assert 'source_event_count_tendsto' in control['lean_declarations']['GoldenSourceCountLimit']
    assert 'aligned_endpoint_layer_control' in control['lean_declarations']['GoldenSourceCountLimit']
    assert control['equal_cell_weights_are_lumped_action_weights'] is False
    assert 'golden_quadrature_closed' in control['lean_declarations']['GoldenSourceAssignment']
    assert control['tent_integrand_lipschitz_coefficient'] == 'K'
    assert 'partition_quadrature_tendsto' in control['lean_declarations']['SourcePopulationQuadrature']
    assert control['observed_postdiction'] is False


def test_protected_population_control_fresh_replay():
    control = ledger._protected_population_control()
    packet = json.loads((ledger.REPO / control['receipt']).read_bytes())
    for key, value in packet['summary'].items():
        assert control['independent_verifier_result'][key] == value
    assert control['mathematical_replay'] is True
    assert control['independent_verifier_result']['source_population_produced'] is False
    assert control['independent_verifier_result']['persistent_memory_supplied'] is True
    assert control['observed_postdiction'] is False


def test_protected_population_control_rejects_resealed_missing_parent(tmp_path):
    packet = json.loads((ledger.CODE / 'source_population/runtime/population_receipt.json').read_bytes())
    next(event for event in packet['events'] if event['op'] == 'pair_mean')['parents'] = []
    previous = '0' * 64
    for event in packet['events']:
        event['previous_hash'] = previous
        raw = (json.dumps({k:v for k,v in event.items() if k != 'event_hash'},
                          sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
        previous = event['event_hash'] = hashlib.sha256(raw).hexdigest()
    packet['final_event_hash'] = previous
    path = tmp_path / 'resealed-population.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='read-from parent'):
        ledger._protected_population_control(path)


def test_local_sm_control_labels_custody_only_and_reads_inventory():
    control = ledger._local_sm_action_control()
    packet = json.loads((ledger.REPO / control['receipt']).read_bytes())
    assert control['mathematical_replay'] is False
    assert control['custody_verifier_result']['mathematical_replay'] is False
    assert control['packet_inventory']['action_sectors'] == len(packet['action_coefficients'])
    assert control['packet_inventory']['matter_currents'] == len(packet['matter_current_coefficients'])
    assert control['observed_postdiction'] is False
    assert 'accepted' not in control


def test_local_sm_control_rejects_quantum_promotion(tmp_path):
    packet = json.loads((ledger.CODE / 'sm_local_action/local_action_receipt.json').read_bytes())
    packet['scope']['renormalized_quantum_theory'] = True
    path = tmp_path / 'promoted-action.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='scientific scope'):
        ledger._local_sm_action_control(path)


def test_local_sm_control_rejects_custody_as_mathematical_replay(monkeypatch):
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (
        {}, {'accepted_custody': True, 'mathematical_replay': True}, {}))
    with pytest.raises(SystemExit, match='masquerade'):
        ledger._local_sm_action_control()


def test_fresh_structural_loader_reloads_same_size_source_and_restores_helper(tmp_path):
    import types
    import sys
    helper = tmp_path / 'isolated_test_helper.py'
    helper.write_text('VALUE=1\n', encoding='utf-8')
    verifier = tmp_path / 'check.py'
    verifier.write_text('from isolated_test_helper import VALUE\n', encoding='utf-8')
    existing = types.ModuleType('isolated_test_helper')
    existing.VALUE = 999
    saved = sys.modules.get('isolated_test_helper')
    sys.modules['isolated_test_helper'] = existing
    try:
        assert ledger._fresh_structural_verifier(verifier, 'isolated_test_helper').VALUE == 1
        stamp = helper.stat()
        helper.write_text('VALUE=2\n', encoding='utf-8')
        import os
        os.utime(helper, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
        assert ledger._fresh_structural_verifier(verifier, 'isolated_test_helper').VALUE == 2
        assert sys.modules['isolated_test_helper'] is existing
    finally:
        if saved is None:
            sys.modules.pop('isolated_test_helper', None)
        else:
            sys.modules['isolated_test_helper'] = saved


def test_spatial_and_local_action_enrich_existing_structural_rows(result):
    rows = {r['id']:r for r in result['sections']['forced_structure']}
    assert len(rows) == 59
    spatial = rows['source_derived_finite_one_three_causal_carrier']
    assert spatial['common_free_scalar_control']['mathematical_replay'] is True
    assert spatial['protected_population_control']['mathematical_replay'] is True
    assert spatial['metric_scalar_continuum']['quantum_continuum_claim'] is False
    assert rows['hypercharge_spectrum']['declared_local_action_control']['mathematical_replay'] is False


def test_scalar_execution_projection_uses_full_exact_replay():
    control = ledger._source_scalar_execution_control()
    assert control['mathematical_replay'] is True
    assert control['all_field_arithmetic_exact'] is True
    summary = control['independent_verifier_result']
    assert summary['events_replayed'] == 5888
    assert summary['dynamic_field_reads_replayed'] == 34944
    assert summary['executed_times'] == 21
    assert ledger.Fraction(summary['full_mass_norm_error_final'][1]) < ledger.Fraction(1, 100)
    for key in ('observed_postdiction', 'spatial_continuum_error_certified',
                'field_history_joined_to_causal_count_clock', 'quantum_covariance_or_probability_queried'):
        assert control[key] is False


def test_scalar_execution_projection_rejects_resealed_wrong_actual_writer(tmp_path):
    packet = json.loads((ledger.CODE / 'source_scalar_execution/source_scalar_execution_receipt.json').read_bytes())
    trace = packet['traces']['ascending_baseline']
    # Both old values are zero, but a real different writer did not write this resource.
    event = trace['events'][128]
    assert event['reads'][0][2] == [0, 0]
    event['reads'][0][2] = [0, 1]
    previous = '0'*64
    for event in trace['events']:
        event['ledger_parent'] = previous
        raw = (json.dumps({k: event[k] for k in ('id', 'reads', 'write', 'ledger_parent')},
                          sort_keys=True, separators=(',', ':'), ensure_ascii=True)+'\n').encode('ascii')
        previous = event['hash'] = hashlib.sha256(raw).hexdigest()
    trace['final_hash'] = previous
    path = tmp_path/'false-actual-writer.json'
    path.write_bytes(json.dumps(packet).encode('ascii'))
    with pytest.raises(SystemExit, match='read/writer/version/value'):
        ledger._source_scalar_execution_control(path)


def test_scalar_execution_projection_rejects_clock_promotion(tmp_path):
    packet = json.loads((ledger.CODE / 'source_scalar_execution/source_scalar_execution_receipt.json').read_bytes())
    packet['scope']['old_count_clock_theorem_applied'] = True
    path = tmp_path/'false-clock.json'
    path.write_bytes(json.dumps(packet).encode('ascii'))
    with pytest.raises(SystemExit, match='scientific scope'):
        ledger._source_scalar_execution_control(path)


def test_cartan_execution_projection_full_replay_is_not_continuous_enclosure():
    control = ledger._cartan_scalar_execution_control()
    assert control['mathematical_replay'] is True
    assert control['exact_symbolic_local_classical_reduction'] is True
    assert control['trajectory_values_numerical'] is True
    assert control['continuous_history_error_enclosure'] is False
    assert control['global_compact_subgroup_or_quantum_truncation'] is False
    assert control['observed_postdiction'] is False
    summary = control['independent_verifier_result']
    assert summary['exact_first_current'] == '2/125'
    assert summary['exact_first_electric_kick'] == '1/25000'
    assert summary['events'] == 5678
    assert summary['producer_reverse_diagnostic_independently_replayed'] is False


def test_cartan_execution_projection_rejects_wrong_current_beyond_custody(tmp_path):
    packet = json.loads((ledger.CODE / 'sm_abelian_reduction/abelian_receipt.json').read_bytes())
    packet['exact_first_edge_control']['current'] = '0'
    path = tmp_path/'false-Cartan-current.json'
    path.write_bytes(json.dumps(packet).encode('ascii'))
    with pytest.raises(SystemExit, match='exact current witness'):
        ledger._cartan_scalar_execution_control(path)


def test_cartan_execution_projection_rejects_custody_in_place_of_replay(monkeypatch):
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (
        {}, {'accepted_custody': True, 'mathematical_replay': False}, {}))
    with pytest.raises(SystemExit, match='distinguish mathematical replay'):
        ledger._cartan_scalar_execution_control()


def test_scalar_and_cartan_executions_preserve_existing_rows(result):
    rows = {r['id']: r for r in result['sections']['forced_structure']}
    assert len(rows) == 59
    spatial = rows['source_derived_finite_one_three_causal_carrier']
    action = rows['hypercharge_spectrum']
    assert spatial['authenticated_scalar_execution_control']['mathematical_replay'] is True
    assert spatial['common_free_scalar_control']['field_history_joined_to_causal_count_clock'] is False
    assert action['cartan_scalar_execution_control']['mathematical_replay'] is True
    assert action['declared_local_action_control']['mathematical_replay'] is False
    assert action['match'] == 'exact'


def test_scalar_quantum_projection_propagates_original_covariance_with_full_parent():
    control = ledger._source_scalar_quantum_control()
    assert control['mathematical_replay'] is True
    assert control['original_finite_action_vacuum_retained'] is True
    assert control['quantum_covariance_and_bounded_effect_probability'] is True
    summary = control['independent_verifier_result']
    assert summary['parent_full_mathematical_replay'] is True
    assert summary['parent_events_replayed'] == 5888
    assert len(summary['resolved_steps']) == 15 and summary['times'] == 21
    row = summary['step16']
    assert row['signal_sign'] == -1
    assert ledger.Fraction(row['probability_error_upper']) < ledger.Fraction(row['split_response_abs_lower'])
    for key in ('observed_postdiction', 'observed_quantum_outcomes',
                'spatial_continuum_error_certified', 'intermediate_time_error_enclosure',
                'field_history_joined_to_causal_count_clock'):
        assert control[key] is False
    original = json.loads((ledger.CODE/'source_scalar_execution/source_scalar_execution_receipt.json').read_bytes())
    assert original['scope']['quantum_covariance_or_probability_queried'] is False


def test_scalar_quantum_projection_rejects_omitted_covariance(tmp_path):
    packet = json.loads((ledger.CODE/'source_scalar_quantum/quantum_probability_receipt.json').read_bytes())
    packet['covariance']['increase_upper'] = '0'
    path = tmp_path/'false-stationary-split-vacuum.json'
    path.write_bytes(json.dumps(packet).encode('ascii'))
    with pytest.raises(SystemExit, match='original-vacuum probability certificate'):
        ledger._source_scalar_quantum_control(path)


def test_scalar_quantum_projection_rejects_continuum_promotion(tmp_path):
    packet = json.loads((ledger.CODE/'source_scalar_quantum/quantum_probability_receipt.json').read_bytes())
    packet['scope']['spatial_continuum_error_transfer'] = True
    path = tmp_path/'false-continuum-probability.json'
    path.write_bytes(json.dumps(packet).encode('ascii'))
    with pytest.raises(SystemExit, match='original-vacuum probability certificate'):
        ledger._source_scalar_quantum_control(path)


def test_scalar_quantum_projection_rejects_arithmetic_only_parent(monkeypatch):
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (
        {}, {'verdict': 'PASS', 'parent_full_mathematical_replay': False}, {}))
    with pytest.raises(SystemExit, match='full authenticated parent replay'):
        ledger._source_scalar_quantum_control()


def test_clock_quantum_projection_retains_all_times_and_historical_scopes(result):
    rows = {r['id']: r for r in result['sections']['forced_structure']}
    assert len(rows) == 59
    spatial = rows['source_derived_finite_one_three_causal_carrier']
    control = spatial['reconstructed_clock_quantum_control']
    summary = control['independent_verifier_result']
    assert summary['times'] == 21
    assert len(summary['resolved_steps']) == 15
    assert sorted(summary['resolved_steps'] + summary['unresolved_steps']) == list(range(1, 22))
    assert control['reference_probability_bound_for_every_time_in_each_interval'] is True
    assert control['clock_only_uncertainty'] is True
    assert control['hidden_common_stationary_history_required'] is True
    row = summary['step16']
    assert ledger.Fraction(row['total_probability_error_upper']) > ledger.Fraction(row['original_probability_error_upper'])
    assert ledger.Fraction(row['continuous_response_abs_lower']) > 0
    for key in ('observed_postdiction', 'record_noise_certifies_quantum_state_error',
                'arbitrary_noisy_history_existence_proved', 'continuous_split_trajectory_enclosed',
                'spatial_continuum_error_certified', 'field_history_joined_to_causal_count_clock'):
        assert control[key] is False
    assert spatial['finite_scalar_quantum_probability_control']['intermediate_time_error_enclosure'] is False
    assert spatial['authenticated_scalar_execution_control']['quantum_covariance_or_probability_queried'] is False


@pytest.mark.parametrize('missing', ['full_clock_replay', 'full_original_quantum_replay'])
def test_clock_quantum_projection_requires_both_parent_proofs(monkeypatch, missing):
    summary = {'verdict': 'PASS', 'full_clock_replay': True,
               'full_original_quantum_replay': True, 'parent_events_replayed': 5888}
    summary[missing] = False
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: ({}, summary, {}))
    with pytest.raises(SystemExit, match='both full parent proofs'):
        ledger._source_scalar_clock_quantum_control()


def test_clock_quantum_projection_rejects_forged_physical_clock(tmp_path):
    packet = json.loads((ledger.CODE/'source_scalar_clock_quantum/clock_quantum_receipt.json').read_bytes())
    packet['scope']['physical_clock_or_count_volume_identified'] = True
    path = tmp_path/'false-physical-clock.json'
    path.write_text(json.dumps(packet), encoding='ascii')
    with pytest.raises(SystemExit, match='independent all-time-interval quantum comparison'):
        ledger._source_scalar_clock_quantum_control(path)


def test_time_refinement_projection_separates_mathematical_limit_from_execution(result):
    rows = {r['id']: r for r in result['sections']['forced_structure']}
    spatial = rows['source_derived_finite_one_three_causal_carrier']
    control = spatial['joint_scalar_time_refinement_control']
    summary = control['independent_verifier_result']
    assert summary['levels'] == 4 and summary['q'] == 233
    assert ledger.Fraction(summary['time_step']) == ledger.Fraction(1, 2**17)
    assert ledger.Fraction(summary['nearby_continuum_error_upper']) < ledger.Fraction(summary['split_response_lower'])
    assert control['joint_scalar_space_time_detector_limit'] is True
    assert control['asymptotic_prearrival_detector_convergence'] is True
    assert control['nearest_in_window_update_readout'] is True
    for key in ('complete_observer_event_log_executed', 'q5_clock_error_or_preparation_transferred',
                'finite_q_prearrival_zero_response_certified', 'raw_ancestry_or_count_volume_identified',
                'full_fock_state_norm_convergence', 'observed_postdiction'):
        assert control[key] is False
    assert spatial['common_free_scalar_control']['field_history_joined_to_causal_count_clock'] is False


def test_time_refinement_projection_rejects_execution_promotion(tmp_path):
    packet = json.loads((ledger.CODE/'source_scalar_time_refinement/source_scalar_time_refinement_receipt.json').read_bytes())
    packet['scope']['observer_event_log_executed'] = True
    path = tmp_path/'false-q233-execution.json'
    path.write_text(json.dumps(packet), encoding='ascii')
    with pytest.raises(SystemExit, match='independent time refinement reconstruction'):
        ledger._source_scalar_time_refinement_control(path)


def test_time_refinement_projection_requires_fresh_spatial_proof(monkeypatch):
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (
        {}, {'verdict': 'PASS', 'full_spatial_parent_replayed': False,
             'observer_event_log_executed': False}, {}))
    with pytest.raises(SystemExit, match='full spatial proof'):
        ledger._source_scalar_time_refinement_control()


def test_regional_frontier_control_retains_full_exterior_and_weyl_scope():
    control = ledger._source_scalar_regional_control()
    report = control['independent_verifier_result']
    assert report['region_count'] == 8 and report['missing_row_counterexamples'] == 62
    central = next(r for r in report['regions'] if r['region'] == 'central_cube')
    assert central == {'region': 'central_cube', 'sites': 8,
                       'coordinate_collar_sites': 24, 'minimal_linear_collar': 8}
    assert control['all_real_weyl_parameters_are_supplied'] is True
    for flag in ('classical_record_trace_provides_weyl_access', 'physical_regional_time_slice',
                 'observed_quantum_outcomes', 'observed_postdiction',
                 'weyl_algebra_and_rank_consequences_formalized'):
        assert control[flag] is False


def test_routing_frontier_control_retains_complete_cost_and_noise_scope():
    control = ledger._source_seam_routing_control()
    report = control['independent_verifier_result']
    assert len(report) == 12 and sum(r['events'] for r in report) == 1449
    last = control['baseline_costs'][-1]
    assert (last['path_edges'], last['seam_mean_operations'], last['all_events'],
            last['protected_sample_words'], last['max_rational_storage_bits']) == (24, 300, 351, 25, 111)
    assert report[-1]['source_noise_gain'] == '1855077841'
    gain = control['uniform_sample_error_gain']
    assert gain['sharp_for_independent_bounded_sample_errors'] is True
    assert gain['all_depths_theorem_formalized'] is False
    assert control['conditioning_analytic_proof']['label'] == 'thm:source-seam-exactconditioning'
    bound = control['supplied_gate_sample_and_decoder_error_bound']
    assert bound['initial_vector_error_bound'] == 'g_d*(T*eta+sigma)+delta'
    assert bound['all_error_budgets_are_supplied'] is True
    assert bound['sharp_gate_error_bound'] is False
    assert bound['attained_physical_precision'] is False
    assert bound['formalized_in_lean'] is False
    for row in report:
        # The Pell recurrence is independently evaluated against the replayed
        # source row of each finite observation inverse.
        a, b = 1, 3
        for _ in range(row['path_edges']):
            a, b = b, 2 * b + a
        assert ledger.Fraction(row['source_noise_gain']) == a
    for flag in ('destination_decoder_reads_remote_baselines',
                 'receipt_noise_bound_includes_dynamical_or_arithmetic_error',
                 'full_metric_neighbor_order_refinement', 'observed_postdiction'):
        assert control[flag] is False


def test_charged_frontier_control_is_sibling_of_immutable_execution(result):
    action = next(r for r in result['sections']['forced_structure'] if r['id'] == 'hypercharge_spectrum')
    old = action['cartan_scalar_execution_control']
    new = action['cartan_scalar_continuous_readout_control']
    assert old['continuous_history_error_enclosure'] is False
    assert new['continuous_classical_current_and_readout_bound'] is True
    assert new['immutable_execution_parent_sha256'] == old['receipt_pin']['sha256']
    assert new['independent_verifier_result']['completed_checkpoints'] == 9
    assert ledger.Fraction(new['neighbor_response']['magnitude_lower_at_end']) == ledger.Fraction(9, 200000)
    for flag in ('first_kick_is_exact_time_sample', 'whole_state_numerical_trajectory_enclosed',
                 'quantum_dynamics_or_born_outcomes', 'observed_postdiction'):
        assert new[flag] is False


@pytest.mark.parametrize('control,relative,mutation', [
    ('_source_scalar_regional_control', 'source_scalar_regional/regional_time_slice_receipt.json', 'collar'),
    ('_source_seam_routing_control', 'source_routing/runtime/path_tomography_receipt.json', 'events'),
    ('_cartan_scalar_continuous_readout_control', 'sm_abelian_readout/readout_receipt.json', 'error'),
])
def test_completion_frontier_receipt_forgery_is_replayed(control, relative, mutation, tmp_path):
    packet = json.loads((ledger.CODE/relative).read_text(encoding='utf-8'))
    if mutation == 'collar':
        packet['regions'][3]['minimal_linear_collar_dimension'] = 0
        packet['regions'][3]['no_collar_recovers_original_regional_algebra'] = True
    elif mutation == 'events':
        packet['episodes'][-1]['events'].pop()
    else:
        packet['checkpoints'][2]['error_upper'] = '0'
    path = tmp_path/'forged.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        getattr(ledger, control)(path)


@pytest.mark.parametrize('control,relative', [
    ('_source_scalar_regional_control', 'source_scalar_regional/verify_regional_time_slice.py'),
    ('_source_seam_routing_control', 'source_routing/verify_routing.py'),
    ('_cartan_scalar_continuous_readout_control', 'sm_abelian_readout/verify_readout.py'),
])
def test_completion_frontier_uses_fresh_verifier_bytes(control, relative, monkeypatch):
    path = ledger.CODE/relative
    original = Path.read_bytes
    def altered(p):
        content = original(p)
        if p == path:
            content += b"\ndef verify(*args, **kwargs):\n    raise ValueError('fresh frontier byte control')\n"
        return content
    monkeypatch.setattr(Path, 'read_bytes', altered)
    with pytest.raises(SystemExit, match='fresh frontier byte control'):
        getattr(ledger, control)()


def test_regional_frontier_refreshes_parent_import_and_restores_environment(monkeypatch):
    sentinel = ledger.types.ModuleType('verify_source_scalar_execution')
    monkeypatch.setitem(ledger.sys.modules, 'verify_source_scalar_execution', sentinel)
    before = ledger.sys.path[:]
    assert ledger._source_scalar_regional_control()['mathematical_replay'] is True
    assert ledger.sys.modules['verify_source_scalar_execution'] is sentinel
    assert ledger.sys.path == before


def test_charged_frontier_rejects_omitted_parent_replay(monkeypatch):
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (
        {}, {'verified': True, 'mathematical_replay': False, 'parent_events': None}, {}))
    with pytest.raises(SystemExit, match='full fresh parent replay'):
        ledger._cartan_scalar_continuous_readout_control()


@pytest.mark.parametrize('name', ['golden_floor', 'cell_mass', 'golden_quadrature',
                                'golden_quadrature_tendsto', 'one_cell_width_counterexample'])
def test_golden_frontier_requires_actual_lean_declaration(name, monkeypatch, tmp_path):
    original = ledger.LEAN_RECEIPTS['GoldenSourceAssignment']
    source = original.read_text(encoding='utf-8').replace('theorem '+name+' ', 'theorem removed_'+name+' ')
    path = tmp_path/'GoldenSourceAssignment.lean'
    path.write_text(source, encoding='utf-8')
    monkeypatch.setitem(ledger.LEAN_RECEIPTS, 'GoldenSourceAssignment', path)
    with pytest.raises(SystemExit, match='Lean declaration missing'):
        ledger._source_population_quadrature_control()


def test_feedback_transport_control_fresh_complete_replay():
    control = ledger._source_feedback_transport_control()
    assert control['independent_verifier_result'] == {'episodes':16,'events':4928,'hops':740}
    assert control['mathematical_replay'] is True
    assert control['semantic_read_from_order_checked_for_all_logical_pairs'] is True
    assert control['baseline_costs'][-1]['protected_scalar_registers'] == 141
    assert control['baseline_costs'][-1]['scratch_scalar_registers'] == 25
    assert 'exact_induced_logical_order' in control['lean_declarations']['SourceFeedbackTransport']
    assert control['copy_reset_factor_two_readout_and_schedule_supplied'] is True
    for key in ('source_selected_feedback_or_quantum_operations',
                'full_q13_q21_routed_histories_executed',
                'physical_causal_order_capacity_or_clock_identified', 'observed_postdiction'):
        assert control[key] is False


@pytest.mark.parametrize('mutation', ['mean', 'capacity', 'physical_scope'])
def test_feedback_transport_control_rejects_forged_evidence(tmp_path, mutation):
    packet = json.loads((ledger.CODE/'source_feedback_transport/transport_receipt.json').read_bytes())
    if mutation == 'mean':
        run = packet['episodes'][0]
        event = next(e for e in run['events'] if e['op'] == 'mean')
        event['writes'][0]['value'] = str(ledger.Fraction(event['writes'][0]['value'])+1)
        previous = '0'*64
        for event in run['events']:
            event['previous_hash'] = previous
            raw = (json.dumps({k:v for k,v in event.items() if k!='event_hash'},
                sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
            previous = event['event_hash'] = hashlib.sha256(raw).hexdigest()
        run['event_root'] = previous
    elif mutation == 'capacity':
        packet['episodes'][0]['costs']['protected_scalar_registers'] -= 1
    else:
        packet['scope'] = 'PHYSICAL_SOURCE_DERIVED_TRANSPORT'
    path = tmp_path/'transport.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        ledger._source_feedback_transport_control(path)


def test_sequential_instrument_control_replays_full_parents_and_preserves_scope():
    control = ledger._source_scalar_sequential_instrument_control()
    result = control['independent_verifier_result']
    assert result['mathematical_replay'] is True and result['parent_events_replayed'] == 5888
    assert result['readout_times'] == 21 and result['prior_pair_brackets'] == 210
    assert len(result['uniformly_resolved_steps']) == 15
    assert len(control['continuous_comparison_rows']) == 21
    assert len(control['regional_gadget']['ghz_cnot_tree']) == 31
    assert control['reported_probabilities_are_unconditional_sequential_marginals'] is True
    assert control['reference_has_same_repeated_instrument'] is True
    for key in ('sampled_quantum_outcome_histories', 'branch_conditioned_clock_reconstructed',
                'pointer_action_graph_identified_with_W12_seams',
                'finite_duration_or_noisy_quantum_gates_certified',
                'instrument_theorems_formalized_in_lean', 'physical_clock_or_observed_outcomes',
                'observed_postdiction'):
        assert control[key] is False


@pytest.mark.parametrize('mutation', ['error', 'commutator', 'gadget', 'outcome_scope'])
def test_sequential_instrument_control_rejects_forged_evidence(tmp_path, mutation):
    packet = json.loads((ledger.CODE/'source_scalar_instruments/sequential_instrument_receipt.json').read_bytes())
    if mutation == 'error':
        packet['rows'][-1]['sequential_probability_error_upper'] = '0'
    elif mutation == 'commutator':
        packet['cross_time_brackets'][0]['commutator_over_i_Qphi'] = ['0','0']
    elif mutation == 'gadget':
        packet['regional_gadget']['ideal_operations_per_readout']['edge_cnots'] -= 1
    else:
        packet['scope']['physical_clock_or_observed_outcomes'] = True
    path = tmp_path/'instrument.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        ledger._source_scalar_sequential_instrument_control(path)


def test_sequential_instrument_control_rejects_arithmetic_only_replay(monkeypatch):
    packet = json.loads((ledger.CODE/'source_scalar_instruments/sequential_instrument_receipt.json').read_bytes())
    summary = {'verified':True,**packet['summary'],'mathematical_replay':False,
               'parent_events_replayed':None,'physical_outcomes':False,'branch_conditioned_clock':False}
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k:(packet,summary,{}))
    with pytest.raises(SystemExit, match='full fresh parent replay'):
        ledger._source_scalar_sequential_instrument_control()


def test_feedback_transport_control_reads_current_verifier_bytes(monkeypatch):
    original = Path.read_bytes
    target = ledger.CODE/'source_feedback_transport/verify_transport.py'
    def altered(path):
        data = original(path)
        if path == target:
            data += b'\ndef verify(packet):\n    raise ValueError("fresh feedback verifier sentinel")\n'
        return data
    monkeypatch.setattr(Path, 'read_bytes', altered)
    with pytest.raises(SystemExit, match='fresh feedback verifier sentinel'):
        ledger._source_feedback_transport_control()


def test_sequential_instrument_control_reads_current_public_verifier_bytes(monkeypatch):
    original = Path.read_bytes
    target = ledger.CODE/'source_scalar_instruments/verify_sequential_instrument.py'
    def altered(path):
        data = original(path)
        if path == target:
            data += b'\ndef verify(packet):\n    raise ValueError("fresh public instrument verifier sentinel")\n'
        return data
    monkeypatch.setattr(Path, 'read_bytes', altered)
    with pytest.raises(SystemExit, match='fresh public instrument verifier sentinel'):
        ledger._source_scalar_sequential_instrument_control()


def test_finite_instrument_control_retains_noise_baseline_and_physical_boundary():
    control = ledger._source_scalar_finite_instrument_control()
    proof = control['independent_verifier_result']
    assert proof['forced_oscillator_identities_verified'] is True
    assert proof['parent_events_replayed'] == 5888
    assert proof['readouts'] == 21 and proof['prior_pairs'] == 210
    assert proof['retained_operations'] == 2709 and len(proof['resolved_steps']) == 15
    row = control['comparison_rows'][15]
    assert Fraction(row['paired_response_noise_error_upper']) == 2*Fraction(row['each_probability_noise_error_upper'])
    assert Fraction(row['noisy_paired_response_abs_lower']) > Fraction('0.0082353')
    assert Fraction(row['noisy_baseline_probability_interval'][0]) < Fraction(1, 2)
    for key in ('physical_clock_or_observed_outcomes', 'native_quantum_controls',
                'noisy_energy_bound', 'branch_conditioned_clock',
                'original_equal_time_detector_region_preserved', 'observed_postdiction'):
        assert control[key] is False


@pytest.mark.parametrize('mutation', ['baseline', 'noise', 'schedule', 'energy'])
def test_finite_instrument_control_rejects_false_precision_or_scope(tmp_path, mutation):
    packet = json.loads((ledger.CODE/'source_scalar_finite_instrument/finite_instrument_receipt.json').read_bytes())
    if mutation == 'baseline':
        packet['rows'][15]['noisy_baseline_probability_interval'] = ['1/2', '1/2']
    elif mutation == 'noise':
        packet['rows'][15]['paired_response_noise_error_upper'] = '0'
    elif mutation == 'schedule':
        packet['schedule'][0]['operations'].pop()
    else:
        packet['scope']['noise_or_total_hardware_energy_bound'] = True
    path = tmp_path/'finite-instrument.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        ledger._source_scalar_finite_instrument_control(path)


@pytest.mark.parametrize('missing', ['mathematical_replay', 'forced_oscillator_identities_verified'])
def test_finite_instrument_control_requires_operator_and_parent_proof(monkeypatch, missing):
    packet = json.loads((ledger.CODE/'source_scalar_finite_instrument/finite_instrument_receipt.json').read_bytes())
    summary = {'verified': True, **packet['summary'], 'mathematical_replay': True,
               'parent_events_replayed': 5888, 'physical_outcomes': False,
               'native_quantum_controls': False, 'forced_oscillator_identities_verified': True}
    summary[missing] = False
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: (packet, summary, {}))
    with pytest.raises(SystemExit, match='forced-oscillator proof and full fresh parent replay'):
        ledger._source_scalar_finite_instrument_control()


def test_source_operator_join_control_replays_with_cp1252_default(monkeypatch):
    expected = ledger._source_operator_join_control()
    original_read_text = Path.read_text

    def windows_read_text(path, encoding=None, errors=None, **kwargs):
        return original_read_text(path, encoding=encoding or "cp1252", errors=errors, **kwargs)

    monkeypatch.setattr(Path, "read_text", windows_read_text)
    with pytest.raises(UnicodeDecodeError):
        (ledger.REPO / "Lean/QFT/SourceOperatorGeneration.lean").read_text()
    assert ledger._source_operator_join_control() == expected


def test_source_operator_join_control_separates_replay_and_universal_proof():
    row = ledger._triple_observer_carrier_row()
    assert row['id'] == 'triple_observer_carrier_coupling'
    control = row['full_operator_join_control']
    proof = control['independent_verifier_result']
    assert proof['transition_carrier_points'] == 73346
    assert proof['chronological_intervals'] == 496
    assert proof['actual_transitions'] == 31
    assert proof['sparse_matrix_unit_witnesses'] == 18
    assert proof['lean_kernel_replay_performed'] is False
    assert proof['universal_operator_theorems_verified_by_python'] is False
    source = control['lean_source_evidence']
    assert source['kernel_replay_performed_by_ledger'] is False
    assert source['universal_proof_requires_separate_Lean_check'] is True
    assert source['sha256'] == hashlib.sha256((ledger.REPO/source['source']).read_bytes()).hexdigest()
    declarations = row['lean_declarations']['TripleCarrierOperatorJoin']
    for name in ('operatorJoin_intersection', 'operatorJoin_generate_top',
                 'operatorJoin88_counted_readout', 'operatorJoin247_counted_readout',
                 'operatorJoin88_evolve', 'operatorJoin247_evolve',
                 'operatorJoinAdmitted_generates', 'operatorJoinAdmitted_evolved'):
        assert name in declarations
    assert 'tripleCarrier_join_receipt' in row['lean_declarations']['TripleCarrierJoin']
    assert control['carrier']['common_complex_vector_dimension'] == 169
    assert control['carrier']['triple_dimension'] == 2366
    assert control['nonstationary_state_witness'] == {
        'transition': 0, 'triple_label': [1, 4, 6], 'denominator': 32,
        'before_count': 0, 'after_count': 1}
    for key in ('python_sparse_witnesses_prove_universal_generation',
                'partial_trace_is_multiplicative',
                'counted_state_identified_with_uniform_tower_state',
                'counted_state_invariant_under_transport', 'source_selected_quantum_controls',
                'tower_or_refinement_morphism', 'physical_time_slice_or_clock',
                'same_carrier_as_scalar64_instrument', 'observed_postdiction'):
        assert control[key] is False


@pytest.mark.parametrize('mutation', ['transition', 'wraparound', 'hinge', 'uniform_state',
                                    'singleton', 'physical_clock', 'scalar_carrier', 'parent'])
def test_source_operator_join_control_rejects_forged_packet(tmp_path, mutation):
    packet = json.loads((ledger.CODE/'source_operator_join/operator_join_packet.json').read_bytes())
    if mutation == 'transition':
        packet['transitions'][0]['permutations'][0] = list(range(13))
    elif mutation == 'wraparound':
        packet['transitions'].append({'row': 31, 'next_row': 0, 'permutations': []})
    elif mutation == 'hinge':
        packet['carrier']['common_complex_vector_dimension'] = 13
    elif mutation == 'uniform_state':
        packet['interpretation']['tower_state'] = 'same_as_empirical'
    elif mutation == 'singleton':
        packet['interpretation']['admitted_generation'] = 'singleton_path_projector'
    elif mutation == 'physical_clock':
        packet['interpretation']['physical_clock'] = True
    elif mutation == 'scalar_carrier':
        packet['interpretation']['scalar64_instrument_same_carrier'] = True
    else:
        first = next(iter(packet['parents']))
        packet['parents'][first] = '0'*64
    path = tmp_path/'operator-join.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(SystemExit, match='structural certificate rejected'):
        ledger._source_operator_join_control(path)


@pytest.mark.parametrize(('key', 'value'), [
    ('transition_carrier_points', 73345), ('chronological_intervals', 495),
    ('actual_transitions', 32), ('lean_kernel_replay_performed', True),
    ('universal_operator_theorems_verified_by_python', True),
    ('physical_adequacy', True), ('observed_quantum_outcomes', True),
    ('physical_adequacy', 0),
])
def test_source_operator_join_control_rejects_incomplete_or_promoted_summary(monkeypatch, key, value):
    summary = {
        'verdict': 'verified_source_operator_join_packet', 'source_rows': 32,
        'actual_transitions': 31, 'transition_carrier_points': 73346,
        'chronological_intervals': 496, 'sparse_matrix_unit_witnesses': 18,
        'lean_kernel_replay_performed': False,
        'universal_operator_theorems_verified_by_python': False,
        'physical_adequacy': False, 'observed_quantum_outcomes': False,
    }
    summary[key] = value
    monkeypatch.setattr(ledger, '_structural_packet', lambda *a, **k: ({}, summary, {}))
    with pytest.raises(SystemExit, match='Python/Lean proof boundary'):
        ledger._source_operator_join_control()


def test_source_operator_join_control_reads_current_verifier(monkeypatch):
    original = Path.read_bytes
    target = ledger.CODE/'source_operator_join/verify.py'
    def altered(path):
        data = original(path)
        if path == target:
            data += b'\ndef verify(packet):\n    raise ValueError("fresh operator join verifier sentinel")\n'
        return data
    monkeypatch.setattr(Path, 'read_bytes', altered)
    with pytest.raises(SystemExit, match='fresh operator join verifier sentinel'):
        ledger._source_operator_join_control()


def test_source_operator_join_control_rejects_changed_lean_with_retained_names(monkeypatch):
    original = Path.read_bytes
    target = ledger.LEAN_RECEIPTS['TripleCarrierOperatorJoin']
    def altered(path):
        data = original(path)
        if path == target:
            data += b'\ntheorem counterfeit_join : False := by sorry\n'
        return data
    monkeypatch.setattr(Path, 'read_bytes', altered)
    with pytest.raises(SystemExit, match='audited Lean source changed'):
        ledger._source_operator_join_control()


def test_source_operator_join_control_requires_all_universal_declarations(monkeypatch):
    original = Path.read_text
    target = ledger.LEAN_RECEIPTS['TripleCarrierOperatorJoin']
    def altered(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        if path == target:
            text = text.replace('theorem operatorJoin_intersection :',
                                'theorem omitted_operatorJoin_intersection :')
        return text
    monkeypatch.setattr(Path, 'read_text', altered)
    with pytest.raises(SystemExit, match='Lean declaration missing'):
        ledger._source_operator_join_control()


def test_source_operator_join_control_rechecks_lean_after_python_replay(monkeypatch):
    original_read = Path.read_bytes
    original_replay = ledger._structural_packet
    target = ledger.LEAN_RECEIPTS['TripleCarrierOperatorJoin']
    finished = False
    def altered(path):
        data = original_read(path)
        return data+b'\n' if finished and path == target else data
    def replay(*args, **kwargs):
        nonlocal finished
        result = original_replay(*args, **kwargs)
        finished = True
        return result
    monkeypatch.setattr(Path, 'read_bytes', altered)
    monkeypatch.setattr(ledger, '_structural_packet', replay)
    with pytest.raises(SystemExit, match='Lean source changed during replay'):
        ledger._source_operator_join_control()


def test_fermion_controls_keep_mean_and_abelian_quantum_constraints_distinct():
    control = ledger._fermion_source_current_control()
    assert control['mathematical_replay'] and not control['observed_postdiction']
    mean = control['controls']['mean_field']['independent_verifier_result']
    quantum = control['controls']['quantum_link']['independent_verifier_result']
    assert mean['operator_Gauss'] is False
    assert mean['initial_local_charge_variance'] == '945/512'
    assert quantum['operator_Gauss'] is True
    assert quantum['gauge_constraint_group'] == 'U1_hypercharge'
    assert quantum['continuous_point_current'] == '3*sqrt(3)'
    coupled = control['controls']['coupled']['independent_verifier_result']
    assert coupled['same_declared_hybrid_action'] is True
    assert coupled['electric_to_link_to_matter_feedback'] is True
    assert coupled['quadratic_Maxwell_kinetic'] is False
    assert coupled['operator_Gauss'] is False


@pytest.mark.parametrize('target,key,value', [
    ('mean_field', 'operator_Gauss', True),
    ('mean_field', 'exact_rational_replay', False),
    ('quantum_link', 'independent_quantum_state_replay', False),
    ('quantum_link', 'exact_continuous_quantum_edge_solution', False),
    ('quantum_link', 'gauge_constraint_group', 'full_Standard_Model'),
    ('quantum_link', 'physical_current_attached', True),
    ('coupled', 'electric_to_link_to_matter_feedback', False),
    ('coupled', 'Higgs_Yukawa_expected_first_variations_zero', False),
    ('coupled', 'quadratic_Maxwell_kinetic', True),
    ('coupled', 'global_Z6_quotient_selected', True),
])
def test_fermion_control_rejects_custody_only_and_physical_upgrades(monkeypatch, target, key, value):
    def changed(package, verifier, receipt, **kwargs):
        which = 'coupled' if 'coupled' in verifier else ('quantum_link' if 'quantum' in verifier else 'mean_field')
        summary = {'verified': True, 'physical_current_attached': False,
                   'operator_Gauss': which == 'quantum_link',
                   'exact_rational_replay': True, 'independent_quantum_state_replay': True,
                   'exact_continuous_quantum_edge_solution': True,
                   'gauge_constraint_group': 'U1_hypercharge',
                   'same_declared_hybrid_action': True, 'electric_to_link_to_matter_feedback': True,
                   'nonabelian_expected_first_variations_zero': True,
                   'Higgs_Yukawa_expected_first_variations_zero': True,
                   'electric_kinetic': 'fixed_lambda_convex_nonquadratic',
                   'quadratic_Maxwell_kinetic': False, 'global_Z6_quotient_selected': False}
        if which == target:
            summary[key] = value
        return {}, summary, {}
    monkeypatch.setattr(ledger, '_structural_packet', changed)
    with pytest.raises(SystemExit):
        ledger._fermion_source_current_control()


def test_pauli_rows_remain_conditional_and_unobserved(result):
    rows = {r['id']: r for r in result['sections']['forced_structure']}
    names = {'spin_exchange_selection_boundary', 'conditional_pauli_scalar_q_selection',
             'conditional_pauli_operation_selection'}
    assert names <= rows.keys()
    for name in names:
        row = rows[name]
        assert row['observed_postdiction'] is False
        assert row['source_selected_statistics'] is False
        assert row['physical_spin_statistics_theorem'] is False
        assert row['independent_verifier_result']['verified'] is True
        raw = (ledger.REPO/row['artifact_refs'][1]).read_bytes()
        assert row['receipt_pin'] == {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
    assert rows['conditional_pauli_operation_selection']['independent_verifier_result']['source_premises_derived'] is False


@pytest.mark.parametrize('package,field,value', [
    ('spin_exchange', 'physical_spin_statistics_theorem', True),
    ('pauli_stability', 'physical_spin_statistics_theorem', True),
    ('pauli_source_selection', 'physical_spin_statistics_theorem', True),
])
def test_pauli_ledger_rejects_false_physical_promotion(monkeypatch, package, field, value):
    original = ledger._structural_packet
    def forged(name, *args, **kwargs):
        packet, summary, pin = original(name, *args, **kwargs)
        if name == package:
            packet = deepcopy(packet)
            packet['scope'][field] = value
        return packet, summary, pin
    monkeypatch.setattr(ledger, '_structural_packet', forged)
    with pytest.raises(SystemExit, match='spin-statistics boundary'):
        ledger._pauli_structure_rows()


def test_pauli_ledger_rejects_false_source_derivation(monkeypatch):
    original = ledger._structural_packet
    def forged(name, *args, **kwargs):
        packet, summary, pin = original(name, *args, **kwargs)
        if name == 'pauli_source_selection':
            summary = {**summary, 'source_premises_derived': True}
        return packet, summary, pin
    monkeypatch.setattr(ledger, '_structural_packet', forged)
    with pytest.raises(SystemExit, match='field/Kraus identification'):
        ledger._pauli_structure_rows()


def test_retrospective_tilt_attaches_pinned_overlapping_likelihood():
    comparison = ledger._retrospective_tilt_reference()
    diagnostic = comparison["act_only_conditional_likelihood"]
    parent = json.loads((ledger.REPO / diagnostic["artifact_ref"]).read_text())
    assert diagnostic["delta_chi2"] == parent["delta_chi2"]
    assert diagnostic["new_independent_data"] is False
    assert diagnostic["FZ15_verdict"] is None
    assert diagnostic["full_cosmological_profile"] is False
    assert diagnostic["native_primordial_generator_identified"] is False


def test_retrospective_tilt_likelihood_rejects_changed_custody(tmp_path):
    path = ledger.CODE / "cosmology/tilt_posterior/likelihood_slice.json"
    parent = json.loads(path.read_text())
    parent["delta_chi2"] += 1
    changed = tmp_path / "changed-likelihood.json"
    changed.write_text(json.dumps(parent))
    with pytest.raises(SystemExit, match="custody mismatch"):
        ledger._retrospective_tilt_likelihood_reference(changed)


def test_shared_parameter_reference_preserves_interval_and_no_joint_verdict():
    row = ledger._retrospective_shared_parameter_reference()
    assert row["tail_weight_interval"] == [3114, 3114]
    assert row["multiplicity"] == 508985
    assert row["coverage"]["joint_likelihood"] is None
    assert row["coverage"]["full_common_physical_construction_established"] is False


def test_shared_parameter_renderer_retains_unequal_tail_endpoints():
    retained = json.loads(ledger.DEFAULT_OUT.read_text(encoding="utf-8"))
    shared = ledger._retrospective_shared_parameter_reference()
    shared["tail_weight_interval"] = [3114, 3115]
    shared["posterior_mass_interval"] = [str(Fraction(n, 508985)) for n in (3114, 3115)]
    retained["retrospective_tilt_comparison"]["shared_forward_parameter"] = shared
    rendered = ledger._render_md(retained)
    assert "posterior mass interval is `[3114,3115]/508985`" in rendered
    assert "both interval endpoints retain ACT tail weight" not in rendered


@pytest.mark.parametrize("mutation", ["claim", "P", "tilt"])
def test_shared_parameter_ledger_rejects_promotion_or_retuning(monkeypatch, mutation):
    original = ledger._public_receipt

    def changed(relative, schema):
        packet, digest = original(relative, schema)
        if schema == "oph.predictions.shared_forward_parameter.v1":
            if mutation == "claim":
                packet["coverage"]["full_common_physical_construction_established"] = True
            elif mutation == "P":
                packet["prediction"]["P_interval"] = ["1.6", "1.6"]
            else:
                packet["prediction"]["scalar_tilt_interval"] = ["0.974", "0.975"]
        return packet, digest

    monkeypatch.setattr(ledger, "_public_receipt", changed)
    with pytest.raises(SystemExit, match="shared-P"):
        ledger._retrospective_shared_parameter_reference()
