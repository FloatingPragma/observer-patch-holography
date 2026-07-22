from __future__ import annotations

import json
from pathlib import Path

import pytest

from correctable_public_record_capacity import evaluate_terminal, public_global_sections
from source_derived_public_checkpoint_packet import (
    PORTS,
    build_source_derived_bundle,
    build_source_derived_packet,
    build_terminal_fiber_manifest,
    certify_source_derived_packet,
    classify_terminal_fiber,
    continuation_actions,
    icosahedral_edges,
    oriented_slots,
    verify_local_marginal_consistency,
    write_runtime_receipts,
)


@pytest.fixture(scope="module")
def packet():
    return build_source_derived_packet()


@pytest.fixture(scope="module")
def manifest():
    return build_terminal_fiber_manifest()


@pytest.fixture(scope="module")
def certificate(packet, manifest):
    return certify_source_derived_packet(packet, terminal_manifest=manifest)


def test_structural_trial_fiber_is_complete_output_blind_singleton(manifest):
    assert manifest["terminal_fiber_complete"] is True
    assert manifest["trial_count"] == 67
    assert manifest["terminal_fiber_kind"] == "SINGLETON"
    assert manifest["terminal_world_ids"] == ["q_echosahedral_oriented_exact"]
    assert manifest["terminal_membership_output_blind"] is True
    assert all(manifest["terminal_membership_source_audit"].values())
    assert sum(trial["terminal"] for trial in manifest["trials"]) == 1
    completeness = manifest["terminal_fiber_completeness_certificate"]
    assert completeness["declared_trials_enumerated"] == completeness["expected_trial_count"] == 67
    assert completeness["missing_trial_ids"] == []
    assert completeness["duplicate_trial_ids"] == []
    assert completeness["all_trial_ids_unique"] is True
    assert completeness["all_candidates_materialized"] is True
    assert completeness["complete"] is True
    assert manifest["terminal_fiber_manifest_hash_scope"].endswith(
        "except terminal_fiber_manifest_sha256"
    )
    assert len(manifest["terminal_fiber_manifest_sha256"]) == 64
    assert all(len(trial["candidate_sha256"]) == 64 for trial in manifest["trials"])


def test_source_record_diagram_has_total_atoms_and_exact_csp_count(packet):
    assert len(packet["observers"]) == 12
    assert len(packet["interfaces"]) == 30
    assert len(icosahedral_edges()) == 30
    assert len(oriented_slots()) == 24
    assert all(len(atoms) == 24 for atoms in packet["observers"].values())
    for interface in packet["interfaces"]:
        assert set(interface["left_readout"]) == set(packet["observers"][interface["left_observer"]])
        assert set(interface["right_readout"]) == set(packet["observers"][interface["right_observer"]])
        assert len(interface["interface_atoms"]) == 24
    sections = public_global_sections(packet["observers"], packet["interfaces"])
    assert len(sections) == 24
    assert len(packet["reachability_witnesses"]) == 24
    assert all(packet["reachability_witnesses"].values())
    assert packet["publicness_policy"] == [list(PORTS)]
    assert len(packet["public_global_sections"]) == 24
    assert packet["packet_hash_scope"].endswith("except packet_sha256")
    assert len(packet["packet_sha256"]) == 64


def test_complete_reversible_continuation_family_and_marginals(packet, certificate):
    assert len(continuation_actions()) == 40
    assert packet["continuation_family_order"] == 40
    assert len(packet["global_checkpoint_kernels"]) == 40
    assert certificate["checkpoint_family"]["support_relation_semigroup_size"] == 40
    assert certificate["checkpoint_family"]["composition"]["composition_entries_checked"] == 1600
    assert verify_local_marginal_consistency(packet)["status"] == "PASS"
    receipts = certificate["checkpoint_family"]["injective_generator_receipts"]
    assert len(receipts) == 40
    assert all(receipt["injective"] for receipt in receipts.values())
    assert all(receipt["worst_input_success"] == 1.0 for receipt in receipts.values())


def test_exact_compound_capacity_decoder_and_carrier_saturation(packet, certificate):
    evaluation = evaluate_terminal(packet)
    assert evaluation["status"] == "PASS"
    assert evaluation["public_global_section_count"] == 24
    assert evaluation["exact_zero_error_capacity"] == 24
    assert all(not neighbours for neighbours in evaluation["confusability_graph"].values())
    assert len(evaluation["independent_set_witness"]) == 24
    assert evaluation["dimension_bound_passed"] is True
    assert evaluation["saturation_rank_one_complete"] is True
    assert certificate["capacity"]["fast_branch_identity"] == "M_0(q)=|X_reach(q)|"
    assert certificate["packet_hash_valid"] is True
    scalarization = certificate["capacity"]["whole_fiber_scalarization_result"]
    assert scalarization["fiber_kind"] == "SINGLETON"
    assert certificate["capacity"]["robust_closure_at_frozen_D"] is True
    count_kernel = certificate["capacity"]["count_kernel_row"]
    assert count_kernel["status"] == "PASS"
    assert count_kernel["exact_reachable_section_count"] == 24
    assert certificate["carrier"]["orthogonal_rank_one_projections"] is True
    assert certificate["carrier"]["rank_sum"] == 24
    assert certificate["carrier"]["rank_one_saturation"] is True


def test_approximate_worst_input_and_tv_robustness(certificate):
    approximate = certificate["approximate_worst_input_branch"]
    assert approximate["exact_family_worst_input_success"] == 1.0
    assert approximate["exact_M_epsilon_at_epsilon_0"] == 24
    assert all(receipt["certified_M_epsilon"] == 24 for receipt in approximate["tv_robustness"])
    noisy = certificate["controls"]["full_support_noise"]
    assert noisy["status"] == "PASS"
    assert noisy["all_rows_full_support"] is True
    assert noisy["exact_zero_error_capacity"] == 1
    assert noisy["approximate_capacity_at_epsilon"] == 24
    assert noisy["tv_identity"] is True
    assert noisy["mixture_weight"] == pytest.approx(1.0e-6)
    assert noisy["epsilon_for_full_code"] < 1.0e-6


def test_all_required_adversarial_controls_pass(certificate):
    controls = certificate["controls"]
    isomorphic = controls["isomorphic_packet"]
    assert isomorphic["status"] == "PASS"
    assert isomorphic["section_bijection_size"] == 24
    assert isomorphic["original_capacity"] == isomorphic["renamed_capacity"] == 24
    assert controls["cyclic_permutation"]["status"] == "PASS"
    alternative = controls["alternative_joint_coupling"]
    assert alternative["same_local_marginals"] is True
    assert alternative["joint_capacity_parity"] == 2
    assert alternative["joint_capacity_independent"] == 1
    assert controls["circular_definition"]["status"] == "PASS"
    assert controls["circular_definition"]["evaluator_status"] == "CIRCULAR_CAPACITY_DEFINITION"
    assert controls["target_taint"]["evaluator_status"] == "TARGET_TAINTED"
    assert set(controls["target_taint"]["statuses_by_taint"].values()) == {"TARGET_TAINTED"}
    order = controls["order_families"]
    assert order["identity_fixed_every_dimension"] is True
    assert order["erasure_fixed_only_singleton"] is True
    fibers = controls["terminal_fibers"]
    assert fibers["empty"]["fiber_kind"] == "EMPTY"
    assert fibers["incomplete"]["fiber_kind"] == "INCOMPLETE"
    assert fibers["ambiguous"]["fiber_kind"] == "AMBIGUOUS"
    assert fibers["singleton"]["fiber_kind"] == "SINGLETON"
    suffix = controls["finite_suffix_nonpromotion"]
    assert suffix["status"] == "PASS"
    assert suffix["equal_suffix"] == [24, 24, 24]
    assert suffix["saturation_or_exhaustion_receipt"] is False
    assert suffix["promoted_to_limit"] is False


def test_separate_extension_and_fixed_D_refinement_injections(certificate):
    receipts = certificate["no_new_confusability_injections"]
    assert receipts["capacity_extension"]["from_dimension"] == 24
    assert receipts["capacity_extension"]["to_dimension"] == 48
    assert receipts["capacity_extension"]["no_new_confusability"] is True
    assert receipts["capacity_extension"]["new_confusability_control_rejected"] is True
    assert receipts["fixed_D_refinement"]["dimension"] == 24
    assert receipts["fixed_D_refinement"]["no_new_confusability"] is True
    assert receipts["fixed_D_refinement"]["new_confusability_control_rejected"] is True


def test_complete_bundle_and_canonical_runtime_receipts(tmp_path: Path, packet, certificate):
    bundle = build_source_derived_bundle()
    assert bundle["certificate"]["status"] == "PASS"
    assert len(bundle["terminal_packets"]) == 1
    paths = write_runtime_receipts(tmp_path)
    emitted_packet = json.loads(Path(paths["packet"]).read_text(encoding="utf-8"))
    emitted_certificate = json.loads(Path(paths["certificate"]).read_text(encoding="utf-8"))
    emitted_manifest = json.loads(Path(paths["terminal_fiber_manifest"]).read_text(encoding="utf-8"))
    assert emitted_packet["packet_sha256"] == packet["packet_sha256"]
    assert emitted_certificate["status"] == certificate["status"] == "PASS"
    assert emitted_manifest["terminal_world_ids"] == ["q_echosahedral_oriented_exact"]


def test_explicit_fiber_classifier(packet):
    assert classify_terminal_fiber([], manifest_complete=True)["fiber_kind"] == "EMPTY"
    assert classify_terminal_fiber([packet], manifest_complete=False)["fiber_kind"] == "INCOMPLETE"
    assert classify_terminal_fiber([packet], manifest_complete=True)["fiber_kind"] == "SINGLETON"
    assert classify_terminal_fiber([packet, packet], manifest_complete=True)["fiber_kind"] == "COMPLETE_SCALAR"

