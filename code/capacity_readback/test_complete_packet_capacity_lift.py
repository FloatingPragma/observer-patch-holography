from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from capacity_indexed_source_family import expected_capacity
from complete_packet_capacity_lift import (
    RECEIPT_PATH,
    SAMPLE_RUNGS,
    a2_naturality_receipt,
    a3_feasible_receipt,
    branch_completion_kernel,
    build_receipt,
    canonical_json_bytes,
    capacity_row,
    erasure_kernel,
    extension_receipt,
    generated_lifted_group,
    has_source_ancestry,
    joint_kernel_consistency_receipt,
    lifted_records,
    manifest_closure_receipt,
    mutation_controls,
    public_sections_receipt,
    reachability_receipt,
    sewing_receipt,
    tagged_sha256,
    terminal_fiber_receipt,
    is_terminal_world_at_rung,
    lifted_world,
)

HERE = Path(__file__).resolve().parent


def _copy_independent_verifier_sandbox(tmp_path: Path) -> tuple[Path, Path]:
    sandbox = tmp_path / "code" / "capacity_readback"
    runtime = sandbox / "runtime"
    runtime.mkdir(parents=True)
    (sandbox / "verify_complete_packet_lift_independent.py").write_bytes(
        (HERE / "verify_complete_packet_lift_independent.py").read_bytes()
    )
    for name in (
        "complete_packet_capacity_lift_receipt.json",
        "complete_packet_capacity_lift_certificate.json",
        "source_derived_public_checkpoint_packet.json",
    ):
        (runtime / name).write_bytes((HERE / "runtime" / name).read_bytes())
    for name in ("capacity_indexed_source_family.py", "F_READBACK_SPEC.md"):
        (sandbox / name).write_bytes((HERE / name).read_bytes())
    return sandbox, runtime


def _rehash_receipt(receipt: dict) -> None:
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))


def test_terminal_fiber_unique_world_and_no_self_read() -> None:
    for k in (1, 2, 5):
        receipt = terminal_fiber_receipt(k)
        assert receipt["unique_terminal_world"] is True
        assert receipt["planted_capacity_metadata_ignored"] is True
        assert receipt["trial_count"] >= 67


def test_terminal_gate_rejects_wrong_generation_count() -> None:
    world = lifted_world(3)
    assert is_terminal_world_at_rung(world, 3) is True
    world["generation_count"] = 4
    assert is_terminal_world_at_rung(world, 3) is False


def test_sections_histories_and_sewing_track_the_register() -> None:
    for k in (1, 3):
        sections = public_sections_receipt(k)
        assert sections["equalizer_section_count"] == 24 * k
        assert sections["connected_components"] == 1
        assert sections["one_section_per_record"] is True
        assert sections["readout_injective"] is True
        histories = reachability_receipt(k)
        assert histories["all_sections_have_histories"] is True
        sewing = sewing_receipt(k)
        assert sewing["regions_connected"] is True
        assert sewing["glued_section_count"] == 24 * k
        assert sewing["fiber_product_count"] == 24 * k
        assert sewing["fiber_product_matches"] is True
    assert (
        sewing_receipt(2, corrupt_seam_readout=True)["fiber_product_matches"] is False
    )


def test_publicness_refinement_and_pinned_reconciliation() -> None:
    from complete_packet_capacity_lift import (
        publicness_receipt,
        refinement_receipt,
        joint_kernel_consistency_receipt,
    )

    publicness = publicness_receipt("reversible_identity", 2)
    assert publicness["publicness_frozen"] is True
    assert publicness["family_matches_pinned_packet"] is True
    assert publicness["kernels_authorized_by_collective_set"] is True

    refinement = refinement_receipt("copy_collapse_erasure", 3)
    assert refinement["capacity_stable_along_refinement"] is True
    assert refinement["capacity_before"] == refinement["capacity_after"] == 24
    broken = refinement_receipt("reversible_identity", 2, corrupt_transport=True)
    assert broken["capacity_stable_along_refinement"] is False

    kernels = joint_kernel_consistency_receipt("reversible_identity", 2)
    reconciliation = kernels["pinned_packet_reconciliation"]
    assert reconciliation["pinned_packet_reconciled"] is True
    assert reconciliation["compared_rows"] == 40 * 12 * 24


def test_lifted_group_order_is_forty_at_every_sampled_rung() -> None:
    for k in (1, 2, 3):
        assert len(generated_lifted_group(k)) == 40


def test_collapse_kernels_lack_source_ancestry() -> None:
    group = generated_lifted_group(2)
    assert has_source_ancestry(erasure_kernel(2), group) is False
    identity_map = {record: record for record in lifted_records(2)}
    assert has_source_ancestry(identity_map, group) is True


def test_manifest_closure_sizes() -> None:
    assert manifest_closure_receipt("reversible_identity", 2)["closure_size"] == 40
    assert manifest_closure_receipt("copy_collapse_erasure", 2)["closure_size"] == 80
    assert manifest_closure_receipt("capped_two_class", 3)["closure_size"] == 80


def test_joint_kernels_marginal_consistency_and_parity_control() -> None:
    receipt = joint_kernel_consistency_receipt("copy_collapse_erasure", 2)
    assert receipt["local_marginals_consistent"] is True
    assert receipt["local_packet_tamper_detected"] is True
    control = receipt["parity_independent_channel_control"]
    assert control["single_observer_marginals_agree"] is True
    assert control["parity_channel_capacity"] == 2
    assert control["independent_channel_capacity"] == 1
    assert control["local_marginals_do_not_determine_joint"] is True


def test_a2_naturality_survivors_and_oscillation_failure() -> None:
    for branch in ("reversible_identity", "copy_collapse_erasure", "capped_two_class"):
        assert a2_naturality_receipt(branch, 2)["a2_natural"] is True
    oscillation = a2_naturality_receipt("parity_oscillation", 3)
    assert oscillation["a2_natural"] is False
    assert oscillation["extension_square_failure_count"] > 0


def test_a3_state_determinacy_spectator_failure() -> None:
    fine = a3_feasible_receipt("reversible_identity", 2)
    assert fine["cover_state_determining"] is True
    spectator = a3_feasible_receipt("hidden_spectator", 1, spectator_multiplicity=2)
    assert spectator["cover_state_determining"] is False
    assert spectator["raw_family_multiplicity_per_public_class"] == 4096


def test_extension_no_new_confusability_for_survivors() -> None:
    for branch in ("reversible_identity", "copy_collapse_erasure", "capped_two_class"):
        for k in (1, 2, 3):
            assert extension_receipt(branch, k)["no_new_confusability"] is True
    assert extension_receipt("parity_oscillation", 3)["no_new_confusability"] is False


def test_capacities_match_bounded_family_formulas() -> None:
    for branch in ("reversible_identity", "copy_collapse_erasure", "capped_two_class"):
        for k in SAMPLE_RUNGS:
            row = capacity_row(branch, k)
            assert row["zero_error_capacity"] == expected_capacity(
                branch, k, spectator_multiplicity=1
            )


def test_receipt_verdict_structure() -> None:
    receipt = build_receipt()
    assert receipt["scientific_verdict"] == (
        "BOUNDED_GENERATION_REGISTER_COUNTERMODEL__UNIVERSAL_MEMBERSHIP_OPEN"
    )
    sampled = receipt["sampled_control_assessment"]
    assert sampled["branches_passing_all_sampled_transported_controls"] == [
        "reversible_identity",
        "copy_collapse_erasure",
        "capped_two_class",
    ]
    assert set(sampled["branches_failing_sampled_transported_controls"]) == {
        "hidden_spectator",
        "parity_oscillation",
    }
    assert sampled["does_not_certify_universal_membership"] is True
    arithmetic = receipt["all_rung_capacity_arithmetic"]
    assert arithmetic["identity_slack_zero_set"] == "every positive rung"
    assert arithmetic["identity_has_unique_positive_slack_zero"] is False
    assert arithmetic["inequivalent_formula_zero_sets"] is True
    source_status = receipt["source_contract_status"]
    assert source_status["universal_all_rung_membership_proved"] is False
    assert source_status["executable_lean_membership_bridge_proved"] is False
    assert (
        source_status["complete_a1_a3_source_class_nonidentifiability_proved"] is False
    )
    assert source_status["direct_n_status"] == (
        "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
    )
    assert receipt["mutation_controls"]["all_mutations_detected"] is True
    assert receipt["bounded_family_cross_check"]["consistent"] is True


def test_committed_receipt_is_byte_stable() -> None:
    committed = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
    rebuilt = build_receipt()
    assert canonical_json_bytes(committed) == canonical_json_bytes(rebuilt)


def test_independent_verifier_passes() -> None:
    result = subprocess.run(
        [sys.executable, str(HERE / "verify_complete_packet_lift_independent.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "BOUNDED_PACKET_LIFT_INDEPENDENT_VALID" in result.stdout


def test_independent_verifier_rejects_tampered_receipt(tmp_path: Path) -> None:
    sandbox, runtime = _copy_independent_verifier_sandbox(tmp_path)
    receipt = json.loads(
        (runtime / "complete_packet_capacity_lift_receipt.json").read_text(
            encoding="utf-8"
        )
    )
    receipt["all_rung_capacity_arithmetic"]["branch_zero_sets"][1][
        "sampled_zero_rungs"
    ] = [1, 2]
    _rehash_receipt(receipt)
    (runtime / "complete_packet_capacity_lift_receipt.json").write_text(
        json.dumps(receipt), encoding="utf-8"
    )
    result = subprocess.run(
        [sys.executable, str(sandbox / "verify_complete_packet_lift_independent.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "ZERO_SET_DRIFT" in result.stdout + result.stderr


def test_independent_verifier_rejects_membership_promotion(tmp_path: Path) -> None:
    sandbox, runtime = _copy_independent_verifier_sandbox(tmp_path)
    receipt_path = runtime / "complete_packet_capacity_lift_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["source_contract_status"]["universal_all_rung_membership_proved"] = True
    _rehash_receipt(receipt)
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(sandbox / "verify_complete_packet_lift_independent.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "UNIVERSAL_MEMBERSHIP_PROMOTED" in result.stdout + result.stderr


def test_independent_verifier_rejects_control_table_promotion(tmp_path: Path) -> None:
    sandbox, runtime = _copy_independent_verifier_sandbox(tmp_path)
    receipt_path = runtime / "complete_packet_capacity_lift_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    parity = receipt["sampled_control_assessment"]["control_table"][-1]
    parity["per_rung_controls"]["3"]["a2_natural"] = True
    _rehash_receipt(receipt)
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(sandbox / "verify_complete_packet_lift_independent.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "CONTROL_VALUE" in result.stdout + result.stderr


def test_independent_verifier_rejects_upstream_pin_drift(tmp_path: Path) -> None:
    sandbox, runtime = _copy_independent_verifier_sandbox(tmp_path)
    receipt_path = runtime / "complete_packet_capacity_lift_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["upstream_pins"]["readback_spec_sha256"] = "sha256:" + "0" * 64
    _rehash_receipt(receipt)
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(sandbox / "verify_complete_packet_lift_independent.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "UPSTREAM_PIN" in result.stdout + result.stderr


def test_mutation_controls_all_detected() -> None:
    controls = mutation_controls()
    assert controls["all_mutations_detected"] is True


def test_branch_kernels_reject_unknown_branch() -> None:
    with pytest.raises(ValueError):
        branch_completion_kernel("unknown_branch", 2)
