from copy import deepcopy

from reversible_public_checkpoint_packet import (
    build_reference_packet,
    certify_reversible_packet,
    icosahedral_edges,
)


def test_icosahedral_reference_packet_closes_exact_reversible_branch():
    packet = build_reference_packet(capacity_dimension=4)
    certificate = certify_reversible_packet(packet)

    assert len(packet["observers"]) == 12
    assert len(icosahedral_edges()) == 30
    assert certificate["status"] == "PASS"
    assert certificate["reachable_public_record_count"] == 4
    assert certificate["exact_zero_error_capacity"] == 4
    assert certificate["robust_terminal_saturation"] is True
    assert certificate["rank_one_complete"] is True
    assert all(item["injective"] for item in certificate["checkpoint_generator_receipts"])


def test_noninjective_generator_fails_the_fast_branch():
    packet = build_reference_packet(capacity_dimension=3)
    records = sorted(packet["reachability_witnesses"])
    packet["global_checkpoint_kernels"][1]["rows"] = {
        record: {records[0]: 1.0} for record in records
    }

    certificate = certify_reversible_packet(packet)
    assert certificate["status"] == "NONINJECTIVE_CHECKPOINT_GENERATOR"


def test_target_taint_remains_fail_closed_in_reference_packet():
    packet = deepcopy(build_reference_packet())
    packet["lambda_used"] = True
    certificate = certify_reversible_packet(packet)
    assert certificate["status"] == "TARGET_TAINTED"
