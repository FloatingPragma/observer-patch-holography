from copy import deepcopy

import pytest

from correctable_public_record_capacity import (
    approximate_public_capacity,
    certify_unique_slack_zero,
    compound_confusability_graph,
    evaluate_terminal,
    evaluate_terminal_fiber,
    greatest_fixed_point,
    no_new_confusability,
    support_relation_semigroup,
    tv_robustness_bound,
)


def binary_packet(terminal_id="q0"):
    records = ["alice=0|bob=0", "alice=1|bob=1"]
    return {
        "terminal_id": terminal_id,
        "capacity_dimension": 2,
        "observers": {"alice": ["0", "1"], "bob": ["0", "1"]},
        "interfaces": [
            {
                "left_observer": "alice",
                "right_observer": "bob",
                "left_readout": {"0": "z", "1": "o"},
                "right_readout": {"0": "z", "1": "o"},
            }
        ],
        "reachability_witnesses": {record: ["write", record] for record in records},
        "publicness_policy": [["alice", "bob"]],
        "continuation_manifest_complete": True,
        "global_checkpoint_kernels": [
            {
                "authorized_observers": ["alice", "bob"],
                "continuation_id": "identity",
                "rows": {
                    records[0]: {"00": 1.0},
                    records[1]: {"11": 1.0},
                },
            }
        ],
        "local_marginal_consistency_passed": True,
        "projection_supports": {records[0]: [0], records[1]: [1]},
    }


def test_exact_packet_saturates_carrier_and_scalarizes_terminal_fiber():
    packet = binary_packet()
    result = evaluate_terminal(packet)
    assert result["status"] == "PASS"
    assert result["exact_zero_error_capacity"] == 2
    assert result["saturation_rank_one_complete"] is True

    fiber = evaluate_terminal_fiber([packet, binary_packet("q1")], manifest_complete=True)
    assert fiber["status"] == "PASS"
    assert fiber["terminal_fiber_capacity_set"] == [2]
    assert fiber["robust_closure"] is True


def test_cyclic_permutation_preserves_capacity_despite_trivial_fixed_algebra():
    packet = binary_packet()
    records = sorted(packet["reachability_witnesses"])
    packet["global_checkpoint_kernels"][0] = {
        "authorized_observers": ["alice", "bob"],
        "continuation_id": "swap",
        "rows": {records[0]: {records[1]: 1.0}, records[1]: {records[0]: 1.0}},
    }
    result = evaluate_terminal(packet)
    assert result["exact_zero_error_capacity"] == 2
    assert all(not neighbors for neighbors in result["confusability_graph"].values())


def test_joint_coupling_changes_capacity_and_cannot_be_replaced_by_marginals():
    packet = binary_packet()
    records = sorted(packet["reachability_witnesses"])
    packet["global_checkpoint_kernels"][0]["rows"] = {
        records[0]: {"00": 0.5, "11": 0.5},
        records[1]: {"00": 0.5, "11": 0.5},
    }
    assert evaluate_terminal(packet)["exact_zero_error_capacity"] == 1

    packet["continuation_manifest_complete"] = False
    assert evaluate_terminal(packet)["status"] == "NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING"


def test_terminal_fiber_ambiguity_is_not_silently_scalarized():
    saturated = binary_packet("saturated")
    erased = deepcopy(saturated)
    erased["terminal_id"] = "erased"
    records = sorted(erased["reachability_witnesses"])
    erased["global_checkpoint_kernels"][0]["rows"] = {
        records[0]: {"same": 1.0},
        records[1]: {"same": 1.0},
    }
    result = evaluate_terminal_fiber([saturated, erased], manifest_complete=True)
    assert result["status"] == "AMBIGUOUS_CAPACITY_READBACK"
    assert result["terminal_fiber_capacity_set"] == [1, 2]


def test_approximate_worst_input_capacity_and_tv_bound():
    records = ["0", "1"]
    noisy = [
        {
            "rows": {
                "0": {"0": 0.9, "1": 0.1},
                "1": {"0": 0.1, "1": 0.9},
            }
        }
    ]
    assert approximate_public_capacity(records, noisy, 0.09)["capacity"] == 1
    assert approximate_public_capacity(records, noisy, 0.10)["capacity"] == 2
    assert tv_robustness_bound(0.10, 0.03) == pytest.approx(0.13)


def test_semigroup_order_and_finite_size_certificates():
    records = ["0", "1"]
    swap = {"rows": {"0": {"1": 1.0}, "1": {"0": 1.0}}}
    relations = support_relation_semigroup(records, [swap])
    assert len(relations) == 2

    graph = compound_confusability_graph(records, [swap])
    assert no_new_confusability(graph, graph, {"0": "0", "1": "1"}) is True

    order = greatest_fixed_point({1: 1, 2: 1, 3: 3, 4: 3})
    assert order["path"] == [4, 3]
    assert order["greatest_fixed_point"] == 3
    assert certify_unique_slack_zero({1: 1, 2: 1, 3: 2}, 1)["status"] == "PASS"
    assert certify_unique_slack_zero({1: 1, 2: 2}, 1)["status"] == "FINITE_SIZE_SELECTOR_FAILED"


def test_target_taint_and_incomplete_carrier_fail_closed():
    packet = binary_packet()
    packet["ew_bridge_used"] = True
    assert evaluate_terminal(packet)["status"] == "TARGET_TAINTED"

    packet = binary_packet()
    packet["projection_supports"].pop(next(iter(packet["projection_supports"])))
    assert evaluate_terminal(packet)["status"] == "NO_CAPACITY_CARRIER_REPRESENTATION"
