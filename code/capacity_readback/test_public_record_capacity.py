#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from public_record_capacity import (
    boundary_basis_trial,
    capacity_extension_naturality,
    evaluate_terminal,
    evaluate_terminal_fiber,
    example_receipt,
    greatest_fixed_point,
    refinement_stabilization,
)


HERE = Path(__file__).resolve().parent


class PublicRecordCapacityTests(unittest.TestCase):
    def test_boundary_basis_identity_gives_F_of_D_equal_D(self) -> None:
        for dimension in range(1, 8):
            terminal = boundary_basis_trial(
                dimension,
                checkpoint="identity",
                terminal_id=f"identity-{dimension}",
            )
            result = evaluate_terminal(dimension, terminal)
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["stable_public_code_size"], dimension)
            self.assertTrue(result["saturation_rank_one_complete_code"])

    def test_checkpoint_ablation_changes_capacity_not_normal_form(self) -> None:
        identity = boundary_basis_trial(
            4, checkpoint="identity", terminal_id="identity"
        )
        erasure = boundary_basis_trial(
            4, checkpoint="erase_to_first", terminal_id="erasure"
        )
        identity_result = evaluate_terminal(4, identity)
        erasure_result = evaluate_terminal(4, erasure)
        self.assertEqual(
            identity_result["normal_form_hash"], erasure_result["normal_form_hash"]
        )
        self.assertEqual(identity_result["stable_public_code_size"], 4)
        self.assertEqual(erasure_result["stable_public_code_size"], 1)
        self.assertFalse(erasure_result["nontrivial_record_capacity"])

    def test_checkpoint_must_descend_to_public_equalizer(self) -> None:
        terminal = boundary_basis_trial(
            3, checkpoint="identity", terminal_id="bad-descent"
        )
        terminal["observers"]["bob"]["atoms"][1]["checkpoint_next"] = "b2"
        result = evaluate_terminal(3, terminal)
        self.assertEqual(result["status"], "NO_PUBLIC_RECORD_DESCENT")

    def test_reachability_must_descend(self) -> None:
        terminal = boundary_basis_trial(
            2, checkpoint="identity", terminal_id="bad-reachability"
        )
        terminal["observers"]["bob"]["atoms"][1]["reachable"] = False
        result = evaluate_terminal(2, terminal)
        self.assertEqual(result["status"], "NO_PUBLIC_REACHABILITY_DESCENT")

    def test_boundary_supports_must_be_orthogonal_and_complete(self) -> None:
        terminal = boundary_basis_trial(
            2, checkpoint="identity", terminal_id="bad-boundary-support"
        )
        terminal["observers"]["alice"]["atoms"][1]["boundary_support"] = [0]
        terminal["observers"]["bob"]["atoms"][1]["boundary_support"] = [0]
        result = evaluate_terminal(2, terminal)
        self.assertEqual(result["status"], "NO_BOUNDARY_REPRESENTATION")

    def test_zero_error_public_reread_is_required(self) -> None:
        terminal = boundary_basis_trial(
            2, checkpoint="identity", terminal_id="bad-reread"
        )
        terminal["public_reread_zero_error"] = False
        result = evaluate_terminal(2, terminal)
        self.assertEqual(result["status"], "NO_PUBLIC_RECORD_CHANNEL")

    def test_fiber_is_set_valued_until_complete_agreement(self) -> None:
        identity = boundary_basis_trial(
            4, checkpoint="identity", terminal_id="identity"
        )
        erasure = boundary_basis_trial(
            4, checkpoint="erase_to_first", terminal_id="erasure"
        )
        ambiguous = evaluate_terminal_fiber(
            4, [identity, erasure], fiber_manifest_complete=True
        )
        incomplete = evaluate_terminal_fiber(
            4, [identity], fiber_manifest_complete=False
        )
        self.assertEqual(ambiguous["status"], "AMBIGUOUS_CAPACITY_READBACK")
        self.assertEqual(ambiguous["set_valued_readback"], [1, 4])
        self.assertEqual(ambiguous["count_kernel_row"], {"1": 1, "4": 1})
        self.assertEqual(incomplete["status"], "INCOMPLETE_TERMINAL_FIBER")

    def test_self_read_predicate_and_external_targets_fail_closed(self) -> None:
        cases = [
            ("self_read_predicate_injected", "CIRCULAR_CAPACITY_DEFINITION"),
            ("target_metadata_read_by_producer", "TARGET_TAINTED"),
            ("ew_bridge_target_used", "EW_BRIDGE_USED_AS_CAPACITY_PRODUCER"),
            ("rho_used_as_capacity_producer", "RHO_USED_AS_CAPACITY_PRODUCER"),
        ]
        for flag, expected in cases:
            with self.subTest(flag=flag):
                terminal = boundary_basis_trial(
                    2, checkpoint="identity", terminal_id=flag
                )
                terminal[flag] = True
                self.assertEqual(evaluate_terminal(2, terminal)["status"], expected)

    def test_hidden_inert_observer_duplication_preserves_capacity(self) -> None:
        terminal = boundary_basis_trial(
            3, checkpoint="identity", terminal_id="base"
        )
        duplicated = copy.deepcopy(terminal)
        duplicated["terminal_id"] = "duplicated"
        duplicated["observers"]["worker-copy"] = copy.deepcopy(
            duplicated["observers"]["alice"]
        )
        duplicated["interfaces"].append(
            {
                "left_observer": "alice",
                "right_observer": "worker-copy",
                "atom_pairs": [
                    {"left_atom": f"b{index}", "right_atom": f"b{index}"}
                    for index in range(3)
                ],
            }
        )
        self.assertEqual(
            evaluate_terminal(3, terminal)["stable_public_code_size"],
            evaluate_terminal(3, duplicated)["stable_public_code_size"],
        )

    def test_hidden_capacity_label_is_not_a_producer_input(self) -> None:
        left = boundary_basis_trial(3, checkpoint="identity", terminal_id="left")
        right = copy.deepcopy(left)
        right["terminal_id"] = "right"
        left["hidden_supplied_capacity_label"] = "label-a"
        right["hidden_supplied_capacity_label"] = "label-b"
        self.assertEqual(
            evaluate_terminal(3, left)["stable_public_code_size"],
            evaluate_terminal(3, right)["stable_public_code_size"],
        )

    def test_greatest_fixed_point_theorem_and_caveats(self) -> None:
        result = greatest_fixed_point({0: 0, 1: 1, 2: 1, 3: 2, 4: 3})
        identity = greatest_fixed_point({value: value for value in range(5)})
        constant_one = greatest_fixed_point(
            {0: 0, **{value: 1 for value in range(1, 5)}}
        )
        self.assertEqual(result["iteration_from_top"], [4, 3, 2, 1])
        self.assertEqual(result["greatest_fixed_point"], 1)
        self.assertEqual(identity["fixed_points"], [0, 1, 2, 3, 4])
        self.assertFalse(identity["fixed_point_unique"])
        self.assertEqual(constant_one["greatest_fixed_point"], 1)
        self.assertFalse(constant_one["positive_nontrivial_fixed_point"])

    def test_greatest_fixed_point_rejects_missing_hypotheses(self) -> None:
        nonmonotone = greatest_fixed_point({0: 0, 1: 1, 2: 2, 3: 1})
        inflationary = greatest_fixed_point({0: 0, 1: 2, 2: 2})
        self.assertEqual(nonmonotone["status"], "NOT_MONOTONE")
        self.assertEqual(inflationary["status"], "NOT_DEFLATIONARY")

    def test_capacity_extension_requires_embeddings_and_monotonicity(self) -> None:
        complete = capacity_extension_naturality(
            {1: 1, 2: 1, 3: 2}, {(1, 2): True, (2, 3): True}
        )
        missing = capacity_extension_naturality(
            {1: 1, 2: 1, 3: 2}, {(1, 2): True}
        )
        destroyed = capacity_extension_naturality(
            {1: 1, 2: 2, 3: 1}, {(1, 2): True, (2, 3): True}
        )
        self.assertEqual(complete["status"], "PASS")
        self.assertEqual(missing["status"], "CAPACITY_EXTENSION_NATURALITY_FAILED")
        self.assertEqual(
            destroyed["status"], "CAPACITY_EXTENSION_NATURALITY_FAILED"
        )

    def test_refinement_stabilization_is_exact_but_not_premature(self) -> None:
        open_tail = refinement_stabilization(
            5, [1, 2, 2], embeddings_certified=True
        )
        saturated = refinement_stabilization(
            5, [1, 3, 5], embeddings_certified=True
        )
        bad = refinement_stabilization(5, [2, 1], embeddings_certified=True)
        self.assertTrue(open_tail["eventual_stabilization_theorem_applies"])
        self.assertFalse(open_tail["permanent_stabilization_certified"])
        self.assertTrue(saturated["permanent_stabilization_certified"])
        self.assertEqual(bad["status"], "NON_NATURAL_REFINEMENT")

    def test_runtime_receipt_is_current_and_nonphysical(self) -> None:
        receipt = example_receipt()
        runtime = json.loads(
            (HERE / "runtime" / "public_record_capacity_example.json").read_text()
        )
        self.assertEqual(runtime, receipt)
        self.assertFalse(receipt["physical_content"])
        self.assertFalse(receipt["moves_cl7"])
        self.assertEqual(receipt["cl7_status"], "open")
        self.assertEqual(
            receipt["rho_op_role"], "independent estimator only; never an input to F"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
