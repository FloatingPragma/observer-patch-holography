"""Fail-closed tests for the bounded issue-506 accounting verdict."""

from __future__ import annotations

import copy
import json
import shutil
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

PACKET_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKET_DIR))

import build_alpha_hvp_verdict as builder  # noqa: E402


class AlphaHvpVerdictTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol = builder.load_protocol()
        cls.snapshot_blobs = builder.verify_snapshot_commit(cls.protocol)
        cls.pinned = builder.verify_exact_pins(cls.protocol)
        cls.endpoint = builder._load_pinned_json(cls.pinned, builder.ENDPOINT_REL)
        cls.bridge = builder._load_pinned_json(cls.pinned, builder.BRIDGE_REL)
        cls.payload_source = builder._load_pinned_json(cls.pinned, builder.PAYLOAD_REL)
        cls.measure = builder._load_pinned_json(cls.pinned, builder.MEASURE_REL)
        cls.payload = builder.build_verdict()

    def evaluate(
        self,
        *,
        endpoint: dict | None = None,
        bridge: dict | None = None,
        payload: dict | None = None,
        measure: dict | None = None,
        protocol: dict | None = None,
    ) -> dict:
        return builder.evaluate_tabulated_accounting_replay(
            endpoint or copy.deepcopy(self.endpoint),
            bridge or copy.deepcopy(self.bridge),
            payload or copy.deepcopy(self.payload_source),
            measure or copy.deepcopy(self.measure),
            protocol or copy.deepcopy(self.protocol),
        )

    def test_stored_verdict_matches_rebuild(self) -> None:
        stored = json.loads(builder.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, self.payload)

    def test_all_declared_source_bytes_are_exactly_pinned(self) -> None:
        self.assertEqual(set(self.pinned), builder.REQUIRED_PIN_PATHS)
        self.assertEqual(set(self.snapshot_blobs), builder.REQUIRED_PIN_PATHS)
        self.assertTrue(self.payload["source_snapshot"]["commit_ancestor_verified"])
        self.assertTrue(self.payload["source_snapshot"]["commit_blob_pins_verified"])
        self.assertTrue(self.payload["source_snapshot"]["all_exact_pins_verified"])

    def test_overall_verdict_is_multi_class_not_evaluable(self) -> None:
        self.assertEqual(
            self.payload["verdict"],
            ("MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"),
        )
        cross = self.payload["cross_class_agreement"]
        self.assertEqual(cross["recorded_accounting_replay_count"], 1)
        self.assertEqual(cross["independently_evaluated_class_count"], 0)
        self.assertEqual(cross["verdict"], "NOT_EVALUABLE_NO_INDEPENDENT_CLASS")

    def test_scope_refuses_independence_freeze_and_prediction(self) -> None:
        scope = self.payload["scope"]
        self.assertEqual(scope["comparison_timing"], "retrospective")
        self.assertFalse(scope["prospective_freeze"])
        self.assertFalse(scope["independent_hvp_implementation_supplied"])
        self.assertFalse(scope["empirical_input_promoted_to_source_output"])
        self.assertFalse(scope["physical_alpha_prediction_emitted"])
        gate = self.payload["serialized_acceptance_gates"][
            "independent_multi_class_evaluation"
        ]
        self.assertEqual(gate["status"], "NOT_ATTAINED")
        self.assertFalse(gate["passed"])

    def test_recorded_accounting_replay_is_compatible(self) -> None:
        row = self.payload["class_matrix"]["tabulated_dispersive"]
        self.assertEqual(
            row["class_verdict"],
            "COMPATIBLE_RECORDED_ACCOUNTING_REPLAY",
        )
        self.assertFalse(row["independent_hvp_implementation"])
        self.assertFalse(row["prospective_comparison"])
        self.assertTrue(row["all_internal_gates_pass"])
        self.assertTrue(all(gate["passed"] for gate in row["gates"].values()))

    def test_unavailable_classes_are_not_evaluable_with_requirements(self) -> None:
        for name in builder.UNAVAILABLE_CLASSES:
            row = self.payload["class_matrix"][name]
            self.assertEqual(
                row["class_verdict"],
                "NOT_EVALUABLE_MISSING_FROZEN_INGEST",
            )
            self.assertFalse(row["independent_hvp_implementation"])
            self.assertTrue(row["ingest_requirement"])
            self.assertTrue(row["fabrication_excluded"])

    def test_protocol_cannot_claim_a_prospective_freeze(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["scope"]["prospective_freeze"] = True
        with self.assertRaisesRegex(
            builder.VerdictInputError, "cannot claim a prospective freeze"
        ):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "protocol.json"
                path.write_text(json.dumps(protocol), encoding="utf-8")
                builder.load_protocol(path)

    def test_protocol_pin_mutation_fails_closed(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["preexisting_source_snapshot"]["sha256"][builder.ENDPOINT_REL] = (
            "sha256:" + "0" * 64
        )
        with self.assertRaisesRegex(builder.VerdictInputError, "pinned source drift"):
            builder.verify_exact_pins(protocol)

    def test_snapshot_commit_must_resolve(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["preexisting_source_snapshot"]["commit"] = "0" * 40
        with self.assertRaisesRegex(
            builder.VerdictInputError, "commit does not resolve"
        ):
            builder.verify_snapshot_commit(protocol)

    def test_snapshot_blob_pin_mutation_fails_closed(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["preexisting_source_snapshot"]["sha256"][builder.ENDPOINT_REL] = (
            "sha256:" + "0" * 64
        )
        with self.assertRaisesRegex(builder.VerdictInputError, "snapshot blob drift"):
            builder.verify_snapshot_commit(protocol)

    def test_pinned_source_byte_mutation_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for relative_path in builder.REQUIRED_PIN_PATHS:
                source = builder.REPO_ROOT / relative_path
                target = root / relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
            endpoint = root / builder.ENDPOINT_REL
            endpoint.write_bytes(endpoint.read_bytes() + b" ")
            with self.assertRaisesRegex(
                builder.VerdictInputError, "pinned source drift"
            ):
                builder.verify_exact_pins(self.protocol, root)

    def test_endpoint_primitive_mutation_is_detected(self) -> None:
        endpoint = copy.deepcopy(self.endpoint)
        endpoint["inputs"]["source_anchor_inv_alpha_MZ"] = "0"
        row = self.evaluate(endpoint=endpoint)
        self.assertEqual(row["class_verdict"], "INTERNAL_INCONSISTENCY")
        self.assertFalse(row["gates"]["endpoint_recomputed_exactly"]["passed"])
        self.assertFalse(row["all_internal_gates_pass"])

    def test_payload_primitive_mutation_is_detected(self) -> None:
        payload = copy.deepcopy(self.payload_source)
        payload["integral"]["value"] = 0.9
        row = self.evaluate(payload=payload)
        self.assertEqual(row["class_verdict"], "INTERNAL_INCONSISTENCY")
        self.assertFalse(row["gates"]["payload_value_matches"]["passed"])

    def test_reference_primitive_mutation_is_detected(self) -> None:
        bridge = copy.deepcopy(self.bridge)
        bridge["reference_decomposition_compare_only"]["Delta_had5"] = 0.9
        row = self.evaluate(bridge=bridge)
        self.assertEqual(row["class_verdict"], "INTERNAL_INCONSISTENCY")
        self.assertFalse(row["gates"]["reference_decomposition_recomputed"]["passed"])

    def test_protocol_gap_mutation_is_detected(self) -> None:
        protocol = copy.deepcopy(self.protocol)
        protocol["accounting_rule"]["same_scheme_gap_interval"] = [
            "0.1",
            "0.2",
        ]
        row = self.evaluate(protocol=protocol)
        self.assertEqual(row["class_verdict"], "INTERNAL_INCONSISTENCY")
        self.assertFalse(row["gates"]["protocol_gap_matches_recomputed_gap"]["passed"])

    def test_valid_outside_band_is_refutation_not_protocol_failure(self) -> None:
        bridge = copy.deepcopy(self.bridge)
        reference = bridge["reference_decomposition_compare_only"]
        reference["Delta_had5"] = 0.04
        alpha_mz = Decimal(str(reference["alpha_inv_0"])) * (
            Decimal(1)
            - Decimal(str(reference["Delta_lep"]))
            - Decimal(str(reference["Delta_had5"]))
            - Decimal(str(reference["Delta_top"]))
        )
        anchor = Decimal(str(bridge["anchor_provenance"]["a0_oph"]))
        reference["alpha_inv_mz_phys_on_shell"] = str(alpha_mz)
        reference["gap_phys_minus_oph"] = str(alpha_mz - anchor)
        bridge["verdict"]["reference_deficit_inside_certified_gap"] = False

        row = self.evaluate(bridge=bridge)
        self.assertTrue(row["all_internal_gates_pass"])
        self.assertFalse(row["reference_deficit_inside_gap"])
        self.assertEqual(row["class_verdict"], "REFUTED_RECORDED_ACCOUNTING_REPLAY")
        self.assertEqual(
            builder._overall_verdict(row["class_verdict"]),
            ("MULTI_CLASS_NOT_EVALUABLE__RECORDED_ACCOUNTING_REPLAY_REFUTED"),
        )

    def test_target_guard_mutation_is_detected(self) -> None:
        endpoint = copy.deepcopy(self.endpoint)
        endpoint["guards"]["measured_alpha_in_solve_path"] = True
        row = self.evaluate(endpoint=endpoint)
        self.assertEqual(row["class_verdict"], "INTERNAL_INCONSISTENCY")
        self.assertFalse(row["gates"]["target_and_promotion_guards"]["passed"])

    def test_verdict_digest_covers_the_serialized_payload(self) -> None:
        expected = builder.sha256_bytes(
            builder.canonical_json(
                {
                    key: value
                    for key, value in self.payload.items()
                    if key != "verdict_sha256"
                }
            ).encode("utf-8")
        )
        self.assertEqual(self.payload["verdict_sha256"], expected)


if __name__ == "__main__":
    unittest.main()
