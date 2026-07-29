"""Fail-closed tests for the issue-506 fixed-verdict builder."""

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

PACKET_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKET_DIR))

import build_alpha_hvp_verdict as builder  # noqa: E402


class AlphaHvpVerdictTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = builder.build_verdict()
        cls.endpoint = json.loads(
            builder.ENDPOINT_PATH.read_text(encoding="utf-8")
        )
        cls.bridge = json.loads(
            builder.BRIDGE_PATH.read_text(encoding="utf-8")
        )

    def test_stored_verdict_matches_rebuild(self) -> None:
        stored = json.loads(builder.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, self.payload)

    def test_tabulated_class_confirmed_within_band(self) -> None:
        row = self.payload["class_matrix"]["tabulated_dispersive"]
        self.assertEqual(row["class_verdict"], "CONFIRMED_WITHIN_FROZEN_BAND")
        self.assertTrue(row["deficit_recomputation_agrees"])
        self.assertTrue(row["recorded_flag_agrees"])
        deficit = Decimal(row["recomputed_reference_deficit"])
        lower = Decimal(row["certified_gap_interval"][0])
        upper = Decimal(row["certified_gap_interval"][1])
        self.assertTrue(lower <= deficit <= upper)

    def test_unavailable_classes_are_not_evaluable_with_requirements(self) -> None:
        for name in ("raw_dispersive", "independent_code", "lattice_hvp"):
            row = self.payload["class_matrix"][name]
            self.assertEqual(
                row["class_verdict"], "NOT_EVALUABLE_MISSING_FROZEN_INGEST"
            )
            self.assertTrue(row["ingest_requirement"])
            self.assertTrue(row["fabrication_excluded"])

    def test_cross_class_not_evaluable_with_single_class(self) -> None:
        cross = self.payload["cross_class_agreement"]
        self.assertEqual(cross["evaluated_class_count"], 1)
        self.assertEqual(cross["verdict"], "NOT_EVALUABLE_SINGLE_CLASS")

    def test_band_mutation_flips_the_class_verdict(self) -> None:
        doctored = json.loads(json.dumps(self.bridge))
        doctored["certified_gap_from_endpoint"][
            "same_scheme_anchor_gap_interval"
        ] = [0.64, 0.65]
        row = builder.evaluate_tabulated_class(self.endpoint, doctored)
        self.assertFalse(row["deficit_inside_certified_gap"])
        self.assertEqual(row["class_verdict"], "REFUTED_OR_INCONSISTENT")

    def test_decomposition_tamper_is_detected(self) -> None:
        doctored = json.loads(json.dumps(self.bridge))
        doctored["reference_decomposition_compare_only"][
            "gap_phys_minus_oph"
        ] = 0.9
        row = builder.evaluate_tabulated_class(self.endpoint, doctored)
        self.assertFalse(row["deficit_recomputation_agrees"])
        self.assertEqual(row["class_verdict"], "REFUTED_OR_INCONSISTENT")

    def test_promotion_guards_are_refused(self) -> None:
        guards = self.payload["guards"]
        self.assertTrue(guards["payload_guards_verified"])
        self.assertFalse(guards["empirical_input_promoted_to_source_output"])
        self.assertFalse(guards["physical_alpha_prediction_emitted"])

    def test_kill_band_was_declared_before_evaluation(self) -> None:
        protocol = json.loads(
            builder.PROTOCOL_PATH.read_text(encoding="utf-8")
        )
        band = protocol["shared_protocol"]["kill_band"]
        self.assertTrue(band["declared_before_this_evaluation"])
        row = self.payload["class_matrix"]["tabulated_dispersive"]
        self.assertEqual(
            [str(Decimal(v)) for v in band["certified_gap_interval"]],
            row["certified_gap_interval"],
        )


if __name__ == "__main__":
    unittest.main()
