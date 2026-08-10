#!/usr/bin/env python3
"""Regression and rehashed semantic-mutation tests for the B14 invariant-metric
phase packet."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "invariant_metric_phase.certificate.json"
PRODUCER = HERE / "invariant_metric_phase.py"
VERIFIER_PATH = HERE / "verify_invariant_metric_phase.py"


def load_verifier():
    spec = importlib.util.spec_from_file_location("invariant_metric_verifier", VERIFIER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VERIFY = load_verifier()


class InvariantMetricPhaseTests(unittest.TestCase):
    def test_committed_certificate_is_current_and_independently_replays(self) -> None:
        subprocess.run(
            ["python3", str(PRODUCER), "--check"], check=True,
            stdout=subprocess.PIPE, text=True,
        )
        summary = VERIFY.verify_certificate(CERTIFICATE)
        self.assertEqual(summary["commutant_dimension"], 4)
        self.assertEqual(summary["balanced_slice_selects"], "G")
        self.assertTrue(summary["P_excluded_everywhere"])
        self.assertEqual(summary["samples_replayed"], 6)

    def test_fresh_production_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "certificate.json"
            subprocess.run(
                ["python3", str(PRODUCER), "--output", str(output)], check=True,
                stdout=subprocess.PIPE, text=True,
            )
            self.assertEqual(output.read_bytes(), CERTIFICATE.read_bytes())

    def assert_semantic_mutation_rejected(self, mutate) -> None:
        body = json.loads(CERTIFICATE.read_text())
        mutate(body)
        body.pop("certificate_sha256", None)
        body["certificate_sha256"] = VERIFY.object_hash(body)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.json"
            path.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
            with self.assertRaises(VERIFY.VerificationError):
                VERIFY.verify_certificate(path)

    def test_mutated_closed_form_coefficient_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["closed_forms"]["d_G_squared"]
                .__setitem__("inv_gamma", [61, 11, -12, 11])
        )

    def test_mutated_sample_value_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["samples"][1]["families"]
                .__setitem__("F", body["samples"][1]["families"]["G"])
        )

    def test_mutated_phase_constant_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["phase_diagram"]["constants"]
                .__setitem__("C_plus", [105, 11, -45, 11])
        )

    def test_mutated_box_endpoint_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["phase_diagram"]["phase_box"]
                .__setitem__("parabola_endpoint_high", [1, 1, 0, 1])
        )

    def test_mutated_commutant_dimension_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["completeness"]
                .__setitem__("commutant_dimension_over_Q", 5)
        )

    def test_mutated_channel_norm_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["channel_norm_monomials"]["t_pf_to_f"]
                .__setitem__("norm", [1200, 1, 240, 1])
        )

    def test_mutated_channel_monomial_label_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["channel_norm_monomials"]["t_pf_to_f"]
                .__setitem__("monomial", "inv_gamma")
        )

    def test_mutated_face_channel_coordinate_is_rejected(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["face_channel_coordinates"]
                .__setitem__("t_mf_to_f", [0, 1, 1, 20])
        )

    def test_mutated_witness_ordering_is_rejected(self) -> None:
        def swap(body):
            witness = body["phase_diagram"]["F_region_witness"]
            witness["d_F_squared"], witness["d_G_squared"] = (
                witness["d_G_squared"], witness["d_F_squared"])
        self.assert_semantic_mutation_rejected(swap)

    def test_mutated_control_values_are_rejected(self) -> None:
        def swap(body):
            control = body["non_induced_control"]
            control["d_G_control"], control["d_F_control"] = (
                control["d_F_control"], control["d_G_control"])
        self.assert_semantic_mutation_rejected(swap)


if __name__ == "__main__":
    unittest.main()
