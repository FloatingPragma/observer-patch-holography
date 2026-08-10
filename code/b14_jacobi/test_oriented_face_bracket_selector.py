#!/usr/bin/env python3
"""Regression and rehashed semantic-mutation tests for the B14 packet."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "oriented_face_bracket_selector.certificate.json"
PRODUCER = HERE / "oriented_face_bracket_selector.py"
VERIFIER_PATH = HERE / "verify_oriented_face_bracket_selector.py"


def load_verifier():
    spec = importlib.util.spec_from_file_location("oriented_face_verifier", VERIFIER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VERIFY = load_verifier()


class OrientedFaceBracketSelectorTests(unittest.TestCase):
    def test_committed_certificate_is_current_and_independently_replays(self) -> None:
        subprocess.run(
            ["python3", str(PRODUCER), "--check"], check=True,
            stdout=subprocess.PIPE, text=True,
        )
        result = VERIFY.verify_certificate(CERTIFICATE)
        self.assertEqual(result["identity"], "B_face=60*R13")
        self.assertEqual(result["jacobi_nonzero"], 240)
        self.assertEqual(result["unique_nearest"], "G")
        self.assertEqual(result["robust_norms"], ["L1", "L2", "Linfinity"])

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
        body["certificate_sha256"] = VERIFY.obj_sha(body)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.json"
            path.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
            with self.assertRaises(VERIFY.VerificationError):
                VERIFY.verify_certificate(path)

    def test_mutated_r13_coefficient_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["source_face_bracket"]["reynolds_coordinates"][0].__setitem__(1, 59)
        )

    def test_mutated_jacobi_count_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["jacobi_failure"].__setitem__("positive_count", 119)
        )

    def test_mutated_distance_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["orthogonal_compact_locus_discriminator"]["families"]["G"]
                .__setitem__("squared_distance", [28, 1, -123, 22])
        )

    def test_false_source_derivation_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["orthogonal_compact_locus_discriminator"]
                .__setitem__("minimum_hs_or_jacobi_repair_is_source_derived", True)
        )

    def test_mutated_nearest_label_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["orthogonal_compact_locus_discriminator"]
                .__setitem__("unique_nearest_family", "F")
        )

    def test_mutated_l1_distance_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["endpoint_norm_robustness"]["l1"]["families"]["G"]
                .__setitem__("distance", [-29, 1, 30, 1])
        )

    def test_mutated_linf_dual_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["endpoint_norm_robustness"]["linfinity"]["families"]["G"]
                ["dual_entries"][0].__setitem__(3, 2)
        )

    def test_false_three_norm_source_selection_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["endpoint_norm_robustness"]
                .__setitem__("repair_norm_or_minimization_rule_is_source_derived", True)
        )


if __name__ == "__main__":
    unittest.main()
