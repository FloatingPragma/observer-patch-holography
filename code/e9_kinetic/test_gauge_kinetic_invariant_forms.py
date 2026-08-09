#!/usr/bin/env python3
"""Regression and semantic-mutation tests for the exact E9 1a packet."""

from __future__ import annotations

import copy
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "gauge_kinetic_invariant_forms.certificate.json"
PRODUCER = HERE / "gauge_kinetic_invariant_forms.py"
VERIFIER = HERE / "verify_gauge_kinetic_invariant_forms.py"


def load_verifier():
    spec = importlib.util.spec_from_file_location("gauge_kinetic_verifier", VERIFIER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


VERIFY = load_verifier()


class GaugeKineticInvariantFormsTests(unittest.TestCase):
    def test_committed_certificate_is_current_and_replays(self) -> None:
        subprocess.run(
            ["python3", str(PRODUCER), "--check"], check=True, cwd=HERE.parents[1],
            stdout=subprocess.PIPE, text=True,
        )
        result = VERIFY.verify_certificate(CERTIFICATE)
        self.assertEqual(result["dimensions"], {"F": 2, "G": 2, "P": 2})

    def test_fresh_production_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "certificate.json"
            subprocess.run(
                ["python3", str(PRODUCER), "--output", str(output)], check=True,
                cwd=HERE.parents[1], stdout=subprocess.PIPE, text=True,
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

    def test_mutated_f_relation_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["families"]["F"]["rref_constraints"][0][2].__setitem__(2, 1)
        )

    def test_mutated_plane_control_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["families"]["P"].__setitem__("ad_invariant_dimension", 1)
        )

    def test_weakened_mirror_boundary_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["mirror_common_control"].__setitem__("not_source_selected", False)
        )

    def test_mutated_upstream_pin_is_rejected_after_rehash(self) -> None:
        self.assert_semantic_mutation_rejected(
            lambda body: body["upstream"].__setitem__("stage2_system_sha256", "sha256:" + "0" * 64)
        )


if __name__ == "__main__":
    unittest.main()
