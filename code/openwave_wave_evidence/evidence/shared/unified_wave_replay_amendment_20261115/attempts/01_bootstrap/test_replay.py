#!/usr/bin/env python3
"""Regression probes for custody, numerical comparison and readout changes."""
import argparse
import copy
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import replay as r

POLICY = r.strict_json((r.HERE / "POLICY.json").read_bytes())
RUN = None


class ComparatorTests(unittest.TestCase):
    def rejects(self, a, b):
        with self.assertRaises(r.ReplayError):
            r.compare_json(a, b, POLICY)

    def test_roundoff_is_allowed(self):
        self.assertEqual(r.compare_json([0.0, 1.0], [1e-13, 1.0 + 1e-13], POLICY)["changed_floats"], 2)

    def test_changed_float_fails(self):
        self.rejects({"value": 0.8059537648223112}, {"value": 0.156520121200273})
        self.rejects([1.0], [1.0 + 1e-10])

    def test_discrete_records_and_types_are_exact(self):
        for a, b in [(1, 2), (1, 1.0), (True, 1), ("pass", "fail"), (None, False)]:
            with self.subTest(a=a, b=b):
                self.rejects([a], [b])

    def test_structure_is_exact(self):
        for a, b in [({"a": 1}, {"b": 1}), ({"a": 1}, {"a": 1, "b": 2}), ([1], [1, 2])]:
            with self.subTest(a=a, b=b):
                self.rejects(a, b)

    def test_nonfinite_values_fail_in_either_operand(self):
        for value in [float("inf"), float("-inf"), float("nan")]:
            self.rejects(value, 0.0)
            self.rejects(0.0, value)
            self.rejects(value, value)

    def test_invalid_json_fails(self):
        for value in ['{"a":1,"a":2}', '[NaN]', '[Infinity]', '[-Infinity]', '[1e999]']:
            with self.subTest(value=value), self.assertRaises(r.ReplayError):
                r.strict_json(value)

    def test_path_escape_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            for value in ["../outside", "/absolute", "C:/absolute", "..\\outside", "."]:
                with self.subTest(value=value), self.assertRaises(r.ReplayError):
                    r.checked_path(Path(directory), value)


class SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="oph_amendment_tests_")
        cls.root = Path(cls.temp.name)
        cls.bundle = r.extract_git(cls.root / "source", POLICY)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_exact_blob_custody(self):
        result = r.custody(self.bundle, POLICY)
        self.assertEqual(result["outer_manifest_entries"], 420)
        self.assertEqual(result["copied_pins"], 409)

    def test_crlf_and_corruption_fail(self):
        path = self.bundle / "evidence/boost/CELL.md"
        original = path.read_bytes()
        try:
            for changed in [original.replace(b"\n", b"\r\n"), original + b"corrupted\n"]:
                path.write_bytes(changed)
                with self.assertRaises(r.ReplayError):
                    r.custody(self.bundle, POLICY)
        finally:
            path.write_bytes(original)

    def test_origin_mutation_fails(self):
        path = self.bundle / "ORIGINS.json"
        original = path.read_bytes()
        try:
            origins = r.strict_json(original)
            origins["free"]["manifest_sha256"] = "0" * 64
            r.save(path, origins)
            with self.assertRaises(r.ReplayError):
                r.custody(self.bundle, POLICY)
            completed = subprocess.run([sys.executable, "-B", str(r.HERE / "replay.py"),
                "--custody-only", "--bundle", str(self.bundle)], capture_output=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn(b"Manifest mismatch", completed.stderr)
        finally:
            path.write_bytes(original)

    def test_optimized_entry_rejected(self):
        for flags, optimize in [(["-O"], None), (["-OO"], None), ([], "1"), ([], "2")]:
            env = os.environ.copy()
            env.pop("PYTHONOPTIMIZE", None)
            if optimize:
                env["PYTHONOPTIMIZE"] = optimize
            completed = subprocess.run([sys.executable, "-B", *flags, str(r.HERE / "replay.py"),
                                        "--custody-only"], capture_output=True, env=env)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn(b"Optimized Python is unsupported", completed.stderr)

    def test_child_assertion_survives_inherited_optimization(self):
        script = self.root / "failing_child.py"
        script.write_text("assert False, 'deliberate scientific failure'\n", encoding="utf-8")
        with patch.dict(os.environ, {"PYTHONOPTIMIZE": "2"}):
            with self.assertRaises(r.ReplayError):
                r.run_child(script, self.root / "failing_child_output")
        receipt = r.strict_json((self.root / "failing_child_output/PROCESS.json").read_bytes())
        self.assertEqual(receipt["environment"]["optimize"], 0)
        self.assertNotEqual(receipt["exit_code"], 0)
        self.assertIn("deliberate scientific failure", (self.root / "failing_child_output/stderr.txt").read_text())

    def test_patch_is_narrow_and_cannot_reapply(self):
        source = (self.bundle / r.ORBITAL / "probe.py").read_bytes()
        path = self.root / "probe.py"
        path.write_bytes(source)
        r.patch_orbital(path)
        changed = path.read_text(encoding="utf-8")
        self.assertNotIn("np.argmax(np.diag(proj))", changed)
        for before, after in reversed(r.PATCHES):
            changed = changed.replace(after, before)
        self.assertEqual(changed.encode("utf-8"), source)
        with self.assertRaises(r.ReplayError):
            r.patch_orbital(path)

    def test_historical_orbital_reference_is_rejected(self):
        original = self.bundle / r.ORBITAL
        with self.assertRaisesRegex(r.ReplayError, "Float mismatch"):
            r.compare_artifacts("orbital", original, original, POLICY)

    def test_canonical_reference_has_no_sign_ambiguity(self):
        source = self.bundle / r.ORBITAL
        receipt = r.strict_json((source / "RECEIPT.json").read_bytes())
        records = r.strict_json(gzip.decompress((source / "RECORDS.json.gz").read_bytes()))
        a, _ = r.amend_orbital_reference(receipt, records)
        for population in receipt["population"]:
            for band in population["bands"]:
                band["ground"]["profile"] = [-x for x in band["ground"]["profile"]]
        b, _ = r.amend_orbital_reference(receipt, records)
        self.assertEqual(a, b)


class GeneratedArtifactTests(unittest.TestCase):
    def setUp(self):
        if RUN is None:
            self.skipTest("Pass --run DIRECTORY for generated-artifact probes")
        self.temp = tempfile.TemporaryDirectory(prefix="oph_amendment_artifact_test_")
        self.target = Path(self.temp.name) / "free"
        shutil.copytree(RUN / "free", self.target)

    def tearDown(self):
        if hasattr(self, "temp"):
            self.temp.cleanup()

    def compare(self):
        return r.compare_artifacts("free", RUN / "free", self.target, POLICY)

    def test_recompressed_archive_passes(self):
        path = self.target / "RECORDS.json.gz"
        raw = gzip.decompress(path.read_bytes())
        path.write_bytes(gzip.compress(raw, compresslevel=1, mtime=123))
        self.assertFalse(self.compare()["artifacts"][path.name]["byte_equal"])

    def test_wrong_embedded_hash_fails(self):
        path = self.target / "RECEIPT.json"
        receipt = r.strict_json(path.read_bytes())
        receipt["raw_uncompressed_sha256"] = "0" * 64
        r.save(path, receipt)
        with self.assertRaisesRegex(r.ReplayError, "embedded raw hash"):
            self.compare()

    def test_changed_record_fails_even_with_updated_hash(self):
        path = self.target / "RECORDS.json.gz"
        data = r.strict_json(gzip.decompress(path.read_bytes()))
        # Branch records are integers; update the digest to test the actual record gate.
        data["branches"][0][0] += 1
        raw = (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()
        path.write_bytes(gzip.compress(raw, mtime=0))
        receipt_path = self.target / "RECEIPT.json"
        receipt = r.strict_json(receipt_path.read_bytes())
        receipt["raw_uncompressed_sha256"] = r.sha_bytes(raw)
        r.save(receipt_path, receipt)
        with self.assertRaisesRegex(r.ReplayError, "Discrete value mismatch"):
            self.compare()

    def test_receipt_mutations_fail(self):
        path = self.target / "RECEIPT.json"
        original = r.strict_json(path.read_bytes())
        for field, value in [("actual_native_calls", 0), ("all_registered_checks", "fail"),
                             ("max_covariance_relative_error", 0.1), ("extra_field", True)]:
            receipt = copy.deepcopy(original)
            receipt[field] = value
            r.save(path, receipt)
            with self.subTest(field=field), self.assertRaises(r.ReplayError):
                self.compare()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path)
    args, remaining = parser.parse_known_args()
    RUN = args.run.resolve() if args.run else None
    unittest.main(argv=[sys.argv[0], *remaining])
