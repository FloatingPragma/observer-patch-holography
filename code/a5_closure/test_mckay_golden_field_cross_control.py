"""Fail-closed tests for the frozen OPH raw-output cross-control."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

SCRIPT_DIR = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "mckay_cross_verifier", SCRIPT_DIR / "verify_mckay_golden_field_cross_control.py")
cross = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cross)

class McKayCrossControlTests(unittest.TestCase):
    def copied_packet(self, temp_root: Path) -> Path:
        for rel in (
            "code/a5_closure/receipts/arithmon_mckay_golden_field_cross_control.json",
            "code/a5_closure/receipts/arithmon_mckay_golden_field_certificate.json",
            "code/a5_closure/receipts/arithmon_mckay_golden_field_reference.json",
            "code/a5_closure/receipts/external/oph_sl2f5_mckay_e8_4ae2148a.json",
        ):
            src = SCRIPT_DIR.parents[1] / rel
            dst = temp_root / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
        return temp_root / "code/a5_closure/receipts/arithmon_mckay_golden_field_cross_control.json"

    def run_in(self, root: Path):
        with patch.object(cross, "ROOT", root), patch.object(
            cross, "FIXTURE_PATH", root / "code/a5_closure/receipts/arithmon_mckay_golden_field_cross_control.json"
        ):
            return cross.verify()

    def test_01_frozen_raw_opH_receipt_agrees(self):
        result = cross.verify()
        self.assertEqual(result["comparison_result"], "INDEPENDENT_EXACT_AGREEMENT")
        self.assertEqual(result["oph_receipt_sha256"], cross.PINNED_OPH_RECEIPT_SHA256)
        self.assertTrue(all(x["match"] for x in result["compared_invariants"].values()))

    def test_02_fixture_does_not_copy_opH_invariant_values(self):
        fixture = json.loads(cross.FIXTURE_PATH.read_text())
        self.assertNotIn("compared_invariants", fixture)
        self.assertIn("raw_receipt_path", fixture["oph"])
        self.assertIn("raw_receipt_sha256", fixture["oph"])

    def test_03_tampered_raw_opH_and_rehashed_fixture_fails_pin(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture_path = self.copied_packet(root)
            raw_path = root / "code/a5_closure/receipts/external/oph_sl2f5_mckay_e8_4ae2148a.json"
            raw = json.loads(raw_path.read_text())
            raw["source"]["order"] = 121
            raw_path.write_text(json.dumps(raw, sort_keys=True) + "\n")
            fixture = json.loads(fixture_path.read_text())
            fixture["oph"]["raw_receipt_sha256"] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
            fixture_path.write_text(json.dumps(fixture, indent=2) + "\n")
            with self.assertRaises(ValueError), patch.object(cross, "ROOT", root), patch.object(
                cross, "FIXTURE_PATH", fixture_path
            ):
                cross.verify()

    def test_04_tampered_arithmon_receipt_and_rehashed_fixture_fails_pin(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture_path = self.copied_packet(root)
            producer_path = root / "code/a5_closure/receipts/arithmon_mckay_golden_field_certificate.json"
            producer_path.write_text(producer_path.read_text().replace('"order": 120', '"order": 121', 1))
            fixture = json.loads(fixture_path.read_text())
            data = producer_path.read_bytes()
            fixture["arithmon_receipts"]["producer"]["sha256"] = hashlib.sha256(data).hexdigest()
            fixture["arithmon_receipts"]["producer"]["git_blob_sha"] = cross.git_blob_sha(data)
            fixture_path.write_text(json.dumps(fixture, indent=2) + "\n")
            with self.assertRaises(ValueError), patch.object(cross, "ROOT", root), patch.object(
                cross, "FIXTURE_PATH", fixture_path
            ):
                cross.verify()

if __name__ == "__main__":
    unittest.main()
