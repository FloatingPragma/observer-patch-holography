from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("collar_gap_certificate", HERE / "collar_gap_certificate.py")
assert SPEC and SPEC.loader
CERTIFICATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CERTIFICATE)


def manifest(name: str) -> dict:
    return json.loads((HERE / "manifests" / name).read_text())


class CollarGapCertificateTests(unittest.TestCase):
    def test_atomic_calibration_has_the_announced_exact_floor(self) -> None:
        receipt = CERTIFICATE.validate(manifest("atomic_4d_ising_calibration.json"))
        self.assertEqual(receipt["active_type_count"], 244)
        self.assertEqual(receipt["c_floor"], "1")
        self.assertEqual(receipt["eta_upper"], "1/2")
        self.assertEqual(receipt["gap_lower"], "1/2")
        self.assertFalse(receipt["physical_clay_receipt"])

    def test_binary_float_input_is_rejected(self) -> None:
        payload = manifest("atomic_4d_ising_calibration.json")
        payload["calibration_type_family"]["template"]["rate_lower"] = 1.0
        with self.assertRaisesRegex(ValueError, "floating"):
            CERTIFICATE.validate(payload)

    def test_nonpositive_rate_is_rejected(self) -> None:
        payload = manifest("atomic_4d_ising_calibration.json")
        payload["calibration_type_family"]["template"]["rate_lower"] = "0"
        with self.assertRaisesRegex(ValueError, "positive"):
            CERTIFICATE.validate(payload)

    def test_noncontractive_influence_is_rejected(self) -> None:
        payload = manifest("atomic_4d_ising_calibration.json")
        payload["calibration_type_family"]["template"]["influences"][0]["upper"] = "1"
        with self.assertRaisesRegex(ValueError, "< 1"):
            CERTIFICATE.validate(payload)

    def test_explicit_production_table_mode_is_equivalent(self) -> None:
        payload = manifest("atomic_4d_ising_calibration.json")
        payload["type_table"] = CERTIFICATE.expand_types(payload)
        del payload["calibration_type_family"]
        self.assertEqual(CERTIFICATE.validate(payload)["active_type_count"], 244)

    def test_tampered_receipt_is_rejected(self) -> None:
        payload = manifest("atomic_4d_ising_calibration.json")
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            manifest_path = directory / "manifest.json"
            receipt_path = directory / "receipt.json"
            manifest_path.write_text(json.dumps(payload))
            CERTIFICATE.certify(manifest_path, receipt_path)
            receipt = json.loads(receipt_path.read_text())
            receipt["gap_lower"] = "9/10"
            receipt_path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "exactly recompute"):
                CERTIFICATE.verify(manifest_path, receipt_path)

    def test_physical_placeholder_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "type_table"):
            CERTIFICATE.validate(manifest("physical_compact_gauge_uninstantiated.json"))
