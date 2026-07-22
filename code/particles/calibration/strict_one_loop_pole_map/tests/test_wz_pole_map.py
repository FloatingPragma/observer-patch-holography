#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "checker"))

from check_receipt import check as check_receipt

from wz_pole_map import (
    CONVENTION_ID,
    EXTERNAL_EVIDENCE_GATES,
    PoleMapError,
    SELF_ENERGY_TRANSLATION,
    build_receipt,
    det2,
    energy_pole_coordinates,
    neutral_inverse_matrix,
    tree_masses,
)


class PoleMapTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture_path = ROOT / "data" / "conditional_smdr_order1_fixture.json"
        self.fixture = json.loads(self.fixture_path.read_text(encoding="utf-8"))

    def test_tree_relations(self) -> None:
        p = self.fixture["renormalized_parameters"]
        got = tree_masses(p["g"], p["gp"], p["v_F_GeV"])
        self.assertAlmostEqual(got.mW0, 79.96948667898089, places=12)
        self.assertAlmostEqual(got.mZ0, 90.8636942619187, places=12)

    def test_conditional_fixture_regression(self) -> None:
        receipt = build_receipt(self.fixture, fixture_path=self.fixture_path)
        self.assertFalse(receipt["physical_promotion_allowed"])
        self.assertEqual(receipt["status"], "CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL")
        w = receipt["complex_poles"]["W"]["energy_pole_readout"]
        z = receipt["complex_poles"]["Z"]["energy_pole_readout"]
        self.assertAlmostEqual(w["M_GeV"], 80.37934572164735, places=12)
        self.assertAlmostEqual(w["Gamma_GeV"], 1.9971890954296638, places=12)
        self.assertAlmostEqual(z["M_GeV"], 90.68783666697097, places=12)
        self.assertAlmostEqual(z["Gamma_GeV"], 2.407078720027981, places=12)

    def test_strict_readout_is_not_sqrt_readout(self) -> None:
        receipt = build_receipt(self.fixture)
        strict = receipt["one_loop_coefficients"]["strict_W"]
        exact = receipt["complex_poles"]["W"]["energy_pole_readout"]
        self.assertGreater(abs(strict["M_strict_1L_GeV"] - exact["M_GeV"]), 1e-4)
        self.assertGreater(abs(strict["Gamma_strict_1L_GeV"] - exact["Gamma_GeV"]), 1e-4)

    def test_energy_pole_roundtrip(self) -> None:
        s = complex(6459.842027569383, -160.5327527730451)
        got = energy_pole_coordinates(s)
        reconstructed = complex(got["M_GeV"], -got["Gamma_GeV"] / 2.0) ** 2
        self.assertLess(abs(reconstructed - s), 2e-12)

    def test_neutral_matrix_keeps_mixing_but_root_excludes_it(self) -> None:
        z = 100.0
        matrix = neutral_inverse_matrix(
            complex(z, 0),
            z=z,
            delta_aa=complex(1, 2),
            delta_az=complex(3, 4),
            delta_za=complex(5, 6),
            delta_zz=complex(7, 8),
        )
        self.assertEqual(det2(matrix), matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0])
        fixture = json.loads(json.dumps(self.fixture))
        fixture["renormalized_one_loop_self_energy"]["Delta_AA_at_z_GeV2"] = [1.0, 2.0]
        fixture["renormalized_one_loop_self_energy"]["Delta_AZ_at_z_GeV2"] = [3.0, 4.0]
        fixture["renormalized_one_loop_self_energy"]["Delta_ZA_at_z_GeV2"] = [5.0, 6.0]
        receipt = build_receipt(fixture)
        self.assertFalse(receipt["neutral_matrix"]["diagnostics"]["mixing_used_in_strict_1L_root"])
        self.assertEqual(receipt["neutral_matrix"]["diagnostics"]["mixing_loop_power"], 2)

    def test_order_guard_rejects_order_2_5(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["perturbative_order"] = "2.5"
        with self.assertRaises(PoleMapError):
            build_receipt(fixture)

    def test_incomplete_mask_rejected(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["contribution_mask"]["one_loop_complete_for_declared_backend"] = False
        with self.assertRaises(PoleMapError):
            build_receipt(fixture)

    def test_sign_convention_is_frozen(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["inverse_propagator_convention"] = CONVENTION_ID + "_MUTATED"
        with self.assertRaises(PoleMapError):
            build_receipt(fixture)

    def test_sign_translation_is_explicit(self) -> None:
        receipt = build_receipt(self.fixture)
        self.assertEqual(receipt["self_energy_sign_translation"], SELF_ENERGY_TRANSLATION)

    def test_caller_evidence_booleans_cannot_self_promote(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["evidence_gates"] = {
            name: True
            for name in EXTERNAL_EVIDENCE_GATES
            if name != "full_neutral_matrix_supplied"
        }
        fixture["renormalized_one_loop_self_energy"].update(
            {
                "Delta_AA_at_z_GeV2": [1.0, -1.0],
                "Delta_AZ_at_z_GeV2": [2.0, -2.0],
                "Delta_ZA_at_z_GeV2": [2.0, -2.0],
            }
        )
        receipt = build_receipt(fixture)
        self.assertFalse(receipt["physical_promotion_allowed"])
        self.assertTrue(receipt["evidence_gates"]["neutral_matrix_point_diagnostic_supplied"])
        for gate in EXTERNAL_EVIDENCE_GATES:
            self.assertFalse(receipt["evidence_gates"][gate])
        self.assertTrue(all(receipt["unverified_evidence_claims"].values()))

    def test_checker_rejects_self_consistent_unrelated_receipt(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["renormalized_parameters"]["g"] = 0.4
        fixture["renormalized_one_loop_self_energy"]["Delta_WW_at_w_GeV2"] = [1.0, -2.0]
        forged = build_receipt(fixture, fixture_path=self.fixture_path)
        with tempfile.TemporaryDirectory() as directory:
            receipt_path = Path(directory) / "forged.json"
            receipt_path.write_text(json.dumps(forged), encoding="utf-8")
            with self.assertRaises((SystemExit, jsonschema.ValidationError)):
                check_receipt(
                    receipt_path,
                    self.fixture_path,
                    ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json",
                )

    def test_checker_rejects_tolerance_inflation(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["tolerances"]["absolute_GeV2"] = 1.0e9
        fixture["tolerances"]["absolute_GeV"] = 1.0e9
        fixture["tolerances"]["relative"] = 1.0
        with tempfile.TemporaryDirectory() as directory:
            fixture_path = Path(directory) / "inflated_fixture.json"
            fixture_path.write_text(json.dumps(fixture), encoding="utf-8")
            receipt = build_receipt(fixture, fixture_path=fixture_path)
            receipt_path = Path(directory) / "inflated_receipt.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises((SystemExit, jsonschema.ValidationError)):
                check_receipt(
                    receipt_path,
                    fixture_path,
                    ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json",
                )

    def test_schema_rejects_promotion_flip(self) -> None:
        receipt = build_receipt(self.fixture, fixture_path=self.fixture_path)
        receipt["physical_promotion_allowed"] = True
        with tempfile.TemporaryDirectory() as directory:
            receipt_path = Path(directory) / "promoted.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises(jsonschema.ValidationError):
                check_receipt(
                    receipt_path,
                    self.fixture_path,
                    ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json",
                )

    def test_checker_rejects_corrupted_derived_fields(self) -> None:
        mutations = {
            "tree_mW0": lambda receipt: receipt["tree"].__setitem__("mW0_GeV", 1.0),
            "strict_M0": lambda receipt: receipt["one_loop_coefficients"]["strict_W"].__setitem__("M0_GeV", 1.0),
            "strict_M": lambda receipt: receipt["one_loop_coefficients"]["strict_W"].__setitem__("M_strict_1L_GeV", 1.0),
            "sqrt_field": lambda receipt: receipt["complex_poles"]["W"]["energy_pole_readout"].__setitem__("sqrt_s_GeV", [1.0, -1.0]),
            "reconstructed_field": lambda receipt: receipt["complex_poles"]["W"]["energy_pole_readout"].__setitem__("s_reconstructed_GeV2", [1.0, -1.0]),
            "coordinate_note": lambda receipt: receipt["complex_poles"]["W"]["energy_pole_readout"].__setitem__("coordinate_transform_note", "mutated"),
            "pole_convention": lambda receipt: receipt.__setitem__("pole_convention", "mutated"),
        }
        schema = ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json"
        with tempfile.TemporaryDirectory() as directory:
            for name, mutate in mutations.items():
                with self.subTest(name=name):
                    receipt = build_receipt(self.fixture, fixture_path=self.fixture_path)
                    mutate(receipt)
                    receipt_path = Path(directory) / f"{name}.json"
                    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
                    with self.assertRaises((SystemExit, jsonschema.ValidationError)):
                        check_receipt(receipt_path, self.fixture_path, schema)

    def test_checker_recomputes_neutral_diagnostics(self) -> None:
        fixture = json.loads(json.dumps(self.fixture))
        fixture["renormalized_one_loop_self_energy"].update(
            {
                "Delta_AA_at_z_GeV2": [1.0, -1.0],
                "Delta_AZ_at_z_GeV2": [2.0, -2.0],
                "Delta_ZA_at_z_GeV2": [3.0, -3.0],
            }
        )
        schema = ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json"
        with tempfile.TemporaryDirectory() as directory:
            fixture_path = Path(directory) / "neutral_fixture.json"
            fixture_path.write_text(json.dumps(fixture), encoding="utf-8")
            receipt = build_receipt(fixture, fixture_path=fixture_path)
            receipt["neutral_matrix"]["diagnostics"]["det_Gamma_N_at_z_GeV4"] = [0.0, 0.0]
            receipt_path = Path(directory) / "corrupt_neutral.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises(SystemExit):
                check_receipt(receipt_path, fixture_path, schema)

    def test_checker_rejects_empty_fixture_substitution(self) -> None:
        schema = ROOT / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json"
        receipt = build_receipt(self.fixture, fixture_path=self.fixture_path)
        with tempfile.TemporaryDirectory() as directory:
            fixture_path = Path(directory) / "empty.json"
            fixture_path.write_text("{}", encoding="utf-8")
            receipt_path = Path(directory) / "receipt.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises(jsonschema.ValidationError):
                check_receipt(receipt_path, fixture_path, schema)


if __name__ == "__main__":
    unittest.main()
