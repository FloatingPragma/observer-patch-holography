#!/usr/bin/env python3
"""Regression and adversarial tests for the #616/#617 multiplicity certificate."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import multiplicity_window_certificate as cert  # noqa: E402


class MultiplicityWindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "multiplicity_window_reference.json"
        cls.stored = cert.load_json(cls.manifest_path)
        cls.payload = cert.certificate_payload()

    # --- schema and determinism ------------------------------------------------

    def test_schema_and_self_hash(self) -> None:
        self.assertEqual(self.stored["schema"], "oph.multiplicity_window_certificate.v1")
        self.assertEqual(self.stored["issues"], [616, 617])
        body = {k: v for k, v in self.stored.items() if k != "manifest_sha256"}
        self.assertEqual(
            self.stored["manifest_sha256"], "sha256:" + cert.sha256_json(body)
        )

    def test_payload_is_deterministic(self) -> None:
        again = cert.certificate_payload()
        self.assertEqual(self.payload, again)
        self.assertEqual(cert.sha256_json(self.payload), cert.sha256_json(again))

    def test_stored_manifest_is_exactly_recomputable(self) -> None:
        verified = cert.verify_manifest(self.manifest_path)
        body = {k: v for k, v in verified.items() if k != "manifest_sha256"}
        self.assertEqual(body, self.payload)

    def test_stored_manifest_round_trips_through_json(self) -> None:
        self.assertEqual(
            json.loads(json.dumps(self.payload, sort_keys=True)),
            {k: v for k, v in self.stored.items() if k != "manifest_sha256"},
        )

    def test_cli_verifies_and_exits_zero(self) -> None:
        result = subprocess.run(
            [sys.executable, str(MODULE_DIR / "multiplicity_window_certificate.py")],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("verified", result.stdout)
        self.assertIn(self.stored["manifest_sha256"], result.stdout)

    # --- Part A: scalar response multiplicity ---------------------------------

    def test_scalar_compatibility_block_matches_pinned_manifest(self) -> None:
        block = self.payload["scalar_response_multiplicity"]["compatibility_block"]
        self.assertEqual(block["admissible_scalar_charges"], [3, -3])
        self.assertEqual(block["representative_scalar_charge"], 3)
        self.assertEqual(block["channel_count"], 3)
        self.assertEqual(
            block["invariant_yukawa_channels"],
            [["L", "Sbar", "e_c"], ["Q", "S", "u_c"], ["Q", "Sbar", "d_c"]],
        )
        matter = cert.load_json(MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json")
        self.assertEqual(
            block["pinned_matter_lift_manifest"]["sha256"], cert.sha256_json(matter)
        )
        self.assertIn(
            "conditional matter fixture",
            block["pinned_matter_lift_manifest"]["scope"],
        )
        declared = sorted(
            sorted([row[0], row[2]])[:1] + [row[1]] + sorted([row[0], row[2]])[1:]
            for row in matter["exterior_matter_contract"]["yukawa_channels"]
        )
        self.assertEqual(block["invariant_yukawa_channels"], declared)

    def test_scalar_countermodels_pass_every_visible_check(self) -> None:
        battery = self.payload["scalar_response_multiplicity"]["countermodel_battery"]
        names = [row["name"] for row in battery["configurations"]]
        self.assertEqual(
            names, ["n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"]
        )
        for configuration in battery["configurations"]:
            self.assertTrue(configuration["all_checks_pass"])
            for check in configuration["checks"]:
                self.assertTrue(check["pass"], msg=str(check))
        duplicate = battery["configurations"][1]
        self.assertEqual(
            [row["invariant_channel_count"] for row in duplicate["scalar_copies"]],
            [3, 3],
        )
        inert = battery["configurations"][2]
        self.assertEqual(
            [row["invariant_channel_count"] for row in inert["scalar_copies"]],
            [3, 0],
        )
        self.assertEqual(battery["configurations"][0]["copy_count"], 0)
        self.assertTrue(battery["every_configuration_passes_every_check"])

    def test_scalar_verdict_and_registry_row(self) -> None:
        verdict = self.payload["scalar_response_multiplicity"]["verdict"]
        self.assertEqual(verdict["scalar_existence"], "not_source_determined")
        self.assertEqual(verdict["scalar_multiplicity"], "independence_limited")
        self.assertEqual(len(verdict["countermodels"]), 3)
        self.assertEqual(
            verdict["registry_row"],
            {
                "id": "scalar_existence_and_multiplicity",
                "class": "conditional_open_interface",
                "area": "matter",
                "owner_issue": 616,
            },
        )

    def test_fermion_anomaly_forms_vanish_exactly(self) -> None:
        forms = cert.fermion_anomaly_forms(Fraction(-1, 3), Fraction(1, 2))
        self.assertEqual(forms["grav"], 0)
        self.assertEqual(forms["su3_sq_u1"], 0)
        self.assertEqual(forms["su2_sq_u1"], 0)
        self.assertEqual(forms["u1_cubed"], 0)
        self.assertEqual(forms["weyl_doublets"], 4)

    def test_incompatible_scalar_charge_fails_closed(self) -> None:
        scan = cert.scalar_pair_scan(Fraction(-1, 3), Fraction(1, 2))
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.scalar_configuration_record(
                "bad_charge",
                [{"label": "X", "charge_q6Y": 1, "yukawa_coupling": "unit"}],
                scan,
                Fraction(-1, 3),
                Fraction(1, 2),
                "control",
            )
        self.assertEqual(ctx.exception.code, "SCALAR_CONFIGURATION")

    # --- Part B: family multiplicity window -----------------------------------

    def test_ckm_count_table_is_exact(self) -> None:
        edge = self.payload["family_multiplicity_window"]["cp_capability_lower_edge"]
        self.assertEqual(
            edge["count_table"],
            {
                "2": {"angles": 1, "phases": 0},
                "3": {"angles": 3, "phases": 1},
                "4": {"angles": 6, "phases": 3},
                "5": {"angles": 10, "phases": 6},
                "6": {"angles": 15, "phases": 10},
            },
        )
        self.assertEqual(edge["excluded_counts"], [1, 2])
        self.assertEqual(edge["lower_edge"], 3)
        self.assertTrue(edge["paper_convention"]["matches_certificate_formula"])

    def test_su2_coefficient_table_is_exact(self) -> None:
        edge = self.payload["family_multiplicity_window"]["su2_ultraviolet_upper_edge"]
        self.assertEqual(edge["b_at_5"], "1/2")
        self.assertEqual(edge["b_at_6"], "-5/6")
        self.assertEqual(
            edge["value_table"],
            {
                "1": "35/6",
                "2": "9/2",
                "3": "19/6",
                "4": "11/6",
                "5": "1/2",
                "6": "-5/6",
                "7": "-13/6",
                "8": "-7/2",
            },
        )
        self.assertEqual(edge["upper_edge"], 5)
        self.assertTrue(edge["paper_convention"]["includes_higgs_term"])
        self.assertTrue(edge["paper_convention"]["matches_certificate_formula"])
        self.assertEqual(cert.b_su2(5), Fraction(1, 2))
        self.assertEqual(cert.b_su2(6), Fraction(-5, 6))

    def test_a5_carrier_arithmetic(self) -> None:
        carrier = self.payload["family_multiplicity_window"]["screen_carrier"]
        self.assertEqual(carrier["irreducible_dimensions"], [1, 3, 3, 4, 5])
        self.assertTrue(carrier["character_table_orthonormal"])
        self.assertTrue(carrier["no_two_dimensional_irreducible"])
        self.assertEqual(carrier["smallest_nontrivial_carrier_dimension"], 3)
        self.assertTrue(carrier["three_slot_carrier_exists"])
        self.assertEqual(
            carrier["two_dimensional_modules"],
            [{"multiplicities": [2, 0, 0, 0, 0], "irreducible": False, "content": "1 + 1"}],
        )

    def test_in_window_non_selection(self) -> None:
        section = self.payload["family_multiplicity_window"]["in_window_non_selection"]
        self.assertEqual(section["window"], [3, 4, 5])
        self.assertEqual([row["n_g"] for row in section["members"]], [3, 4, 5])
        self.assertEqual([row["cp_phases"] for row in section["members"]], [1, 3, 6])
        self.assertEqual(
            [row["b_su2"] for row in section["members"]], ["19/6", "11/6", "1/2"]
        )
        for row in section["members"]:
            self.assertTrue(row["cp_clause_pass"])
            self.assertTrue(row["uv_clause_pass"])
            self.assertTrue(row["k_copies_grammar_admissible"])
            self.assertEqual(row["k_copies_witten_parity"], 0)
        self.assertFalse(section["copy_counting_clause_exists"])
        self.assertFalse(section["count_inside_window_source_selected"])

    def test_family_verdict_and_registry_row(self) -> None:
        verdict = self.payload["family_multiplicity_window"]["verdict"]
        self.assertEqual(verdict["window"], [3, 4, 5])
        self.assertTrue(verdict["window_is_exact"])
        self.assertEqual(verdict["count_inside_window"], "not_source_selected")
        self.assertEqual(
            verdict["registry_row"]["id"], "family_attachment_and_multiplicity"
        )
        self.assertEqual(verdict["registry_row"]["class"], "conditional_open_interface")

    # --- fail-closed controls --------------------------------------------------

    def test_controls_reject_with_typed_codes(self) -> None:
        controls = self.payload["fail_closed_controls"]
        self.assertEqual(
            [row["error_code"] for row in controls],
            ["CP_COUNT_REJECTION", "UV_SIGN_REJECTION", "A5_IRREP_REJECTION"],
        )
        self.assertTrue(all(row["rejected"] for row in controls))

    def test_control_functions_fail_closed_directly(self) -> None:
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.check_cp_phase_claim(2, 1)
        self.assertEqual(ctx.exception.code, "CP_COUNT_REJECTION")
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.check_asymptotic_freedom_claim(6)
        self.assertEqual(ctx.exception.code, "UV_SIGN_REJECTION")
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.check_a5_irreducible_dimension_claim(2)
        self.assertEqual(ctx.exception.code, "A5_IRREP_REJECTION")

    def test_control_functions_accept_true_values(self) -> None:
        cert.check_cp_phase_claim(3, 1)
        cert.check_asymptotic_freedom_claim(5)
        cert.check_a5_irreducible_dimension_claim(3)
        cert.check_a5_irreducible_dimension_claim(4)

    # --- mutation gates ---------------------------------------------------------

    def _write_mutated(self, mutate) -> Path:
        mutated = copy.deepcopy(self.stored)
        mutate(mutated)
        path = self.manifest_path.parent / "_mutated_multiplicity_window.json"
        path.write_text(json.dumps(mutated, indent=2, sort_keys=True) + "\n")
        self.addCleanup(path.unlink)
        return path

    def test_mutated_verdict_is_rejected(self) -> None:
        def mutate(manifest) -> None:
            manifest["scalar_response_multiplicity"]["verdict"][
                "scalar_existence"
            ] = "source_determined"

        path = self._write_mutated(mutate)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_manifest(path)
        self.assertEqual(ctx.exception.code, "MANIFEST_HASH")

    def test_mutated_window_with_recomputed_hash_is_rejected(self) -> None:
        def mutate(manifest) -> None:
            manifest["family_multiplicity_window"]["verdict"]["window"] = [3]
            body = {k: v for k, v in manifest.items() if k != "manifest_sha256"}
            manifest["manifest_sha256"] = "sha256:" + cert.sha256_json(body)

        path = self._write_mutated(mutate)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_manifest(path)
        self.assertEqual(ctx.exception.code, "MANIFEST_DRIFT")

    def test_mutated_b0_value_with_recomputed_hash_is_rejected(self) -> None:
        def mutate(manifest) -> None:
            manifest["family_multiplicity_window"]["su2_ultraviolet_upper_edge"][
                "b_at_6"
            ] = "5/6"
            body = {k: v for k, v in manifest.items() if k != "manifest_sha256"}
            manifest["manifest_sha256"] = "sha256:" + cert.sha256_json(body)

        path = self._write_mutated(mutate)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_manifest(path)
        self.assertEqual(ctx.exception.code, "MANIFEST_DRIFT")

    def test_mutated_schema_is_rejected(self) -> None:
        def mutate(manifest) -> None:
            manifest["schema"] = "oph.multiplicity_window_certificate.v0"

        path = self._write_mutated(mutate)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_manifest(path)
        self.assertEqual(ctx.exception.code, "SCHEMA")

    def test_mutated_countermodel_check_is_rejected(self) -> None:
        def mutate(manifest) -> None:
            battery = manifest["scalar_response_multiplicity"]["countermodel_battery"]
            battery["configurations"][0]["checks"][0]["pass"] = False
            body = {k: v for k, v in manifest.items() if k != "manifest_sha256"}
            manifest["manifest_sha256"] = "sha256:" + cert.sha256_json(body)

        path = self._write_mutated(mutate)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_manifest(path)
        self.assertEqual(ctx.exception.code, "MANIFEST_DRIFT")

    # --- exact arithmetic guards -----------------------------------------------

    def test_orthonormality_uses_exact_field_arithmetic(self) -> None:
        table = cert._a5_character_table()
        phi = cert.F5(Fraction(1, 2), Fraction(1, 2))
        psi = cert.F5(Fraction(1, 2), Fraction(-1, 2))
        self.assertEqual(table["3"][3], phi)
        self.assertEqual(table["3"][4], psi)
        self.assertEqual(phi + psi, cert.F5(1))
        self.assertEqual(phi * psi, cert.F5(-1))

    def test_scan_matches_matter_lift_certificate(self) -> None:
        scan = cert.scalar_pair_scan(Fraction(-1, 3), Fraction(1, 2))
        self.assertEqual(sorted(scan), [-3, 3])
        self.assertEqual(len(scan[3]), 3)
        self.assertEqual(len(scan[-3]), 3)
        conjugated = sorted(
            [left, "Sbar" if scalar == "S" else "S", right]
            for left, scalar, right in scan[3]
        )
        self.assertEqual(conjugated, scan[-3])


if __name__ == "__main__":
    unittest.main()
