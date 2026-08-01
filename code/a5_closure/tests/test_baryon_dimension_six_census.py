#!/usr/bin/env python3
"""Regression and adversarial tests for the conditional #641 census."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import baryon_dimension_six_census as producer  # noqa: E402
import verify_baryon_dimension_six_census_independent as independent  # noqa: E402


def resign(document: dict) -> dict:
    document = copy.deepcopy(document)
    body = {key: value for key, value in document.items() if key != "receipt_sha256"}
    document["receipt_sha256"] = independent.digest_bytes(independent.canonical_bytes(body))
    return document


class BaryonDimensionSixCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = producer.build_receipt()

    def assert_independent_rejects(self, mutation) -> None:
        doctored = copy.deepcopy(self.receipt)
        mutation(doctored)
        with self.assertRaises(independent.IndependentVerificationError):
            independent.verify_document(resign(doctored), verify_files=False)

    def test_exact_operator_census(self) -> None:
        census = self.receipt["census"]
        self.assertEqual(census["two_fermion_multisets_checked"], 55)
        self.assertEqual(census["two_fermion_baryon_violating_gauge_lorentz_survivors"], 0)
        self.assertEqual(census["four_fermion_multisets_checked"], 715)
        self.assertEqual(census["oriented_monomial_count"], 8)
        self.assertEqual(census["hermitian_conjugacy_class_count"], 4)
        self.assertEqual(
            set(census["hermitian_conjugacy_classes"]),
            {"QQQL", "QQUE", "DUQL", "DUUE"},
        )

    def test_every_class_and_conjugate_preserves_b_minus_l(self) -> None:
        for row in self.receipt["census"]["oriented_monomials"]:
            self.assertEqual(row["q6"], 0)
            self.assertEqual(abs(row["b3"]), 3)
            self.assertEqual(abs(row["lepton_number"]), 1)
            self.assertEqual(row["three_times_delta_B_minus_L"], 0)

    def test_all_epsilon_pairings_are_quotiented_exactly(self) -> None:
        contractions = {
            row["class"]: row
            for row in self.receipt["explicit_nonzero_contractions"]
        }
        self.assertEqual(contractions["QQQL"]["epsilon_pairing_patterns_enumerated"], 9)
        self.assertEqual(contractions["QQQL"]["nonzero_pairing_patterns"], 9)
        self.assertEqual(contractions["DUUE"]["epsilon_pairing_patterns_enumerated"], 3)
        self.assertEqual(contractions["DUUE"]["nonzero_pairing_patterns"], 2)
        for name, row in contractions.items():
            self.assertEqual(
                row["grassmann_span_rank_after_schouten_fierz_and_pauli_relations"],
                1,
                name,
            )
            self.assertTrue(row["one_generation_representative_survives"], name)
            self.assertGreater(row["nonzero_grassmann_monomials"], 0, name)
            self.assertTrue(row["displayed_direct_contraction_lies_in_span"], name)

    def test_qqql_and_duue_survive_identical_field_antisymmetry(self) -> None:
        contractions = {
            row["class"]: row
            for row in self.receipt["explicit_nonzero_contractions"]
        }
        self.assertGreater(contractions["QQQL"]["nonzero_grassmann_monomials"], 0)
        self.assertGreater(contractions["DUUE"]["nonzero_grassmann_monomials"], 0)
        self.assertLess(
            contractions["DUUE"]["nonzero_pairing_patterns"],
            contractions["DUUE"]["epsilon_pairing_patterns_enumerated"],
            "the test must see the Pauli-vanishing u_R-u_R pairing",
        )

    def test_four_contraction_classes_are_globally_independent(self) -> None:
        basis = self.receipt["global_contraction_basis"]
        self.assertEqual(basis["exact_grassmann_span_rank"], 4)
        self.assertEqual(
            set(basis["representative_classes"]),
            {"QQQL", "QQUE", "DUQL", "DUUE"},
        )

    def test_diagonal_z6_does_not_remove_any_class(self) -> None:
        self.assertEqual(
            set(self.receipt["diagonal_z6"]["field_phases_sixths"].values()),
            {0},
        )
        self.assertIn("cannot remove", self.receipt["diagonal_z6"]["conclusion"])

    def test_physical_boundary_is_fail_closed(self) -> None:
        self.assertTrue(self.receipt["physical_boundary"])
        self.assertTrue(all(value is False for value in self.receipt["physical_boundary"].values()))
        self.assertIn(
            "declared accidental-charge labels",
            self.receipt["grammar"]["baryon_and_lepton_number_status"],
        )

    def test_stored_receipt_replays_and_is_byte_exact(self) -> None:
        self.assertEqual(producer.verify_stored()["status"], "PASS")
        with tempfile.TemporaryDirectory() as scratch:
            target = Path(scratch) / "receipt.json"
            producer.write_json(target, producer.build_receipt())
            self.assertEqual(target.read_bytes(), producer.RECEIPT_PATH.read_bytes())

    def test_independent_verifier_passes(self) -> None:
        verdict = independent.audit()
        self.assertEqual(verdict["status"], "PASS")
        self.assertEqual(verdict["independent_classes"], 4)
        self.assertEqual(verdict["global_contraction_rank"], 4)
        self.assertTrue(verdict["qqql_survives"])
        self.assertTrue(verdict["duue_survives"])

    def test_resigned_charge_mutation_is_rejected(self) -> None:
        self.assert_independent_rejects(lambda d: d["fields"]["u_R"].__setitem__("q6", 5))

    def test_resigned_omitted_operator_is_rejected(self) -> None:
        self.assert_independent_rejects(lambda d: d["census"]["oriented_monomials"].pop())

    def test_resigned_contraction_rank_mutation_is_rejected(self) -> None:
        self.assert_independent_rejects(
            lambda d: d["explicit_nonzero_contractions"][0].__setitem__(
                "grassmann_span_rank_after_schouten_fierz_and_pauli_relations", 0
            )
        )

    def test_resigned_grassmann_witness_mutation_is_rejected(self) -> None:
        self.assert_independent_rejects(
            lambda d: d["explicit_nonzero_contractions"][0]["first_nonzero_witness"].__setitem__(
                "coefficient", 0
            )
        )

    def test_resigned_scalar_baryon_charge_is_rejected(self) -> None:
        self.assert_independent_rejects(lambda d: d["grammar"]["declared_scalar"].__setitem__("b3", 3))

    def test_resigned_family_scope_widening_is_rejected(self) -> None:
        self.assert_independent_rejects(lambda d: d["grammar"].__setitem__("family_scope", "three_generations"))

    def test_resigned_physical_promotion_is_rejected(self) -> None:
        self.assert_independent_rejects(
            lambda d: d["physical_boundary"].__setitem__("proton_decay_predicted", True)
        )

    def test_resigned_right_handed_neutrino_boundary_deletion_is_rejected(self) -> None:
        self.assert_independent_rejects(
            lambda d: d["grammar"]["excluded_extensions"].remove("right-handed neutrino")
        )

    def test_raw_source_pin_mutation_is_rejected(self) -> None:
        doctored = copy.deepcopy(self.receipt)
        doctored["source_pins"][0]["sha256"] = "sha256:" + "0" * 64
        with self.assertRaises(independent.IndependentVerificationError):
            independent.verify_document(resign(doctored), verify_files=True)

    def test_json_round_trip_is_stable(self) -> None:
        encoded = producer.canonical_json_bytes(self.receipt)
        self.assertEqual(json.loads(encoded), self.receipt)
        self.assertTrue(encoded.endswith(b"\n"))


if __name__ == "__main__":
    unittest.main()
