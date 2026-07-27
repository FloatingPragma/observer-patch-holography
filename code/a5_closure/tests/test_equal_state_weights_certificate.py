#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #610 equal-state-weights certificate."""

from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import equal_state_weights_certificate as cert  # noqa: E402


class EqualStateWeightsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "equal_state_weights_reference.json"
        cls.expected = cert.build_payload()

    def test_reference_manifest_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(self.manifest_path)
        self.assertEqual(stored, self.expected)
        body = {key: value for key, value in stored.items() if key != "payload_sha256"}
        self.assertEqual(stored["payload_sha256"], "sha256:" + cert.sha256_json(body))

    def test_rerun_writes_byte_identical_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            first_bytes = first.read_bytes()
            self.assertEqual(first_bytes, second.read_bytes())
            self.assertEqual(first_bytes, self.manifest_path.read_bytes())

    def test_schema_and_identity_fields(self) -> None:
        self.assertEqual(self.expected["schema"], "oph.equal_state_weights_certificate.v1")
        self.assertEqual(self.expected["issue"], 610)
        self.assertEqual(
            self.expected["generated_by"],
            "code/a5_closure/equal_state_weights_certificate.py",
        )
        spine = self.expected["lean_spine"]
        self.assertEqual(spine["module"], "Lean/Screen/EqualStateWeights.lean")
        self.assertEqual(spine["action_module"], "Lean/Screen/A5PortAction.lean")
        self.assertEqual(
            spine["composition_theorem"], "OPH.EqualStateWeights.equal_state_weights"
        )

    def test_deck_action_report(self) -> None:
        deck = self.expected["deck_action"]
        self.assertEqual(deck["full_incidence_automorphism_order"], 120)
        self.assertEqual(deck["rotation_group_order"], 60)
        self.assertTrue(deck["identity_present"])
        self.assertTrue(deck["closure_verified"])
        self.assertTrue(deck["inverses_verified"])
        self.assertTrue(deck["adjacency_preserved"])
        self.assertTrue(
            deck["oriented_faces_preserved_with_consistent_cyclic_orientation"]
        )
        self.assertTrue(deck["antipode_i_to_11_minus_i_commutes"])
        self.assertEqual(
            deck["pairwise_transitivity"],
            {"ordered_port_pairs_checked": 144, "all_pairs_connected": True},
        )
        self.assertTrue(deck["matches_lean_listed_rotations"])

    def test_positive_lane_is_exactly_one_twelfth(self) -> None:
        lane = self.expected["positive_lane"]
        self.assertEqual(lane["projection"], ["1/12"] * 12)
        self.assertEqual(lane["result_per_port"], "1/12")
        stationarity = lane["kkt_stationarity"]
        self.assertEqual(stationarity["density_ratios_rho_over_tau"], ["1"] * 12)
        self.assertTrue(stationarity["gradient_components_all_equal"])
        self.assertFalse(stationarity["numerics_used"])
        constraint = lane["invariant_linear_constraint_case"]
        self.assertTrue(constraint["deck_invariant"])
        self.assertTrue(constraint["uniform_point_feasible"])
        self.assertTrue(constraint["uniform_point_stationary"])
        self.assertTrue(lane["reference_invariant_under_deck_action"])

    def test_projection_witness_is_exact(self) -> None:
        witness = cert.projection_witness(cert.WRONG_REFERENCE_TAU, Fraction(1))
        self.assertEqual(witness["projection_exact"], cert.WRONG_REFERENCE_TAU)
        self.assertTrue(
            all(isinstance(v, Fraction) for v in witness["projection_exact"])
        )
        with self.assertRaises(cert.CertificateError) as caught:
            cert.projection_witness((Fraction(1, 2),) * 12, Fraction(1))
        self.assertEqual(caught.exception.code, "REFERENCE_NORMALIZATION")
        with self.assertRaises(cert.CertificateError) as caught:
            cert.projection_witness(cert.UNIFORM, Fraction(0))
        self.assertEqual(caught.exception.code, "WEIGHTS_POSITIVE")

    def test_stabilizer_of_antipodal_pair_has_order_ten(self) -> None:
        rotations = cert.deck_rotations()
        subgroup = cert.stabilizer_of_antipodal_pair(rotations, 0)
        self.assertEqual(len(subgroup), 10)
        orbit_sizes = sorted(len(orbit) for orbit in cert.orbits_under(subgroup))
        self.assertEqual(orbit_sizes, [2, 10])

    def test_controls_record_required_failures(self) -> None:
        controls = self.expected["controls"]
        self.assertEqual(
            set(controls),
            {"wrong_reference", "non_transitive", "non_unique", "incomplete_cover"},
        )
        for name, verdict in controls.items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        wrong = controls["wrong_reference"]["witness"]
        self.assertEqual(wrong["tau_prime"][0], "2/13")
        self.assertEqual(wrong["projection"][0], "2/13")
        self.assertTrue(wrong["projection_equals_tau_prime"])
        non_transitive = controls["non_transitive"]["witness"]
        self.assertEqual(non_transitive["stabilizer_order"], 10)
        self.assertEqual(non_transitive["stabilizer_orbit_sizes"], [2, 10])
        self.assertEqual(non_transitive["projection_value_on_pair_orbit"], "1/4")
        self.assertEqual(non_transitive["projection_value_on_ring_orbit"], "1/20")
        self.assertFalse(non_transitive["projection_is_uniform"])
        non_unique = controls["non_unique"]["witness"]
        self.assertNotEqual(non_unique["minimizer_one"], non_unique["minimizer_two"])
        self.assertEqual(non_unique["objective_values"], ["1", "1"])
        cover = controls["incomplete_cover"]["witness"]
        self.assertEqual(cover["omitted_port"], 0)
        self.assertEqual(cover["state_one_at_omitted_port"], "1/12")
        self.assertEqual(cover["state_two_at_omitted_port"], "5/12")
        self.assertTrue(cover["restrictions_to_cover_equal"])
        self.assertFalse(cover["restriction_map_injective"])

    def test_block_trace_and_state_weight_are_distinct_fields(self) -> None:
        record = self.expected["state_weight_vs_block_trace"]
        self.assertEqual(record["measured_central_block_trace"], "1/12")
        self.assertEqual(record["derived_state_weight"], "1/12")
        self.assertTrue(record["distinct_objects_that_agree_at_one_twelfth"])
        carrier_manifest = cert.load_json(
            MODULE_DIR / "manifests" / "echosahedral_federation_reference.json"
        )
        self.assertEqual(
            record["carrier_manifest_sha256"], cert.sha256_json(carrier_manifest)
        )

    def test_tampered_rotation_breaks_closure(self) -> None:
        rotations = cert.deck_rotations()
        # Composing one rotation with the antipode keeps adjacency and
        # antipode commutation and flips orientation, so the first exact
        # failure is group closure.
        tampered = list(rotations)
        victim = next(g for g in tampered if g != tuple(range(12)))
        index = tampered.index(victim)
        tampered[index] = tuple(11 - v for v in victim)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_rotation_group(tampered)
        self.assertEqual(caught.exception.code, "GROUP_CLOSURE")

    def test_tampered_rotation_breaks_adjacency(self) -> None:
        rotations = list(cert.deck_rotations())
        victim = list(rotations[1])
        victim[0], victim[1] = victim[1], victim[0]
        rotations[1] = tuple(victim)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_rotation_group(rotations)
        self.assertEqual(caught.exception.code, "ADJACENCY_PRESERVATION")

    def test_dropped_rotation_breaks_group_order(self) -> None:
        rotations = cert.deck_rotations()
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_rotation_group(rotations[:59])
        self.assertEqual(caught.exception.code, "GROUP_ORDER")

    def test_forced_uniform_wrong_reference_fails_closed(self) -> None:
        original = cert.WRONG_REFERENCE_TAU
        cert.WRONG_REFERENCE_TAU = cert.UNIFORM
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.build_payload()
            self.assertEqual(caught.exception.code, "CONTROL_NOT_FAILED")
            with tempfile.TemporaryDirectory() as scratch:
                output = Path(scratch) / "forced.json"
                self.assertEqual(cert.main(["--output", str(output)]), 1)
                self.assertFalse(output.exists())
        finally:
            cert.WRONG_REFERENCE_TAU = original

    def test_no_floats_anywhere_in_payload(self) -> None:
        cert.require_no_floats(self.expected)
        poisoned = copy.deepcopy(self.expected)
        poisoned["positive_lane"]["poison"] = 1 / 12
        with self.assertRaises(cert.CertificateError) as caught:
            cert.require_no_floats(poisoned)
        self.assertEqual(caught.exception.code, "FLOAT_FORBIDDEN")

    def test_computed_rotations_match_lean_listing(self) -> None:
        self.assertEqual(
            set(cert.deck_rotations()), set(cert.lean_listed_rotations())
        )


if __name__ == "__main__":
    unittest.main()
