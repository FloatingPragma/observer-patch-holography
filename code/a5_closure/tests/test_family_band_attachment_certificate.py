#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #569 family-band certificate."""

from __future__ import annotations

import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import family_band_attachment_certificate as cert  # noqa: E402

F5 = cert.F5


class FamilyBandAttachmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()
        cls.carrier, cls.rotations, _ = cert.load_carrier()
        cls.adjacency = cert.lift(cert.adjacency_int(cls.carrier))
        cls.projectors = cert.spectral_projectors(cls.adjacency)

    def test_reference_manifest_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(cert.MANIFEST_PATH)
        body = {k: v for k, v in stored.items() if k != "manifest_sha256"}
        self.assertEqual(body, self.payload)
        self.assertEqual(stored["manifest_sha256"], "sha256:" + cert.sha256_json(body))
        self.assertEqual(cert.verify_stored()["status"], "PASS")

    def test_rerun_writes_byte_identical_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_schema_issue_and_bounded_exit(self) -> None:
        self.assertEqual(self.payload["schema"], cert.SCHEMA)
        self.assertEqual(self.payload["issue"], 569)
        self.assertEqual(self.payload["bounded_exit"], "exact_named_realization")
        self.assertTrue(self.payload["invisibility_preserved"])
        interface = self.payload["named_interface"]
        self.assertEqual(interface["class"], "conditional_open_interface")
        self.assertEqual(
            sorted(interface["clauses"]), ["R_realization", "S_selection"]
        )
        self.assertEqual(
            interface["clause_controls"],
            {"R_realization": "external_copy_reduct", "S_selection": "excluded_cone"},
        )

    def test_selection_is_the_three_band_with_rank_45(self) -> None:
        self.assertEqual(self.payload["selection"]["minimizer"], "3")
        self.assertTrue(self.payload["selection"]["strict"])
        order = [row["object"] for row in self.payload["selection"]["order"]]
        self.assertEqual(order, ["3", "5", "3p"])
        self.assertEqual(self.payload["attachment"]["family_dimension"], 3)
        self.assertEqual(self.payload["attachment"]["complex_rank"], 45)
        self.assertEqual(self.payload["band_action_kernels"], {"1": 60, "3": 1, "3p": 1, "5": 1})

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        self.assertEqual(
            self.payload["controls"]["excluded_cone"]["excluded_readback_minimizer"],
            "3p",
        )
        self.assertEqual(
            self.payload["controls"]["dropped_faithfulness"]["unguarded_minimizer"],
            "1",
        )
        self.assertEqual(
            self.payload["controls"]["measured_channel_swap"]["swapped_minimizer"],
            "3p",
        )

    def test_measured_receipt_realizes_clause_s(self) -> None:
        measured = self.payload["measured_receipt"]
        self.assertEqual(
            measured["clause_S"], "simulator_realized_on_declared_channel"
        )
        self.assertEqual(
            measured["clause_R"],
            "simulator_realized_for_declared_response_resolvent__matter_pole_identification_open",
        )
        self.assertEqual(measured["measured_minimizer"], "frame_band (the 3 band)")
        order = [row["object"] for row in measured["measured_cost_order"]]
        self.assertEqual(order, ["3", "5", "3p"])
        self.assertEqual(
            self.payload["named_interface"]["clause_status"],
            {
                "R_realization": (
                    "simulator_realized_for_declared_response_resolvent__"
                    "matter_pole_identification_open"
                ),
                "S_selection": "simulator_realized_on_declared_channel",
            },
        )
        self.assertEqual(
            self.payload["upstream_pins"]["measured_response_artifact"]["issue"], 599
        )

    def test_artifact_carrier_pin_is_enforced(self) -> None:
        _, _, carrier_pin = cert.load_carrier()
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.pin_response_artifact({"sha256": "sha256:doctored"})
        self.assertIn("ARTIFACT_CARRIER_MISMATCH", str(ctx.exception))
        artifact, pin = cert.pin_response_artifact(carrier_pin)
        self.assertEqual(pin["issue"], 599)

    def test_channel_parser_is_exact(self) -> None:
        self.assertEqual(cert.parse_channel("5"), cert.F5(5, 0))
        self.assertEqual(cert.parse_channel("-1"), cert.F5(-1, 0))
        self.assertEqual(cert.parse_channel("0 + 1*sqrt(5)"), cert.F5(0, 1))
        self.assertEqual(cert.parse_channel("0 + -1*sqrt(5)"), cert.F5(0, -1))

    def test_doctored_measured_channels_fail(self) -> None:
        _, rotations, carrier_pin = cert.load_carrier()
        artifact, _ = cert.pin_response_artifact(carrier_pin)
        doctored = {k: (dict(v) if isinstance(v, dict) else v) for k, v in artifact.items()}
        doctored["response_basis"] = dict(artifact["response_basis"])
        doctored["response_basis"]["adjacency_channel_values"] = dict(
            artifact["response_basis"]["adjacency_channel_values"]
        )
        doctored["response_basis"]["adjacency_channel_values"]["frame_band"] = "5"
        kernels = {"1": 60, "3": 1, "3p": 1, "5": 1}
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.measured_band_receipt(
                self.carrier, self.adjacency, self.projectors, doctored, (5, -1), kernels
            )
        self.assertIn("MEASURED_CHANNEL_MISMATCH", str(ctx.exception))

    def test_spectral_gates_reject_a_doctored_adjacency(self) -> None:
        neighbor = sorted(self.carrier.adjacency[0])[0]
        self.assertEqual(self.adjacency[0][neighbor], F5(1, 0))
        doctored = [row[:] for row in self.adjacency]
        doctored[0][neighbor] = F5(0, 0)
        doctored[neighbor][0] = F5(0, 0)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_spectral_resolution(
                doctored, cert.spectral_projectors(doctored)
            )
        self.assertIn("PROJ", str(ctx.exception))

    def test_equivariance_gate_rejects_a_doctored_projector(self) -> None:
        doctored = {k: [row[:] for row in v] for k, v in self.projectors.items()}
        doctored["3"][0][0] = doctored["3"][0][0] + F5(1, 0)
        with self.assertRaises(cert.CertificateError):
            cert.verify_equivariance(doctored, self.rotations)

    def test_strict_minimizer_rejects_a_tied_cost_table(self) -> None:
        tied = {"3": F5(6, 0), "5": F5(6, 0), "3p": F5(5, 1)}
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.strict_minimizer(["3", "5", "3p"], tied)
        self.assertIn("COST_ORDER_NOT_STRICT", str(ctx.exception))

    def test_f5_sort_is_numeric_not_lexicographic(self) -> None:
        # 6 - sqrt5 < 5 < 6 + sqrt5; a lexicographic pair sort would put 5 first.
        costs = {"a": F5(6, 1), "b": F5(6, -1), "c": F5(5, 0)}
        ordered = cert.f5_sorted(["a", "b", "c"], lambda n: costs[n])
        self.assertEqual(ordered, ["b", "c", "a"])

    def test_candidate_enumeration_requires_the_exact_survivor_set(self) -> None:
        costs = cert.band_costs((5, -1))
        kernels = {"1": 60, "3": 1, "3p": 1, "5": 1}
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.enumerate_candidates((3, 4), kernels, costs)
        self.assertIn("CANDIDATE_SET", str(ctx.exception))

    def test_generation_certificate_rejects_drifted_charges(self) -> None:
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.generation_certificate(
                {"color_block": Fraction(-1, 3), "weak_block": Fraction(1, 3)}
            )
        self.assertIn("GENERATION_ANOMALY", str(ctx.exception))

    def test_tampered_stored_manifest_fails_verify(self) -> None:
        stored = cert.load_json(cert.MANIFEST_PATH)
        original = cert.MANIFEST_PATH.read_bytes()
        try:
            tampered = dict(stored)
            tampered["attachment"] = dict(stored["attachment"])
            tampered["attachment"]["complex_rank"] = 60
            cert.write_json(cert.MANIFEST_PATH, tampered)
            with self.assertRaises(cert.CertificateError):
                cert.verify_stored()
        finally:
            cert.MANIFEST_PATH.write_bytes(original)

    def test_schema_is_v4(self) -> None:
        self.assertEqual(cert.SCHEMA, "oph.family_band_attachment_certificate.v4")
        self.assertEqual(self.payload["schema"], cert.SCHEMA)

    def test_pole_residue_receipt_realizes_clause_r(self) -> None:
        receipt = self.payload["pole_residue_receipt"]
        self.assertEqual(receipt["family_band_residue"]["band"], "frame")
        self.assertEqual(receipt["family_band_residue"]["measured_rank"], 3)
        self.assertTrue(receipt["family_band_residue"]["equals_exact_frame_projector"])
        self.assertEqual(receipt["attachment_rank"]["complex_rank"], 45)
        self.assertEqual(
            self.payload["named_interface"]["clause_status"]["R_realization"],
            "simulator_realized_for_declared_response_resolvent__matter_pole_identification_open",
        )
        self.assertEqual(
            self.payload["upstream_pins"]["measured_pole_residue_artifact"]["issue"],
            569,
        )

    def test_pole_artifact_chain_pin_is_enforced(self) -> None:
        _, _, carrier_pin = cert.load_carrier()
        response, _ = cert.pin_response_artifact(carrier_pin)
        doctored = dict(response)
        doctored["artifact_sha256"] = "sha256:doctored"
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.pin_pole_residue_artifact(carrier_pin, doctored)
        self.assertIn("POLE_ARTIFACT_CHAIN", str(ctx.exception))

    def test_pole_table_mutation_control_fired(self) -> None:
        verdict = self.payload["controls"]["pole_table_mutation"]
        self.assertTrue(verdict["expected_failure"])
        self.assertTrue(verdict["failed"])

    def test_window_pin_matches_the_stored_receipt(self) -> None:
        lower, upper, pin = cert.pin_window()
        self.assertEqual((lower, upper), (3, 5))
        self.assertEqual(pin["issue"], 617)
        _, excluded, cone_pin = cert.pin_cost_cone()
        self.assertEqual(tuple(excluded), (6, 1))
        self.assertEqual(cone_pin["issue"], 625)
        charges, matter_pin = cert.pin_matter()
        self.assertEqual(charges["color_block"], Fraction(-1, 3))
        self.assertEqual(matter_pin["issue"], 314)


if __name__ == "__main__":
    unittest.main()
