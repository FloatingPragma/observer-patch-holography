"""Tests for the directed-interval certified contour verifier."""

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "producers"))

from mpmath import iv  # noqa: E402

import ball_arithmetic as ba  # noqa: E402
import certified_wz_contours as cc  # noqa: E402
from complex_interval import CInterval, SheetError  # noqa: E402


class ComplexIntervalTests(unittest.TestCase):
    def setUp(self) -> None:
        iv.prec = 128
        self.gate = (iv.mpf(157) / iv.mpf(100)).a

    def test_square_of_zero_straddling_interval_is_nonnegative(self) -> None:
        z = CInterval.box(
            Fraction(-1, 10), Fraction(1, 10), Fraction(1, 2), Fraction(1, 2)
        )
        self.assertGreaterEqual(float(z.abs2().a), 0.0)
        self.assertFalse(z.contains_zero())

    def test_division_by_zero_interval_refused(self) -> None:
        z = CInterval.from_fraction(1)
        origin = CInterval.box(
            Fraction(-1, 10), Fraction(1, 10), Fraction(-1, 10), Fraction(1, 10)
        )
        with self.assertRaises(SheetError):
            _ = z / origin

    def test_log_refuses_the_principal_cut(self) -> None:
        on_cut = CInterval.box(
            Fraction(-2), Fraction(-1), Fraction(-1, 10), Fraction(1, 10)
        )
        with self.assertRaises(SheetError):
            on_cut.log(self.gate)
        off_cut = CInterval.box(
            Fraction(-2), Fraction(-1), Fraction(1, 10), Fraction(2, 10)
        )
        off_cut.log(self.gate)

    def test_sqrt_matches_squaring(self) -> None:
        z = CInterval.box(
            Fraction(3, 2), Fraction(3, 2), Fraction(1, 3), Fraction(1, 3)
        )
        root = z.sqrt(self.gate)
        squared = root * root
        self.assertLessEqual(float(squared.re.a), 1.5)
        self.assertGreaterEqual(float(squared.re.b), 1.5)
        self.assertLessEqual(float(squared.im.a), 1.0 / 3)
        self.assertGreaterEqual(float(squared.im.b), 1.0 / 3)


class IntervalLoopTests(unittest.TestCase):
    def setUp(self) -> None:
        iv.prec = 128
        self.gate = (iv.mpf(157) / iv.mpf(100)).a

    def _contains(self, interval: CInterval, ball) -> bool:
        return bool(
            float(interval.re.a) - 1e-9 <= float(ball.mid_re)
            and float(ball.mid_re) <= float(interval.re.b) + 1e-9
            and float(interval.im.a) - 1e-9 <= float(ball.mid_im)
            and float(ball.mid_im) <= float(interval.im.b) + 1e-9
        )

    def test_b0_contains_sampled_convention_in_all_mass_cases(self) -> None:
        point = complex(0.111, 0.001)
        s = CInterval.box(
            Fraction(111, 1000), Fraction(111, 1000),
            Fraction(1, 1000), Fraction(1, 1000),
        )
        for m1, m2 in (
            (Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1, 900)),
            (Fraction(1, 9), Fraction(0)),
            (Fraction(1, 2500), Fraction(1, 400)),
        ):
            interval = cc.b0_interval(s, m1, m2, Fraction(1), self.gate)
            sampled = ba.b0_fin(point, float(m1), float(m2), 1.0, precision=192)
            self.assertTrue(self._contains(interval, sampled), (m1, m2))


class SyntheticWindingTests(unittest.TestCase):
    class _LinearStub:
        """Duck-typed evaluator for f(s) = s - root."""

        def __init__(self, root: CInterval) -> None:
            self.root = root
            self.loop_cache: dict = {}

        def inverse_propagator(self, s, s_key, tree_mass):
            return s - self.root

        def inverse_propagator_derivative(self, s, s_key, tree_mass):
            return CInterval.from_fraction(1)

    def setUp(self) -> None:
        iv.prec = 128
        self.gate = (iv.mpf(157) / iv.mpf(100)).a
        self.box = {
            "re": (Fraction(0), Fraction(1)),
            "im": (Fraction(0), Fraction(1)),
        }

    def test_enclosed_root_gives_winding_one(self) -> None:
        stub = self._LinearStub(
            CInterval.from_fraction(Fraction(1, 2), Fraction(1, 2))
        )
        verdict = cc.certify_winding(stub, self.box, Fraction(0), self.gate)
        self.assertTrue(verdict["certified"])
        self.assertEqual(verdict["winding"], 1)

    def test_external_root_gives_winding_zero(self) -> None:
        stub = self._LinearStub(
            CInterval.from_fraction(Fraction(3), Fraction(3))
        )
        verdict = cc.certify_winding(stub, self.box, Fraction(0), self.gate)
        self.assertTrue(verdict["certified"])
        self.assertEqual(verdict["winding"], 0)

    def test_root_on_boundary_fails_closed(self) -> None:
        stub = self._LinearStub(
            CInterval.from_fraction(Fraction(1, 2), Fraction(0))
        )
        verdict = cc.certify_winding(stub, self.box, Fraction(0), self.gate)
        self.assertFalse(verdict["certified"])


class FrozenReceiptTests(unittest.TestCase):
    def test_frozen_receipt_status_and_flags(self) -> None:
        receipt = json.loads(
            (PACKAGE_ROOT / "outputs" / "certified_wz_contours.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            receipt["status"], "PRINCIPAL_SHEET_ZERO_EXCLUSION_CERTIFIED"
        )
        promotion = receipt["promotion"]
        self.assertTrue(promotion["complex_ball_certified"])
        self.assertTrue(promotion["principal_sheet_zero_exclusion_certified"])
        self.assertTrue(promotion["root_count_certified_on_declared_boxes"])
        self.assertFalse(promotion["pole_enclosure_certified"])
        self.assertFalse(promotion["second_sheet_certified"])
        self.assertFalse(promotion["laurent_residue_certified"])
        self.assertFalse(promotion["bmhv_restoration_certified"])
        self.assertFalse(promotion["physical_current_claim"])
        self.assertFalse(promotion["oph_native"])
        self.assertFalse(promotion["unit_claim"])
        self.assertIn("artifact_finding", receipt)
        for name in ("W", "Z"):
            for row in receipt["results"][name].values():
                self.assertTrue(row["zero_exclusion_certified"])
                self.assertEqual(row["boundary_winding"]["winding"], 0)
                self.assertTrue(
                    row["interior_holomorphy"]["holomorphic_on_box"]
                )
            self.assertTrue(
                receipt["precision_nesting"][name][
                    "enclosures_nested_with_precision"
                ]
            )
            per_quantity = receipt["precision_nesting"][name][
                "per_quantity_probe_nesting"
            ]
            self.assertIn("inverse_propagator", per_quantity)
            self.assertIn("inverse_propagator_derivative", per_quantity)
            self.assertTrue(all(per_quantity.values()))
            self.assertIn(
                name, receipt["dimensional_prefactor_finite_correction"]
            )
            self.assertNotEqual(
                receipt["dimensional_prefactor_finite_correction"][name], "0"
            )
        gates = receipt["serialized_gates"]
        for key in (
            "precisions_bits",
            "initial_segments_per_edge",
            "max_subdivision_depth",
            "holomorphy_method",
            "boundary_method",
            "arg_width_gate",
            "winding_tolerance",
        ):
            self.assertIn(key, gates)


class SecondSheetFrozenReceiptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.receipt = json.loads(
            (
                PACKAGE_ROOT / "outputs" / "certified_second_sheet_poles.json"
            ).read_text(encoding="utf-8")
        )

    def test_status_and_windings(self) -> None:
        self.assertEqual(
            self.receipt["status"],
            "SECOND_SHEET_POLE_CERTIFIED_ON_DECLARED_CONTINUATION",
        )
        for name in ("W", "Z"):
            for row in self.receipt["results"][name].values():
                self.assertTrue(row["pole_certified"])
                self.assertTrue(row["simple_root_certified"])
                self.assertEqual(row["boundary_winding"]["winding"], 1)
                null_vectors = row["interval_newton"]["null_vectors"]
                self.assertTrue(null_vectors["left_residual_contains_zero"])
                self.assertTrue(null_vectors["right_residual_contains_zero"])
                self.assertTrue(
                    null_vectors["laurent_denominator_excludes_zero"]
                )

    def test_claim_boundary_flags_stay_false(self) -> None:
        promotion = self.receipt["promotion"]
        self.assertFalse(promotion["matrix_rank_laurent_certified"])
        self.assertFalse(promotion["bmhv_restoration_certified"])
        self.assertFalse(promotion["physical_current_claim"])
        self.assertFalse(promotion["oph_native"])
        self.assertFalse(promotion["unit_claim"])

    def test_continuation_probes_and_nesting(self) -> None:
        for name in ("W", "Z"):
            continuation = self.receipt["declared_continuation"][name]
            self.assertTrue(continuation["consistency_probes"]["all_passed"])
            nesting = self.receipt["precision_nesting"][name]
            self.assertTrue(nesting["enclosures_nested_with_precision"])
            self.assertTrue(
                nesting["per_segment_enclosure_nesting"]["all_nested"]
            )
            self.assertTrue(all(nesting["newton_ball_nesting"].values()))

    def test_zero_exclusion_receipt_is_pinned(self) -> None:
        import hashlib

        stored = self.receipt["pins"]["zero_exclusion_receipt_sha256"]
        actual = "sha256:" + hashlib.sha256(
            (
                PACKAGE_ROOT / "outputs" / "certified_wz_contours.json"
            ).read_bytes()
        ).hexdigest()
        self.assertEqual(stored, actual)


if __name__ == "__main__":
    unittest.main()
