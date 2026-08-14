"""Tests for the directed-interval certified contour verifier."""

import copy
import importlib.util
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

principal_checker_spec = importlib.util.spec_from_file_location(
    "principal_contour_evidence_checker",
    PACKAGE_ROOT / "checkers" / "check_certified_wz_contours.py",
)
principal_checker = importlib.util.module_from_spec(principal_checker_spec)
principal_checker_spec.loader.exec_module(principal_checker)

second_checker_spec = importlib.util.spec_from_file_location(
    "declared_chart_pole_evidence_checker",
    PACKAGE_ROOT / "checkers" / "check_certified_second_sheet_poles.py",
)
second_checker = importlib.util.module_from_spec(second_checker_spec)
second_checker_spec.loader.exec_module(second_checker)


def _rehash(receipt, checker_module) -> None:
    body = {
        key: value
        for key, value in receipt.items()
        if key != "receipt_sha256"
    }
    if hasattr(checker_module, "canonical_digest"):
        receipt["receipt_sha256"] = checker_module.canonical_digest(body)
    else:
        receipt["receipt_sha256"] = checker_module.sha256_bytes(
            checker_module.canonical_json(body).encode("utf-8")
        )


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

    def test_one_mass_b0_is_mass_exchange_symmetric_off_axis(self) -> None:
        mass = Fraction(1, 100)
        for imaginary_part in (Fraction(1, 1000), Fraction(-1, 1000)):
            s = CInterval.box(
                Fraction(1, 10),
                Fraction(1, 10),
                imaginary_part,
                imaginary_part,
            )
            b0_left = cc.b0_interval(s, 0, mass, 1, self.gate)
            b0_right = cc.b0_interval(s, mass, 0, 1, self.gate)
            b0p_left = cc.b0p_interval(s, 0, mass, 1, self.gate)
            b0p_right = cc.b0p_interval(s, mass, 0, 1, self.gate)
            self.assertEqual(
                cc.serialize_cinterval(b0_left),
                cc.serialize_cinterval(b0_right),
            )
            self.assertEqual(
                cc.serialize_cinterval(b0p_left),
                cc.serialize_cinterval(b0p_right),
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
            self.assertIn("probe:center:inverse_propagator", per_quantity)
            self.assertIn(
                "probe:center:inverse_propagator_derivative", per_quantity
            )
            self.assertTrue(all(per_quantity.values()))
            self.assertTrue(
                receipt["precision_nesting"][name][
                    "coefficient_denominator_nesting"
                ]["all_nested"]
            )
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
        scope = receipt["acceptance_scope"]
        self.assertTrue(scope["auxiliary_principal_sheet_zero_exclusion_only"])
        self.assertFalse(scope["independent_numerical_replay_certified"])
        self.assertEqual(
            scope["engine_inverse_propagator_convention"],
            "G(s)=s-m_tree^2-Pi_engine(s)",
        )
        self.assertFalse(scope["theorem_self_energy_sign_bridge_certified"])
        self.assertFalse(scope["precision_ladder_receipt_supplied"])
        self.assertFalse(
            scope["independent_numerical_reevaluation_supplied"]
        )
        self.assertFalse(scope["root_laurent_evidence_supplied"])
        self.assertFalse(scope["full_physical_pole_evidence_supplied"])


class PrincipalCheckerMutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = json.loads(
            (
                PACKAGE_ROOT / "outputs" / "certified_wz_contours.json"
            ).read_text(encoding="utf-8")
        )

    def assert_rehashed_mutation_rejected(self, mutate):
        forged = copy.deepcopy(self.receipt)
        mutate(forged)
        _rehash(forged, principal_checker)
        with self.assertRaises(
            principal_checker.ReceiptValidationError
        ) as caught:
            principal_checker.validate_receipt(forged)
        return caught.exception

    def test_canonical_v3_receipt_passes_evidence_checker(self) -> None:
        verdict = principal_checker.validate_receipt(self.receipt)
        self.assertTrue(verdict["evidence_contract_valid"])
        self.assertFalse(
            verdict["independent_physics_reevaluation_performed"]
        )

    def test_fixture_mutation_with_valid_self_hash_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["fixture"].__setitem__("g1", "999")
        )

    def test_deleted_loop_chart_with_valid_self_hash_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"][
                "interior_holomorphy"
            ]["loop_charts"].pop()
        )

    def test_empty_quantity_map_with_valid_self_hash_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"].__setitem__(
                "quantity_enclosures", {}
            )
        )

    def test_non_nested_quantity_with_valid_self_hash_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            key = "probe:center:inverse_propagator"
            outer_lo = Fraction(
                receipt["results"]["W"]["128"]["quantity_enclosures"][key][
                    "re"
                ]["lo"]
            )
            receipt["results"]["W"]["192"]["quantity_enclosures"][key]["re"][
                "lo"
            ] = str(outer_lo - 1)

        self.assert_rehashed_mutation_rejected(mutate)

    def test_forged_argument_increment_with_valid_hash_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            increment = receipt["results"]["W"]["128"]["boundary_winding"][
                "segment_evidence"
            ][0]["endpoint_increment"]
            increment["lo"] = "0"
            increment["hi"] = "0"

        self.assert_rehashed_mutation_rejected(mutate)

    def test_forged_subdivision_depth_with_valid_hash_is_rejected(self) -> None:
        error = self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"][
                "boundary_winding"
            ].__setitem__("max_depth_used", 1)
        )
        self.assertTrue(
            any("max_depth_used" in problem for problem in error.problems)
        )

    def test_double_ccw_traversal_with_valid_hash_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            box = principal_checker.BOXES["W"]
            re_lo, re_hi = box["re"]
            im_lo, im_hi = box["im"]
            corners = (
                (re_lo, im_lo),
                (re_hi, im_lo),
                (re_hi, im_hi),
                (re_lo, im_hi),
            )

            def loop(edge_counts):
                segments = []
                for edge, count in enumerate(edge_counts):
                    edge_start = corners[edge]
                    edge_end = corners[(edge + 1) % 4]
                    for step in range(count):
                        t0 = Fraction(step, count)
                        t1 = Fraction(step + 1, count)
                        start = (
                            edge_start[0]
                            + (edge_end[0] - edge_start[0]) * t0,
                            edge_start[1]
                            + (edge_end[1] - edge_start[1]) * t0,
                        )
                        end = (
                            edge_start[0]
                            + (edge_end[0] - edge_start[0]) * t1,
                            edge_start[1]
                            + (edge_end[1] - edge_start[1]) * t1,
                        )
                        segments.append((start, end))
                return segments

            one = {
                "re": {"lo": "1", "hi": "1"},
                "im": {"lo": "0", "hi": "0"},
            }
            zero = {
                "re": {"lo": "0", "hi": "0"},
                "im": {"lo": "0", "hi": "0"},
            }
            real_zero = {"lo": "0", "hi": "0"}

            def segment_record(start, end):
                start_text = [str(start[0]), str(start[1])]
                end_text = [str(end[0]), str(end[1])]
                midpoint = (
                    (start[0] + end[0]) / 2,
                    (start[1] + end[1]) / 2,
                )
                offset = {
                    "re": {
                        "lo": str(
                            min(start[0], end[0]) - midpoint[0]
                        ),
                        "hi": str(
                            max(start[0], end[0]) - midpoint[0]
                        ),
                    },
                    "im": {
                        "lo": str(
                            min(start[1], end[1]) - midpoint[1]
                        ),
                        "hi": str(
                            max(start[1], end[1]) - midpoint[1]
                        ),
                    },
                }
                return {
                    "segment_id": principal_checker.canonical_digest(
                        {"start": start_text, "end": end_text}
                    ),
                    "start": start_text,
                    "end": end_text,
                    "midpoint": [str(midpoint[0]), str(midpoint[1])],
                    "center_value": one,
                    "derivative_hull": zero,
                    "offset": offset,
                    "image": one,
                    "rotated_image": one,
                    "image_zero_exclusion_abs2_lower": "1",
                    "rotated_zero_exclusion_abs2_lower": "1",
                    "rotated_argument": real_zero,
                    "rotated_argument_width": "0",
                    "argument_width_slack": "157/100",
                    "start_value": one,
                    "end_value": one,
                    "endpoint_ratio": one,
                    "endpoint_increment": real_zero,
                    "endpoint_increment_slack": "157/100",
                }

            # Two complete laps with the frozen 51-record count.  The first
            # lap deliberately reaches edge 3 before the second returns to
            # edge 0, which the old orientation-only check accepted.
            forged_segments = loop((1, 1, 1, 8)) + loop((8, 8, 8, 16))
            evidence = [
                segment_record(start, end)
                for start, end in forged_segments
            ]
            segment_ids = [
                segment["segment_id"] for segment in evidence
            ]
            partition_digest = principal_checker.canonical_digest(
                [
                    {
                        "segment_id": segment["segment_id"],
                        "start": segment["start"],
                        "end": segment["end"],
                    }
                    for segment in evidence
                ]
            )
            quantity_keys = principal_checker._expected_quantity_keys(
                "W", segment_ids
            )
            for precision in ("128", "192", "256"):
                row = receipt["results"]["W"][precision]
                boundary = row["boundary_winding"]
                boundary["segment_evidence"] = copy.deepcopy(evidence)
                boundary["partition_sha256"] = partition_digest
                boundary["max_depth_used"] = (
                    1 if precision == "128" else 0
                )
                boundary["total_variation_interval"] = real_zero
                boundary["winding_residual"] = real_zero
                boundary["winding_tolerance_slack"] = "157/100"
                row["quantity_enclosures"] = {
                    key: zero for key in quantity_keys
                }
            receipt["precision_nesting"]["W"][
                "per_quantity_probe_nesting"
            ] = {key: True for key in quantity_keys}

        error = self.assert_rehashed_mutation_rejected(mutate)
        self.assertTrue(
            any(
                "not exactly one CCW traversal" in problem
                for problem in error.problems
            )
        )

    def test_zero_denominator_witness_with_valid_hash_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            witness = receipt["results"]["W"]["128"][
                "interior_holomorphy"
            ]["coefficient_denominator_witnesses"][0]
            witness["enclosure"] = {
                "re": {"lo": "0", "hi": "0"},
                "im": {"lo": "0", "hi": "0"},
            }
            witness["zero_exclusion_abs2_lower"] = "0"
            witness["excludes_zero"] = False

        self.assert_rehashed_mutation_rejected(mutate)

    def test_issue_row_overclaim_with_valid_self_hash_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["acceptance_scope"].__setitem__(
                "root_laurent_evidence_supplied", True
            )
        )

    def test_nested_prose_overclaim_with_valid_hash_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"][
                "interior_holomorphy"
            ]["loop_charts"][0].__setitem__(
                "certificate",
                "FULL SECOND-SHEET POLE AND LAURENT CERTIFIED",
            )
        )

    def test_consistently_downgraded_incomplete_receipt_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            segment = receipt["results"]["W"]["128"]["boundary_winding"][
                "segment_evidence"
            ][0]
            segment["image_zero_exclusion_abs2_lower"] = "0"
            boundary = receipt["results"]["W"]["128"]["boundary_winding"]
            boundary["certified"] = False
            boundary["reason"] = "forged incomplete row"
            receipt["results"]["W"]["128"][
                "zero_exclusion_certified"
            ] = False
            for key in (
                "complex_ball_certified",
                "sheet_certified_on_declared_boxes",
                "principal_sheet_zero_exclusion_certified",
                "root_count_certified_on_declared_boxes",
            ):
                receipt["promotion"][key] = False
            receipt["status"] = "CERTIFICATION_INCOMPLETE"

        self.assert_rehashed_mutation_rejected(mutate)


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
            "SCALAR_POLE_CERTIFIED_ON_DECLARED_ALGEBRAIC_CHART",
        )
        for name in ("W", "Z"):
            for row in self.receipt["results"][name].values():
                self.assertTrue(row["scalar_pole_certified"])
                self.assertTrue(row["simple_scalar_root_certified"])
                self.assertEqual(row["boundary_winding"]["winding"], 1)
                null_vectors = row["interval_newton"]["null_vectors"]
                self.assertTrue(null_vectors["left_residual_contains_zero"])
                self.assertTrue(null_vectors["right_residual_contains_zero"])
                self.assertTrue(
                    null_vectors["laurent_denominator_excludes_zero"]
                )

    def test_claim_boundary_flags_stay_false(self) -> None:
        scope = self.receipt["claim_scope"]
        self.assertTrue(scope["declared_chart_only"])
        self.assertFalse(scope["continuation_identity_independently_certified"])
        self.assertFalse(scope["standard_second_sheet_identification_certified"])
        self.assertFalse(scope["full_matrix_rank_n_minus_1_laurent_certified"])
        self.assertFalse(scope["matrix_rank_laurent_evidence_supplied"])
        self.assertFalse(scope["independent_numerical_replay_certified"])
        self.assertEqual(
            scope["engine_inverse_propagator_convention"],
            "G(s)=s-m_tree^2-Pi_engine(s)",
        )
        self.assertFalse(scope["theorem_self_energy_sign_bridge_certified"])
        self.assertFalse(scope["precision_ladder_receipt_supplied"])
        self.assertFalse(
            scope["independent_numerical_reevaluation_supplied"]
        )
        self.assertFalse(scope["full_physical_pole_evidence_supplied"])
        promotion = self.receipt["promotion"]
        self.assertFalse(
            promotion["continuation_identity_independently_certified"]
        )
        self.assertFalse(
            promotion["standard_second_sheet_identification_certified"]
        )
        self.assertFalse(
            promotion["theorem_self_energy_sign_bridge_certified"]
        )
        self.assertFalse(
            promotion["full_matrix_rank_n_minus_1_laurent_certified"]
        )
        self.assertFalse(promotion["matrix_rank_laurent_evidence_supplied"])
        self.assertFalse(
            promotion["precision_ladder_receipt_supplied"]
        )
        self.assertFalse(
            promotion["independent_numerical_reevaluation_supplied"]
        )
        self.assertFalse(promotion["full_physical_pole_evidence_supplied"])
        self.assertFalse(promotion["bmhv_restoration_certified"])
        self.assertFalse(promotion["physical_current_claim"])
        self.assertFalse(promotion["oph_native"])
        self.assertFalse(promotion["unit_claim"])
        self.assertFalse(promotion["unitarity_claim"])

    def test_declared_mixed_chart_and_nesting(self) -> None:
        for name in ("W", "Z"):
            continuation = self.receipt["declared_continuation"][name]
            diagnostics = continuation["consistency_diagnostics"]
            self.assertEqual(
                diagnostics["role"],
                "non_certifying_finite_delta_diagnostic",
            )
            self.assertFalse(diagnostics["gates_scalar_certificate"])
            nesting = self.receipt["precision_nesting"][name]
            self.assertTrue(nesting["enclosures_nested_with_precision"])
            self.assertTrue(
                nesting["per_segment_enclosure_nesting"]["all_nested"]
            )
            self.assertTrue(all(nesting["newton_ball_nesting"].values()))
            self.assertTrue(
                nesting["coefficient_denominator_nesting"]["all_nested"]
            )

        w_channels = {
            (channel["m1"], channel["m2"]): channel
            for channel in self.receipt["declared_continuation"]["W"][
                "channels"
            ]
        }
        self.assertEqual(w_channels[("0", "1/9")]["sheet_action"], "principal")
        self.assertEqual(
            w_channels[("0", "1/6400")]["crossing_correction"],
            "2*pi*i*(1-(1/6400)/s)",
        )
        self.assertEqual(
            w_channels[("1/9", "0")]["crossing_correction"],
            "2*pi*i*(1-(1/9)/s)",
        )

    def test_zero_exclusion_receipt_is_pinned(self) -> None:
        import hashlib

        stored = self.receipt["pins"][
            "principal_zero_exclusion_receipt_sha256"
        ]
        actual = "sha256:" + hashlib.sha256(
            (
                PACKAGE_ROOT / "outputs" / "certified_wz_contours.json"
            ).read_bytes()
        ).hexdigest()
        self.assertEqual(stored, actual)


class SecondSheetCheckerMutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = json.loads(
            (
                PACKAGE_ROOT / "outputs" / "certified_second_sheet_poles.json"
            ).read_text(encoding="utf-8")
        )

    def assert_rehashed_mutation_rejected(self, mutate) -> None:
        forged = copy.deepcopy(self.receipt)
        mutate(forged)
        _rehash(forged, second_checker)
        with self.assertRaises(second_checker.ReceiptValidationError):
            second_checker.validate_receipt(forged)

    def test_canonical_v2_receipt_passes_evidence_checker(self) -> None:
        self.assertIsNone(second_checker.validate_receipt(self.receipt))

    def test_physical_sheet_overclaim_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["claim_scope"].__setitem__(
                "standard_second_sheet_identification_certified", True
            )
        )

    def test_channel_sheet_action_mutation_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["declared_continuation"]["W"]["channels"][
                0
            ].__setitem__("sheet_action", "principal")
        )

    def test_finite_delta_diagnostic_promotion_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["declared_continuation"]["W"][
                "consistency_diagnostics"
            ].__setitem__("gates_scalar_certificate", True)
        )

    def test_forged_subdivision_depth_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"][
                "boundary_winding"
            ].__setitem__("max_depth_used", 1)
        )

    def test_deleted_quantity_enclosure_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            receipt["results"]["W"]["128"]["quantity_enclosures"].popitem()

        self.assert_rehashed_mutation_rejected(mutate)

    def test_forged_newton_margin_is_rejected(self) -> None:
        def mutate(receipt) -> None:
            receipt["results"]["W"]["128"]["interval_newton"][
                "strict_inclusion_margin"
            ]["minimum_margin"] = "999"

        self.assert_rehashed_mutation_rejected(mutate)

    def test_claim_bearing_reading_overclaim_is_rejected(self) -> None:
        self.assert_rehashed_mutation_rejected(
            lambda receipt: receipt["results"]["W"]["128"].__setitem__(
                "reading", "ISSUE 593 FULLY COMPLETE"
            )
        )


if __name__ == "__main__":
    unittest.main()
