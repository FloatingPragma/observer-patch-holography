"""Independent controls for native operators, normalization and full record bounds."""
from fractions import Fraction as F
import unittest
import sympy as s

from .operators import operator_controls, reference_controls, reduced_EOM
from .bounds import atomic, response, interval_mul, interval_square
from .alignment import certificate as admitted_states
from model.native_quark_response.response import positive_response
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data


class JoinTests(unittest.TestCase):
    def test_full_original_and_reference_current_actions(self):
        controls = operator_controls()
        self.assertEqual(sum(r['actual_EOM_components'] for r in controls['family_controls']), 24)
        self.assertEqual(sum(r['gauge_component_checks'] for r in controls['family_controls']), 576)
        ref = reference_controls()
        self.assertEqual(ref['full_joint_CAR_modes'], 60)
        self.assertEqual(sum(r['all_joint_current_commutators'] for r in ref['families']), 102)

    def test_complete_actual_target_EOM_images(self):
        expected = {'uc': 26, 'dc': 26, 'Q': 52}
        for name, size in expected.items():
            row = reduced_EOM(name)
            self.assertEqual(row['retained_target_dimension'], size)
            self.assertEqual(row['full_weight_dimension'], 236)
            self.assertTrue(row['literal_original_CAR_Higgs_image_equals_Ht_A_minus_A_Hs'])

    def test_exact_source_moments_and_raw_normalization(self):
        _, _, field = field_data()
        e = field.unit
        for name in ('uc', 'dc', 'Q'):
            row = atomic(name)
            m0, m1, m2 = positive_response(name)['moments'][:3]
            self.assertEqual(m0, field.convert(s.Rational(1, 2)))
            self.assertEqual(row['totals'][1]-row['totals'][-1], (m1-e*m0)/2)
            self.assertEqual(row['totals'][1]+row['totals'][-1], (m0+(m2-2*e*m1+e*e*m0)/4)/2)
            self.assertEqual(sum(r['exactly_zero_residue'] for r in row['records']), 12 if name == 'Q' else 0)
        self.assertEqual(atomic('uc')['cross'], atomic('dc')['cross'])

    def test_all_six_complete_original_attempt_records_are_positive(self):
        for name in ('uc', 'dc', 'Q'):
            row = response(name)
            self.assertGreater(F(row['original_attempt_joint_contrast_lower']), F(96, 1000))
            self.assertGreater(F(row['exact_ground_normalized_composite_bin_probability_lower']), F(75, 1000))
            self.assertEqual(row['reference_success_cost'], '4')
            self.assertEqual(row['complete_successful_preparation_labels'], 8)
            for setting in row['settings']:
                self.assertEqual(setting['failed_pump_energy_hit_probability'], '0')
                for bounds in setting['original_attempt_joint_setting_records'].values():
                    lo, hi = map(F, bounds)
                    self.assertGreater(lo, 0)
                    self.assertLessEqual(lo, hi)
                    self.assertLessEqual(hi, F(1, 2))
                prepared = [tuple(map(F, x)) for x in setting['prepared_conditional_records'].values()]
                self.assertLessEqual(sum(p[0] for p in prepared), 1)
                self.assertGreaterEqual(sum(p[1] for p in prepared), 1)

    def test_interference_and_full_EOM_in_noncommuting_control(self):
        h = s.diag(0, 3)
        v = s.Matrix([[0, s.Rational(1, 7)], [s.Rational(1, 7), 0]])
        a = s.Matrix([[0, 0], [1, 0]])
        h0 = h-v
        c = h0*a-a*h0
        r = v*a-a*v
        self.assertEqual(h*a-a*h, c+r)
        g = s.Matrix([1, 0])
        p = s.diag(0, 1)
        self.assertEqual(p*c*g, h*p*a*g-p*r*g)
        plus, minus = (a+c/2)/2, (a-c/2)/2
        cross = (g.T*a.T*p*c*g)[0]
        self.assertEqual((g.T*(plus.T*p*plus-minus.T*p*minus)*g)[0], cross/2)
        # Fair-setting contrast is a half-difference of effects. Its range
        # can attain one, so trace distance delta is the sharp universal error.
        effect_difference = (s.diag(1, 0)-s.diag(0, 1))/2
        self.assertEqual(effect_difference[0, 0]-effect_difference[1, 1], 1)

    def test_interval_sign_and_zero_boundaries(self):
        self.assertEqual(interval_square((F(-2), F(3))), (F(0), F(9)))
        self.assertEqual(interval_square((F(-3), F(-2))), (F(4), F(9)))
        self.assertEqual(interval_mul((F(-2), F(1)), (F(-3), F(4))), (F(-8), F(6)))
        for name in ('uc', 'Q'):
            for row in atomic(name)['records']:
                if row['exactly_zero_residue']:
                    self.assertEqual(row['cross'], (0, 0))
                    self.assertEqual(row['composite'], (0, 0))
                    self.assertEqual(row['fringe'], {1: (0, 0), -1: (0, 0)})

    def test_actual_preparation_normalized_history_and_unequal_rates(self):
        value = admitted_states()
        self.assertTrue(value['original_two_source_preparation'])
        self.assertFalse(value['additional_P0'])
        for row in value['rows'].values():
            self.assertLess(F(row['common_history_total_variation_upper']), F(5193, 10**9))
            self.assertFalse(row['same_history_policy_may_read_pump_setting'])
            self.assertFalse(row['different_classical_setting_label_in_compared_state'])
            self.assertTrue(row['distinct_pump_and_bin_success_probabilities_certified'])
            self.assertGreater(F(row['original_attempt_elementary_in_bin_probability_lower']), F(49, 100))
            self.assertGreater(F(row['original_attempt_composite_B_over_2_in_bin_probability_lower']), F(75, 1000))
            for record in row['complete_setting_records'].values():
                law = [tuple(map(F, pair)) for pair in record['prepared_conditional'].values()]
                self.assertLessEqual(sum(x[0] for x in law), 1)
                self.assertGreaterEqual(sum(x[1] for x in law), 1)
                self.assertTrue(all(x[0] > 0 for x in law))

    def test_sharp_geometric_and_flagged_normalization_controls(self):
        a = s.Matrix([1, 0])
        b = s.Matrix([s.Rational(16, 5), s.Rational(12, 5)])
        self.assertEqual((b-5*a).dot(b-5*a), 9)
        self.assertEqual(1-(a.dot(b))**2/b.dot(b), s.Rational(3, 5)**2)
        for aa, bb in (((F(1, 5), F(1, 10)), (F(1, 10), F(1, 5))),
                       ((F(1, 5), F(1, 10)), (F(1, 10), F(1, 10)))):
            p, q = sum(aa), sum(bb)
            trace = sum(abs(x-y) for x, y in zip(aa, bb, strict=True))
            delta = (trace+abs(p-q))/2
            normalized = sum(abs(x/p-y/q) for x, y in zip(aa, bb, strict=True))/2
            self.assertLessEqual(normalized, delta/q)


if __name__ == '__main__':
    unittest.main()
