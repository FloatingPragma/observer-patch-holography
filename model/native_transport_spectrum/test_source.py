from fractions import Fraction as F
import math
import unittest
import numpy as np
from . import source as m


class NativeSpectrumTest(unittest.TestCase):
    def test_full_matrix_spectral_values_lie_in_exact_isolating_intervals(self):
        for k in (1,2,3):
            lo,hi,mult,diagnostic=m.smallest_positive_interval(k)
            self.assertGreater(lo,0)
            self.assertLessEqual(float(lo)-2e-12,diagnostic)
            self.assertGreaterEqual(float(hi)+2e-12,diagnostic)
            self.assertEqual(mult,3 if k==3 else 1)

    def test_direct_integer_parents_match_the_complete_native_compression(self):
        g=m.motion.carrier()
        for k in (0,3):
            actual=g.events*(np.eye(g.states)-m.motion.apply_T(g,np.eye(g.states),k))
            self.assertLess(np.max(abs(actual-m.integer_laplacian(k))),1e-13)

    def test_independent_prime_field_characteristic_polynomials(self):
        self.assertEqual(len(m.finite_field_control()),18)

    def test_phase_interval_and_zero_branch(self):
        self.assertEqual(m.arccos_interval(F(1)),(F(0),F(0)))
        lo,hi=m.arccos_interval(F(99,100))
        self.assertLessEqual(float(lo),math.acos(.99)+2e-15)
        self.assertGreaterEqual(float(hi),math.acos(.99)-2e-15)
        with self.assertRaises(ValueError):m.arccos_interval(F(-1))

    def test_fixed_reference_and_gap_do_not_gain_a_uniform_limit(self):
        result=m.enlargement_controls(F(1))
        first,last=result['rows'][0],result['rows'][-1]
        self.assertLess(F(last['prepared_cyclic_frequency_squared_upper_exact']),F(first['prepared_cyclic_frequency_squared_upper_exact'])/200)
        self.assertFalse(result['four_dimensional_continuum_or_Yang_Mills_limit_proved'])
