import unittest
import numpy as np
from fractions import Fraction as F
from .response import operators,clocks
class StressTests(unittest.TestCase):
 def test_source_derivatives_and_finite_boundary(self):
  r=operators()
  self.assertLess(max(r['sylvester_residuals']),1e-10)
  self.assertLess(r['velocity_residual'],1e-10)
  for a in range(2):
   residual=r['actual_momentum_rate'][a]-r['spring'][a]+r['native_force'][a]+r['field_force'][a]-r['finite_boundary_correction'][a]
   self.assertLess(np.linalg.norm(residual),1e-10)
   self.assertGreater(np.linalg.norm(r['finite_boundary_correction'][a]),.01)
 def test_certified_clock_nonuniversality_and_domain(self):
  c=clocks();self.assertEqual(len(c['all64_native_vacuum_block_certificates']),64)
  self.assertGreater(F(c['strict_clock_ratio_spread_interval'][0]),0)
  self.assertTrue(all(r['same_source_clock_shift_ppm_interval'][1]<0 for r in c['all7_same_source_clock_ratios']))
  for row in c['complete_dressed_clock_records']:self.assertAlmostEqual(np.sum(row['all_X_clock_geometry_number_probabilities']),1,places=12)
  self.assertLess(F(c['exact_test_clock_theorem']['at_probe_coefficient_one_thousandth_error_upper']),F(26,10**6))
if __name__=='__main__':unittest.main()
