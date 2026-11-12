import unittest
from fractions import Fraction as F
import numpy as np
from .clock import phase_records,phase_bounds,instrument_certificate

class NativeClockTest(unittest.TestCase):
 def test_complete_clock_instrument_retains_all_native_levels(self):
  r=instrument_certificate();self.assertEqual(r['setting_and_native_outcomes'],16);self.assertEqual(r['exact_effect_sum'],'I_8')
  v=[np.array([1.,0.],complex),np.array([1j,0.],complex)];p,total,z=phase_records(v);self.assertAlmostEqual(total,1);self.assertAlmostEqual(p.sum(),1)
  inferred=2*((p[0].sum()-p[2].sum())+1j*(p[1].sum()-p[3].sum()));self.assertAlmostEqual(abs(inferred-z),0)
 def test_phase_disk_contains_independently_perturbed_phases(self):
  z=.3+.4j;e=.01;lo,hi=map(float,phase_bounds(z,F(e)))
  for angle in np.linspace(0,2*np.pi,41):
   phase=np.angle(z+e*np.exp(1j*angle));self.assertLessEqual(lo,phase);self.assertGreaterEqual(hi,phase)
 def test_common_clock_origin_cancels_from_relative_phase(self):
  z0=.08+.02j;z1=.07+.03j
  for phase in(.4,1.7,2.8):
   w=np.exp(1j*phase);self.assertAlmostEqual(abs((w*z1)*(w*z0).conjugate()-z1*z0.conjugate()),0)
if __name__=='__main__':unittest.main()
