import unittest
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from .source import compression,initial,exact_geometry,native
from .propagation import advance,step_error,coefficients
from .spectrum import wave_schur

class GeometryRefinementTest(unittest.TestCase):
 def test_complete_chebyshev_action_against_independent_eigensystem(self):
  source=compression(1,0,1);v=initial(1);e,V=eigh(source['H']);reference=V@(np.exp(-.5j*e)*(V.T@v));actual=advance(source,v)
  self.assertLess(np.linalg.norm(actual-reference),float(step_error(source)))
  self.assertLess(coefficients()['tail'],F(1,10**30))
 def test_schur_isolates_the_embedded_trial(self):
  r=wave_schur(1,2);lo,hi=map(F,r['master_ground_interval_rational']);self.assertLess(lo,hi);self.assertLess(hi,F(r['master_second_eigenvalue_lower']))
  self.assertGreater(F(r['exact_dyadic_residual_squared']),0)
 def test_reexponentiation_has_an_actual_missing_walk(self):
  x1=np.array([[0.,1.],[1.,0.]]);x2=np.array([[0.,1.,0.],[1.,0.,np.sqrt(2)],[0.,np.sqrt(2),0.]])
  self.assertAlmostEqual(np.linalg.matrix_power(x1,4)[0,0],1)
  self.assertAlmostEqual(np.linalg.matrix_power(x2,4)[0,0],3)
  self.assertGreater(native()[2][1].lo,0)
if __name__=='__main__':unittest.main()
