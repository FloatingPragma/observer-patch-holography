import unittest
import numpy as np
from scipy.linalg import eigh
from .waves import matrix,apply_function,tail_bound,rays
class WaveTests(unittest.TestCase):
 def test_complete_function_against_dense(self):
  n=7;q=np.zeros((n,n));q[3,3]=-.0001;K=matrix(q);E,O=eigh(K.toarray());v=np.arange(n*n)+1j*np.arange(n*n)[::-1];v=v/np.linalg.norm(v)
  expected=O@(np.exp(-8j*np.sqrt(E))*(O.T@v));actual=apply_function(K,v,n,8)
  self.assertLess(np.linalg.norm(expected-actual),2e-12)
  self.assertLess(float(tail_bound(31,24)),2e-21)
 def test_gamma_and_velocity_instrument(self):
  n=7;q=np.zeros((n,n));self.assertEqual((matrix(q,-1)!=matrix(q,1)).nnz,0)
  q[3,3]=-.0001;K=matrix(q);E,O=eigh(K.toarray());H=(O*np.sqrt(E))@O.T;Y=np.diag(np.tile(np.arange(n),n));V=1j*(H@Y-Y@H);B=3*(n-1)
  for sign in(-1,1):self.assertGreater(np.linalg.eigvalsh((np.eye(n*n)+sign*V/B)/2)[0],0)
  self.assertLess(E[-1],9)
 def test_principal_ray_control(self):
  q=np.zeros((31,31));q[15,15]=-.0001
  for x in rays(q,-1,16):self.assertEqual(x['first_order_transverse_velocity_shift'],0);self.assertEqual(x['first_order_transverse_position_shift'],0)
if __name__=='__main__':unittest.main()
