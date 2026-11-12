import unittest
import numpy as np
from fractions import Fraction as F
from .source import shape,geometry,field,block
from .certificate import spectral_certificate,observations,wave_square_roots
class DynamicWaveTests(unittest.TestCase):
 def test_harmonic_shape_and_complete_positive_source(self):
  s=shape();g=geometry();self.assertEqual(g['dim'],4)
  self.assertLess(np.linalg.norm(g['Q'][0]@g['Q'][1]-g['Q'][1]@g['Q'][0]),1e-17)
  self.assertTrue(np.allclose(s['F'][s['source_sites']],np.eye(2)))
  self.assertGreater(np.linalg.eigvalsh(field()['H'])[0],.25)
  self.assertEqual(np.linalg.eigvalsh(block(left=0)['vacuum'])[0],0)
 def test_exact_positive_roots_and_complete_arrival(self):
  roots,_=wave_square_roots()
  self.assertTrue(all(F(r['positive_root_Gershgorin_lower'])>0 and F(r['exact_positive_root_operator_error_upper'])<F(1,10**9)for r in roots))
  c=spectral_certificate(1,0);r=observations(c)
  for t in r['single_time_records']:
   self.assertAlmostEqual(sum(t['all_vacuum_geometry_number_probabilities'])+np.sum(t['all_wave_site_geometry_number_probabilities']),1,places=8)
   self.assertAlmostEqual(t['source_energy_before'],t['source_energy_after'],places=8)
  a=r['arrival_instrument'];total=sum(x['unconditional_hit_probability']for x in a['all_arrival_records'])+sum(a['all_vacuum_timeout_geometry_probabilities'])+np.sum(a['complete_wave_timeout_geometry_site_probabilities'])
  self.assertAlmostEqual(total,1,places=8)
if __name__=='__main__':unittest.main()
