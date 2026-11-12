import unittest
import numpy as np
from fractions import Fraction as F
from .source import geometry,clock_state,partial_transpose,independent_geometry_control,complete_records,block
from .certificate import block_certificate,witness
class QuantumGeometryTests(unittest.TestCase):
 def test_complete_boundary_and_commuting_geometry(self):
  for cap in(1,2,4):
   g=geometry(cap);self.assertEqual(len(g['basis']),(cap+1)**2)
   self.assertLess(np.linalg.norm(g['Q'][0]@g['Q'][1]-g['Q'][1]@g['Q'][0]),1e-17)
   self.assertEqual(np.linalg.eigvalsh(block(cap,0,0))[0],0.)
   self.assertGreater(np.linalg.eigvalsh(block(cap,0,1))[0],0.)
 def test_certificate_and_complete_records(self):
  rows=[block_certificate(1,i,j)for i in range(8)for j in range(8)];w=witness(1,10000,rows)
  self.assertLess(F(w['exact_source_partial_transpose_witness_interval'][1]),0)
  rec=complete_records(1,10000);self.assertAlmostEqual(rec['source_energy_before'],rec['source_energy_after'],places=12)
  for r,c in zip(rec['all9_complete_joint_instruments'],w['complete_certified_joint_instruments']):
   self.assertEqual(r['complete_outcome_count'],256)
   for row,intervalrow in zip(r['active_four_clock_lines_times_all_geometry_occupations'],c['all_active_clock_geometry_probability_intervals']):
    for x,(lo,hi)in zip(row,intervalrow):self.assertLessEqual(lo,x);self.assertLessEqual(x,hi)
 def test_no_spatial_edge_and_initial_product(self):
  for cap in(1,4):
   self.assertGreater(min(independent_geometry_control(cap,10000)['partial_transpose_eigenvalues']),-1e-14)
   self.assertGreater(np.linalg.eigvalsh(partial_transpose(clock_state(cap,0)[1]))[0],-1e-14)
if __name__=='__main__':unittest.main()
