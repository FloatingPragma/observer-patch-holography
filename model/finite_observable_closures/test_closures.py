"""Independent complete-law, no-signalling and arithmetic controls."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy import linalg
import sympy as s

from . import coherence,counting,measurement
from model.observer_time_fields import federation
from model.common_source import finite_signal,scalar_coupling


class ObservableClosures(unittest.TestCase):
    def test_joint_interference_and_untouched_reference(self):
        r=coherence.calculate()
        self.assertEqual(r['completed_probability_difference'],'101/968')
        self.assertFalse(r['reference_alone_distinguishes_C_and_E'])
        for domain in ('completed_joint_records','capped_joint_records'):
            for rows in r[domain].values():
                self.assertEqual(len(rows),42)
                self.assertEqual(sum(F(x['probability']) for x in rows),1)
        for run in r['exact_small_native_executions']:
            self.assertEqual(run['reference_marginal'],[['1/2','0'],['0','1/2']])

    def test_exact_full_sector_law_dominates_every_rounded_native_entry(self):
        for feedback in (False,True):
            exact,record=counting.exact_law(feedback)
            actual=federation.execute(cap=1,ticks=4,grid_bits=51,feedback=feedback)
            scale=1 << 51; denominator=int(record['common_denominator'])
            self.assertEqual(record['total_count'],record['common_denominator'])
            self.assertEqual(record['record_count'],124)
            for column,key in enumerate(actual.records):
                for row in np.flatnonzero(actual.counts[:,column]):
                    self.assertLessEqual(int(actual.counts[row,column])*denominator,int(exact[key][row])*scale)
            self.assertEqual(sum(F(x['probability']) for x in record['records']),1)
            self.assertEqual(F(record['success_probability']),F(65,630)**4)
            altered=sum(int(x['count']) for x in record['records'])-1
            self.assertNotEqual(altered,denominator)

    def test_exact_finite_generator_and_both_pointer_completeness(self):
        exact=finite_signal.hamiltonian_matrix()
        actual=np.array([[float(a)+1j*float(b) for a,b in row] for row in exact])
        np.testing.assert_allclose(actual,scalar_coupling.ideal_hamiltonian(2,.25),atol=2e-15,rtol=2e-15)
        x0=s.zeros(16);x2=s.zeros(16)
        for i in range(16):x0[i^1,i]=1;x2[i^4,i]=1
        feedback=(3*s.eye(16)-4*s.I*x0)/5
        self.assertEqual(feedback.H*feedback,s.eye(16))
        first=[s.diag(*(int((i&1)==b) for i in range(16))) for b in (0,1)]
        second=[(s.eye(16)+sign*x2)/2 for sign in (-1,1)]
        for projectors in (first,second):
            self.assertEqual(sum(projectors,s.zeros(16)),s.eye(16))
            self.assertTrue(all(p*p==p and p.H==p for p in projectors))
        r=measurement.calculate()
        self.assertLess(F(r['continuous_contrast_from_fresh_execution_interval'][1]),0)
        self.assertEqual(len(r['fresh_source_executions']),2)
        self.assertTrue(all(len(x['records'])==4 for x in r['fresh_source_executions']))

    def test_binary_density_certificate_rejects_corruption_and_ignores_phase(self):
        polynomial=((F(1),F(0)),(F(0),F(0)))
        error,mass=measurement.density_distance_upper([[0.,1.],[0.,0.]],polynomial,F(0))
        self.assertEqual(error,0);self.assertEqual(mass,1)
        error,_=measurement.density_distance_upper([[1.,0.],[.1,0.]],polynomial,F(0))
        self.assertGreater(error,F(1,20))
        for x in (F(0),F(1,7),F(2),F(10**-6)):
            bound=measurement.sqrt_upper(x)
            self.assertGreaterEqual(bound*bound,x)


if __name__=='__main__':unittest.main()
