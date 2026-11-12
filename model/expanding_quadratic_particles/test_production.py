from fractions import Fraction as F
import math
import unittest
import numpy as np
from scipy.linalg import expm
from . import production as p


class ExpansionTests(unittest.TestCase):
    def test_literal_source_phases_and_all_modes(self):
        d=p.source();self.assertEqual(len(d['pairs']),90)
        self.assertEqual(set(sum((x['mode_indices'] for x in d['pairs']),[])),set(range(180)))
        self.assertEqual({x['canonical_pairing_sign'] for x in d['pairs']},{-1,1})
        self.assertLess(d['controls'][0]['canonical_named_frame_phase_error'],1e-12)

    def test_sudden_quench_and_massless_control(self):
        m=.3;e=p.pair_evolution(['3/10','3/10'],0,0)
        wanted=(1-(m*m+.25)/(math.hypot(1,m)*math.hypot(.25,m)))/2
        self.assertAlmostEqual(p.mid(e['probability']),wanted,14)
        self.assertEqual(p.pair_evolution(['0','0'],F(8),64)['probability'],p.ZERO)

    def test_full_even_state_against_literal_two_mode_fock(self):
        m=.3;duration=F(1);steps=64;e=p.pair_evolution(['3/10','3/10'],duration,steps)
        for sign in (-1,1):
            def frame(momentum):
                en=math.hypot(momentum,m);c=math.sqrt((en+momentum)/(2*en));s=m/math.sqrt(2*en*(en+momentum))
                return np.array([[c,sign*s],[-sign*s,c]])
            u=np.eye(2,dtype=complex)
            for k in range(steps):
                mom=float(1/p.ramp_scale(F(2*k+1,2*steps)))
                u=expm(-1j*float(duration/steps)*np.array([[-mom,sign*m],[sign*m,mom]]))@u
            g=frame(.25).T@u@frame(1)
            x=np.array([[0,1],[1,0]]);d=np.diag([1,sign]);expected=d@x@e['numeric']@x@d
            self.assertLess(np.linalg.norm(g-expected),2e-13)
            for temperature in (None,F(1)):
                state=p.pair_state(e,temperature,sign);f=0 if temperature is None else 1/(1+math.exp(math.hypot(1,m)))
                density=g@np.diag([(1-f)**2,f*f])@g.conj().T
                c=state['even_vacuum_pair_coherence'];coherence=complex(sum(map(F,c['real']))/2,sum(map(F,c['imaginary']))/2)
                self.assertLess(abs(density[0,1]-coherence),2e-13)

    def test_complete_number_law_and_ramp_scope(self):
        result=p.run(0,0)
        for prep in result['preparations']:
            law=prep['complete_number_readout'];lo,hi=map(F,law['normalization_interval'])
            self.assertLessEqual(lo,1);self.assertGreaterEqual(hi,1)
            self.assertEqual(len(law['all181_total_number_probability_intervals']),181)
            self.assertEqual(len(law['complete84_group_count_marginals']),84)
            self.assertEqual(len(law['complete12_species_number_laws']),12)
        coarse=F(p.smooth_discretization_bound('1',F(8),64))
        fine=F(p.smooth_discretization_bound('1',F(8),256))
        self.assertLess(fine,coarse/15)

    def test_high_resolution_no_interval_wrapping(self):
        result=p.run(F(8),1024)
        top=next(x for x in result['species_production'] if x['species']=='top')
        self.assertLess(float(F(top['certified_Bogoliubov_matrix_operator_error'])),1e-20)
        for prep in result['preparations']:
            law=prep['complete_number_readout']
            self.assertGreaterEqual(law['expected_total_number'],0)
            self.assertLessEqual(law['expected_total_number'],180)
            self.assertLess(law['point_law_normalization_error'],1e-12)
            self.assertLess(law['maximum_total_number_probability_interval_width'],1e-20)
            bound=prep['continuous_ramp_mean_number_bounds']
            self.assertEqual(len(bound['species']),12)
            self.assertLess(float(F(bound['total_absolute_error_bound'])),.2)
        self.assertLess(float(F(result['preparations'][0]['continuous_ramp_mean_number_bounds']['total_exact_mean_interval'][1])),.23)

if __name__=='__main__':unittest.main()
