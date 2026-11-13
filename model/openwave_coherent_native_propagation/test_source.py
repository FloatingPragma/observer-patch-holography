import unittest
from fractions import Fraction as F
from . import source as s

class PropagationTests(unittest.TestCase):
    def test_complete_arc_control_and_wrong_continuation(self):
        x=s.arc_control()
        self.assertNotEqual(x['reflection_two_step_origin_probability'],'1')
        self.assertEqual(x['identity_continuation_two_step_origin_probability'],'1')
        self.assertEqual(x['matrix_unit_terminal_port_effect_checks'],8640)
        self.assertEqual(x['full_spectrum_census']['full_W_minus_one_multiplicity'],678)
    def test_one_step_is_native_population(self):
        g=s.small_graph();D=2*len(g.seams);a=[int(i==0)for i in range(12)];b=[0]*12
        aa,bb=s.cyclic_step(g,a,b);num=s.port_numerators(g,aa,bb)
        self.assertEqual([F(x,D**3)for x in num],[F(x,D)for x in s.integer_transfer(g,a)])
    def test_compressed_amplitude_is_not_complete_probability(self):
        g=s.small_graph();D=2*len(g.seams);v=[int(i==0)for i in range(12)]
        t=s.integer_transfer(g,v)
        self.assertLess(sum(F(x*x,D*D)for x in t),1)
        a,b=s.cyclic_step(g,v,[0]*12)
        self.assertEqual(sum(s.port_numerators(g,a,b)),D**3)
    def test_uniform_mode_and_interval_direction(self):
        g=s.small_graph();a=[1]*12;b=[0]*12;D=2*len(g.seams);p=1
        for _ in range(5):a,b=s.cyclic_step(g,a,b);p*=D
        self.assertEqual([F(x,D*p*p*12)for x in s.port_numerators(g,a,b)],[F(1,12)]*12)
        for x in (F(1,3),F(-1,3),F(1),F()):
            lo,hi=map(F,s.interval(x));self.assertLessEqual(lo,x);self.assertGreaterEqual(hi,x)
    def test_held_out_width_all_outcomes(self):
        x=s.packet(s.small_graph(),2,(0,1,4))
        self.assertTrue(x['held_out_preparation_width'])
        self.assertEqual(len(x['records']),3)
        self.assertTrue(all(z['total_probability_exact']=='1'for z in x['records']))
