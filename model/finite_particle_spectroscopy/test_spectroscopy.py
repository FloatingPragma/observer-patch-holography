import unittest
from fractions import Fraction as F
import mpmath as mp
from . import spectroscopy as s

class SpectroscopyTests(unittest.TestCase):
    def test_grid_pointer_and_non_grid_completeness(self):
        with mp.workdps(100):
            self.assertEqual(s.kernel(mp.mpf(1)/4),[mp.mpc(int(y==8)) for y in range(32)])
            a=s.kernel(mp.mpf('0.321'))
            self.assertLess(abs(sum(abs(x)**2 for x in a)-1),mp.mpf('1e-90'))
            tail=sum(abs(x)**2 for y,x in enumerate(a) if abs(mp.mpf(y)/32-mp.mpf('.321'))>mp.mpf(6)/32)
            self.assertLess(tail,s.number(s.TAIL))

    def test_complete_adaptive_operator_words(self):
        with mp.workdps(100):
            for bits in (1,2,3,5):
                for phase in (mp.mpf(1)/4,mp.mpf('.321'),mp.mpf(1)/4+mp.mpf('1e-24')):
                    direct=s.kernel(phase,1<<bits);adaptive=s.adaptive_kernel(phase,bits)
                    self.assertLess(max(abs(x-y) for x,y in zip(direct,adaptive)),mp.mpf('1e-90'))
            # Complete complex-amplitude equality preserves arbitrary reference coherences.
            left=s.adaptive_kernel(mp.mpf('.321'));right=s.adaptive_kernel(mp.mpf('.407'))
            direct_left=s.kernel(mp.mpf('.321'));direct_right=s.kernel(mp.mpf('.407'))
            self.assertLess(max(abs(a*mp.conj(b)-c*mp.conj(d)) for a,b,c,d in zip(left,right,direct_left,direct_right)),mp.mpf('1e-90'))
        with self.assertRaises(ValueError):s.adaptive_kernel(mp.mpf('.3'),0)

    def test_inconclusive_bins_are_retained(self):
        self.assertIsNone(s.confidence_interval(0,F(8)))
        v=s.confidence_interval(8,F(8));self.assertEqual(F(v['energy_GeV_rational'][0]),1)
        self.assertEqual(F(v['mass_GeV_rational'][0]),0)
        with self.assertRaises(ValueError):s.confidence_interval(1,F(8),16)

    def test_neutrino_stable_mass_frame(self):
        with mp.workdps(100):
            square=F(1,10**24);phase=s.phase_for(square,F(1))
            self.assertGreater(phase,mp.mpf(1)/4)
            self.assertAlmostEqual(float((phase-mp.mpf(1)/4)/s.number(square)),.5)
            self.assertEqual(float(mp.sqrt(1+s.number(square))),1.)
            a=s.kernel(phase);self.assertGreater(abs(a[9]),0)

    def test_confidence_mass_interval_contains_grid_source(self):
        # phi=3/8 at R=8 gives E=2 and m=sqrt(3).
        v=s.confidence_interval(12,F(8));lo,hi=map(F,v['mass_GeV_rational'])
        self.assertLess(lo*lo,3);self.assertGreater(hi*hi,3)

if __name__=='__main__':unittest.main()
