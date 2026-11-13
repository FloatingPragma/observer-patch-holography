from fractions import Fraction as F
import unittest
import numpy as np
from . import source as m


class MovingClockTest(unittest.TestCase):
    def test_inverse_preserves_unprepared_event_and_reference_amplitudes(self):
        g=m.geometry()[0]
        rng=np.random.default_rng(771)
        state=rng.normal(size=(g.states,g.events,2))+1j*rng.normal(size=(g.states,g.events,2))
        state/=np.linalg.norm(state)
        for k in (0,1,3,5):
            returned=m.transport(m.transport(state,k,3),k,-3)
            self.assertLess(np.linalg.norm(returned-state),2e-12)
            self.assertLess(np.max(abs(returned.reshape(-1,2).conj().T@returned.reshape(-1,2)-state.reshape(-1,2).conj().T@state.reshape(-1,2))),2e-12)

    def test_fixed_code_really_loses_weight_without_tracking(self):
        g,_,logical,_,_=m.geometry()
        wave=m.frame(1,1)@logical[:,0]
        projected=m.motion.matching_code_projection(g,wave)
        self.assertAlmostEqual(float(np.sum(abs(projected)**2)),float(F(67,90)),places=12)
        self.assertAlmostEqual(float(np.sum(abs(wave)**2)),1.,places=12)

    def test_complete_physical_quadrature_includes_the_complement(self):
        for axis in ('X','Y'):
            actual=m.physical_effects(1,4,0,axis)
            expected=m.clock.effects(0,axis)
            self.assertLess(max(np.max(abs(a-b))for a,b in zip(actual,expected)),2e-12)
            self.assertLess(np.max(abs(sum(actual)-np.eye(36))),2e-12)
            self.assertGreater(np.trace(actual[2]).real,33.9)

    def test_invalid_character_is_rejected(self):
        with self.assertRaises(ValueError):m.step(m.frame(0,0),True)
        with self.assertRaises(ValueError):m.step(m.frame(0,0),6)
