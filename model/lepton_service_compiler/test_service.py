"""Independent phase, reference, padding, archive and service controls."""
import unittest
import numpy as np

from . import service as s


class ServiceTests(unittest.TestCase):
    def test_blank_phase_complete_word_and_global_phase(self):
        for bits in (1,3,6):
            word=s.zero_phase_word(bits,-.37)
            actual=s.source.execute(word,bits)
            target=np.eye(1<<bits,dtype=complex);target[0,0]=np.exp(-.37j)
            self.assertLess(np.linalg.norm(actual-target),4e-11)

    def test_rank_one_identity_with_untouched_reference(self):
        rng=np.random.default_rng(193)
        raw=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16));v=np.linalg.qr(raw)[0]
        state=rng.normal(size=(16,7))+1j*rng.normal(size=(16,7))
        z=np.eye(16,dtype=complex);z[0,0]=np.exp(-.4j)
        self.assertLess(np.linalg.norm(s.rank_kick(state,v[:,0],.4)-v@z@v.conj().T@state),3e-14)

    def test_archive_factors_retain_within_return_coherence(self):
        local={r:np.zeros((4096,4,3),complex) for r in s.instrument.LABELS}
        local['yes_yes'][0,0,0]=1/np.sqrt(2);local['yes_yes'][1,0,0]=1j/np.sqrt(2)
        local['x_no'][7,1,1]=1;local['y_no'][9,2,2]=1
        rows,factors=s.history_products(local,3)
        self.assertEqual(len(rows),64)
        self.assertAlmostEqual(sum(r['probability_maximally_entangled_register_reference'] for r in rows),1)
        self.assertAlmostEqual(abs(factors['yes_yes'][0,0]*factors['yes_yes'][1,0].conjugate()),.5)
        with self.assertRaises(ValueError):s.history_products(local,4)

    def test_complete_coherent_service_countercontrol(self):
        values=s.service_controls()
        self.assertAlmostEqual(values['2']['dephasing_negative_control_unhalved_distance'],1)
        self.assertAlmostEqual(values['4']['dephasing_negative_control_unhalved_distance'],1.5)


if __name__=='__main__':unittest.main()
