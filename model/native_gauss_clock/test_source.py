from fractions import Fraction as F
import unittest
import numpy as np
from . import source as c


class GaussClockTest(unittest.TestCase):
    def test_all_reference_columns_return_without_link_reset(self):
        for n in (0,1,2):
            keys,a,b,den,matrix,inverse=c.transported_columns(n)
            self.assertEqual(len(inverse),12)
            self.assertTrue(all(inverse))
            self.assertLess(np.max(abs(matrix.conj().T@matrix-np.eye(12))),2e-12)
            self.assertTrue(all(c.joint.electric.gauss(c.joint.state(key[0],1))==(0,)*12 for key in keys))

    def test_untracked_canonical_decoder_retains_large_complement(self):
        wave=c.joint.prepare()
        for _ in range(2):wave=c.joint.step(wave,electric_phase=True)
        decoded=c.joint.decode_canonical_dressing(wave)
        self.assertEqual(F(decoded['canonical_dressing_probability']),F(48179,81000))
        self.assertGreater(F(decoded['complete_outside_dressing_probability']),F(2,5))
        self.assertFalse(decoded['reset_performed'])

    def test_complete_records_accept_mixed_source_states_without_reweighting(self):
        rho=np.diag(np.arange(1,13)/78)
        result=c.complete_joint_observation(2,rho)
        self.assertAlmostEqual(result['full_probability_sum_diagnostic'],1,places=12)
        self.assertAlmostEqual(result['complete_Gauss_zero_probability_diagnostic'],1,places=12)
        self.assertFalse(result['normalization_performed'])

    def test_frame_domain_and_source_operator_shape(self):
        with self.assertRaises(ValueError):c.transported_columns(3)
        with self.assertRaises(ValueError):c.complete_joint_observation(1,np.eye(6))
