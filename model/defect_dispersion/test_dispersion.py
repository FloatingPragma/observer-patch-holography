import unittest
import numpy as np
from scipy import linalg as la
from .dispersion import bloch, homogeneous_spectrum, two_cell_fourier_control, compare

class DispersionTests(unittest.TestCase):
    def test_independent_complete_fourier_graph(self):
        for axis in range(3):
            self.assertLess(two_cell_fourier_control(axis=axis)['all_row_intertwiner_error'], 1e-10)

    def test_homogeneous_all_rows_and_twist(self):
        t=np.array([.2,-.3,.4]); h=bloch(6,8,t,uniform=True)
        self.assertLess(la.norm(h-h.conj().T),1e-12)
        np.testing.assert_allclose(la.eigvalsh(h),homogeneous_spectrum(6,8,t),atol=1e-11,rtol=0)
        self.assertGreater(la.norm(h-bloch(6,8,[0,0,0],uniform=True)), .1)

    def test_band_count_and_negative_secant(self):
        e0=la.eigvalsh(bloch(6,8,[0,0,0])); e1=la.eigvalsh(bloch(6,8,[0,0,np.pi/2]))
        for e in [e0,e1]:
            self.assertEqual(sum(abs(e)<3),4)
            self.assertGreaterEqual(min(abs(e[abs(e)>=3])),5)
        self.assertTrue(np.all(e1[abs(e1)<3][2:] < e0[abs(e0)<3][2:]-1e-5))

    def test_inputs_and_receipt_mutation(self):
        for side,m,t in [(5,8,[0,0,0]),(6,7,[0,0,0]),(6,8,[1j,0,0]),(6,8,[np.nan,0,0])]:
            with self.assertRaises(ValueError):bloch(side,m,t)
        self.assertFalse(compare({'mass':False,'roots':[1.]},{'mass':True,'roots':[1.]}))
        self.assertFalse(compare({'roots':[1.]},{'roots':[1.01]}))

if __name__=='__main__':unittest.main()
