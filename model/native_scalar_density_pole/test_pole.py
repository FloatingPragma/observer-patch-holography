import unittest
from fractions import Fraction as F
import numpy as np
from scipy import linalg as la
from model.defect_scalar_response.response import ground_pair
from .pole import build

class ScalarPoleTests(unittest.TestCase):
    def test_full_generalized_spectrum_and_both_outcomes(self):
        d=ground_pair();h=np.array(d['block']['H'],complex);g=np.array(d['block']['gram'],complex)
        n=np.array(d['n'],complex);energies,v=la.eigh(g@h,g)
        delta=1-2*energies[0];psi=v[:,0];population=float(np.real(psi.conj()@g@n@psi))
        total=0
        for operator,p in [(n,population),(np.eye(len(n))-n,1-population)]:
            weights=abs(v.conj().T@g@operator@psi)**2
            ground=np.abs(energies-energies[0])<1e-9
            excited=np.abs(energies-energies[0]-delta)<1e-9
            self.assertAlmostEqual(float(sum(weights[ground])),p*p,places=11)
            self.assertAlmostEqual(float(sum(weights[excited])),population*(1-population),places=11)
            self.assertLess(float(sum(weights[~(ground|excited)])),1e-20)
            total+=sum(weights)
        self.assertAlmostEqual(float(total),1,places=12)
        r=build();lo,hi=map(lambda x:float(F(x)),r['positive_pole_source_energy']['interval'])
        self.assertAlmostEqual(delta,(lo+hi)/2,places=10)

    def test_exact_outcome_intervals_and_scope(self):
        r=build();rows=r['complete_occupation_energy_records'];self.assertEqual(len(rows),6)
        bounds=[tuple(map(F,x['original_probability']['interval'])) for x in rows]
        self.assertLessEqual(sum(x[0] for x in bounds),1)
        self.assertGreaterEqual(sum(x[1] for x in bounds),1)
        for x in rows:
            if x['energy_outcome']=='outside':
                self.assertEqual(tuple(map(F,x['original_probability']['interval'])),(0,0))
                self.assertIsNone(x['source_energy_transfer'])
        self.assertIsNone(r['scope']['physical_rest_mass_GeV'])
        self.assertFalse(r['scope']['elementary_scalar_phi_poles_identical'])

if __name__=='__main__':unittest.main()
