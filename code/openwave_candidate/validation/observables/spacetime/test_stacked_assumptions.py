import importlib.util
from pathlib import Path
import unittest
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stacked_assumption_audit',HERE/'audit_stacked_assumptions.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

class StackedAssumptionTests(unittest.TestCase):
    def test_actual_zero_vertical_link_control(self):
        result=a.zero_links_control()
        self.assertEqual(result['zero_links_q3_energy_difference'],0)
        self.assertLess(result['zero_links_z_marginal_change'],1e-12)
        self.assertGreater(result['active_links_q3_energy_difference'],.1)

    def test_explicit_formula_input_domain(self):
        for args in ((True,.23,.2,[0,0,0],[0]),(4,0,.2,[0,0,0],[0]),
                     (4,-.23,.2,[0,0,0],[0]),(4,.23,0,[0,0,0],[0]),
                     (4,.23,.2,[0,0],[0]),(4,.23,.2,[0,0,np.nan],[0]),
                     (4,.23,.2,[0,0,0],[]),(4,.23,.2,[0,0,0],[np.inf])):
            with self.assertRaises(ValueError):a.declared_input_contract(*args)
        self.assertEqual(a.declared_input_contract(4,.23,.2,[0,0,0],[0])['amplitude_components'],64*111)

    def test_producer_rejects_unsupported_inputs_and_anchors_nonzero_time(self):
        probes=a.input_probes()
        self.assertTrue(all(probe['status']=='input_rejected' for probe in probes))
        result,_,_=a.m.simulate(n=4,times=[2.,2.37])
        samples=result['observables']['samples'];initial=np.array(samples[0]['centroid'])
        np.testing.assert_allclose(samples[0]['weighted_group_velocity_ray'],initial,atol=1e-14)
        np.testing.assert_allclose(samples[1]['weighted_group_velocity_ray'],initial+.37*np.array(result['observables']['weighted_group_velocity']),atol=1e-13)
        self.assertIn('finite-grid',result['supplied']['band_preparation'])
        self.assertIn('stored_spatial_apparatus',result)

    def test_static_offsets_do_not_create_single_band_centroid_motion(self):
        grid=a.m.band_grid(4,.23,True)
        rng=np.random.default_rng(91);amplitude=rng.normal(size=(4,4,4))+1j*rng.normal(size=(4,4,4))
        amplitude/=np.linalg.norm(amplitude.ravel())
        populations=[]
        for time in (0.,.37):
            psi=np.fft.ifftn(amplitude[...,None]*grid['spinor']*np.exp(-1j*time*grid['energy'])[...,None],axes=(0,1,2),norm='ortho')
            populations.append(np.sum(abs(psi)**2,axis=(0,1,2)))
        np.testing.assert_allclose(populations[0],populations[1],atol=1e-14)
        offsets=a.m.internal_offsets()
        np.testing.assert_allclose(populations[0]@offsets,populations[1]@offsets,atol=1e-14)
        self.assertEqual(len(np.unique(offsets,axis=0)),7)

if __name__=='__main__':unittest.main()
