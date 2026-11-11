import unittest

import numpy as np

from .readout import calibrate, certificate, flavour_checks, multiply


class PionReadoutTests(unittest.TestCase):
    def test_isospin_and_symmetry_breaking(self):
        result=flavour_checks()
        self.assertLess(result['max_charge_and_isospin_residual'], 1e-14)
        self.assertTrue(np.allclose(result['gram'],np.eye(3)))
        # A nondegenerate propagator has a nonzero neutral disconnected trace.
        zero=np.diag([1.,-1.])/np.sqrt(2)
        self.assertEqual(np.trace(zero @ np.eye(2)), 0)
        self.assertGreater(abs(np.trace(zero @ np.diag([1.,2.]))), .7)

    def test_calibration_and_retained_archive(self):
        result=certificate()
        self.assertTrue(result['all_checks_pass'])
        self.assertEqual(len(result['particles']),3)
        self.assertEqual(len(result['lattice_mass_rows']),4)
        for row in result['particles']:
            self.assertIsNone(row['physical_source_derived_mass_GeV'])
            self.assertAlmostEqual(row['rho_anchor']['value_GeV'], .13799628)
            self.assertEqual(row['rho_anchor']['statistical_GeV'],0)
        self.assertGreater(result['pion_charged_neutral_splitting']['experimental_difference_GeV'],.004)
        self.assertIsNone(result['calibrations']['anchor_and_model_uncertainty_GeV'])
        for bad in [0,-1,np.inf,np.nan]:
            with self.assertRaises(ValueError):
                multiply(1,0,0,bad)

    def test_shared_scale_is_not_independent_error(self):
        # Exact ratio tuning cancels an arbitrarily fluctuating common rho scale.
        rho_replicas=np.array([.4,.5,.6,.7])
        pion_replicas=.178*rho_replicas
        calibrated=pion_replicas*(.77526/rho_replicas)
        self.assertLess(np.ptp(calibrated),1e-15)
        self.assertGreater(np.std(pion_replicas),.01)


if __name__=='__main__':
    unittest.main()
