import unittest
import numpy as np
from . import scale

class ScaleTests(unittest.TestCase):
    def test_bare_critical_zero_and_tuned_cancellation(self):
        c={'kappa_c':.17,'kappa_phys':.169,'mrho_phys':.6,'mN_phys':.9,'sqrt_sigma':.4,'r0':3.,'mpi2_phys':.01}
        row=scale.coordinates(c,[.16,.165],{'rho':.8,'string':.44})
        self.assertEqual(row['a_mq/critical'],0.)
        self.assertAlmostEqual(row['bare_Wilson_GeV/rho/input_0.16'],(.8/.6)*.5*(1/.16-1/.17))
        self.assertAlmostEqual(row['inverse_a_GeV/rho']/row['inverse_a_GeV/string'],row['scale_ratio/rho/string'])

    def test_full_paired_chain_and_parent_values(self):
        r=scale.calculate();v=r['conditional_coordinates']['primary']
        self.assertEqual(len(v['replicas']),32)
        self.assertEqual(len(r['all_input_kappas']),4)
        self.assertGreater(r['tuned_mass_control']['fixed_kappa_c_countercontrol_statistical'],r['tuned_mass_control']['paired_statistical'])
        for row in r['paired_covariance_controls']:self.assertLess(row['covariance_identity_error'],1e-13)
        self.assertFalse(r['scope']['total_uncertainty_or_significance_claimed'])

    def test_common_calibration_rescaling(self):
        c={'kappa_c':.17,'kappa_phys':.169,'mrho_phys':.6,'mN_phys':.9,'sqrt_sigma':.4,'r0':3.,'mpi2_phys':.01}
        a=scale.coordinates(c,[.16],{'rho':.8,'string':.44});b=scale.coordinates(c,[.16],{'rho':1.6,'string':.88})
        self.assertEqual(a['scale_ratio/rho/string'],b['scale_ratio/rho/string'])
        self.assertEqual(2*a['bare_Wilson_GeV/rho/input_0.16'],b['bare_Wilson_GeV/rho/input_0.16'])
        self.assertEqual(a['a_mq/input_0.16'],b['a_mq/input_0.16'])

if __name__=='__main__':unittest.main()
