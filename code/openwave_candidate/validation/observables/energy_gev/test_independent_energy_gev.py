import unittest
import review_energy_gev as audit

class IndependentEnergyAudit(unittest.TestCase):
    def test_exact_spectrum_and_complete_instruments(self):
        result=audit.run()
        self.assertEqual(result['outcome'],'PASS_WITH_SCOPE_NOTES')
        self.assertEqual(result['native_H_negative_m_multiplicity'],40)
        self.assertLess(max(result['errors'].values()),1e-10)
        self.assertFalse(result['full_physical_criterion_met'])

if __name__=='__main__':unittest.main()
