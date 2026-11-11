"""Independent dense spectral controls and complete weight tests."""
from fractions import Fraction as F
import unittest
import numpy as np
from scipy.linalg import eigh

from .resolution import calculate, polarization, interval, midpoint


class ResolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = calculate()

    def test_literal_memory_spectra_and_blank_weights(self):
        for q in (F(2, 3), F(-1, 3), F(-1, 6)):
            matrix = np.array([[0., float(q)], [float(q), 1.]])
            energies, vectors = eigh(matrix)
            exact = polarization(q)
            for i, row in enumerate(exact):
                self.assertAlmostEqual(energies[i], midpoint(row['shift']), delta=1e-14)
                self.assertAlmostEqual(abs(vectors[0, i])**2, midpoint(row['weight']), delta=1e-14)
            self.assertAlmostEqual(np.dot(abs(vectors[0])**2, energies), 0, delta=1e-14)
            self.assertAlmostEqual(np.dot(abs(vectors[0])**2, energies**2), float(q*q), delta=1e-14)

    def test_complete_native_tensor_spectra_and_schur_response(self):
        channels = self.receipt['mechanisms']['neutral_polarization_memory']['channels']
        for family, row in channels.items():
            records = row['complete_records']
            lows = [r for r in records if r['branch'] == 'lower']
            bare = np.array([midpoint(interval(r['original_bare_source_excitation_interval'])) for r in lows])
            q = float(F(row['total_native_target_charge']))
            memory = np.array([[0., q], [q, 1.]])
            full = np.kron(np.diag(bare), np.eye(2))+np.kron(np.eye(len(bare)), memory)
            energies = eigh(full, eigvals_only=True)
            expected = sorted(midpoint(interval(r['source_or_relational_excitation_interval'])) for r in records)
            np.testing.assert_allclose(energies, expected, atol=1e-12, rtol=0)
            z = .3+1.1j
            direct = np.linalg.inv(z*np.eye(len(full))-full)[::2, ::2]
            schur = np.diag(1/(z-bare-q*q/(z-bare-1)))
            np.testing.assert_allclose(direct, schur, atol=1e-12, rtol=0)

    def test_all_outcomes_and_charge_sensitive_second_moment(self):
        channels = self.receipt['mechanisms']['neutral_polarization_memory']['channels']
        for family, row in channels.items():
            self.assertEqual(len(row['complete_records']), 2*row['complete_native_energy_count'])
            self.assertEqual(row['failure_probability'], '1/2')
            weights = [interval(r['original_probability_interval']) for r in row['complete_records']]
            self.assertLessEqual(sum(w[0] for w in weights), F(1, 2))
            self.assertGreaterEqual(sum(w[1] for w in weights), F(1, 2))
            self.assertTrue(all(interval(r['source_or_relational_excitation_interval'])[0] > 0
                                for r in row['complete_records'] if not r['zero_weight']))
            detached = row['detached_coupling_zero_records']
            self.assertTrue(all(r['zero_weight'] for r in detached if r['branch'] == 'upper'))
        self.assertEqual(F(channels['uc']['complete_second_energy_moment_change'])
                         -F(channels['dc']['complete_second_energy_moment_change']), F(1, 6))

    def test_ordering_and_no_fitted_common_scale(self):
        mechanisms = self.receipt['mechanisms']
        self.assertGreater(mechanisms['positive_source_current_square']['one_scale_test']['ratio_numeric'], 1)
        for name in ('compensated_reference_attraction', 'neutral_polarization_memory'):
            self.assertLess(mechanisms[name]['one_scale_test']['ratio_numeric'], 1)
        self.assertTrue(all(not m['one_scale_test']['one_scale_matches_frozen_two_targets'] for m in mechanisms.values()))
        self.assertFalse(self.receipt['scope']['physical_quark_rest_masses_derived'])

    def test_complete_spatial_symmetry_and_unchanged_contact(self):
        symmetry = self.receipt['native_equal_coupling_domain']['microscopic_full_spatial_symmetry']
        self.assertEqual(symmetry['full_compact_rows_checked'], 686)
        self.assertTrue(symmetry['all_native_local_interactions_spin_rotation_invariant'])
        self.assertFalse(symmetry['unrestricted_spatial_global_ground_B_zero_proved'])
        controls = self.receipt['exact_operator_checks']
        self.assertEqual(controls['all_anomalous_contacts_unchanged'], 900)
        self.assertFalse(controls['all_scalar_poles_unchanged_claimed'])


if __name__ == '__main__':
    unittest.main()
