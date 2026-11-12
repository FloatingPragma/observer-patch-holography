from fractions import Fraction as F
import unittest

from .compression import native_spans,exact_leakage,complete_dilation_control,resource_bound


class CompressionTests(unittest.TestCase):
    def test_actual_native_projector_on_every_column(self):
        rows=native_spans()
        self.assertEqual(sum(r['all_projector_columns_checked'] for r in rows),10800)
        self.assertTrue(rows[1]['effect']<rows[2]['effect']<rows[0]['effect'])
        self.assertTrue(all(max(r['errors'].values())<3e-12 for r in rows))

    def test_exact_leakage_retains_occupied_sector_dependence(self):
        value=exact_leakage()
        self.assertEqual(value['all_physical_occupations'],4)
        self.assertEqual(value['leakage_eigenvalues'],['1-e^2','1-e','1-e','(1-e)^2'])
        # A full physical pair still leaks to an ancillary pair: an occupancy
        # cutoff or a Pauli-blocking-only calculation would incorrectly give0.
        for e in (F(1,5),F(1,2),F(9,10)):
            self.assertGreater((1-e)**2,0)
            self.assertGreater(1-e*e,1-e)

    def test_full_quantum_archive_isometry_and_wrong_frame_control(self):
        for row in native_spans():
            c=complete_dilation_control(row['effect'])
            self.assertEqual(c['raw_quantum_archive_coordinates'],4096)
            self.assertEqual(len(c['complete_fine_diagonal_records']),64)
            self.assertEqual(len(c['complete_coarse_history_probabilities']),8)
            self.assertAlmostEqual(sum(c['complete_coarse_history_probabilities'].values()),1,places=12)
            self.assertGreater(c['maximum_retained_offdiagonal_archive_norm'],1e-8)
            self.assertGreater(c['wrong_equal_leg_frame_Hamiltonian_error'],.1)
            self.assertLessEqual(c['complete_output_trace_error'],c['sharper_unhalved_bound']+2e-12)
            self.assertLess(max(c['errors'].values()),3e-12)

    def test_volume_bound_and_all_archive_resources(self):
        a=resource_bound(10**15,sites=1);b=resource_bound(10**15,sites=64)
        self.assertEqual(F(b['complete_unhalved_trace_error_upper_exact']),64*F(a['complete_unhalved_trace_error_upper_exact']))
        self.assertEqual(b['total_CAR_modes'],6528)
        self.assertEqual(b['quantum_archive_modes'],12*64*10**15)
        self.assertLess(F(a['complete_unhalved_trace_error_upper_exact']),F(12,10**15))
        for args in ((0,F(1),1),(1,-F(1),1),(1,F(1),0)):
            with self.assertRaises(ValueError):resource_bound(*args)


if __name__=='__main__':unittest.main()
