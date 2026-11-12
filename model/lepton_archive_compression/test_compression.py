import unittest
import numpy as np
from . import compression as c


class ArchiveTests(unittest.TestCase):
    def test_exact_minors_exclude_zero(self):
        rows=c.exact_rank_certificate()['all_six_nonzero_pair_minors']
        self.assertEqual(len(rows),6)
        for row in rows:
            lo,hi=map(c.F,row['rational_minor_interval'])
            self.assertTrue(hi<0 or lo>0)

    def test_all_histories_keep_rare_quantum_directions(self):
        expected=({1:2,2:2,3:0},{1:4,2:8,3:4},{1:16,2:64,3:176})
        for n,counts in zip((1,2,4),expected):
            run=c.history_run(n)
            self.assertEqual(run['rank_counts'],{str(k):v for k,v in counts.items()})
            self.assertEqual(len(run['all_original_histories']),4**n)
            self.assertAlmostEqual(run['complete_probability_sum'],1,places=12)
            self.assertGreater(run['wrong_single_environment_vector_unhalved_reference_output_distance'],1e-7)
        self.assertGreater(run['relative_eigenvalue_threshold_misclassified_histories'],100)

    def test_streaming_inverse_preserves_entangled_family(self):
        data=c.streaming_sequence(('x_no','both_no','x_no','both_no'))
        compressed=data['columns'][:3]/np.sqrt(3)
        expected=data['original']/np.sqrt(3)
        self.assertLess(np.linalg.norm(data['embedding']@compressed-expected),1e-12)
        self.assertLess(np.linalg.norm(data['embedding'].T@data['embedding']-np.eye(3)),1e-12)
        with self.assertRaises(ArithmeticError):c.complete_frame(data['original'],2)

    def test_native_small_compressor_all_input_columns(self):
        data=c.streaming_sequence(('both_no','x_no'))
        tape=c.native.compile_real_orthogonal(data['unitary'])
        actual=c.native.execute(tape,2)
        self.assertLess(np.linalg.norm(actual-data['unitary']),1e-11)
        self.assertLess(np.linalg.norm(actual.conj().T@actual-np.eye(4)),1e-11)


if __name__=='__main__':unittest.main()
