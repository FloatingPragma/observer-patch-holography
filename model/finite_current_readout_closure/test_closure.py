import unittest
from fractions import Fraction as F
import numpy as np
from model.defect_integer_lift import model as lift
from .closure import axial, integer_lift


class CurrentClosureTests(unittest.TestCase):
    def test_complete_native_record_law(self):
        r=axial();counts=np.asarray(r['complete_two_attempt_word_counts'])
        self.assertEqual(counts.shape,(180,180))
        self.assertEqual(counts.sum(),66*60**2)
        self.assertGreater(r['retained_zero_word_count'],0)
        self.assertGreater(F(r['exact_four_step_uniform_minorization']),0)
        self.assertFalse(r['physical_position_or_lab_time_identified'])

    def test_integer_source_attachment_and_nonidentifiability(self):
        r=integer_lift()
        self.assertTrue(all(r['reference_checks'].values()))
        self.assertEqual(r['complete_load_coin_class_contexts'],2592)
        self.assertGreater(r['noop_contexts'],0)
        self.assertFalse(r['archive_census_reexecuted'])
        fed=lift.Federation(0)
        base=[0]*fed.seams;other=base.copy();other[0]=6
        a=lift.curl(fed,base);b=lift.curl(fed,other)
        self.assertNotEqual(a,b)
        self.assertEqual([v%6 for v in a],[v%6 for v in b])
        self.assertEqual(sum(b),0)


if __name__=='__main__':unittest.main()
