import importlib.util
from pathlib import Path
import unittest
import sympy as sp

from oph_openwave import general_cut as m

class GeneralCutTests(unittest.TestCase):
    def test_symbolic_identity_and_winding_alias(self):
        r=m.run()
        self.assertTrue(r['symbolic_region_identity'])
        self.assertTrue(r['cycle_alias']['same_integer_gauss_sector'])
        self.assertTrue(r['cycle_alias']['same_modular_record'])
        self.assertTrue(any(r['cycle_alias']['integer_half_record_difference']))
        self.assertTrue(all(value==0 for value in r['cycle_alias']['closed_boundary_current_differences']))
        self.assertTrue(r['sector_alias']['different_integer_gauss'])

    def test_wrong_orientation_rejected_by_internal_cut(self):
        _,columns,a=m.objects();source,target=columns[0]
        region={source,target}
        self.assertEqual(m.cut_coefficients(a,region)[0],0)
        broken=a.copy();broken[target,0]=1
        self.assertNotEqual(m.cut_coefficients(broken,region)[0],0)

    def test_arbitrary_record_cut_identity_not_bridge_specific(self):
        _,_,a=m.objects();region=set(range(0,a.rows,3))
        records=[(i%7)-3 for i in range(a.cols)]
        coeff=m.cut_coefficients(a,region)
        for position in (0,7,a.rows-1):
            g=m.raw_gauss(a,position,records)
            charge=int(position in region)
            current=sum(x*y for x,y in zip(coeff,records))
            self.assertEqual(sum(g[i] for i in region)-current,charge)
            self.assertNotEqual(sum(g[i] for i in region)-(current+1),charge)

if __name__=='__main__':unittest.main()
