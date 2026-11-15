import unittest
from fractions import Fraction as F
from . import source as s

class ArchitectureTests(unittest.TestCase):
    def test_complete_operator_instrument_and_padding(self):
        p=s.event_permutation();cols=s.coherent_columns()
        self.assertEqual(len(set(p)),256)
        self.assertEqual(len(set(sum((list(x) for x in cols),[]))),32)
        for x in range(8):
            self.assertEqual(sum(w for _,w in s.committed_column(x,x).values()),1)
        for x in (6,7):self.assertTrue(all((z&7)==x for z in cols[x]))
    def test_raw_commit_differs_from_recoded_coin_trace(self):
        self.assertEqual(s.committed_column(2,1),{})
        # Re-encoded swap coins retain |10><01| and its swapped image.
        r=s.record_recoder();cols=s.coherent_columns();out={}
        for a,b in zip(cols[2],cols[1],strict=True):
            a,b=r[a],r[b]
            if a>>3==b>>3:out[a&7,b&7]=out.get((a&7,b&7),F())+F(1,4)
        # The basis-dependent coin relabelling changes which columns share
        # an environment label; sum all coherent terms within each law.
        out={}
        for a in map(lambda z:r[z],cols[2]):
            for b in map(lambda z:r[z],cols[1]):
                if a>>3==b>>3:out[a&7,b&7]=out.get((a&7,b&7),F())+F(1,4)
        self.assertEqual(out,{(2,1):F(1,2),(1,2):F(1,2)})
    def test_generic_reversible_compiler_non_involution(self):
        p=list(range(16));p[1],p[5],p[14]=5,14,1
        program=s.compile_permutation(p);full=s.execute_basis(program)
        self.assertEqual(list(full[:16]),p)
        self.assertEqual(len(set(full)),len(full))
        with self.assertRaises(ValueError):s.compile_permutation([0]*8)
    def test_inverse_recoding_restores_every_native_branch(self):
        r=s.record_recoder();inverse={y:x for x,y in enumerate(r)}
        for col in s.coherent_columns():
            for x in col:self.assertEqual(inverse[r[x]],x)
        self.assertNotEqual(r,tuple(range(256)))
    def test_actual_source_primitive_errors(self):
        _,errors,_=s.execute_primitives()
        self.assertLess(max(errors.values()),1e-10)
