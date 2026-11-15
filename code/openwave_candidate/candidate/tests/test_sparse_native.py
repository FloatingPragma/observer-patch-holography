import importlib.util
from pathlib import Path
import unittest
from fractions import Fraction
import numpy as np

from oph_openwave import sparse_native as m

class SparseWindowTests(unittest.TestCase):
    def triangle(self,radius=2):return m.SparseCoarseWindow(m.Construction(3,((0,1),(1,2),(2,0)),.23),6,radius)

    def test_rejections_restore_trace_and_are_not_native_outcomes(self):
        window=self.triangle();size=len(window.states);complete=np.zeros((size,size));dropped=np.zeros_like(complete)
        for key,src,dst in window.maps:
            k=np.zeros((size,size));k[dst,src]=1/np.sqrt(4*len(window.construction.seams))
            complete+=k.T@k
            if key[-1]!='boundary_rejected':dropped+=k.T@k
        np.testing.assert_allclose(complete,np.eye(size),atol=1e-14)
        self.assertGreater(np.linalg.norm(complete-dropped),.1)
        result,rho=window.interval(.37,.29,True)
        self.assertGreater(result['boundary_rejection_probability'],0)
        self.assertTrue(all(r['native_record'] is None and r['applied_current']==0 for r in result['ledger'] if r['status']=='boundary_rejected'))
        self.assertAlmostEqual(np.trace(rho).real,1)
        self.assertGreater(np.linalg.eigvalsh(rho).min(),-1e-12)

    def test_local_null_outcome_preserves_coherence(self):
        window=self.triangle(4)
        key,src,dst=next(item for item in window.maps if item[0][1]=='C' and item[0][3]==0 and item[0][-1]=='native' and len(item[1])>=2)
        source=np.zeros(len(window.states),complex);source[src[:2]]=1/np.sqrt(2)
        branch=np.zeros_like(source);branch[dst]=source[src]
        rho=np.outer(branch,branch.conj())
        self.assertGreater(abs(rho[dst[0],dst[1]]),.4)
        self.assertNotEqual(rho[dst[0],dst[1]],np.diag(np.diag(rho))[dst[0],dst[1]])

    def test_exact_full_comparison_and_actual_ring(self):
        result=m.run(include_ledger=True)
        for row in result['triangle_full_comparison']:
            self.assertGreater(row['actual_trace_distance'],0)
            self.assertLess(row['actual_trace_distance'],row['bound']['trace_distance_upper_bound'])
        ring=result['actual_ring']
        self.assertEqual(ring['geometry'],'flat_ring(1,1)')
        self.assertLess(ring['states'],100)
        self.assertLess(ring['regulator_trace_distance_bound'],.006)
        self.assertAlmostEqual(sum(ring['complete_window_readout']['position_probabilities']),1)
        self.assertGreater(ring['omitted_coherent_derivative_norm'],0)

    def test_rational_numerical_certificate_preserves_norm_defect(self):
        window=self.triangle(4)
        certificate=window.certified_interval(.1,.29)
        real,imag,tail,order=m._rational_taylor_action(window.h,.23,2,.1)
        norm=sum(a*a+b*b for a,b in zip(real,imag))
        self.assertNotEqual(norm,Fraction(1))
        self.assertGreaterEqual(certificate['polynomial_norm_defect_abs_upper'],abs(float(norm-1)))
        self.assertLess(certificate['discarded_state_norm_tail_upper'],1e-14)
        self.assertGreaterEqual(certificate['numerical_effect_error_upper'],certificate['discarded_state_norm_tail_upper'])
        self.assertEqual(len(certificate['positions']),6)
        self.assertTrue(all(row['full_modular_probability_interval'][0] <= row['polynomial_probability'] <= row['full_modular_probability_interval'][1] for row in certificate['positions']))
        # Independent floating propagation is a numerical check, not the proof.
        result,_=window.interval(.1,.29)
        differences=[abs(a-b['polynomial_probability']) for a,b in zip(result['position_probabilities'],certificate['positions'])]
        self.assertLess(max(differences),certificate['numerical_effect_error_upper']+1e-15)

    def test_exact_outward_conversion_and_subnormalized_effect_control(self):
        for q in (Fraction(1,3),Fraction(2,7),Fraction(1,10**100)):
            self.assertLessEqual(Fraction.from_float(m.outward(q,False)),q)
            self.assertGreaterEqual(Fraction.from_float(m.outward(q,True)),q)
        # A scalar polynomial amplitude 1+eps has probability error2eps+eps².
        # Half the trace norm alone is insufficient for this trace-changing A.
        epsilon=Fraction(1,100)
        probability_error=(1+epsilon)**2-1
        half_trace=epsilon*(2+epsilon)/2
        self.assertGreater(probability_error,half_trace)
        self.assertEqual(half_trace+probability_error/2,probability_error)

    def test_bound_inputs_and_zero_interval(self):
        c=m.Construction(3,((0,1),(1,2),(2,0)),.23)
        self.assertEqual(m.trace_bound(c,4,0.,0.)['trace_distance_upper_bound'],0)
        for time,probability in ((float('nan'),.2),(.1,float('nan')),(.1,-.1)):
            with self.assertRaises(ValueError):m.trace_bound(c,4,time,probability)
        with self.assertRaises(ValueError):m.SparseCoarseWindow(c,radius=True)

if __name__=='__main__':unittest.main()


def test_strict_runner_schema_and_compact_default():
    from oph_openwave.runner import validate_result
    result=m.run()
    validate_result('sparse_native',result)
    assert 'complete_window_readout' not in result['actual_ring']
    assert 'complete_certificate' not in result['actual_ring']
    assert set(result['source_sha256'])=={'sparse_native.py','modular_native.py','construction.py','geometry.py'}
    complete=m.run(include_ledger=True)
    for row in complete['actual_ring']['complete_certificate']['ledger']:
        if row['status']=='boundary_rejected':
            assert row['native_record'] is None and row['applied_current']==0
            assert 'proposed_current' in row
