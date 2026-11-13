"""Constraint, exact-flow and independent full-equation controls."""
import unittest
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from model.openwave_capacity_transfer.source import background,rhs,components
from .source import I,ZERO,ONE,mid,start_data,transport,reduction,step_state,exp_point

class TransferTests(unittest.TestCase):
    def test_exact_comoving_reduction(self):
        self.assertEqual(reduction()['symbolic_residuals'],['0','0'])
    def test_interval_step_encloses_known_stable_affine_flow(self):
        c,r,_,_,_=step_state([F(0)],F(0),[[I.point(-2)]],[[I.point(-2)]],[[ZERO]],[ONE],F(1,10))
        exact=(ONE-exp_point(F(-1,5)))/2
        self.assertLessEqual(c[0]-r,exact.lo);self.assertGreaterEqual(c[0]+r,exact.hi)
    def test_contracting_cell_retains_start_radius(self):
        _,_,_,_,tube=step_state([F(0)],F(1),[[I.point(-2)]],[[I.point(-2)]],[[ZERO]],[ZERO],F(1,10))
        self.assertLessEqual(tube[0].lo,-1);self.assertGreaterEqual(tube[0].hi,1)
    def test_original_full_seven_equations_inside_certificate(self):
        cert=transport();T=float(mid(start_data()[4]));N=-1.
        m,r,s,b,A,M,h,H2=background()(N);phi=1.5-T
        for row in cert['finite_wavelength_records']:
            k=float(F(row['k_over_H0_exact']));Q=k*k/(np.exp(2*N)*H2)
            state=np.array([-2*Q*phi/(3*m),T,0.,0.,s*T,b*T+s*(phi+h*T),phi])
            self.assertLess(abs(components(N,state,k,'capacity')[10]),2e-14)
            sol=solve_ivp(lambda t,v:rhs(t,v,k,'capacity'),(N,0),state,method='DOP853',rtol=2e-12,atol=1e-14,max_step=.01)
            self.assertTrue(sol.success)
            for value,key in [(sol.y[6,-1],'metric_transfer_interval'),(sol.y[0,-1],'comoving_dust_interval')]:
                lo,hi=map(lambda v:float(F(v)),row[key]);self.assertLess(lo,value);self.assertLess(value,hi)
            self.assertLess(abs(components(0,sol.y[:,-1],k,'capacity')[10]),1e-10)
            self.assertLess(sol.y[0,-1],0)
    def test_invalid_interval_resolution_rejected(self):
        for n in (0,64,1.5,4097):
            with self.assertRaises(ValueError):transport(n)
if __name__=='__main__':unittest.main()
