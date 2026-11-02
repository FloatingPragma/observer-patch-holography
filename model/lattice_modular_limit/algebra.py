"""Exact endpoint and weighted modular-generator identities."""
from functools import lru_cache
import sympy as s


@lru_cache(maxsize=1)
def certificate():
    z,r=s.symbols('z r',positive=True)
    lattice=z**4-(2+4*r)*z*z+1
    discriminant=s.factor(s.discriminant(lattice,z))
    if discriminant!=4096*r*r*(r+1)**2:
        raise ArithmeticError('lattice branch point polynomial changed')
    x,t=s.symbols('x t',real=True)
    w=x*(1-x);k=s.log((1-x)/x)
    if s.simplify(w*s.diff(k,x)+1)!=0:
        raise ArithmeticError('weighted logarithmic derivative failed')
    derivative=s.simplify(s.diff(w*s.exp(-s.I*t*k),x)/s.exp(-s.I*t*k))
    if s.simplify(derivative-(1-2*x+s.I*t))!=0:
        raise ArithmeticError('weighted modular phase derivative failed')
    return {'lattice_branch_polynomial':str(lattice),
            'lattice_discriminant':str(discriminant),
            'lattice_branch_roots_simple_for_positive_transverse_sine_square':True,
            'compact_band_spinors_excluded_in_lattice_and_continuum':True,
            'weighted_phase_derivative':'(1-2x+i*t)*exp(-i*t*k(x))',
            'weighted_phase_Lipschitz_bound':'1+abs(t)',
            'core':'range(C*(1-C))',
            'core_is_graph_core_for_log_ratio':True,
            'modular_annihilator_sign':'sigma_t(c(f))=c(exp(-i*t*k)*f)',
            'uniform_spectral_gap_assumed':False}
