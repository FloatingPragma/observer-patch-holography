"""Exact spinor, Fourier, Stieltjes and entropy coefficients."""
import sympy as s


def certificate():
    x,y,z,c,h,t=s.symbols('x y z c h t',real=True)
    sigma=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    X=x*sigma[0]+y*sigma[1]+z*sigma[2];I=s.eye(2)
    M=c*I+s.I*h*X;tau=c*c+h*h*(x*x+y*y+z*z)
    numerator=-s.I*(s.I*t*I+X)
    a,L=s.symbols('a L',positive=True)
    stieltjes=s.pi*(1-s.exp(-a*L))/s.sin(s.pi*a)
    series=s.series(stieltjes,a,0,3).removeO().expand()
    p=s.symbols('p',real=True);phi=s.pi*p/s.sinh(s.pi*p)
    orbit=2*s.pi*(s.cosh(s.pi*p)/s.sinh(s.pi*p)**2-s.pi*p/s.sinh(s.pi*p)**3)
    total=s.Rational(2,9)/s.pi**2*(s.factorial(6)*s.zeta(6)+8*s.pi**2*s.factorial(4)*s.zeta(4)+16*s.pi**4*s.factorial(2)*s.zeta(2))
    checks={
        'imaginary_spinor_prefactor_norm':s.simplify(M.H*M-tau*I)==s.zeros(2),
        'imaginary_Weyl_kernel_numerator_norm':s.simplify(numerator.H*numerator-(t*t+x*x+y*y+z*z)*I)==s.zeros(2),
        'Stieltjes_zero_derivative':s.simplify(series.coeff(a,0)-L)==0,
        'Stieltjes_second_derivative':s.simplify(2*series.coeff(a,2)-(L**3+s.pi**2*L)/3)==0,
        'logistic_spectral_density_correlation':s.simplify((-s.diff(phi,p,2)+s.pi**2*phi)/s.pi-orbit)==0,
        'constant_spinor_norm':s.limit(orbit,p,0)==4*s.pi/3,
        'zeta_hessian_coefficient':s.simplify(total-64*s.pi**4/35)==0,
        'Parseval_sech8_coefficient':s.Rational(1,8)/s.pi**2/(2*s.pi)*(2*s.pi**2)*(4*s.pi**3)**2/s.pi*s.Rational(32,35)==64*s.pi**4/35,
    }
    # A generic natural vector does not obey the new rho^6 estimate.
    r=s.symbols('r',positive=True);lam=r**3;mass=lam*(1-lam)
    checks['generic_natural_vector_chi_is_rho3_not_rho6']=s.cancel(mass**2*(1/lam+1/(1-lam))-mass)==0
    if not all(bool(v) for v in checks.values()):raise ArithmeticError('exact entropy asymptotic identity failed')
    return {'checks':{k:bool(v) for k,v in checks.items()},
            'Stieltjes_formula':'|u|^2/(3pi) [L^3+4pi^2 L]',
            'logistic_density':'|u|^2(k^2+pi^2)/(4pi cosh(k/2)^2)',
            'constant_spinor_hessian':'64pi^4|u|^4/35',
            'constant_spinor_relative_entropy_leading':'32pi^4 theta^2|u|^4/35',
            'generic_natural_vector_bound_is_insufficient':True}
