"""Independent floating diagnostics; only verified facts enter the receipt."""
import math
import numpy as np
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss


def finite(value, label):
    if not np.all(np.isfinite(value)):
        raise ArithmeticError('nonfinite '+label)
    return value


def diagnostics():
    pi=math.pi
    def density(k):return (k*k+pi*pi)/(4*pi*math.cosh(k/2)**2)
    norm,error=quad(density,-40,40,epsabs=1e-11)
    finite([norm,error],'density integral')
    if abs(norm-4*pi/3)>1e-8:raise ArithmeticError('spectral density normalization failed')
    orbit_errors=[]
    for t in (0.125,0.25,0.5):
        val,err=quad(lambda k:math.cos(t*k)*density(k),-40,40,epsabs=1e-11)
        finite([val,err],'spectral orbit integral')
        a=pi*t;expected=2*pi*(math.cosh(a)/math.sinh(a)**2-a/math.sinh(a)**3)
        orbit_errors.append(abs(val-expected))
    if max(orbit_errors)>1e-8:raise ArithmeticError('spectral orbit identity failed')
    stieltjes_errors=[]
    for t in (0.125,1.0,8.0):
        val,err=quad(lambda k:density(k)/(1/(1+math.exp(k))+t),-40,40,epsabs=1e-10)
        finite([val,err],'Stieltjes integral')
        L=math.log1p(1/t);expected=(L**3+4*pi*pi*L)/(3*pi)
        stieltjes_errors.append(abs(val-expected))
    if max(stieltjes_errors)>1e-8:raise ArithmeticError('Stieltjes identity failed')
    x,w=leggauss(240);x=x*30;w=w*30
    g=(x*x+pi*pi)/np.cosh(x/2);difference=x[:,None]-x[None,:]
    h=np.divide(difference,np.sinh(difference/2),out=np.full_like(difference,2),where=difference!=0)
    finite([*x,*w,*g],'spectral quadrature nodes')
    finite(h,'divided difference kernel')
    chi=float(np.sum(w[:,None]*w[None,:]*h*g[:,None]*g[None,:])/(8*pi*pi))
    finite(chi,'independent spectral Hessian')
    expected=64*pi**4/35
    if abs(chi/expected-1)>1e-6:raise ArithmeticError('independent spectral Hessian failed')
    return {'spectral_mass_error':abs(norm-4*pi/3),'orbit_max_error':max(orbit_errors),
            'Stieltjes_max_error':max(stieltjes_errors),'hessian_relative_error':abs(chi/expected-1)}


def certificate():
    diagnostics()
    return {'spectral_density_and_real_orbit_checked':True,
            'independent_Stieltjes_integral_checked':True,
            'independent_two_variable_Hessian_checked':True,
            'nonfinite_intermediates_rejected':True,
            'mass_orbit_Stieltjes_absolute_tolerance':'1e-8',
            'Hessian_relative_tolerance':'1e-6',
            'raw_floating_diagnostics_serialized':False,
            'diagnostics_replace_analytic_proof':False}
