"""Complete complex eigensystem certificates with exact integer residuals."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from .operators import binomial_half, SPEED_SQUARED
from model.openwave_transverse_completion.source import root_upper, exact_phase


def polynomial_rounding_bound(mass, degree):
    """IEEE binary64 bound at theta0; Pauli-square powers are exact dyadics."""
    m = F(mass)
    u = F(1,2**53)
    rounded = total = F(0)
    for k in range(1,degree+1):
        c = binomial_half(k)*m**(1-2*k)
        total += abs(c)*SPEED_SQUARED**k
        rounded += abs(c-F(float(c)))*SPEED_SQUARED**k
    # Every power has Gaussian-integer entries /512**k with integer magnitude
    # <=4**k. For k<=6 these powers and their sparse additions are exact.
    assert degree <= 6
    gamma = 4*(degree+1)*u/(1-4*(degree+1)*u)
    return rounded+gamma*(total+rounded)


def complete_eigensystem(H, mass, degree):
    lam,Q = eigh(H,driver='evr')
    n = len(H); dv=1<<26; de=dh=1<<24
    X=np.rint(Q.real*dv).astype(np.int64)
    Y=np.rint(Q.imag*dv).astype(np.int64)
    A=np.rint(H.real*dh).astype(np.int64)
    B=np.rint(H.imag*dh).astype(np.int64)
    E=np.rint(lam*de).astype(np.int64)
    assert np.array_equal(A,A.T) and np.array_equal(B,-B.T)
    qmax=max(int(abs(X).max()),int(abs(Y).max()))
    assert 2*n*qmax*qmax+dv*dv < 2**63
    rowsum=max(sum(abs(int(a))+abs(int(b)) for a,b in zip(ar,br,strict=True)) for ar,br in zip(A,B,strict=True))
    assert rowsum*qmax+qmax*int(abs(E).max()) < 2**63
    gr=X.T@X+Y.T@Y-np.eye(n,dtype=np.int64)*dv*dv
    gi=X.T@Y-Y.T@X
    rr=A@X-B@Y-X*E
    ri=A@Y+B@X-Y*E
    g=root_upper(F(sum(int(x)**2 for x in gr.flat)+sum(int(x)**2 for x in gi.flat),dv**4))
    r=root_upper(F(sum(int(x)**2 for x in rr.flat)+sum(int(x)**2 for x in ri.flat),(dh*dv)**2))
    difference=F(0)
    for ar,br,hr in zip(A,B,H,strict=True):
        for a,b,h in zip(ar,br,hr,strict=True):
            difference += (F(int(a),dh)-F(float(h.real)))**2+(F(int(b),dh)-F(float(h.imag)))**2
    source_error=root_upper(difference)+polynomial_rounding_bound(mass,degree)
    # Fixed outward budgets are checked against the complete live residuals.
    # The result does not depend on a LAPACK basis inside a degenerate space.
    assert g<=F(n,10**8) and r<=F(n,10**7) and source_error<=F(n,10**8)
    g,r,source_error=F(n,10**8),F(n,10**7),F(n,10**8)
    assert g < F(1,100)
    K=F(max(abs(int(x)) for x in E),de)
    spectral_error=(2*K*g+(1+g)*r)/(1-g)+source_error
    return {'Q':Q,'eigenvalues':lam,'X':X,'Y':Y,'E':E,
            'vector_denominator':dv,'energy_denominator':de,
            'Gram_error_upper_exact':str(g),'eigensystem_residual_upper_exact':str(r),
            'exact_polynomial_to_dyadic_matrix_norm_upper_exact':str(source_error),
            'uniform_eigenvalue_error_upper_exact':str(spectral_error),
            'complete_dimension':n}


def actual_vector_certificate(cert,psi,initial_index,time):
    X,Y,E=cert['X'],cert['Y'],cert['E'];dv=cert['vector_denominator'];de=cert['energy_denominator']
    dp=1<<80
    phases=[exact_phase(F(int(e)*time,de)) for e in E]
    pr=np.array([int(p[0]*dp) for p in phases],dtype=object)
    pi=np.array([int(p[1]*dp) for p in phases],dtype=object)
    # Q*|initial> has numerator X_initial-iY_initial.
    cr=pr*X[initial_index]+pi*Y[initial_index]
    ci=pi*X[initial_index]-pr*Y[initial_index]
    re=X.astype(object)@cr-Y.astype(object)@ci
    im=X.astype(object)@ci+Y.astype(object)@cr
    den=dv*dv*dp
    difference=sum(((F(int(a),den)-F(float(z.real)))**2+(F(int(b),den)-F(float(z.imag)))**2
                    for a,b,z in zip(re,im,psi,strict=True)),F(0))
    g=F(cert['Gram_error_upper_exact']);r=F(cert['eigensystem_residual_upper_exact'])
    source=F(cert['exact_polynomial_to_dyadic_matrix_norm_upper_exact'])
    error=g+abs(time)*r*root_upper(1+g)+root_upper(difference)+(1+g)*max(p[2] for p in phases)+abs(time)*source
    budget=F(len(psi),1000000)
    assert error <= budget, ('complete numerical budget exceeded',float(error),float(budget))
    return {'conditional_vector_error_upper_exact':str(budget),
            'budget_checked_against_actual_exact_calculation':True,
            'method':'complete complex dyadic eigensystem, exact integer Gram/residual, rational phases and actual amplitude comparison'}
