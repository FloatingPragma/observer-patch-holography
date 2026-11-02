"""Exact radial Sobolev bounds and the actual symmetric detector form."""
from functools import lru_cache
from fractions import Fraction as F
from math import isqrt
import sympy as s
from model.weyl_gravity_response.algebra import response_matrix
from model.weyl_geometric_modular.algebra import pauli, source_frames


def fraction(x):
    return F(int(s.numer(x)), int(s.denom(x)))


def laurent_integral(expression, r):
    """Exact rational integral, rejecting an uncancelled logarithm."""
    total = s.Rational(0)
    for term in s.Add.make_args(s.expand(expression)):
        power = term.as_powers_dict().get(r, s.Integer(0))
        if power.is_Integer is not True:
            raise ValueError('integer Laurent powers required')
        coefficient = s.cancel(term/r**power)
        if coefficient.is_Rational is not True:
            raise ValueError('rational Laurent coefficients required')
        if power == -1:
            raise ArithmeticError('uncancelled logarithmic moment')
        integral = (s.Integer(2)**(power+1)-1)/(power+1)
        total += coefficient*integral
    return fraction(total)


def ceiling_sqrt(value):
    value = F(value)
    if value < 0:
        raise ValueError('nonnegative rational required')
    n = isqrt(value.numerator//value.denominator)
    return n if n*n == value else n+1


@lru_cache(maxsize=1)
def sobolev():
    r = s.symbols('r', positive=True)
    phi = r*r*(r-1)**4*(2-r)**4
    if any(s.diff(phi,r,m).subs(r,x) != 0 for m in range(4) for x in (1,2)):
        raise ArithmeticError('shell is not weakly H4 across its boundaries')
    # Angular L2 weights after removing their common pi^3.
    weights = {0:s.Integer(396900), 1:s.Integer(176400), 2:s.Rational(99225,2)}
    norm = s.Rational(150088295,4249388)
    if sum(weights.values())*s.integrate(r*r*phi**2,(r,1,2)) != norm:
        raise ArithmeticError('spherical harmonic weights do not reproduce shell norm')
    moments=[]
    for m in range(1,5):
        integrand=s.Integer(0)
        for ell,weight in weights.items():
            value=phi
            for _ in range(m//2):
                value=s.expand(s.diff(value,r,2)+2*s.diff(value,r)/r-ell*(ell+1)*value/r**2)
            density=(r*r*s.diff(value,r)**2+ell*(ell+1)*value**2) if m%2 else r*r*value**2
            integrand += (weight/norm)*density
        moment=laurent_integral(integrand,r)
        if moment<=0:
            raise ArithmeticError('nonpositive Sobolev moment')
        moments.append(moment)
    return {'angular_weights_without_pi_cubed':{str(k):str(v) for k,v in weights.items()},
            'exact_normalized_momentum_moments':[str(v) for v in moments],
            'H1_H2_H3_H4_integer_majorants':[ceiling_sqrt(v) for v in moments],
            'weak_H4_boundary_test':True,
            'H5_membership_claimed':False}


@lru_cache(maxsize=1)
def certificate():
    sigma=pauli();p=s.symbols('p1:4',real=True);u=s.symbols('u0:4',real=True)
    h=sum((p[j]*sigma[j] for j in range(3)),s.zeros(2))
    A=u[0]*s.eye(2)-sum((u[j+1]*sigma[j] for j in range(3)),s.zeros(2))
    D=u[0]*h-sum(u[j+1]*p[j] for j in range(3))*s.eye(2)
    if s.simplify((A*D+D*A)/2-(A*D+D*A).adjoint()/2) != s.zeros(2):
        raise ArithmeticError('detector principal form is not Hermitian')
    # Off-diagonal momenta test the weighted derivative ordering too.
    q=s.symbols('q1:4',real=True)
    hq=sum((q[j]*sigma[j] for j in range(3)),s.zeros(2))
    avec=(s.eye(2),*(-x for x in sigma))
    dp=(h,*(-x*s.eye(2) for x in p));dq=(hq,*(-x*s.eye(2) for x in q))
    contracted=s.zeros(2)
    for mu in range(4):
        for nu in range(4):
            stress=(avec[mu]*dq[nu]+dp[nu]*avec[mu]
                    +avec[nu]*dq[mu]+dp[mu]*avec[nu])/4
            contracted += u[mu]*u[nu]*stress
    Dq=u[0]*hq-sum(u[j+1]*q[j] for j in range(3))*s.eye(2)
    if any(s.expand(x)!=0 for x in contracted-(A*Dq+D*A)/2):
        raise ArithmeticError('weighted detector differs from symmetric covariant stress')
    matrix=response_matrix()
    inverse=max(sum(abs(v) for v in row) for row in matrix.inv().tolist())
    if inverse != s.Rational(68,5):
        raise ArithmeticError('nine-direction response inverse changed')
    return {'sobolev':sobolev(),'frames':source_frames(),
            'detector':'Q_u(w)=1/2 {w A_u,D_u}',
            'A_u':'u0 I - sum u_j sigma_j','D_u':'u0 h - sum u_j p_j',
            'rest_detector':'Q_(1,0)(w)=1/2 {w,h}',
            'response_determinant':str(matrix.det()),'response_inverse_norm':str(inverse),
            'selected_valley':[1,1,1], 'all_eight_valleys_retained':True,
            'off_diagonal_momentum_symmetric_stress_contraction_verified':True}
