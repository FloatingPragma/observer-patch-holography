"""Full finite-CAR checks of modular sign and relative entropy."""
from functools import lru_cache
import sympy as s


@lru_cache(maxsize=1)
def certificate():
    annihilators=[]
    for mode in range(2):
        matrix=s.zeros(4)
        for n in range(4):
            if n&(1<<mode):matrix[n^(1<<mode),n]=(-1)**((n&((1<<mode)-1)).bit_count())
        annihilators.append(matrix)
    rho=s.diag(s.Rational(4,25),s.Rational(1,25),s.Rational(16,25),s.Rational(4,25))
    logarithm=s.diag(*(s.log(x) for x in rho.diagonal()))
    f=s.Matrix([1,1])/s.sqrt(2)
    a=sum((s.conjugate(f[j])*annihilators[j] for j in range(2)),s.zeros(4))
    n=a.adjoint()*a
    unitary=s.eye(4)+(s.I-1)*n
    changed=unitary*rho*unitary.adjoint()
    direct=s.simplify(s.trace(rho*logarithm)-s.trace(changed*logarithm))
    expected=s.Rational(3,5)*s.log(4)
    if s.simplify(direct-expected)!=0:
        raise ArithmeticError('finite-CAR relative entropy formula failed')
    C=s.diag(s.Rational(1,5),s.Rational(4,5));k=s.diag(s.log(4),-s.log(4))
    cov=(f.adjoint()*C*f)[0]*(f.adjoint()*k*f)[0]-(f.adjoint()*C*k*f)[0]
    if s.simplify(2*cov-expected)!=0:
        raise ArithmeticError('one-particle entropy expression failed')
    # Derivative of rho^(it) a rho^(-it), compared with c(exp(-it k)f).
    derivative=s.I*(logarithm*a-a*logarithm)
    predicted=sum((s.conjugate((-s.I*k*f)[j])*annihilators[j] for j in range(2)),s.zeros(4))
    if any(s.simplify(x)!=0 for x in derivative-predicted):
        raise ArithmeticError('modular annihilator sign failed')
    if all(s.simplify(x)==0 for x in derivative+predicted):
        raise ArithmeticError('opposite modular sign was not excluded')
    return {'full_Fock_dimension':4,'density_positive_and_trace_one':True,
            'relative_entropy_at_theta_pi_over_two':'3*log(4)/5',
            'modular_annihilator_sign_checked':True,
            'opposite_sign_negative_control_rejected':True,
            'entropy_formula':'2*(1-cos(theta))*(<C>*<k>-<Ck>)',
            'absolute_continuum_entropy_claimed':False}
