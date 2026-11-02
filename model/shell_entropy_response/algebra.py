"""Exact noncommuting, full-CAR controls for the rank-one entropy family."""
import sympy as s
from model.shrinking_modular_entropy.algebra import finite_example


def certificate():
    e=finite_example();A=e['A'];f=e['F'];X=f*f.T;I=s.eye(2)
    t,z=s.symbols('t z', positive=True)
    q=(f.T*(A+t*I).inv()*f)[0]
    r=(f.T*(I-A+t*I).inv()*f)[0]
    forward=(f.T*(A-z*X+t*I).inv()*f)[0]
    backward=(f.T*(I-A+z*X+t*I).inv()*f)[0]
    m=(f.T*f)[0]
    theta=e['theta'];rho=e['rho']
    checks={
        'source_projection':e['P']**2==e['P'],
        'source_exterior_mode_unit':(e['g'].T*e['g'])[0]==1,
        'full_Fock_rank_one_Wick':s.simplify(e['actual_covariance']-e['expected_covariance'])==s.zeros(2),
        'full_Fock_state_affine':rho.diff(theta,2)==s.zeros(4),
        'off_diagonal_block_factorization':e['B']*e['B'].T==A*(I-A),
        'noncommuting_covariance_perturbation':A*X!=X*A,
        'particle_resolvent_negative_sign':s.cancel(forward-q/(1-z*q))==0,
        'hole_resolvent_positive_sign':s.cancel(backward-r/(1+z*r))==0,
        'natural_endpoint_sum_at_most_one':s.cancel(q.subs(t,0)+r.subs(t,0))<=1,
        'positive_natural_mass':0<m<=s.Rational(1,4),
    }
    # Compression defect is a positive Gram matrix; this survives finite CAR
    # exhaustion without a commuting projection assumption.
    C=s.Matrix([[s.Rational(1,2),s.Rational(1,8),0],
                [s.Rational(1,8),s.Rational(1,2),s.Rational(1,8)],
                [0,s.Rational(1,8),s.Rational(1,2)]])
    P=s.diag(1,1,0)
    checks['compression_defect_Gram_identity']=(P*C*P-(P*C*P)**2-P*(C-C*C)*P)==P*C*(s.eye(3)-P)*C*P
    # One-mode natural family saturates chi(0)<=m exactly.
    x=s.symbols('x',positive=True)
    checks['one_mode_hessian_upper_saturated']=s.cancel((x*(1-x))**2*(1/x+1/(1-x))-x*(1-x))==0
    if not all(bool(v) for v in checks.values()):
        raise ArithmeticError('exact entropy-response algebra failed')
    return {'checks':{k:bool(v) for k,v in checks.items()},
            'example_mass':str(m),'example_endpoint_sum':str(s.cancel(q.subs(t,0)+r.subs(t,0))),
            'full_source_Fock_dimension':16,'regional_Fock_dimension':4,
            'matrix_commutativity_assumed':False}
