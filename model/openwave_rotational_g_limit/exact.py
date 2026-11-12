"""Exact complete-shell identities underlying rational magnetic error bounds."""
from functools import lru_cache
import sympy as s
from .limit import oscillator


@lru_cache(maxsize=1)
def certificate():
    basis=oscillator()['basis'];n=len(basis);where={v:i for i,v in enumerate(basis)};p=[];x=[]
    for axis in range(3):
        a=s.zeros(n)
        for col,v in enumerate(basis):
            if v[axis]:
                w=list(v);w[axis]-=1;a[where[tuple(w)],col]=s.sqrt(v[axis])
        p.append(-s.I*(a-a.H)/s.sqrt(2));x.append((a+a.H)/s.sqrt(2))
    sig=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    gamma=[s.kronecker_product(sig[0],z)for z in sig]
    A=sum((s.kronecker_product(g,pj)for g,pj in zip(gamma,p)),s.zeros(4*n))
    E=(-s.kronecker_product(gamma[0],x[1])+s.kronecker_product(gamma[1],x[0]))/2
    T=A*A;Q=E*E;W=A*E+E*A
    assert (T*T-s.Rational(5,2)*T).applyfunc(s.simplify)==s.zeros(4*n)
    assert (Q*(Q-s.eye(4*n)/4)*(Q-s.eye(4*n)/2)).applyfunc(s.simplify)==s.zeros(4*n)
    assert W[0,0]==1 and T[0,0]==s.Rational(3,2)
    assert Q[0,0]==s.Rational(1,4)
    assert T==T.H and Q==Q.H
    return {'rows':4*n,'squared_kinetic_polynomial':'T(T-5/2 I)=0',
        'squared_probe_polynomial':'Q(Q-1/4 I)(Q-1/2 I)=0',
        'positive_squares':True,'kinetic_norm_squared_upper':'5/2','probe_norm_squared_upper':'1/2',
        'squared_kinetic_probe_norm_upper':'sqrt(5)<3',
        'vacuum_T':'3/2','vacuum_dT_dB':1,'vacuum_quadratic_probe':'1/4',
        'magnetic_spin_coefficient':'the vacuum derivative is native Sigma_z, not an inserted Pauli term'}
