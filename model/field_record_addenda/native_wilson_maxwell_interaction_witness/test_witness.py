from pathlib import Path
from fractions import Fraction as F
import sys
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
from . import model as m


def test_source_symbol_and_complete_current_interpolation():
    assert m.symbolic_controls()['current_interpolation_weight_bound']==515
    alpha,beta=m.family.matrices();B=-(beta+s.I*alpha[0])/2
    def car(mode):
        return s.MutableSparseMatrix(256,256,{(word^(1<<mode),word):(-1)**((word&((1<<mode)-1)).bit_count())
            for word in range(256) if word>>mode&1})
    a=[car(j) for j in range(8)]
    hop=s.zeros(256,256)
    for x in range(4):
        for y in range(4):
            if B.H[x,y]:hop+=B.H[x,y]*a[4+x].H*a[y]
    J=s.SparseMatrix(s.I*(hop-hop.H));I=s.eye(256)
    assert J.H==J and J*(J*J-I)*(J*J-4*I)==s.zeros(256)


def test_independent_bulk_quadrature_and_analytic_identity():
    # Numerical diagnostic only; the authority is the symbolic integration
    # identity and E<=7 bound, not this quadrature.
    q=(np.arange(48)+.5)*2*np.pi/48
    c=np.cos(q);sn=np.sin(q)
    A=4-c[None,:,None]-c[None,None,:]
    C=c[:,None,None];S2=sn[None,:,None]**2+sn[None,None,:]**2
    E=np.sqrt(A*A+1+S2-2*A*C)
    first=np.mean((A*C-1)/E)
    second=-np.mean(((A*C-1)**2+S2)/E**3)
    assert abs(first-second)<1e-10 and first<-1/343
    assert E.max()<=7 and E.min()>=1


def test_exact_transverse_plus_longitudinal_current_identity():
    # The literal cube projector identity checks the sign/orientation seam;
    # its full box proof is incidence algebra, not a fitted covariance.
    data=m.family.field.cube();D=data['D'];G=data['G'];PT=data['PT']
    assert PT+D*G*D.T==s.eye(len(data['edges']))
    assert data['C']*D==s.zeros(data['C'].rows,D.cols)
    assert D[0,0]==-1 and sum(D.row(0))==0


def test_complete_records_and_source_only_margins():
    value=m.scientific();b=value['bounds_exact']
    assert F(b['finite_complete_law_TV_separation_lower'])>F(17,10**11)
    assert F(b['finite_interior_dressed_hopping_magnitude_lower'])>F(1,345)
    # Algebraic fairness/TV lower: no record gets omitted, and comparing
    # the different stored input bits themselves is not the statistic.
    q=F(9,10);pk=F(3,5);pn=F(2,5);p0=F(9,20)
    actual=[q*pk/2,q*(1-pk)/2,q*pn/2,q*(1-pn)/2,(1-q)/2,(1-q)/2]
    free=[q*p0/2,q*(1-p0)/2,q*p0/2,q*(1-p0)/2,(1-q)/2,(1-q)/2]
    assert sum(actual)==sum(free)==1
    assert sum(abs(x-y) for x,y in zip(actual,free))/2>=q*abs(pk-pn)/2


def test_capacity_interval_encloses_independent_factorial_formula():
    import mpmath as mp
    from math import comb
    mp.mp.dps=100
    interval=m.capacity(1,3,2)['log_photon_dimension_nats_interval']
    lower,upper=[mp.mpf(F(x).numerator)/F(x).denominator for x in interval]
    assert lower<mp.log(comb(8,5))<upper
    n=8200;V=(n+1)**3;E=3*n*(n+1)**2;dim=E-V+1;K=363595359
    interval=m.capacity(n,K,7)['log_photon_dimension_nats_interval']
    lower,upper=[mp.mpf(F(x).numerator)/F(x).denominator for x in interval]
    direct=mp.loggamma(K+dim+1)-mp.loggamma(K+1)-mp.loggamma(dim+1)
    assert lower<direct<upper
