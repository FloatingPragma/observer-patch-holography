"""Full70928-state one-charge trial solve; certification is a separate step.

No photon subspace is discarded. The finite Peierls exponential is its exact
degree-eight interpolation polynomial on all total-occupation<=3 states.
The real basis uses U=exp(-i*pi*Nph/2): U(i A_e)U†=a†(f_e)-a(f_e).
The number-conserving photon Hamiltonian and occupation/position projectors
are unchanged by this exact unitary change of basis.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import itertools, json, math, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh
import sympy as s
import mpmath as mp

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from model.affine_gauss_mobile import model as parent
from model.affine_gauss_spectrum import model as compression
HERE=Path(__file__).resolve().parent


def ivalue(x):
    x=s.sympify(x)
    if x.is_Pow and x.args[1].is_Rational and x.args[1].q==4:
        return mp.iv.sqrt(mp.iv.sqrt(ivalue(x.args[0])))**int(x.args[1].p)
    if x.is_Add:return sum(map(ivalue,x.args),mp.iv.mpf(0))
    if x.is_Mul:return math.prod(map(ivalue,x.args))
    return parent.ivalue(x)


def rounded_float(x):
    lo,hi=parent.endpoint(x),parent.endpoint(x,True)
    value=float((lo+hi)/2)
    error=max(abs(F.from_float(value)-lo),abs(F.from_float(value)-hi))
    return value,np.nextafter(float(error),np.inf)


def iv_solve(A,b):
    A=[list(a)+[v] for a,v in zip(A,b)]
    n=len(b)
    for i in range(n):
        p=A[i][i]
        assert not parent.endpoint(p)<=0<=parent.endpoint(p,True)
        A[i]=[x/p for x in A[i]]
        for j in range(n):
            if j!=i:
                c=A[j][i];A[j]=[x-c*y for x,y in zip(A[j],A[i])]
    return [A[j][-1] for j in range(n)]


def peierls_coefficients(variance):
    iv=mp.iv
    # D has imaginary eigenvalues0,±i sigma,±i sqrt3 sigma,
    # ±i sqrt(3−sqrt6) sigma and ±i sqrt(3+sqrt6) sigma.
    xs=[variance*x for x in (iv.mpf(1),iv.mpf(3),3-iv.sqrt(6),3+iv.sqrt(6))]
    even=iv_solve([[(-x)**k for k in range(1,5)] for x in xs],
                  [iv.cos(iv.sqrt(x))-1 for x in xs])
    odd=iv_solve([[(-x)**k for k in range(4)] for x in xs],
                 [iv.sin(iv.sqrt(x))/iv.sqrt(x) for x in xs])
    coeff=[iv.mpf(1)]+[iv.mpf(0) for _ in range(8)]
    for k,v in enumerate(even,1):coeff[2*k]=v
    for k,v in enumerate(odd):coeff[2*k+1]=v
    return coeff


def algebraic_frame(*, verbose=False):
    if verbose: print('Constructing exact rational transverse frame',flush=True)
    mp.iv.dps=80
    geom=parent.parent.geometry();K,PT,powers,polys=geom[6:10]
    from sympy.polys.matrices import DomainMatrix
    _,pivots=DomainMatrix.from_Matrix(PT).to_field().rref()
    columns=[PT[:,j] for j in pivots];qs=[];norms=[]
    def primitive(q):
        denominator=s.ilcm(*[v.q for v in q])
        numerators=[int(v*denominator) for v in q]
        common=math.gcd(*numerators)
        assert common
        return s.Matrix([v//common for v in numerators])
    for col in columns:
        q=primitive(col)
        for old,norm in zip(qs,norms):q=primitive(q-old*((old.T*q)[0]/norm))
        norm=(q.T*q)[0];assert norm>0
        qs.append(q);norms.append(norm)
        if verbose: print('Exact frame direction',len(qs),flush=True)
    Q=s.Matrix.hstack(*qs)
    assert len(qs)==30 and Q.T*Q==s.diag(*norms)
    assert PT*Q==Q
    roots=parent.parent.ROOTS
    omega_coefs=[];field_coefs=[]
    for degree in range(9):
        omega_coefs.append(sum((ivalue(polys[j].nth(degree))*ivalue(s.sqrt(roots[j])) for j in range(1,9)),mp.iv.mpf(0)))
        field_coefs.append(sum((ivalue(polys[j].nth(degree))*ivalue(roots[j]**s.Rational(-1,4))/mp.iv.sqrt(2) for j in range(1,9)),mp.iv.mpf(0)))
    QTpower=[Q.T*p for p in powers[:9]]
    QTKQ=[p*Q for p in QTpower]
    if verbose: print('Exact projected powers ready',flush=True)
    om=np.zeros((30,30));ome=np.zeros_like(om)
    ff=np.zeros((42,30));ffe=np.zeros_like(ff)
    for i in range(30):
        ni=ivalue(norms[i]);normal=mp.iv.sqrt(ni)
        for e in range(42):
            value=sum((c*ivalue(p[i,e]) for c,p in zip(field_coefs,QTpower)),mp.iv.mpf(0))/normal
            ff[e,i],ffe[e,i]=rounded_float(value)
        for j in range(i,30):
            value=sum((c*ivalue(p[i,j]) for c,p in zip(omega_coefs,QTKQ)),mp.iv.mpf(0))/mp.iv.sqrt(ni*ivalue(norms[j]))
            om[i,j],ome[i,j]=rounded_float(value);om[j,i]=om[i,j];ome[j,i]=ome[i,j]
        if verbose: print('Outward coefficient row',i+1,flush=True)
    coefs=[];coeferr=[]
    for edge in (0,12):
        variance=ivalue(parent.exact_peierls_data(edge)[2])
        rounded=[rounded_float(x) for x in peierls_coefficients(variance)]
        coefs.append([v for v,e in rounded]);coeferr.append([e for v,e in rounded])
    charges,_=parent.parent.charge_geometry()
    energy=np.array([float(row['total_energy_above_actual_neutral_ground']) for row in charges])
    energye=np.array([np.nextafter(float(abs(F.from_float(value)-F(row['total_energy_above_actual_neutral_ground']))),np.inf) for value,row in zip(energy,charges)])
    result={'omega':om,'omega_error':ome,'field':ff,'field_error':ffe,
            'peierls':np.array(coefs),'peierls_error':np.array(coeferr),
            'static_energy':energy,'static_energy_error':energye}
    return result


def occupations(modes=30,cutoff=3):
    return [word for n in range(cutoff+1) for word in itertools.combinations_with_replacement(range(modes),n)]


def matrices(frame):
    words=occupations();index={w:i for i,w in enumerate(words)};n=len(words)
    assert n==5456
    rows=[];cols=[];modes=[];values=[]
    for j,word in enumerate(words):
        if len(word)==3:continue
        for mode in range(30):
            rows.append(index[tuple(sorted(word+(mode,)))]);cols.append(j)
            modes.append(mode);values.append(math.sqrt(word.count(mode)+1))
    rows,cols,modes,values=map(np.array,(rows,cols,modes,values))
    Ds=[]
    for f in frame['field']:
        creation=sp.coo_matrix((values*f[modes],(rows,cols)),shape=(n,n)).tocsr()
        Ds.append(creation-creation.T)
    rr=[];cc=[];vv=[]
    for col,word in enumerate(words):
        for j in set(word):
            left=list(word);left.remove(j);nj=word.count(j)
            for i in range(30):
                row=index[tuple(sorted(left+[i]))]
                rr.append(row);cc.append(col)
                vv.append(frame['omega'][i,j]*math.sqrt(nj*(left.count(i)+1)))
    photons=sp.coo_matrix((vv,(rr,cc)),shape=(n,n)).tocsr()
    # Symmetry is a check, never a post-hoc replacement of the source matrix.
    assert np.max(np.abs((photons-photons.T).data),initial=0)==0
    return Ds,photons


class ActualOperator:
    def __init__(self,frame):
        self.frame=frame;self.Ds,self.photons=matrices(frame);self.calls=0
    def apply(self,value):
        shape=np.shape(value);x=np.asarray(value).reshape(13,5456)
        y=(self.photons@x.T).T+self.frame['static_energy'][:,None]*x
        for edge,(u,v) in enumerate(parent.actual_edges()):
            D=self.Ds[edge];c=self.frame['peierls'][int(edge>=12)]
            incoming=np.column_stack((x[u],x[v]));out=c[8]*incoming
            for k in range(7,-1,-1):
                out=D@out;out[:,1]*=-1
                out+=c[k]*incoming
            y[v]-=out[:,0]/24;y[u]-=out[:,1]/24
        self.calls+=1
        return y.reshape(shape)
