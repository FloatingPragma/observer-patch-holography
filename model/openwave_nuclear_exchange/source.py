"""Exact integer components and inherited complete Casimir reduction."""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_six_quark import source as native
from model.openwave_calibrated_hadrons.source import parameters,weights


@lru_cache(None)
def components():
    old=native.reduction();src=old['basis'];n=len(src);keep=old['keep']
    C=native.colour48(src,(0,))+native.colour48(src,(1,));S=sp.csr_matrix((n,n),dtype=np.int64);I=S.copy()
    for x in range(2):
        for spin in (True,False):
            z=native.local_number(src,x,1 if spin else 2);plus=native.raising(src,(x,),spin)
            term=sp.diags(z*z+2*z,dtype=np.int64)+4*(plus.T@plus)
            if spin:S+=term
            else:I+=term
    hop=[]
    for i,(x,ss,f,c)in enumerate(native.MODES):hop.extend(((i,i,1),(i,native.INDEX[1-x,ss,f,c],-1)))
    K=native.op(hop,src,src);n0=native.local_number(src,0);B9=sp.diags(n0*n0+(6-n0)**2,dtype=np.int64)
    Ig=native.raising(src,(0,1),False);Ig4=4*(Ig.T@Ig)
    full={'N':6*sp.eye(n,dtype=np.int64,format='csr'),'D':sp.csr_matrix((n,n),dtype=np.int64),
          'K':K,'C48':C,'S4':S,'I4':I,'B9':B9,'Ig4':Ig4}
    parts={k:a.tocsr()[keep][:,keep].tocsr()for k,a in full.items()}
    P=old['Pnum'];B=old['Binteger'];reduced={}
    for key,A in parts.items():
        assert (native.exact_product(A,P)-native.exact_product(P,A)).nnz==0
        reduced[key]=native.exact_product(B.T,native.exact_product(A,B)).toarray()
    assert (old['Hnum']-(48*parts['N']+6*K[keep][:,keep]+C[keep][:,keep]+3*S[keep][:,keep]+3*I[keep][:,keep])).nnz==0
    return {**old,'parts':parts,'full_parts':full,'integer_reduced_parts':reduced}


@lru_cache(None)
def reduction():
    old=components();p=parameters();w=weights(p);A=sy.zeros(27)
    for key,value in w.items():A+=sy.Rational(value.numerator,value.denominator)*sy.Matrix(old['integer_reduced_parts'][key].tolist())
    L,D=old['Gram'].LDLdecomposition(hermitian=True);Li=L.inv();T=Li*A*Li.T;g=list(D.diagonal())
    H=sy.Matrix(27,27,lambda i,j:T[i,j]/sy.sqrt(g[i]*g[j]))
    h=np.asarray(H.evalf(18),float)
    return {**old,'H_exact':H,'h':h,'parameters':p,
            'full_H':sum(float(v)*old['full_parts'][k]for k,v in w.items()),
            'cartan_H':sum(float(v)*old['parts'][k]for k,v in w.items())}


def exact_matvec(A,v):
    A=A.tocsr();return [sum(int(A.data[k])*int(v[int(A.indices[k])])for k in range(A.indptr[i],A.indptr[i+1]))for i in range(A.shape[0])]


def rational_action(d,v):
    result=[F(0)for _ in v]
    for key,value in weights(parameters()).items():
        if not value:continue
        for i,x in enumerate(exact_matvec(d['parts'][key],v)):result[i]+=value*x
    return result
