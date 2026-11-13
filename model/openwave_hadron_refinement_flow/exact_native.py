"""Exact native rational-frame check of the compact symbol and coupling square."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import sympy as sy
from .bloch import compact_frame


@lru_cache(None)
def certificate():
    z=compact_frame(4);d=z['d'];p=d['parameters']
    ids=np.flatnonzero((d['nu']==1)&(d['ns']==2));g=d['gram'][ids]
    P=z['P'][ids];R=sy.zeros(len(ids),P.shape[1])
    for j in range(P.shape[1]):
        nz=np.flatnonzero(abs(P[:,j])>1e-10);base=P[nz[0],j]/np.sqrt(g[nz[0]])
        for i in nz:
            candidate=F(float(P[i,j]/np.sqrt(g[i])/base)).limit_denominator(100)
            R[int(i),j]=sy.Rational(candidate)
    G=sy.diag(*map(int,g));norm=R.T*G*R
    assert norm.is_diagonal() and all(norm[j,j]>0 for j in range(norm.rows))
    for j in range(P.shape[1]):
        reconstructed=np.sqrt(g)*np.array(R[:,j]).astype(float).ravel()/np.sqrt(float(norm[j,j]))
        assert np.max(abs(reconstructed-P[:,j]))<2e-13
    def integer_part(part):
        a=part.tocsr()[ids][:,ids].toarray()*np.sqrt(np.outer(g,g))
        integer=np.rint(a).astype(np.int64)
        assert np.max(abs(a-integer))<3e-10
        return sy.MutableSparseMatrix(integer.tolist())
    H=sy.zeros(len(ids))
    for name,part in d['parts'].items():H+=sy.Rational(d['weights'][name])*integer_part(part)
    for x in (1,2):
        H-=sy.Rational(p['spin_gradient']/4)*integer_part(4*d['angular']((x,),True))
        H-=sy.Rational(p['eta']/4)*integer_part(4*d['angular']((x,),False))
        H-=sy.Rational(p['t'])*sy.diag(*[int(g[i]*d['n'][ids[i],x]) for i in range(len(ids))])
    for spin in (True,False):
        S4=integer_part(4*d['angular']((0,1,2,3),spin))
        assert (S4-3*G)*R==sy.zeros(*R.shape)
    h=R.T*H*R
    W=H*R
    square=W.T*G.inv()*W-h*norm.inv()*h
    scale=sy.diag(*[1/sy.sqrt(norm[j,j]) for j in range(norm.rows)])
    hn=sy.simplify(scale*h*scale);T=sy.simplify(scale*square*scale)
    groups=z['groups'];edges={x:sum([inds for kind,xx,n,inds in groups if kind=='edge' and xx==x],[]) for x in range(3)}
    t=sy.Rational(p['t']);T0=T.extract(edges[1],edges[1]);T1=T.extract(edges[1],edges[2])
    assert T0==3*t*t*sy.eye(4)
    expected=sy.zeros(4);expected[2,0]=expected[3,1]=t*t/2;expected[2,1]=expected[3,0]=3*t*t/2
    assert T1==expected
    assert T.extract(edges[0],edges[2])==sy.zeros(4)
    assert T.extract(range(4),range(T.cols))==sy.zeros(4,T.cols)
    return {'exact_frame_rank':R.cols,'complete_fixed_projection_native_dimension':R.rows,
            'exact_frame_Gram':list(map(str,norm.diagonal())),
            'exact_native_spin_and_isospin_intertwiners':True,
            'exact_coupling_square_identity':True,'off_diagonal_hopping_uses_all_span_two_outputs':True,
            'compact_H_float_check':float(np.max(abs(np.array(hn).astype(float)-z['H'])))}
