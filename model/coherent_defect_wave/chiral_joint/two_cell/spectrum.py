"""Complete physical Gauss blocks and certified finite energy intervals."""
from functools import lru_cache
from math import gcd,isqrt
from collections import defaultdict
import numpy as np
import scipy.linalg as la
import sympy as s
from .basis import diagonal_basis,physical_basis,coordinate_rows,block_labels
from .hamiltonian import matrix_on


@lru_cache(maxsize=None)
def physical_block(boundary,L,twoSz):
    keys=diagonal_basis(boundary,L,twoSz)
    b,g=physical_basis(boundary,L,twoSz)
    if b.cols==0:return keys,b,g,s.zeros(0),s.zeros(0)
    h=matrix_on(keys)
    # This gauge-invariant diagonal unitary makes the edge coefficients real.
    # The physical embedding is D B, never B alone after this change.
    phases=[s.I**((key[0]&63).bit_count()) for key in keys]
    hr=s.ImmutableSparseMatrix(h.rows,h.cols,{(i,j):s.expand(s.conjugate(phases[i])*v*phases[j])
                                             for (i,j),v in h.todok().items()})
    assert all(v.is_Rational for v in hr.values())
    hb=hr*b
    generator=s.ImmutableSparseMatrix(hb.extract(coordinate_rows(boundary,L,twoSz),range(b.cols)))
    assert b*generator==hb
    k=s.ImmutableSparseMatrix(g*generator)
    assert k==k.T and generator.T*g==g*generator
    embedding=s.ImmutableSparseMatrix(b.rows,b.cols,{(i,j):phases[i]*v for (i,j),v in b.todok().items()})
    assert embedding.H*embedding==g
    return keys,embedding,g,k,generator


def components(k,g):
    parent=list(range(k.rows))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for m in (k,g):
        for i,j in m.todok():parent[find(i)]=find(j)
    out=defaultdict(list)
    for i in range(k.rows):out[find(i)].append(i)
    return tuple(tuple(v) for v in out.values())


@lru_cache(maxsize=None)
def numerical_block(boundary,L,twoSz):
    keys,b,g,k,h=physical_block(boundary,L,twoSz)
    if not b.cols:return np.array([]),np.zeros((0,0))
    gf=np.asarray(g.tolist(),dtype=float);kf=np.asarray(k.tolist(),dtype=float)
    values,vectors=la.eigh(kf,gf,driver='gvd')
    return values,vectors


def diagnostic():
    rows=[]
    for boundary in (0,-1,1):
        for L,twoSz in block_labels(boundary):
            keys,b,g,k,h=physical_block(boundary,L,twoSz)
            if not b.cols:continue
            values,_=numerical_block(boundary,L,twoSz)
            row={'boundary':boundary,'L':L,'twoSz':twoSz,'dim':b.cols,
                 'minimum':float(values[0]),'components':[len(c) for c in components(k,g)]}
            rows.append(row)
            print(row,flush=True)
    return rows


if __name__=='__main__':diagnostic()


def _integer_matrix(m):
    denominator=s.ilcm(1,*(v.q for v in m.values())) if m.values() else s.Integer(1)
    entries={(i,j):int(v*denominator) for (i,j),v in m.todok().items()}
    return int(denominator),entries


def _sparse_integer_product(entries,z):
    result=np.zeros(z.shape,dtype=object)
    for (i,j),coefficient in entries.items():result[i,:]+=coefficient*z[j,:]
    return result


def exact_eigenbasis_error(k,g,values,vectors,bits=44):
    """All-eigenvalue certificate using exact integer arithmetic.

    The floating solver proposes a complete basis.  Exact Gram and residual
    estimates certify it; no floating rank or residual is trusted.
    """
    n=k.rows
    if n==0:return s.Integer(0),s.Integer(0)
    scale=1<<bits
    z=np.rint(vectors*scale).astype(np.int64).astype(object)
    d=np.rint(values*scale).astype(np.int64).astype(object)
    kd,ki=_integer_matrix(k);gd,gi=_integer_matrix(g)
    kz=_sparse_integer_product(ki,z);gz=_sparse_integer_product(gi,z)
    gram=z.T@gz
    gram_denominator=gd*scale**2
    for j in range(n):gram[j,j]-=gram_denominator
    delta=s.Rational(max(sum(abs(v) for v in row) for row in gram),gram_denominator)
    assert delta<1
    # R=K Q-G Q D; Q=Z/scale, D=diag(d)/scale.
    residual=gd*scale*kz-kd*gz*d[None,:]
    residual_denominator=kd*gd*scale**2
    squared=int(sum(v*v for v in z.flat))*int(sum(v*v for v in residual.flat))
    upper=isqrt(squared)
    if upper*upper<squared:upper+=1
    product_bound=s.Rational(upper,scale*residual_denominator)
    maximum=s.Rational(max(abs(v) for v in d),scale)
    error=s.factor((product_bound+2*maximum*delta)/(1-delta))
    return error,delta


@lru_cache(maxsize=None)
def certified_block(boundary,L,twoSz):
    keys,b,g,k,generator=physical_block(boundary,L,twoSz)
    if not b.cols:return {'boundary':boundary,'L':L,'twoSz':twoSz,'physical_dimension':0,'energy_intervals':[]}
    values,vectors=numerical_block(boundary,L,twoSz)
    error,gram_error=exact_eigenbasis_error(k,g,values,vectors)
    # A fixed outward envelope avoids tying the receipt to last-bit solver
    # choices. Every one of these rational intervals is checked exactly.
    scale=1<<32;radius=s.Rational(1,1<<28)
    fine=1<<44
    centers=np.rint(values*scale).astype(np.int64)
    finecenters=np.rint(values*fine).astype(np.int64)
    for c,d in zip(centers,finecenters,strict=True):
        assert abs(s.Rational(int(c),scale)-s.Rational(int(d),fine))+error<=radius
    assert gram_error<radius
    return {'boundary':boundary,'L':L,'twoSz':twoSz,'physical_dimension':b.cols,
            'diagonal_dimension':len(keys),'exact_physical_invariance':True,
            'energy_center_denominator':scale,'energy_centers_numerator':[int(c) for c in centers],
            'energy_interval_radius':str(radius),'energy_intervals':[
                [str(s.Rational(int(c),scale)-radius),str(s.Rational(int(c),scale)+radius)] for c in centers],
            'full_eigenbasis_Gram_error_upper_bound':str(radius),
            'certificate':'exact integer full-basis Gram and residual; Hermitian minmax',
            'numerical_rank_or_vacuum_Krylov_assumption':False}


def sector_certificate(boundary):
    blocks=[certified_block(boundary,L,spin) for L,spin in block_labels(boundary)]
    nonempty=[b for b in blocks if b['physical_dimension']]
    intervals=[(s.Rational(pair[0]),s.Rational(pair[1]),b['L'],b['twoSz'],j)
               for b in nonempty for j,pair in enumerate(b['energy_intervals'])]
    best=min(intervals,key=lambda r:r[1])
    competitors=[r for r in intervals if r!=best]
    ground_lower=min(r[0] for r in intervals);ground_upper=best[1]
    candidates=[r for r in intervals if r[0]<=ground_upper]
    separated=[r[0] for r in intervals if r[0]>ground_upper]
    gap_lower=min(separated)-ground_upper if separated else None
    return {'boundary':boundary,'full_physical_dimension':sum(b['physical_dimension'] for b in blocks),
            'complete_conserved_blocks':blocks,'global_ground_interval':[str(ground_lower),str(ground_upper)],
            'ground_window_candidate_labels':[{'L':r[2],'twoSz':r[3],'eigen_index':r[4]} for r in candidates],
            'ground_window_rank':len(candidates),
            'ground_window_to_rest_gap_lower_bound':str(gap_lower),
            'ground_is_certifiably_simple':len(candidates)==1,
            'all_occupation_sectors_in_Gauss_domain_included':True}
