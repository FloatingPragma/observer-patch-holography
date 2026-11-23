"""Exact cut-invisible stationary family and sharp input memory at half filling."""
import argparse
from itertools import combinations
from fractions import Fraction
import json
from math import gcd, lcm
from functools import reduce
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csr_matrix, diags

from model.executable_model import cell, federation
from .measure import SEED
from ..equilibrium_refinement.measure import seam

HERE = Path(__file__).resolve().parent
SHOTS = 2**32


def source():
    host = federation(0)
    geometry = cell()
    cut = tuple(np.flatnonzero(host.degree[:12] == 6))
    words = list(combinations(range(12),6))
    index = {w:i for i,w in enumerate(words)}
    f = seam()
    assert round(2*f[1,1])==65 and round(2*f[1,2])==-1
    rows, cols, diag = [], [], []
    for i,w in enumerate(words):
        ends=set(w)
        active=0
        for a,b in geometry.seams:
            if (a in ends) ^ (b in ends):
                rows.append(index[tuple(sorted(ends ^ {a,b}))]); cols.append(i)
                active+=1
        diag.append(65*active+64*sum(6-int(host.degree[p]) for p in w))
    shift=round((max(diag)+min(diag))/2)
    g = (coo_matrix((np.full(len(rows),-1,dtype=np.int64),(rows,cols)),shape=(len(words),)*2)
         + diags(np.array(diag,dtype=np.int64)-shift,dtype=np.int64)).tocsr()
    graph=nx.Graph(geometry.seams)
    aut=[s for s in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter()
         if set(s[p] for p in cut)==set(cut)]
    reflections=sorted([s for s in aut if any(s[p]!=p for p in cut)
                        and all(s[s[p]]==p for p in range(12))],key=lambda s:tuple(s[p] for p in cut))
    s,t=reflections[:2]
    bridge=next(a for a in aut if all(a[s[p]]==t[a[p]] for p in range(12)))
    def perm(q):
        target=np.array([index[tuple(sorted(q[p] for p in w))] for w in words])
        return coo_matrix((np.ones(len(words),dtype=np.int64),(target,np.arange(len(words)))),
                          shape=g.shape).tocsr()
    sm,tm,r=map(perm,(s,t,bridge))
    assert (r@sm@r.T!=tm).nnz==0 and all((g@q!=q@g).nnz==0 for q in (sm,tm,r))
    masks=np.array([sum(1<<p for p in w) for w in words]); mask=sum(1<<p for p in cut)
    _,inside=np.unique(masks & ~mask,return_inverse=True)
    _,edge=np.unique(masks & mask,return_inverse=True)
    return words,cut,g,shift,sm-tm,r,inside,edge


def partial_cut(matrix,inside,edge):
    entries=matrix.tocoo()
    keep=inside[entries.row]==inside[entries.col]
    out=np.zeros((edge.max()+1,edge.max()+1),dtype=np.int64)
    np.add.at(out,(edge[entries.row[keep]],edge[entries.col[keep]]),entries.data[keep])
    return out


def witness(g, y, inside, edge):
    powers=[y]
    couplings=[]
    for j in range(6):
        if j:
            powers.append(g@powers[-1])
        bound=max(abs(powers[-1]).sum(axis=1).A.ravel())
        if bound>2**48:
            raise ArithmeticError('Use a higher-precision integer bound before extending the search')
        couplings.append(partial_cut(powers[-1],inside,edge).ravel())
        system=sp.Matrix(np.column_stack(couplings))
        for sol in system.nullspace():
            denom=lcm(*(int(s.q) for s in sol))
            coeff=[int(v*denom) for v in sol]
            divisor=reduce(gcd,coeff)
            coeff=[c//divisor for c in coeff]
            assert system @ sp.Matrix(coeff) == sp.zeros(system.rows,1)
            # Accumulate the solution using Python integers even when the nullspace
            # coefficients exceed machine-word size. Only normalize after the audit.
            integer=np.zeros(g.shape,dtype=object)
            for c,p in zip(coeff,powers):
                entries=p.tocoo()
                integer[entries.row,entries.col] += c*entries.data.astype(object)
            diagonal=integer.diagonal()
            group_sums=[sum(diagonal[inside==i]) for i in range(inside.max()+1)]
            if any(v != 0 for v in group_sums):
                assert sum(diagonal)==0 and np.array_equal(integer,integer.T)
                scale=int(max(abs(integer).sum(axis=1)))
                x=csr_matrix(integer.astype(float)/scale)
                assert max(abs((x@g-g@x).data),default=0)<1e-11
                group_sizes=np.bincount(inside)
                pooled=[Fraction(int(z),int(n)) for z,n in zip(group_sums,group_sizes,strict=True)]
                pooled_scale=max(abs(v) for v in pooled)
                overlap=sum(z*w for z,w in zip(group_sums,pooled,strict=True))/pooled_scale
                frobenius=sum(int(n)**2 for n in integer.ravel() if n)
                lower=Fraction(3,2*len(inside))*overlap**2/frobenius
                return x,coeff,scale,lower
    raise ArithmeticError('No cut-invisible observable witness in the enrolled search')


def run():
    words,cut,g,shift,y,r,inside,edge=source()
    d=len(words)
    x,coeff,scale,lower_exact=witness(g,y,inside,edge)
    base=np.bincount(inside)/d
    dr=np.bincount(inside,weights=x.diagonal())/d
    transferred=dr[dr>0].sum()
    direction=dr/transferred
    states=[]
    rng=np.random.default_rng([SEED,601])
    for a in (-.75,0.,.75):
        probs=base+a*dr
        counts=rng.multinomial(SHOTS, probs/probs.sum())
        freq=counts/SHOTS
        states.append(dict(alpha=a, coarse_counts=counts.tolist(), moved_probability=abs(a)*transferred,
                           normalized_signed_bulk_slope=float(-direction@np.log(probs)),
                           measured_slope=float(-direction@np.log(freq)),
                           slope_se=float(np.sqrt(np.sum(direction**2/freq)/SHOTS))))
    o=dr/base
    o/=max(abs(o))
    launch=o[inside]
    # O is equivariant under port relabeling and is a difference of conjugate diagonal scores.
    assert abs(launch.sum())<1e-10 and np.max(abs(np.bincount(edge,weights=launch)))<1e-10
    conserved_overlap=float(x.diagonal()@launch)
    lower=1.5*conserved_overlap**2/(d*float(np.sum(x.data.astype(float)**2)))
    np.testing.assert_allclose(lower,float(lower_exact),rtol=1e-12,atol=1e-18)
    energies,v=eigh(g.toarray()/2,driver='evd')
    observations=[]
    used_spectrum=float(0.)
    groups=np.split(np.arange(d),np.flatnonzero(np.diff(energies)>1e-8)+1)
    for ix in groups:
        used_spectrum=max(used_spectrum,abs(float(np.sum(launch[:,None]*v[:,ix]**2)))/d)
    for time in (0,1,4,16):
        u=(v*np.exp(-1j*time*energies)) @ v.T
        diag_bias=abs(u)**2 @ launch
        samples=[]
        for a in (-.75,.75):
            sharp=(1+a*diag_bias)/d
            probs=np.bincount(inside,weights=sharp)
            counts=rng.multinomial(SHOTS,probs/probs.sum())
            freq=counts/SHOTS
            score=freq@o
            samples.append(dict(alpha=a,counts=counts.tolist(),observed_score=float(score),
                                se=float(np.sqrt((freq@(o-score)**2)/SHOTS)),
                                expected_score=float(probs@o)))
        observations.append(dict(exposure=time,preparations=samples))
    return dict(shots=SHOTS,seed=SEED,n_words=d,cut_addresses=list(map(int,cut)),center_shift_twice_H=shift,
                integer_polynomial=list(coeff),exact_bound_L=scale,cut_partial_trace_exactly_zero=True,
                conjugating_permutation_verified=True,stationary_states=states,
                sharp_launch_energy_projector_difference_audit=used_spectrum,
                guaranteed_time_averaged_terminal_score_difference=float(lower_exact),
                exact_score_lower_bound=str(lower_exact),sharp_release=observations)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    data=run()
    with parser.parse_args().output.open('x') as stream:
        json.dump(data,stream,separators=(',',':'),allow_nan=False)
        stream.write('\n')
