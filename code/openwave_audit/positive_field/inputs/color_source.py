"""Exact finite three-quark SU(3) Gauss-chain source, compressed without cutoff.

Real algebraic entries are finite sums q sqrt(d), with rational q and squarefree d.
The integer label order is the increasing occupation triple (cell,spin,weak).
No producer eigensolver or ideal spectral preparation is imported.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from functools import lru_cache
from math import isqrt

def radical(n):
    n=int(n)
    if n<0:raise ValueError('nonnegative radicand required')
    if n==0:return 0,1
    a=1;d=n;p=2
    while p*p<=d:
        while d%(p*p)==0:a*=p;d//=p*p
        p+=1
    return a,d

def add(a,b,scale=F(1)):
    out=dict(a)
    for d,q in b.items():out[d]=out.get(d,F())+scale*q
    return {d:q for d,q in out.items() if q}

def mul(a,b):
    out={}
    for d,x in a.items():
        for e,y in b.items():
            n,r=radical(d*e);out=add(out,{r:n*x*y})
    return out

def root(n):
    a,d=radical(n)
    return {d:F(a)} if a else {}

def interval(a,bits):
    lo=hi=F()
    for d,q in a.items():
        n=isqrt(d<<(2*bits));l=F(n,1<<bits);h=l if n*n==d<<(2*bits) else F(n+1,1<<bits)
        if q<0:l,h=h,l
        lo+=q*l;hi+=q*h
    return lo,hi

def value(a):return sum(float(q)*d**.5 for d,q in a.items())

def act(n,i,j):
    if not n[j]:return None,{}
    out=list(n);out[j]-=1;a=root(n[j]*(out[i]+1));out[i]+=1
    return tuple(out),a

@lru_cache(None)
def data(L,kappa=F(1),weakup=1,U=F(1),J=F(1,4)):
    if L not in [2,3,4,5] or kappa not in [F(0),F(1,2),F(1),F(2)] or weakup not in [1,2] or U not in [F(0),F(1)] or J not in [F(0),F(1,4)] or kappa>2*U:
        raise ValueError('frozen finite chain family only')
    labels=[];counts=[]
    for lab in combinations_with_replacement(range(4*L),3):
        if sum(i%2==0 for i in lab)!=weakup or sum(i%4<2 for i in lab)!=2:continue
        labels.append(lab);counts.append(tuple(lab.count(i) for i in range(4*L)))
    ix={n:j for j,n in enumerate(counts)};n=len(labels);bits=(n-1).bit_length();D=1<<bits
    H={};patterns=[];spans=[];dip=[]
    def put(i,j,a,scale=F(1)):
        H[i,j]=add(H.get((i,j),{}),a,scale)
    for col,occ in enumerate(counts):
        nn=[sum(occ[4*x:4*x+4]) for x in range(L)]
        pref=[sum(nn[:x+1]) for x in range(L-1)]
        c=F(2,3)*sum(k*(3-k) for k in nn)
        electric=F(2,3)*sum(k*(3-k) for k in pref)
        span=labels[col][-1]//4-labels[col][0]//4
        assert electric==F(4,3)*span
        patterns.append(tuple(nn));spans.append(span)
        dip.append(sum(F(x,3)*(2*(occ[4*x]+occ[4*x+2])-occ[4*x+1]-occ[4*x+3]) for x in range(L)))
        put(col,col,{1:F(3)+(U-kappa/2)*c+kappa*electric})
        # Positive graph hopping, every spin/weak orbital.
        for x in range(L-1):
            for sw in range(4):
                a,b=4*x+sw,4*(x+1)+sw
                put(col,col,{1:F(occ[a]+occ[b],8)})
                for i,j in [(a,b),(b,a)]:
                    target,z=act(occ,i,j)
                    if target is not None:put(ix[target],col,z,-F(1,8))
        # S_x²=Sz²+(S+S-+S-S+)/2; intermediate states are not projected away.
        for x in range(L):
            z=F(occ[4*x]+occ[4*x+1]-occ[4*x+2]-occ[4*x+3],2)
            put(col,col,{1:z*z},J)
            plus=[(4*x+w,4*x+2+w) for w in range(2)]
            minus=[(b,a) for a,b in plus]
            for outer,inner in [(plus,minus),(minus,plus)]:
                for a,b in inner:
                    mid,c1=act(occ,a,b)
                    if mid is None:continue
                    for c,d in outer:
                        target,c2=act(mid,c,d)
                        if target is not None:put(ix[target],col,mul(c1,c2),J/2)
    H={ij:z for ij,z in H.items() if z}
    assert all(z==H[j,i] for (i,j),z in H.items())
    assert all(set(H.get((i,i),{}))<={1} for i in range(n))
    diagonal=[H.get((i,i),{}).get(1,F()) for i in range(n)]+[F(16)]*(D-n)
    pairs=[(i,j,z) for (i,j),z in sorted(H.items()) if i<j]
    return dict(L=L,kappa=kappa,weakup=weakup,U=U,J=J,labels=labels,counts=counts,dimension=n,bits=bits,code_dimension=D,
                entries=H,diagonal=diagonal,pairs=pairs,patterns=patterns,span=spans,dipole=dip+[F()]*(D-n))

def dense(d,field=F()):
    import numpy as np
    M=np.diag([float(x+field*y) for x,y in zip(d['diagonal'],d['dipole'])])
    for i,j,a in d['pairs']:M[i,j]=M[j,i]=value(a)
    return M

def rounding(d,bits,field=F()):
    diagonal=[x+field*y for x,y in zip(d['diagonal'],d['dipole'])]
    pairs=[];rows=[F()] * d['code_dimension'];normsum=max(map(abs,diagonal))
    for i,j,a in d['pairs']:
        lo,hi=interval(a,bits);mid=(lo+hi)/2;radius=(hi-lo)/2
        pairs.append((i,j,mid));rows[i]+=radius;rows[j]+=radius
        normsum+=max(abs(lo),abs(hi))
    return dict(diagonal=diagonal,pairs=pairs,operator_error=max(rows),normsum_upper=normsum,
                bits=bits,group_count=1+len(pairs))

def exact_gap():
    # eta=1/8 completes the square in the P/Q off-diagonal blocks.
    eta=F(1,8);a2=F(3,32);TD=F(281,48)
    assert TD-a2/eta>F(61,16)
    return dict(cluster_rank=2,cluster_lower=F(55,16),cluster_upper=F(57,16),
                next_level_lower=F(61,16),filter_cut=F(59,16),filter_energy_margin=F(1,8),
                total_source_bound=F(965,48),phase_period=F(32),cluster_colocation_lower=1-a2/(TD-F(57,16))**2,
                ground_trial=F(57,16),separated_floor=TD,coupling_squared=a2)
