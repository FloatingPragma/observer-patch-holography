from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
import numpy as np
from model.cube_eta_protection import cube,source
from model.certified_reference import interval as iv
I=iv.Interval;Z=I.point(0);O=I.point(1)
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def scale(a,s):return(a[0]*s,a[1]*s)
from functools import lru_cache

@lru_cache(None)
def targets():
    zero=(Z,Z);terms=tuple(t for group in cube.GEOMETRY.groups for t in group)
    masks=cube.weights(2).tolist();lookup={x:i for i,x in enumerate(masks)}
    parent=[[zero,zero]for _ in masks]
    for site,coord in enumerate(cube.GEOMETRY.sites):
     i=lookup[3<<(2*site)];sgn=(-1)**sum(coord)
     parent[i]=[(I.point(sgn),Z),(I.point(sgn-(6 if site==0 else 0)),Z)]
    T,_,v=cube.parent_lift()
    R=np.zeros(65536,dtype=np.int64)
    for mask,val in enumerate(source.reference_values()):R[sum(3<<(2*j)for j in range(8)if mask>>j&1)]=val
    vp=v.copy();vp[lookup[3]]-=6
    assert np.array_equal(T@v,6*R[cube.weights(6)])
    assert np.array_equal(T@vp,6*(R*np.array([(-1)**((m&3)==3)for m in range(65536)]))[cube.weights(6)])
    s,c=iv.sine_cosine(F(1,2048));ss,cc=iv.sine_cosine(-F(1,64));contact=(cc,ss)
    pairs={}
    for p,q,z in terms:
     rows=[]
     for mask in masks:
      if mask>>q&1 and not mask>>p&1:
       clear=mask^(1<<q);dest=clear|(1<<p)
       sign=(-1)**((mask&((1<<q)-1)).bit_count()+(clear&((1<<p)-1)).bit_count())
       w=z/abs(z)*sign
       rows.append((lookup[dest],lookup[mask],(s*int(w.imag),-s*int(w.real)),(-s*int(w.imag),-s*int(w.real))))
     pairs[(p,q,z)]=rows
    for step in range(64):
     for schedule in (terms,tuple(reversed(terms))):
      if schedule is not terms:
       for site in range(8):
        j=lookup[3<<(2*site)];parent[j]=[cmul(contact,x)for x in parent[j]]
      for term in schedule:
       for j,i,alpha,beta in pairs[term]:
        a,b=parent[j],parent[i]
        parent[j]=[add(scale(x,c),cmul(alpha,y))for x,y in zip(a,b)]
        parent[i]=[add(scale(y,c),cmul(beta,x))for x,y in zip(a,b)]
    ss,cc=iv.sine_cosine(F(-2));scalar=(cc,ss);normal=iv.sqrt_fraction(F(1,112))/6
    parent=[[scale(cmul(scalar,x),normal)for x in row]for row in parent]
    # Lift both arms into every N6 row, with the exact eta² integer map.
    T=T.tocsr();arms={}
    for i,mask in enumerate(cube.weights(6).tolist()):
     a,b=zero,zero
     for k in range(T.indptr[i],T.indptr[i+1]):
      j=int(T.indices[k]);n=int(T.data[k]);a=add(a,scale(parent[j][0],n));b=add(b,scale(parent[j][1],n))
     arms[mask]=(a,b)
    root=iv.sqrt_fraction(F(1,2));targets=[]
    for phase in (1,1j):
     target={}
     for mask,(a,b)in arms.items():
      w=phase*((-1j)**((mask&3).bit_count()));b=cmul((I.point(int(w.real)),I.point(int(w.imag))),b)
      target[mask]=scale(add(a,b),root)
      target[65536+mask]=scale(add(scale(a,-1),b),root)
     targets.append(target)
    return targets

def compare_arrays(arrays):
    results=[]
    for k,target in enumerate(targets()):
        state=arrays[f'analyzer_{k}_full_success_state'][:,0]
        if state.shape!=(131072,) or not np.isfinite(state).all():
            raise ValueError('complete finite full-state array required')
        error2=F(); p0lo=p0hi=F()
        for index,z in enumerate(state):
            a,b=target.get(index,(Z,Z))
            er=(I.point(F(float(z.real)))-a).absolute_upper()
            ei=(I.point(F(float(z.imag)))-b).absolute_upper()
            error2+=er*er+ei*ei
            if index<65536:
                for cell in (a,b):
                    lo,hi=F(cell.lo,iv.SCALE),F(cell.hi,iv.SCALE)
                    p0lo+=0 if lo<=0<=hi else min(lo*lo,hi*hi)
                    p0hi+=max(lo*lo,hi*hi)
        results.append({'setting':k,'complete_state_error_upper':str(F(iv.sqrt_fraction(error2).hi,iv.SCALE)),
                        'exact_template_pointer_zero_interval':[str(p0lo),str(p0hi)]})
    return results
