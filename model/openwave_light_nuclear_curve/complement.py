"""Exact local exterior characters and colour-singlet multiplicities."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb
import sys
from pathlib import Path
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from model.openwave_common_hyperfine_spectrum.gauss import character, casimir

ROOTS = ((2, -1), (-1, 2), (1, 1))
DEN = tuple((sum(x for x,y in z), sum(y for x,y in z), (-1)**n)
            for n in range(4) for z in combinations(ROOTS,n))


@lru_cache(None)
def local_rows():
    weights = [(a,b,s,i) for a,b in ((1,0),(-1,1),(0,-1))
               for s,i in product((-1,1),repeat=2)]
    char = Counter({(0,0,0,0,0):1})
    for w in weights:
        out=char.copy()
        for (n,a,b,s,i),m in char.items():
            out[n+1,a+w[0],b+w[1],s+w[2],i+w[3]] += m
        char=out
    rows=[]
    for n in range(13):
        dim=0
        for a,b,s,i in product(range(5),range(5),range(7),range(7)):
            mult=sum(sign*(-1)**(u+v)*char[n,a+x,b+y,s+2*u,i+2*v]
                     for x,y,sign in DEN for u,v in product((0,1),repeat=2))
            assert mult>=0
            if not mult:continue
            d=(a+1)*(b+1)*(a+b+2)//2*(s+1)*(i+1)
            rows.append((n,a,b,s,i,mult))
            dim+=d*mult
        assert dim==comb(12,n),(n,dim)
    assert sum((a+1)*(b+1)*(a+b+2)//2*(s+1)*(i+1)*m
               for n,a,b,s,i,m in rows)==4096
    return tuple(rows)


@lru_cache(None)
def singlet_multiplicity(a,b,c):
    tensor=Counter()
    for x,u in character(*a).items():
        for y,v in character(*b).items():
            tensor[x[0]+y[0],x[1]+y[1]]+=u*v
    # a tensor b must contain the conjugate of c.
    m=sum(sign*tensor[c[1]+x,c[0]+y] for x,y,sign in DEN)
    assert m>=0
    return m



def three_cell_singlet_dimension():
    byn={n:[r for r in local_rows() if r[0]==n] for n in range(13)}
    dimension=0
    for n0 in range(10):
        for n1 in range(10-n0):
            for triple in product(byn[n0],byn[n1],byn[9-n0-n1]):
                fusion=singlet_multiplicity(*[r[1:3] for r in triple])
                dimension+=fusion*__import__('math').prod((r[3]+1)*(r[4]+1)*r[5] for r in triple)
    return dimension
