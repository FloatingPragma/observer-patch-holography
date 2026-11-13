"""Complete physical Gauss-chain complement by finite irrep dynamic programming."""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
import sys
from pathlib import Path
from .complement import local_rows,character,casimir,DEN
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients


@lru_cache(None)
def fusion(a,b):
    weights=Counter()
    for x,u in character(*a).items():
        for y,v in character(*b).items():
            weights[x[0]+y[0],x[1]+y[1]]+=u*v
    out=[]
    for c in sorted(k for k in weights if min(k)>=0):
        m=sum(sign*weights[c[0]+x,c[1]+y] for x,y,sign in DEN)
        assert m>=0
        if m:out.append((c,m))
    dim=lambda z:(z[0]+1)*(z[1]+1)*(sum(z)+2)//2
    assert sum(dim(c)*m for c,m in out)==dim(a)*dim(b)
    return tuple(out)


def full_complement(cells,three_body=F(0),vacuum_compatible=True):
    p,_,unit=inputs();kappa=declared_kappa();target=3*cells
    dU,dm=degree_coefficients()
    # key=(total number, outgoing edge irrep, any nonquartet cell, last two N)
    states={(0,(0,0),False,()):(F(0),())}
    checked_fusions=set()
    for x in range(cells):
        out={}
        for (total,edge,excited,last),(value,path) in states.items():
            for row in local_rows():
                n,a,b,s,i,m=row
                if total+n>target:continue
                energy=p['A']*n*n/9+(p['U']-kappa/2)*casimir(a,b)
                energy+=p['Js']*F(s*(s+2),4)+p['Jw']*F(i*(i+2),4)
                if vacuum_compatible and 0<x<cells-1:
                    energy-=p['spin_gradient']*F(s*(s+2),4)+p['eta']*F(i*(i+2),4)+p['t']*n
                    energy+=dU*casimir(a,b)+dm*n
                if len(last)==2 and (*last,n)==(3,3,3):energy-=three_body
                checked_fusions.add((edge,(a,b)))
                for newedge,mult in fusion(edge,(a,b)):
                    if x==cells-1 and (newedge!=(0,0) or total+n!=target):continue
                    candidate=value+energy+(kappa*casimir(*newedge) if x<cells-1 else 0)
                    flag=excited or row[:5]!=(3,0,0,1,1)
                    key=(total+n,newedge,flag,(*last,n)[-2:])
                    if key not in out or candidate<out[key][0]:out[key]=(candidate,(*path,row))
        states=out
    candidates=[v for k,v in states.items() if k[2]]
    minimum,path=min(candidates,key=lambda z:z[0])
    # Fixed Iz is added separately; this floor is uniform before -delta Iz.
    floor=3*cells*p['mbar']+minimum
    return {'floor_relative_to_minus_delta_Iz':floor,
            'floor_GeV_relative_to_minus_delta_Iz':floor*unit,'witness':path,
            'cells':cells,'three_body':three_body,'vacuum_compatible':vacuum_compatible,
            'full_local_irrep_rows':len(local_rows()),'fusion_products_checked':len(checked_fusions)}

