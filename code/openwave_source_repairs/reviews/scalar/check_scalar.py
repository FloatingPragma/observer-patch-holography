#!/usr/bin/env python3
"""Reviewer-built scalar/CP check from graph construction, no producer imports."""
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import json

PACKET=Path(__file__).resolve().parents[2]/'work/2026-11-15_codex_symmetric_source_audit'
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def dot(x,y):
    z=(0,0)
    for a,b in zip(x,y):z=add(z,mul(a,b))
    return z
z=(0,0);p=(0,1);n=(0,-1);o=(1,0);m=(-1,0)
points=[(m,p,z),(o,p,z),(m,n,z),(o,n,z),(z,m,p),(z,o,p),
        (z,m,n),(z,o,n),(p,z,m),(p,z,o),(n,z,m),(n,z,o)]
cell=[(a,b) for a in range(12) for b in range(a+1,12) if dot(points[a],points[b])==p]
edges=[]
for module in range(2):
    off=37*module
    for shell in range(3):edges.extend((off+12*shell+a,off+12*shell+b) for a,b in cell)
    for k in range(12):edges.extend(((off+k,off+24+k),(off+12+k,off+24+k),(off+24+k,off+36)))
edges.append((36,73));edges.sort()
assert len(edges)==253
f={k:1 for k in range(12)}|{k:-1 for k in range(12,24)}
g={k+37:v for k,v in f.items()}
def laplacian(v):
    out=defaultdict(int)
    for a,b in edges:
        difference=v.get(a,0)-v.get(b,0);out[a]+=difference;out[b]-=difference
    return {k:v for k,v in out.items() if v}
assert laplacian(f)==f and laplacian(g)==g
fidelity=Q();escape=Q()
for a,b in edges:
    for coin in (0,1):
        branches=defaultdict(lambda:defaultdict(int))
        for source,amp in f.items():
            x=int(source==a);y=int(source==b);u=(x+y+coin)//2;v=x+y-u
            dest=source if source not in (a,b) else (a if u else b)
            current=x-u
            assert dest not in branches[current]
            branches[current][dest]=amp
        for branch in branches.values():
            fidelity+=Q(sum(f.get(k,0)*v for k,v in branch.items())**2,506*24**2)
            escape+=Q(sum(v*v for k,v in branch.items() if k not in f),506*24)
assert fidelity==Q(11867,12144) and escape==Q(1,506)
lower={k:(1 if k<37 else -1) for k in range(74)}
rayleigh=Q(sum((lower[a]-lower[b])**2 for a,b in edges),sum(v*v for v in lower.values()))
assert rayleigh==Q(2,37)
# With x=1/506, (1-x)^2 > 1-2x proves the shared-event log pair excess.
assert (1-Q(1,506))**2>1-Q(2,506)
r=json.loads((PACKET/'RECEIPT.json').read_text())
assert r['unread_ray_fidelity']==str(fidelity) and r['unread_escape_probability']==str(escape)
assert r['full_row_solved'] is False and r['physical_spin_generated'] is False
print(json.dumps({'passed':True,'producer_imports':[],'native_unread_fidelity':str(fidelity),
                  'native_escape':str(escape),'extended_Rayleigh':str(rayleigh),
                  'compact_eigenvalue':1,'graph_edges':len(edges)},sort_keys=True))
