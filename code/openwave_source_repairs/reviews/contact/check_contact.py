#!/usr/bin/env python3
"""Reviewer-built exact sparse branch check; no producer code imports."""
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import json

PACKET=Path(__file__).resolve().parents[2]/'work/2026-11-15_codex_compact_contact_repair'
edges=list(map(tuple,json.loads((PACKET/'inputs/graph.json').read_text())['seams']))
f={4:1,5:-1,6:-1,7:1}; g={5:-1,6:-1,8:1,11:1}
old_g={1:-1,2:-1,9:1,10:1}


def laplacian(v):
    out=defaultdict(int)
    for a,b in edges:
        delta=v.get(a,0)-v.get(b,0);out[a]+=delta;out[b]-=delta
    return out


def residual(v):
    w=laplacian(v)
    return max(abs(w.get(p,0)-6*v.get(p,0)) for p in range(240))


def move(p,r,a,b,coin,law):
    state={p:1,r:-1}; x,y=state.get(a,0),state.get(b,0)
    if law=='hold' or (law=='swap_hold' and coin==0): u,v=x,y
    elif law=='swap_hold':u,v=y,x
    else:u=(x+y+coin)//2;v=x+y-u
    state[a],state[b]=u,v;state={k:n for k,n in state.items() if n}
    if not state:return None,x-u
    return (next(k for k,n in state.items() if n==1),next(k for k,n in state.items() if n==-1)),x-u


assert residual(f)==residual(g)==0 and residual(old_g)>0
psi=defaultdict(int)
for p,x in f.items():
    for r,y in g.items():psi[p,r]+=x*y;psi[r,p]-=x*y
psi={k:v for k,v in psi.items() if v}
assert len(psi)==sum(v*v for v in psi.values())==24
results={}
for law in ('native','hold','swap_hold'):
    vacuum=Q();fidelity=Q();probability=Q();transfer=defaultdict(int)
    for a,b in edges:
        for coin in (0,1):
            branches=defaultdict(lambda:defaultdict(int))
            for (p,r),amp in psi.items():
                out,j=move(p,r,a,b,coin,law)
                assert out not in branches[j]
                branches[j][out]+=amp;transfer[out]+=amp
            for branch in branches.values():
                probability+=Q(sum(v*v for v in branch.values()),24*1260)
                vacuum+=Q(branch.get(None,0)**2,24*1260)
                fidelity+=Q(sum(psi.get(k,0)*v for k,v in branch.items())**2,24**2*1260)
    assert probability==1
    expected=1260 if law=='hold' else 1248
    assert {k:v for k,v in transfer.items() if v}=={k:expected*v for k,v in psi.items()}
    results[law]={'vacuum':str(vacuum),'fidelity':str(fidelity),'probability':str(probability)}
assert results['native']['vacuum']=='1/1260' and results['native']['fidelity']=='932/945'
assert results['hold']['vacuum']==results['swap_hold']['vacuum']=='0'
receipt=json.loads((PACKET/'RECEIPT.json').read_text())
for law,values in results.items():
    for key,value in values.items():
        assert receipt['outcome_totals'][law]['psi_click' if key=='fidelity' else key]==value
print(json.dumps({'passed':True,'producer_imports':[],'corrected_full_graph_residual':0,
                  'old_preparation_residual':residual(old_g),'results':results},sort_keys=True))
