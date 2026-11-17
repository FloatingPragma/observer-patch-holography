#!/usr/bin/env python3
"""Separate floor/ceiling and graph audit. Imports no producer/source code."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np

HERE=Path(__file__).resolve().parent


def endpoint_pair(p,r,a,b,coin,law='native'):
    state={p:1,r:-1} if p>=0 else {}
    x,y=state.get(a,0),state.get(b,0)
    if law=='hold' or (law=='swap_hold' and coin==0):left,right=x,y
    elif law=='swap_hold':left,right=y,x
    else:
        lo=(x+y)//2; hi=x+y-lo
        left,right=(lo,hi) if coin==0 else (hi,lo)
    state[a],state[b]=left,right
    state={v:n for v,n in state.items() if n}
    if not state:return 0,x-left
    pp=next(v for v,n in state.items() if n==1)
    rr=next(v for v,n in state.items() if n==-1)
    assert len(state)==2
    return 256*pp+rr,x-left


def main():
    start=time.perf_counter()
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    graph=json.loads((HERE/'inputs/graph.json').read_text())
    inputs=json.loads((HERE/'INPUTS.json').read_text())
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in prereg['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    for name,digest in receipt['artifacts'].items():assert hashlib.sha256((HERE/'evidence'/name).read_bytes()).hexdigest()==digest,name
    edges=[tuple(e) for e in graph['seams']]; neighbors=[set() for _ in range(240)]
    for a,b in edges:neighbors[a].add(b);neighbors[b].add(a)
    f={int(k):v for k,v in inputs['f_numerator'].items()}
    g={int(k):v for k,v in inputs['g_numerator'].items()}
    def laplace(vector):
        return {v:len(neighbors[v])*vector.get(v,0)-sum(vector.get(u,0) for u in neighbors[v]) for v in range(240)}
    assert all(laplace(f)[v]==6*f.get(v,0) and laplace(g)[v]==6*g.get(v,0) for v in range(240))
    psi=defaultdict(int)
    for a,x in f.items():
        for b,y in g.items():psi[256*a+b]+=x*y;psi[256*b+a]-=x*y
    psi={z:n for z,n in psi.items() if n}
    assert len(psi)==24 and sum(n*n for n in psi.values())==24
    assert sorted([z,n] for z,n in psi.items())==receipt['state']['support']
    residual=defaultdict(int)
    for z,n in psi.items():
        a,b=divmod(z,256)
        residual[z]+=(len(neighbors[a])+len(neighbors[b])-12)*n
        for u in neighbors[a]:residual[256*u+b]-=n
        for u in neighbors[b]:residual[256*a+u]-=n
    assert not any(residual.values())

    data=np.load(HERE/'evidence/complete_actions.npz',allow_pickle=False)
    assert data['offsets'].shape==(1261,) and data['branch_counts'].shape==(1260,3)
    invalid=[z for z in range(65536) if z and not(z//256<240 and z%256<240 and z//256!=z%256)]
    assert np.array_equal(data['invalid'],invalid) and len(invalid)==8175
    expanded_action=hashlib.sha256()
    aggregates={law:Counter() for law in ('native','hold','swap_hold')}
    transfers={law:defaultdict(int) for law in aggregates}
    expected_records=[]; touched_checks=0
    total_origins=0
    for ei,(a,b) in enumerate(edges):
        touched={(p,r) for v in (a,b) for w in range(240) if v!=w for p,r in ((v,w),(w,v))}
        for c in (0,1):
            changed=[]
            for p,r in sorted(touched):
                out,j=endpoint_pair(p,r,a,b,c)
                if j:changed.append((256*p+r,out,j))
                else:assert out==256*p+r
                touched_checks+=1
            z=2*ei+c; lo,hi=map(int,data['offsets'][z:z+2])
            observed=list(zip(map(int,data['inputs'][lo:hi]),map(int,data['outputs'][lo:hi]),map(int,data['currents'][lo:hi])))
            assert observed==changed,(ei,c,'source action')
            cnt=Counter(j for _,_,j in changed);cnt[0]=57361-len(changed)
            assert list(data['branch_counts'][z])==[cnt[-1],cnt[0],cnt[1]]
            assert sum(cnt.values())==57361
            for current in (-1,1):
                dest=[out for _,out,j in changed if j==current]
                assert len(set(dest))==len(dest)
            out=np.arange(65536,dtype='<u2');j=np.zeros(65536,dtype=np.int8)
            for src,dst,cur in changed:out[src]=dst;j[src]=cur
            expanded_action.update(out.tobytes());expanded_action.update(j.tobytes())
            total_origins+=len(changed)
            for law in aggregates:
                branches=defaultdict(dict)
                for source,n in psi.items():
                    p,r=divmod(source,256);dest,cur=endpoint_pair(p,r,a,b,c,law)
                    assert dest not in branches[cur]
                    branches[cur][dest]=n;transfers[law][dest]+=n
                    if law=='native':expected_records.append([ei,c,cur,source,dest,n])
                for cur,branch in branches.items():
                    norm=sum(n*n for n in branch.values())
                    overlap=sum(n*psi.get(dest,0) for dest,n in branch.items())
                    aggregates[law]['probability']+=F(norm,30240)
                    aggregates[law]['vacuum']+=F(branch.get(0,0)**2,30240)
                    aggregates[law]['psi_click']+=F(overlap*overlap,725760)
    assert expanded_action.hexdigest()==receipt['source']['action_sha256']
    assert expanded_action.hexdigest()==json.loads((HERE/'inputs/archived_action_receipt.json').read_text())['complete_action_sha256']
    assert total_origins==receipt['source']['sparse_changed_columns']
    assert touched_checks==1260*954
    for law,values in aggregates.items():
        for name,value in values.items():assert str(value)==receipt['outcome_totals'][law][name],(law,name)
        assert values['probability']==1
    def clean(d):return {k:v for k,v in d.items() if v}
    assert clean(transfers['native'])=={z:1248*n for z,n in psi.items()}
    assert clean(transfers['hold'])=={z:1260*n for z,n in psi.items()}
    swap_same=clean(transfers['swap_hold'])=={z:1248*n for z,n in psi.items()}
    assert swap_same==receipt['controls']['swap_hold_retains_E12_eigenray']
    with gzip.open(HERE/'evidence/prepared_branches.jsonl.gz','rt') as fobj:observed=[json.loads(line) for line in fobj]
    assert sorted(expected_records)==sorted(observed)
    assert len(observed)==30240
    all_outcomes=json.loads((HERE/'evidence/complete_outcomes.json').read_text())
    assert len(all_outcomes)==3*1260*5
    for row in all_outcomes:
        assert F(row['probability'])==F(row['vacuum'])+F(row['psi_click'])+F(row['valid_complement'])+F(row['invalid'])
        assert all(F(row[k])>=0 for k in ('probability','vacuum','psi_click','valid_complement','invalid'))
    for law in aggregates:
        for name in ('probability','vacuum','psi_click'):
            value=sum((F(row[name]) for row in all_outcomes if row['law']==law),F())
            assert value==aggregates[law][name]
    contact_edges=[i for i,(a,b) in enumerate(edges) if psi.get(256*a+b,0)]
    assert contact_edges==receipt['contact']['support_edge_indices']
    assert aggregates['native']['vacuum']==F(len(contact_edges),630*12)
    assert aggregates['native']['vacuum']>0 and aggregates['hold']['vacuum']==aggregates['swap_hold']['vacuum']==0
    for lam,key in ((6,'E6_interval'),(12,'E12_interval')):
        x=F(lam,1260); s=F()
        for k in range(1,33):s+=1260*x**k/k
        bounds=(s,s+1260*x**33/(33*(1-x)))
        assert list(map(str,bounds))==receipt['energy'][key]
    assert receipt['scope']['native_compiled_pair_preparation'] is False
    assert receipt['scope']['radiation'] is False
    result={'all_pass':True,'independent_floor_ceiling_columns':touched_checks,'untouched_full_domain_proved_by_identity':True,
            'expanded_published_action_sha256':expanded_action.hexdigest(),'prepared_fine_records':len(observed),
            'native_vacuum_probability':str(aggregates['native']['vacuum']),
            'native_psi_fidelity':str(aggregates['native']['psi_click']),
            'swap_hold_keeps_spectral_ray':swap_same,'elapsed_seconds':time.perf_counter()-start,
            'review_group':'Codex same-group; no G6 or award'}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
