"""Separate floor/ceiling and graph route; imports no experiment/native code."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'RECEIPT.json').read_text())
records=json.loads((HERE/'RECORDS.json').read_text())['columns']
graph=json.loads((HERE/'inputs/graph.json').read_text())
counts=defaultdict(lambda: [[0,0],[0,0]])
work=defaultdict(lambda: defaultdict(F))
for phase,n,p,sign,edge,coin,j,y,meta in records:
    a,b=edge
    qa=sign if p==a else 0
    qb=sign if p==b else 0
    total=qa+qb
    low=total//2
    high=total-low
    oa,ob=(low,high) if coin==0 else (high,low)
    assert j==qa-oa and qb+j==ob
    out={p:sign} if p is not None else {}
    out[a]=oa
    out[b]=ob
    support=[k for k,v in out.items() if v]
    assert (support[0] if support else None)==y
    if phase=='hold_active':
        edgeid,h,t,selector=meta
        counts[(edgeid,h,t,sign)][int(p==t)][int(y==t)]+=1
    elif phase=='hold_idle':
        edgeid,h,t,s,position,selector=meta
        counts[(edgeid,h,t,s)][int(position==t)][int(position==t)]+=1
    elif phase=='work':
        edgeid,h,t=meta
        before=int(p==t)
        after=int(y==t)
        work[(edgeid,h,t,sign)][after-before]+=(F(2,3) if before==0 else F(1,3))/2
assert len(records)==r['direct_native_calls']==r['record_columns']
assert len(counts)==120 and all(x==[[3,1],[2,2]] for x in counts.values())
assert len(work)==120
for w in work.values():
    assert dict(w)=={-1:F(1,6),0:F(1,2),1:F(1,3)}
    assert sum(m*F(2)**(-v) for v,m in w.items())==1

L=np.zeros((240,240))
for a,b in graph['seams']:
    L[a,a]+=1
    L[b,b]+=1
    L[a,b]-=1
    L[b,a]-=1
ev=np.linalg.eigvalsh(L)
assert sum(ev<6-1e-8)==80
f=np.zeros(240)
f[[4,5,6,7]]=[.5,-.5,-.5,.5]
assert np.max(abs(L@f-6*f))==0
# Rebuild the complete one-event channel from the supplied record table.
groups=defaultdict(lambda: np.zeros(240))
for phase,n,p,sign,edge,coin,j,y,meta in records:
    if phase=='compact_channel':
        groups[(tuple(edge),coin,j)][y]+=f[p]
rho=sum(np.outer(v,v) for v in groups.values())/1260
assert abs(f@rho@f-float(F(r['baseline']['compact_fidelity_after_one_unread_event'])))<1e-13
assert abs(np.trace(rho)-1)<1e-13
supp=[4,5,6,7]
assert abs(1-sum(rho[i,i] for i in supp)-float(F(r['baseline']['compact_probability_outside_after_one_event'])))<1e-13
pair=511*math.log(1-12/1260)/math.log(1-6/1260)
assert abs(pair-r['baseline']['pair_keV'])<1e-9
assert pair>1022

for row in r['selection']['return_contexts']:
    # Return takes dist successful native half-probability moves.
    for d, observed in zip(row['distance'],row['timeout_after_32']):
        exact=F(sum(math.comb(32,k) for k in range(d)),2**32)
        assert F(observed)==exact
assert F(r['work']['reset_exponential_work_mean'])==F(171,256)
assert F(r['work']['reset_work_probabilities']['0'])+F(r['work']['reset_work_probabilities']['1'])==1
assert sum(F(p)*F(2)**(-int(w)) for w,p in r['work']['reset_work_probabilities'].items())==F(171,256)
K=np.array([[3,1],[2,2]])/4
P=np.array([[float(F(v)) for v in row] for row in r['separation']['independent_epoch_kernel']])
assert np.max(abs(P-np.kron(K,K)))==0
assert np.max(abs(np.array([2/3,1/3])@K-np.array([2/3,1/3])))<1e-15
assert r['separation']['joint_stationary_transfer_energy']==0
assert not r['physical_rest_energy_established'] and not any(r['physical_gates'].values())
print(json.dumps({'independent_record_check':True,'records':len(records),'selected_contexts':len(counts),
                  'physical_rest_energy_established':False},sort_keys=True))
