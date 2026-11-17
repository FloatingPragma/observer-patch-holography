"""Separate record reader for the scheduler-tunability counterexample."""
from pathlib import Path
from collections import defaultdict,Counter
from fractions import Fraction as F
import json
HERE=Path(__file__).resolve().parent
records=json.loads((HERE/'TUNABILITY_RECORDS.json').read_text())
r=json.loads((HERE/'TUNABILITY_RECEIPT.json').read_text())
kernels=defaultdict(lambda:[[F(0),F(0)],[F(0),F(0)]])
work=defaultdict(lambda:defaultdict(F))
for phase,h,t,sign,p,bits,c,j,y in records:
    assert type(c)is int and c in (0,1)
    key=(h,t,sign)
    if phase in ('gate','idle'):
        assert j==0
        if phase=='idle':
            assert p==h and bits[:2]!=[0,0] and y==p
            kernels[key][0][0]+=F(1,8)
    else:
        a,b=sorted((h,t))
        qa=sign if p==a else 0
        qb=sign if p==b else 0
        low=(qa+qb)//2;high=qa+qb-low
        out_a,out_b=(low,high) if c==0 else (high,low)
        assert j==qa-out_a and y==(a if out_a else b)
        if phase=='active':
            assert p==t or bits[:2]==[0,0]
            kernels[key][int(p==t)][int(y==t)]+=F(1,8)
        elif phase=='work':work[key][int(y==t)-int(p==t)]+=(F(4,5) if p==h else F(1,5))/2
        else:raise AssertionError(phase)
assert Counter(x[0] for x in records)=={'gate':3840,'active':1200,'idle':720,'work':480}
assert r['native_calls']==len(records)==6240
assert len(kernels)==len(work)==r['contexts']==120
assert all(k==[[F(7,8),F(1,8)],[F(1,2),F(1,2)]] for k in kernels.values())
assert all(dict(w)=={-1:F(1,10),0:F(1,2),1:F(2,5)} for w in work.values())
assert r['stationary_reference']==['4/5','1/5'] and r['modular_gap']=='log(4)'
assert not r['physical_rest_energy_established']
print(json.dumps({'independent_tunability_check':True,'native_calls':len(records),'contexts':len(kernels)}))
