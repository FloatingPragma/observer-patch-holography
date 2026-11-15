#!/usr/bin/env python3
"""Independent first-hit tape, scheduler and rational population review."""
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import itertools
import json
import math
import random

HERE=Path(__file__).resolve().parent
PACKET=HERE.parents[1]/"work/2026-11-14_codex_first_hit_clock"


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def power(a,n):
    out=[[F(i==j) for j in range(len(a))] for i in range(len(a))]
    while n:
        if n%2:out=matmul(out,a)
        a=matmul(a,a);n//=2
    return out


def solve(a,b):
    rows=[list(row)+[v] for row,v in zip(a,b)];n=len(rows)
    for i in range(n):
        k=next(k for k in range(i,n) if rows[k][i]);rows[i],rows[k]=rows[k],rows[i]
        divisor=rows[i][i];rows[i]=[v/divisor for v in rows[i]]
        for k in range(n):
            if k!=i:
                scale=rows[k][i];rows[k]=[v-scale*w for v,w in zip(rows[k],rows[i])]
    return [row[-1] for row in rows]


def main():
    receipt=json.loads((PACKET/"RECEIPT.json").read_text())
    prereg=json.loads((PACKET/"PREREG.json").read_text())
    assert all(hashlib.sha256((PACKET/f).read_bytes()).hexdigest()==h for f,h in prereg['files'].items())
    raw=gzip.decompress((PACKET/"RECORDS.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest()==receipt['raw_uncompressed_sha256']
    records=json.loads(raw);cells=records['graph']['cells'];edges=records['graph']['edges']
    assert cells==[list(x) for x in itertools.product(range(3),range(3),range(3),range(-2,2))]
    assert len(edges)==8100 and len({tuple(x) for x in edges})==8100
    adj=[[] for _ in range(1296)];vert=[[] for _ in range(1296)]
    for i,(a,b,axis) in enumerate(edges):
        assert 0<=a<1296 and 0<=b<1296 and a!=b and axis in (-1,0,1,2,3)
        adj[a].append(i);adj[b].append(i)
        if axis==3:vert[a].append(i);vert[b].append(i)
    assert all(len(adj[i])==(12 if cells[i//12][3] in (-2,1) else 13) for i in range(1296))
    assert len(records['branches'])==48600
    for index,row in enumerate(records['branches']):
        e,slot,c,j=row
        assert (e,slot,c)==(index//6,(index%6)//2,index%2)
        a,b=int(slot==0),int(slot==1)
        assert j==a-(a+b+c)//2
    rng=random.Random(2026111430);summaries=[];total=48600
    momenta=((1,0,0),(1,1,0),(1,1,1))
    for law in ('incident','vertical_only_control'):
        runs=records['runs'][law];assert len(runs)==4096
        sums=[F(0)]*3;events=0;censored=0;lengths=[]
        for trial,run in enumerate(runs):
            loc=run['start'];assert cells[loc//12]==[0,0,0,0] and loc%12==trial%12
            assert 0<len(run['tape'])<=4096
            for e,c,j in run['tape']:
                assert cells[loc//12][3]<1
                eligible=adj[loc] if law=='incident' else vert[loc]
                assert e==eligible[rng.randrange(len(eligible))] and c==rng.randrange(2)
                a,b,_=edges[e];xa,xb=int(loc==a),int(loc==b)
                assert xa+xb==1 and j==xa-(xa+xb+c)//2
                if j==1:loc=b
                elif j==-1:loc=a
                else:assert j==0
            complete=cells[loc//12][3]==1
            assert complete==run['completed']
            if not complete:assert len(run['tape'])==4096
            censored+=not complete;events+=len(run['tape']);lengths.append(len(run['tape']))
            for k,p in enumerate(momenta):
                phase=sum(a*b for a,b in zip(p,cells[loc//12][:3]))%3
                sums[k]+=(F(1) if phase==0 else F(-1,2)) if complete else 0
        means=[x/4096 for x in sums];r=receipt['summaries'][law]
        assert all(float(x)==y for x,y in zip(means,r['sample_means']))
        assert (events,censored,max(lengths))==(r['actual_attempts'],r['censored'],r['max_attempts'])
        total+=events
        summaries.append({'law':law,'actual_events_verified':events,'all_scheduler_choices_reproduced':True,
                          'censored':censored,'exact_record_means':[str(x) for x in means]})
    assert total==receipt['actual_native_calls']==405755
    rational=[];identity=[[F(i==j) for j in range(3)] for i in range(3)]
    for k in range(4):
        a=3*k
        q=[[F(23-a,24),F(1,24),F(0)],
           [F(1,26),F(24-a,26),F(1,26)],
           [F(0),F(1,26),F(24-a,26)]]
        b=[F(0),F(0),F(1,26)]
        lhs=[[x-y for x,y in zip(r,s)] for r,s in zip(identity,q)]
        h=solve(lhs,b)
        if k==0:
            assert h==[F(1)]*3
            costs=solve(lhs,[F(1)]*3);assert costs[2]==76
            tail=sum(power(q,4096)[2]);assert F(0)<tail<F(4,10**15)
            tail_float=float(tail)
        else:
            r=receipt['population'][k-1]
            assert abs(float(h[2])-r['hit_transform'])<1e-14
            assert all(abs(float(v)-w)<1e-14 for row,row2 in zip(q,r['kernel']) for v,w in zip(row,row2))
            omega=2*math.asinh(math.sqrt(a)/2)
            assert abs(float(h[2])-math.cosh(2.5*omega)/math.cosh(3.5*omega))<1e-14
            rational.append({'spatial_symbol':a,'exact_uncapped_transfer':str(h[2]),
                             'numeric_match_absolute_error':abs(float(h[2])-r['hit_transform'])})
    eps=1.5*math.sqrt(math.log(1200)/(2*4096))
    assert all(abs(x['simultaneous_hoeffding_halfwidth']-eps)<1e-15 for x in receipt['summaries'].values())
    out={'review_type':'same-group independent full record and rational population check',
         'producer_or_native_imported':False,'prereg_hashes_match':True,'actual_native_calls_verified':total,
         'cases':summaries,'population':rational,'exact_mean_attempts':76,
         'exact_rational_tail_less_than':'4/10^15','tail_float_diagnostic':tail_float,
         'simultaneous_six_mean_hoeffding_coverage':0.99,
         'theorem_scope':'Repeated identical first-hit transfer requires infinite-depth vertical translation symmetry; fixed finite bottom gives changing effective depths.'}
    (HERE/'RECEIPT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out))


if __name__=='__main__':main()
