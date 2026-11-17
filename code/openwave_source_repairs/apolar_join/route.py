#!/usr/bin/env python3
"""Continue sealed selected units through native routing, without source reset."""
from collections import deque
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import random
import time
from native_event import native

HERE=Path(__file__).resolve().parent
CAP=32


def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n'
def sha(data):return hashlib.sha256(data).hexdigest()


def freeze():
    r=json.loads((HERE/'PREREG_ROUTE.json').read_text())
    for name,digest in r['files'].items():assert sha((HERE/name).read_bytes())==digest,name


def update(x,row):
    a,b,c,j=row;y=list(x);total=x[a]+x[b]
    y[a]=(total+c)//2;y[b]=total-y[a]
    assert j==x[a]-y[a] and c in (0,1)
    return y


def verify_checkpoint(trial):
    x=trial['initial_audit_only'][:]
    assert len(x)==240 and sum(x)==1 and not any(x[12:])
    audited=0
    for block in trial['selection']:
        for row in block['reference']+[block['body']]:
            assert row is not None
            x=update(x,row);audited+=1
    for row in trial['census']:x=update(x,row);audited+=1
    handoff=trial['basis_handoff']
    assert x==handoff['observed_source']==handoff['actual_source_audit_only']
    prefix={'selection':trial['selection'],'census':trial['census']}
    assert sha(canonical(prefix).encode())==handoff['prefix_sha256']
    assert handoff['occupied']==x.index(1) and handoff['sign']==1
    assert sum(abs(v) for v in x)==1
    return tuple(x),handoff['observed_source'][:],audited


def shortest(start,target,edges):
    queue=deque([[start]]);seen={start}
    while queue:
        path=queue.popleft()
        if path[-1]==target:return path
        neighbours=sorted(b if a==path[-1] else a for a,b in edges if path[-1] in (a,b))
        for nxt in neighbours:
            if nxt not in seen:seen.add(nxt);queue.append(path+[nxt])
    raise AssertionError('No native route')


def locate(observed):
    nz=[i for i,v in enumerate(observed) if v]
    assert len(nz)==1 and observed[nz[0]]==1 and nz[0]<12
    return nz[0]


def route_source(actual,observed,plan,seed,force_hold=False):
    rng=random.Random(seed);legs=[];calls=0
    for leg,(departure,target) in enumerate(zip(plan,plan[1:])):
        assert locate(observed)==departure
        edge=tuple(sorted((departure,target)));records=[]
        moved=False
        for attempt in range(CAP):
            coin=(1 if departure==edge[0] else 0) if force_hold else rng.randrange(2)
            actual,(_,j)=native(actual,edge,coin);calls+=1
            a,b=edge;observed[a]-=j;observed[b]+=j
            assert tuple(observed)==actual
            records.append([a,b,coin,j])
            if j:
                assert locate(observed)==target;moved=True;break
        legs.append({'leg':leg,'departure':departure,'target':target,
                     'records':records,'success':moved})
        if not moved:return actual,observed,legs,calls,'route_timeout'
    return actual,observed,legs,calls,'routed'


def execute():
    freeze()
    graph=json.loads((HERE/'inputs/graph.json').read_text())
    edges=[tuple(e) for e in graph['seams'] if max(e)<12]
    assert len(edges)==30 and (0,1) in edges
    raw=(HERE/'inputs/bridge_RECORDS.json.gz').read_bytes()
    origin=json.loads((HERE/'inputs/bridge_RECEIPT.json').read_text())
    assert sha(raw)==origin['evidence_sha256']
    data=json.loads(gzip.decompress(raw))
    trials=[];native_calls=0;reconstructed=0
    h=1-F(1,2**CAP)
    for r in data['trials']:
        if r['selection_law']!='native':continue
        actual,observed,count=verify_checkpoint(r);reconstructed+=count
        initial=list(actual);start=locate(observed)
        to_zero=shortest(start,0,edges);plan=to_zero+[1,0]
        assert len(plan)-1<=5
        source,observed,legs,calls,status=route_source(actual,observed,plan,832711+r['index'])
        native_calls+=calls
        assert status=='routed' and locate(observed)==0
        frame_leg=legs[len(to_zero)-1]
        assert (frame_leg['departure'],frame_leg['target'])==(0,1)
        assert frame_leg['records'][-1][3]==1
        probability=h**(len(plan)-1)
        failures=[h**k*(1-h) for k in range(len(plan)-1)]
        assert sum(failures,F())+probability==1
        prefix={'bridge_prefix':r['basis_handoff']['prefix_sha256'],'route':legs}
        trials.append({'bridge_trial':r['index'],'input':initial,'plan':plan,
                       'seed':832711+r['index'],'legs':legs,'status':status,
                       'actual_final_source_audit_only':list(source),'observed_final_source':observed,
                       'frame':[0,1],'frame_native_event':frame_leg['records'][-1],
                       'prefix_sha256':sha(canonical(prefix).encode()),
                       'success_probability':str(probability),
                       'first_failed_leg_probabilities':list(map(str,failures))})
    assert len(trials)==36
    forced=[1]+[0]*239
    actual,obs,legs,calls,status=route_source(tuple(forced),forced[:],[0,1],0,True)
    native_calls+=calls
    assert status=='route_timeout' and len(legs[0]['records'])==CAP and list(actual)==forced
    result={'candidate':'A-old240-compact-energy-apolar-join-v1','stage':'actual_native_route',
            'successful_checkpoint_continuations':len(trials),'fresh_native_calls':native_calls,
            'predecessor_events_independently_reconstructed':reconstructed,
            'predecessor_reconstruction_is_not_fresh_native_execution':True,
            'maximum_move_count':max(len(t['plan'])-1 for t in trials),
            'route_failure_upper':str(1-h**5),'checkpoint_conditional_scope':True,
            'compact_state_prepared':False,'all_source_aborts_retained':True,
            'forced_timeout':{'status':status,'legs':legs,'final_source':list(actual)},
            'trials':trials}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    args=parser.parse_args();start=time.perf_counter();result=execute()
    target=HERE/'ROUTE_RECORDS.json'
    if args.record:
        with target.open('x')as f:f.write(canonical(result))
    else:assert json.loads(target.read_text())==result
    print(json.dumps({k:v for k,v in result.items() if k not in ('trials','forced_timeout')}
                     |{'seconds':time.perf_counter()-start},sort_keys=True))
