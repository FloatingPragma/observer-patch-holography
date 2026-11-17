#!/usr/bin/env python3
"""Second implementation, same author/group; imports no producer modules."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import time

sys.set_int_max_str_digits(0)  # All input lengths are bounded by this packet.
HERE = Path(__file__).resolve().parent


def q(a=0,b=0): return F(a),F(b)
def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def mul(x,y):
    a,b=x; c,d=y
    return a*c+b*d,a*d+b*c+b*d
def div(x,y):
    a,b=y; denominator=a*a+a*b-b*b
    return mul(x,((a+b)/denominator,-b/denominator))
def total(items):
    value=q()
    for item in items: value=add(value,item)
    return value
def dot(x,y): return total(mul(a,b) for a,b in zip(x,y))
def scale(v,a): return tuple(mul(x,a) for x in v)
def subtract(x,y): return tuple(add(a,neg(b)) for a,b in zip(x,y))


def own_event(sign,port,edge,coin,law):
    a,b=edge
    qa=sign*int(port==a); qb=sign*int(port==b)
    if law=='native':
        na=(qa+qb+coin)//2
        nb=qa+qb-na
    elif law=='swap': na,nb=qb,qa
    elif law=='hold': na,nb=qa,qb
    else: raise AssertionError(law)
    destination=port if port not in edge else (a if na else b)
    return destination,qa-na


def check():
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    raw=gzip.decompress((HERE/'RECORDS.json.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest()==receipt['records_sha256']
    records=json.loads(raw)
    points=[(q(-1),q(0,1),q()),(q(1),q(0,1),q()),
            (q(-1),q(0,-1),q()),(q(1),q(0,-1),q()),
            (q(),q(-1),q(0,1)),(q(),q(1),q(0,1)),
            (q(),q(-1),q(0,-1)),(q(),q(1),q(0,-1)),
            (q(0,1),q(),q(-1)),(q(0,1),q(),q(1)),
            (q(0,-1),q(),q(-1)),(q(0,-1),q(),q(1))]
    seams=[(a,b) for a in range(12) for b in range(a+1,12)
           if dot(subtract(points[a],points[b]),subtract(points[a],points[b]))==q(4)]
    assert seams==list(map(tuple,json.loads((HERE/'seams.json').read_text())))
    for sign,port,edge,coin,current,out in records['native_columns']:
        if sign:
            expected,j=own_event(sign,port,seams[edge],coin,'native')
            assert (out,current)==(expected,j)
        else: assert out==-1 and current==0
    assert len(records['native_columns'])==24*60+62
    for sign,origin,arrival,law,start,coin,j,out in records['suffix_transition_table']:
        expected,current=own_event(sign,start,tuple(sorted((origin,arrival))),coin,law)
        assert (out,j)==(expected,current)
    histograms={}
    for sign,origin,arrival,law,start,bits,currents,ports in records['suffix_histories']:
        edge=tuple(sorted((origin,arrival)))
        p=start
        assert len(bits)==len(currents)==len(ports)==8
        for t,(coin,j,out) in enumerate(zip(map(int,bits),currents,ports)):
            p,current=own_event(sign,p,edge,coin,law)
            assert (out,j)==(p,current)
            histograms.setdefault((sign,origin,arrival,law,start,t),Counter())[p]+=1
    assert len(records['suffix_histories'])==120*3*2*256
    assert len(histograms)==120*3*2*8
    for (sign,origin,arrival,law,start,t),h in histograms.items():
        assert sum(h.values())==256
        if law=='native': assert h=={origin:128,arrival:128}
        else:
            opposite=arrival if start==origin else origin
            expected=opposite if law=='swap' and (t+1)%2 else start
            assert h=={expected:256}
    for record in receipt['spatial']['rotations']:
        a,b=record['edge']
        m=scale(tuple(add(x,y) for x,y in zip(points[a],points[b])),q(F(1,2)))
        norm2=dot(m,m)
        assert norm2==q(1,1)
        rotated=[subtract(scale(m,div(mul(q(2),dot(m,x)),norm2)),x) for x in points]
        permutation=[points.index(x) for x in rotated]
        assert permutation==record['half_turn_port_permutation']
        assert permutation[a]==b and permutation[b]==a
        assert set(tuple(sorted((permutation[u],permutation[v]))) for u,v in seams)==set(seams)
        # Exact formula has eigenvalues +1 on m and -1 on its orthogonal plane.
        d=subtract(points[b],points[a])
        assert dot(m,d)==q() and dot(d,d)==q(4)
        assert record['native_field_half_turn_TV']=='0'
        assert record['polar_control_half_turn_TV']=='1'
        assert record['asymmetric_plant_half_turn_TV']=='1/2'
    discovery=receipt['discovery']
    qfail=F(1,16)**16; accepted=1-qfail; continued=accepted*F(11,12)
    total_discovery=F(0); total_selector=F(0)
    for (t,move),(u,fail) in zip(discovery['per_directed_first_move_by_time'],
                               discovery['selector_failure_by_time']):
        assert t==u and F(move)==continued**(t-1)*accepted/60
        assert F(fail)==continued**(t-1)*qfail
        total_discovery+=5*F(move); total_selector+=F(fail)
    assert total_discovery+total_selector+continued**256==1
    assert F(discovery['discovery_mass'])==total_discovery
    assert F(discovery['selector_failure_mass'])==total_selector
    assert F(discovery['search_timeout_mass'])==continued**256
    trajectory_events=0
    for record in records['controller_trajectories']:
        sign=record['sign']; port=record['start']; mode='search'; frame=None
        for event in record['records']:
            trajectory_events+=1
            if event[0]=='selector':
                assert event[-1]==0 and event[-2] in (0,1)
                continue
            phase,t,edge,coin,j,out=event
            expected,current=own_event(sign,port,seams[edge],coin,record['law'])
            assert (out,j)==(expected,current)
            if phase=='search' and j:
                assert frame is None
                frame=(port,out); mode='holding'
            if phase=='holding': assert mode=='holding' and port in frame and out in frame
            port=out
        assert (list(frame) if frame else None)==record['frame']
    return {'producer_imports':[], 'independence':'Separate implementation; same Codex author/group',
            'source_table_rows':len(records['native_columns']),
            'histories_recomputed':len(records['suffix_histories']),
            'histogram_checks':len(histograms),'exact_Qphi_rotations':30,
            'controller_events_replayed':trajectory_events,
            'receipt_sha256':hashlib.sha256((HERE/'RECEIPT.json').read_bytes()).hexdigest(),
            'records_sha256':receipt['records_sha256'],'passed':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--record',action='store_true')
    args=parser.parse_args(); started=time.perf_counter(); result=check()
    if args.record:
        with (HERE/'INDEPENDENT_CHECK.json').open('x')as f:
            f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else: assert json.loads((HERE/'INDEPENDENT_CHECK.json').read_text())==result
    print(json.dumps(dict(result,seconds=time.perf_counter()-started),sort_keys=True))
