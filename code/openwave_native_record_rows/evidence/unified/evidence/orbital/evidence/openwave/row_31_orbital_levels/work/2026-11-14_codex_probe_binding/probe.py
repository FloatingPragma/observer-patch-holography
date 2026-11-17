#!/usr/bin/env python3
"""Native-probe feedback calibration; no supplied potential table or gamma."""
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import random
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
TRIALS=4096
PROBE_CAP=256
ROUTE_CAP=128
SEED=2026111431
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
       (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
       (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))


def encoded(x):
    return (json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def solve(matrix,rhs):
    a=[list(row)+[v] for row,v in zip(matrix,rhs)]
    for j in range(len(a)):
        pivot=next(k for k in range(j,len(a)) if a[k][j])
        a[j],a[pivot]=a[pivot],a[j]
        v=a[j][j];a[j]=[x/v for x in a[j]]
        for k in range(len(a)):
            if k!=j:
                v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[j])]
    return [row[-1] for row in a]


def main():
    path=HERE/'inputs/native.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()== \
        'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('probe_binding_native',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module);native=module.native
    cells=list(itertools.product(range(-2,3),repeat=3));ci={x:i for i,x in enumerate(cells)}
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})
    edges=[]
    for x in cells:
        c=ci[x];edges.extend((12*c+a,12*c+b) for a,b in internal)
        for axis in range(3):
            if x[axis]<2:
                xx=list(x);xx[axis]+=1
                edges.extend((12*c+p,12*ci[tuple(xx)]+p) for p in range(12))
    assert len(edges)==7350
    adj=[[] for _ in range(1500)];lookup={}
    for eid,(a,b) in enumerate(edges):
        adj[a].append(eid);adj[b].append(eid);lookup[frozenset((a,b))]=eid
    branches=[];current={}
    for eid in range(len(edges)):
        for slot in range(3):
            loads=[0,0,0];loads[slot]=1
            for coin in (0,1):
                out,rec=native(loads,(0,1),coin)
                rebuilt=list(loads);rebuilt[0]-=rec[1];rebuilt[1]+=rec[1]
                assert tuple(rebuilt)==out
                branches.append([eid,slot,coin,rec[1]])
                current[eid,slot,coin]=rec[1]
    calls=len(branches)

    def destination(pos,eid,coin):
        a,b=edges[eid];j=current[eid,0 if pos==a else 1,coin]
        return b if j==1 else a if j==-1 else pos

    q=[[F(0) for _ in range(3)] for _ in range(3)];hit=[F(0)]*3;escape=[F(0)]*3
    orbit_rows={}
    for x in itertools.product(range(-1,2),repeat=3):
        if x==(0,0,0):continue
        orbit=sum(v!=0 for v in x)-1
        for port in range(12):
            pos=12*ci[x]+port;assert len(adj[pos])==11
            row=[F(0)]*3;bh=F(0);be=F(0)
            for eid in adj[pos]:
                for coin in (0,1):
                    y=cells[destination(pos,eid,coin)//12]
                    if y==(0,0,0):bh+=F(1,22)
                    elif max(map(abs,y))==2:be+=F(1,22)
                    else:row[sum(v!=0 for v in y)-1]+=F(1,22)
            old=orbit_rows.setdefault(orbit,(row,bh,be))
            assert old==(row,bh,be)
    for i in range(3):q[i],hit[i],escape[i]=orbit_rows[i]
    imq=[[F(i==j)-q[i][j] for j in range(3)] for i in range(3)]
    h=solve(imq,hit);mean=solve(imq,[F(1)]*3)
    assert h==[F(5,22),F(1,11),F(1,22)]
    assert mean==[F(31,3),F(10),F(26,3)]
    qf=np.asarray(q,dtype=float);power=np.linalg.matrix_power(qf,PROBE_CAP)
    capped_hit=np.linalg.solve(np.eye(3)-qf,(np.eye(3)-power)@np.asarray(hit,dtype=float))
    capped_escape=np.linalg.solve(np.eye(3)-qf,(np.eye(3)-power)@np.asarray(escape,dtype=float))
    tail=power@np.ones(3);tail_bound=(10/11)**PROBE_CAP
    assert max(tail)<=tail_bound and tail_bound<3e-11

    rng=random.Random(SEED);source={'i':12*ci[(0,0,0)],'o':12*ci[(0,0,0)]}
    visible=source.copy();all_events=[]

    def step(which,eid,events):
        nonlocal calls
        coin=rng.randrange(2)
        if which=='z':
            out,rec=native((0,0),(0,1),coin);assert out==(0,0) and rec[1]==0
        else:
            a,b=edges[eid];pos=source[which]
            out,rec=native((int(pos==a),int(pos==b)),(0,1),coin)
            if rec[1]:
                assert visible[which]==source[which]==(a if rec[1]==1 else b)
                visible[which]=source[which]=b if rec[1]==1 else a
            assert out==(int(source[which]==a),int(source[which]==b))
        calls+=1;record=[which,eid,rec[0],rec[1]];events.append(record);all_events.append(record)
        return rec[1]

    def route(which,target,events,port_zero=False):
        if port_zero and visible[which]%12:
            assert visible[which]%12==1
            cell=visible[which]//12;eid=lookup[frozenset((12*cell,12*cell+1))]
            for _ in range(ROUTE_CAP):
                step(which,eid,events)
                if visible[which]%12==0:break
            else:return False
        for axis in range(3):
            while cells[visible[which]//12][axis]!=target[axis]:
                pos=visible[which];xx=list(cells[pos//12])
                xx[axis]+=1 if xx[axis]<target[axis] else -1
                dest=12*ci[tuple(xx)]+pos%12;eid=lookup[frozenset((pos,dest))]
                for _ in range(ROUTE_CAP):
                    step(which,eid,events)
                    if visible[which]==dest:break
                else:return False
        return True

    cases=[((0,0,0),True),((1,0,0),True),((1,1,0),True),((1,1,1),True),((1,0,0),False)]
    trial_records=[];summaries=[];epsilon=math.sqrt(math.log(4000)/(2*TRIALS))
    preparation_bias_bound=25/2**ROUTE_CAP
    for x,feedback in cases:
        observed_hit=[];observed_move=[];censors=0;case_events=0
        for _ in range(TRIALS):
            records=[];ok=route('o',x,records,True) and route('i',x,records)
            prep_count=len(records)
            outcome='preparation_timeout';probe_steps=0
            if ok:
                for k in range(PROBE_CAP+1):
                    y=cells[visible['i']//12]
                    if y==(0,0,0):outcome='hit';break
                    if max(map(abs,y))==2:outcome='escape';break
                    if k==PROBE_CAP:outcome='probe_timeout';break
                    eligible=adj[visible['i']]
                    step('i',eligible[rng.randrange(len(eligible))],records);probe_steps+=1
            censored=outcome.endswith('timeout');censors+=censored
            held=(not ok) or (feedback and outcome!='escape')
            if held:move=step('z',-1,records)
            else:
                base=12*ci[x];eid=lookup[frozenset((base,base+1))]
                assert visible['o']==base
                move=step('o',eid,records);assert move in (0,1)
            observed_hit.append(int(outcome=='hit'));observed_move.append(move)
            case_events+=len(records)
            trial_records.append(dict(x=x,feedback=feedback,preparation_events=prep_count,
                probe_events=probe_steps,outcome=outcome,held=held,records=records))
        orbit=sum(v!=0 for v in x)-1
        hit_target=1. if orbit<0 else float(capped_hit[orbit])
        escape_target=0. if orbit<0 else float(capped_escape[orbit])
        move_target=escape_target/2 if feedback else .5
        hit_mean=float(np.mean(observed_hit));move_mean=float(np.mean(observed_move))
        assert abs(hit_mean-hit_target)<=epsilon+preparation_bias_bound
        assert abs(move_mean-move_target)<=epsilon+preparation_bias_bound
        assert censors<=2
        summaries.append(dict(x=x,feedback=feedback,hit_mean=hit_mean,hit_target=hit_target,
            move_mean=move_mean,move_target=move_target,simultaneous_halfwidth=epsilon,
            censors=censors,actual_native_calls=case_events))
    assert summaries[-1]['move_mean']-summaries[1]['move_mean']>.05
    # Independent current-only continuation through every preparation, query and hold.
    replay={'i':12*ci[(0,0,0)],'o':12*ci[(0,0,0)]}
    for trial in trial_records:
        for which,eid,coin,j in trial['records']:
            if which=='z':assert j==0
            elif j:
                a,b=edges[eid];assert replay[which]==(a if j==1 else b)
                replay[which]=b if j==1 else a
        yy=cells[replay['i']//12]
        if trial['outcome']=='hit':assert yy==(0,0,0)
        if trial['outcome']=='escape':assert max(map(abs,yy))==2
        assert trial['held']==(trial['outcome']=='preparation_timeout' or
                              trial['feedback'] and trial['outcome']!='escape')
    assert replay==visible

    # Explicit finite operator consequence, not a hydrogen spectrum claim.
    outer_cells=list(itertools.product(range(-1,2),repeat=3));oi={x:i for i,x in enumerate(outer_cells)}
    ls=np.zeros((27,27))
    for x in outer_cells:
        for axis in range(3):
            xx=list(x);xx[axis]=(xx[axis]+2)%3-1
            i,j=oi[x],oi[tuple(xx)];ls[i,i]+=1;ls[j,j]+=1;ls[i,j]-=1;ls[j,i]-=1
    acceptance=np.asarray([0. if x==(0,0,0) else capped_escape[sum(v!=0 for v in x)-1]
                           for x in outer_cells])
    lint=np.zeros((12,12))
    for eid,(a,b) in enumerate(edges[:30]):
        # First cell's internal branches determine its exchange Laplacian.
        for slot,pos in ((0,a),(1,b)):
            for coin in (0,1):
                target=destination(pos,eid,coin)
                if target!=pos:lint[pos%12,pos%12]+=1;lint[pos%12,target%12]-=1
    source_bands=np.linalg.eigvalsh(lint);lam_values=[]
    for lam in source_bands:
        if lam>1e-8 and not any(abs(lam-v)<1e-8 for v in lam_values):lam_values.append(float(lam))
    wells=[]
    for lam in lam_values:
        matrix=ls+lam*np.diag(acceptance);eigs,vecs=np.linalg.eigh(matrix)
        control=float(np.linalg.eigvalsh(ls+lam*np.eye(27))[0])
        assert 0<eigs[0]<85*lam/99+1e-10 and abs(control-lam)<1e-10
        wells.append(dict(native_lambda=lam,lowest_squared_energy=float(eigs[0]),
            uniform_trial_bound=85*lam/99,feedback_disabled_threshold=control,
            ground_profile=vecs[:,0].tolist(),spatial_cells=outer_cells,
            hydrogen_ladder_established=False))
    raw=encoded(dict(graph=dict(cells=cells,edges=edges),branches=branches,trials=trial_records))
    (HERE/'RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    receipt=dict(schema='oph.native_probe_binding.v1',candidate='MFR-native-probe-binding',
        class_='N+',official='untested',actual_native_calls=calls,census_calls=len(branches),
        source_hit_probabilities=[str(x) for x in h],probe_mean_calls=[str(x) for x in mean],
        transient_Q=[[str(x) for x in row] for row in q],capped_hit=capped_hit.tolist(),
        capped_escape=capped_escape.tolist(),probe_censor_bound=tail_bound,
        preparation_censor_bound=preparation_bias_bound,confidence_family_coverage=.99,
        summaries=summaries,finite_wells=wells,record_only_replay=True,
        supplied_numeric_potential=False,supplied_gamma=False,all_registered_checks='pass',
        hydrogen_levels_established=False,physical_charge_assigned=False,
        raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest())
    (HERE/'RECEIPT.json').write_bytes(encoded(receipt))
    print(json.dumps({k:v for k,v in receipt.items() if k!='finite_wells'},indent=2,sort_keys=True))


if __name__=='__main__':main()
