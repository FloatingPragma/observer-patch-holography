#!/usr/bin/env python3
"""Native first-hit record clock; retains censoring and changed schedules."""
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
R=2
RUNS=4096
CAP=4096
SEED=2026111430
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
       (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
       (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))


def encoded(x):
    return (json.dumps(x,sort_keys=True,separators=(",",":"),
                       allow_nan=False)+"\n").encode()


def main():
    path=HERE/'inputs/native.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()== \
        'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('first_hit_native',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod
    spec.loader.exec_module(mod);native=mod.native
    cells=list(itertools.product(range(3),range(3),range(3),range(-R,2)))
    ci={z:i for i,z in enumerate(cells)}
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})
    edges=[]
    for z in cells:
        cell=ci[z]
        edges.extend((12*cell+a,12*cell+b,-1) for a,b in internal)
        for axis in range(3):
            zz=list(z);zz[axis]=(zz[axis]+1)%3
            edges.extend((12*cell+p,12*ci[tuple(zz)]+p,axis) for p in range(12))
        if z[3]<1:
            zz=z[:3]+(z[3]+1,)
            edges.extend((12*cell+p,12*ci[zz]+p,3) for p in range(12))
    nv=12*len(cells);adj=[[] for _ in range(nv)];vertical=[[] for _ in range(nv)]
    for eid,(a,b,axis) in enumerate(edges):
        adj[a].append(eid);adj[b].append(eid)
        if axis==3:vertical[a].append(eid);vertical[b].append(eid)
    assert (len(cells),nv,len(edges))==(108,1296,8100)
    branches=[];records={}
    for eid in range(len(edges)):
        for slot in range(3):
            loads=[0,0,0];loads[slot]=1
            for coin in (0,1):
                out,rec=native(loads,(0,1),coin)
                rebuilt=list(loads);rebuilt[0]-=rec[1];rebuilt[1]+=rec[1]
                assert tuple(rebuilt)==out
                branches.append([eid,slot,coin,rec[1]])
                records[eid,slot,coin]=rec[1]
    calls=len(branches)

    def endpoint(pos,eid,coin,call=False):
        nonlocal calls
        a,b,_=edges[eid]
        if call:
            _,rec=native((int(pos==a),int(pos==b)),(0,1),coin)
            calls+=1;j=rec[1]
        else:
            j=records[eid,0 if pos==a else 1,coin]
        return (b if j==1 else a if j==-1 else pos),j

    indices=((0,0,0),(1,0,0),(1,1,0),(1,1,1))
    population=[]
    for index in indices:
        p=np.asarray(index)*2*np.pi/3
        per_port=[]
        for port in range(12):
            q=np.zeros((R+1,R+1),dtype=complex);bvec=np.zeros(R+1,dtype=complex)
            for height in range(-R,1):
                source=12*ci[(0,0,0,height)]+port
                for eid in adj[source]:
                    for coin in (0,1):
                        target,_=endpoint(source,eid,coin)
                        zz=cells[target//12]
                        value=np.exp(1j*np.dot(p,zz[:3]))/(2*len(adj[source]))
                        if zz[3]==1:bvec[height+R]+=value
                        else:q[height+R,zz[3]+R]+=value
            per_port.append((q,bvec))
        q,bvec=per_port[0]
        assert max(np.max(np.abs(q-x))+np.max(np.abs(bvec-y)) for x,y in per_port)<1e-12
        assert np.max(np.abs(q.imag))<1e-12 and np.max(np.abs(bvec.imag))<1e-12
        q=q.real;bvec=bvec.real
        hit=np.linalg.solve(np.eye(R+1)-q,bvec)
        censored=np.linalg.solve(np.eye(R+1)-q,(np.eye(R+1)-np.linalg.matrix_power(q,CAP))@bvec)
        if index==(0,0,0):
            tail=float((np.linalg.matrix_power(q,CAP)@np.ones(R+1))[R])
            cost=float(np.linalg.solve(np.eye(R+1)-q,np.ones(R+1))[R])
            assert abs(hit[R]-1)<1e-12 and abs(cost-76)<1e-10 and tail<1e-4
        else:
            a=float(np.sum((2*np.sin(p/2))**2));omega=2*math.asinh(math.sqrt(a)/2)
            target=math.cosh(omega*(R+.5))/math.cosh(omega*(R+1.5))
            error=hit[R]/math.exp(-omega)-1
            assert abs(hit[R]-target)<1e-12
            assert -1e-12<=error<=math.exp(-2*omega*(R+.5))+1e-12
            population.append(dict(index=list(index),kernel=q.tolist(),absorption=bvec.tolist(),
                hit_transform=float(hit[R]),capped_transform=float(censored[R]),
                target_reflecting_strip=target,infinite_depth_transfer=math.exp(-omega),
                relative_depth_error=float(error),slice_energy=omega))
    vertical_q=np.zeros((R+1,R+1))
    for height in range(-R,1):
        source=12*ci[(0,0,0,height)]
        for eid in vertical[source]:
            for coin in (0,1):
                target,_=endpoint(source,eid,coin)
                zz=cells[target//12]
                if zz[3]<1:
                    vertical_q[height+R,zz[3]+R]+=1/(2*len(vertical[source]))
    vertical_tail=float((np.linalg.matrix_power(vertical_q,CAP)@np.ones(R+1))[R])
    rng=random.Random(SEED);runs={};summaries={}
    epsilon=1.5*math.sqrt(math.log(12/.01)/(2*RUNS))
    for law in ('incident','vertical_only_control'):
        series=[];outcomes=[];lengths=[];censored_count=0
        for trial in range(RUNS):
            start=12*ci[(0,0,0,0)]+trial%12
            source=reader=start;tape=[]
            for step in range(CAP):
                eligible=adj[reader] if law=='incident' else vertical[reader]
                eid=eligible[rng.randrange(len(eligible))];coin=rng.randrange(2)
                source,j=endpoint(source,eid,coin,call=True)
                a,b,_=edges[eid]
                if j:
                    assert reader==(a if j==1 else b)
                    reader=b if j==1 else a
                assert source==reader
                tape.append([eid,coin,j])
                if cells[reader//12][3]==1:break
            success=cells[reader//12][3]==1
            censored_count+=not success
            outcome=[]
            for item in population:
                phase=int(np.dot(item['index'],cells[reader//12][:3]))%3
                outcome.append((1. if phase==0 else -.5) if success else 0.)
            outcomes.append(outcome);lengths.append(len(tape))
            series.append(dict(start=start,tape=tape,completed=success))
        means=np.mean(outcomes,axis=0)
        target=np.asarray([x['capped_transform'] for x in population]) if law=='incident' else np.full(3,1-vertical_tail)
        assert max(abs(means-target))<=epsilon
        assert censored_count<=4
        runs[law]=series
        summaries[law]=dict(sample_means=means.tolist(),population_targets=target.tolist(),
            simultaneous_hoeffding_halfwidth=epsilon,simultaneous_coverage_lower_bound=.99,
            censored=censored_count,actual_attempts=sum(lengths),mean_attempts=float(np.mean(lengths)),
            max_attempts=max(lengths))
    assert min(np.asarray(summaries['vertical_only_control']['sample_means'])-
               np.asarray(summaries['incident']['sample_means']))>.5
    # Positive reader control: point endpoint has transform1; uniform spatial
    # endpoints have transform0 at every registered nonzero momentum.
    uniform=[sum(math.cos(2*math.pi/3*np.dot(item['index'],z))
                 for z in itertools.product(range(3),repeat=3))/27 for item in population]
    assert max(abs(x) for x in uniform)<1e-12
    # Replay reads only declared start and retained current letters.
    for series in runs.values():
        for run in series:
            location=run['start']
            for eid,coin,j in run['tape']:
                assert cells[location//12][3]<1
                a,b,_=edges[eid]
                if j:
                    assert location==(a if j==1 else b)
                    location=b if j==1 else a
            assert (cells[location//12][3]==1)==run['completed']
    raw=encoded(dict(graph=dict(cells=cells,edges=edges),branches=branches,runs=runs))
    (HERE/'RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    receipt=dict(schema='oph.native_first_hit_clock.v1',candidate='MFR-first-hit-clock',
        class_='N+',official='untested',branch_calls=len(branches),actual_native_calls=calls,
        population=population,source_tail_at_cap=tail,control_tail_at_cap=vertical_tail,
        source_exact_mean_attempts=cost,
        summaries=summaries,record_only_replay=True,raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest(),
        all_registered_checks='pass',clock='native first-hit Euclidean transfer',
        lorentzian_description='effective reconstruction, not primitive attempt evolution')
    (HERE/'RECEIPT.json').write_bytes(encoded(receipt))
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__':main()
