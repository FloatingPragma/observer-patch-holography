#!/usr/bin/env python3
"""Independent full control replay; prior corrected source custody retained."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import gzip,hashlib,json,math,random,re
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
P=HERE.parent.parent/'work/2026-11-14_codex_unified_source_join'
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))

def make_graph(inner):
    cells=list(product(range(-2,3),repeat=3)) if inner else list(product(range(-1,2),range(-1,2),range(-1,2),range(3)))
    ci={v:i for i,v in enumerate(cells)};edges=[]
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})
    for x in cells:
        c=ci[x];edges += [(12*c+a,12*c+b,-1) for a,b in internal]
        for axis in range(len(x)):
            y=list(x);y[axis]+=1
            if inner:
                if y[axis]>2:continue
            elif axis==3:y[axis]%=3
            elif y[axis]>1:y[axis]=-1
            edges += [(12*c+p,12*ci[tuple(y)]+p,axis) for p in range(12)]
    adj=[[] for _ in range(12*len(cells))];pair={}
    for eid,(a,b,_) in enumerate(edges):adj[a].append(eid);adj[b].append(eid);pair[frozenset((a,b))]=eid
    return cells,ci,edges,adj,pair,internal

def main():
    oldhash='f36f7e8e7b1c179632fb70fdc17de0571431ee423842e4e5f9705fc49121a001'
    rawhash='6d1f37566679caef2493a85e9cd7de12f769f3c68d9165f1e59f6c8c242b0beb'
    assert hashlib.sha256((P/'inputs/PRODUCTION_RECORDS.json.gz').read_bytes()).hexdigest()==oldhash
    assert hashlib.sha256((P/'CONTROL_RECORDS.json.gz').read_bytes()).hexdigest()==rawhash
    receipt=json.loads((P/'RECEIPT.json').read_text());raw=gzip.decompress((P/'CONTROL_RECORDS.json.gz').read_bytes());r=json.loads(raw)
    assert hashlib.sha256(raw).hexdigest()==receipt['control_raw_uncompressed_sha256']
    oc,oi,oe,oa,op,internal=make_graph(False);ic,ii,ie,ia,ip,_=make_graph(True);W=len(oe);NV=972;n=32768
    rng=random.Random(2026111435);outer=inner=internal_moves=downward=0;counts=[];stats=[]
    for ch in r['channels']:
        reader=rng.randrange(NV);assert reader==ch['start'];ipos=12*ii[(0,0,0)];origin=oc[reader//12][3]
        ticks=0;clock=[];occ=[0]*NV;calls=0;censors=[]
        for step,entry in enumerate(ch['trace'],1):
            eid=rng.randrange(W);rng.randrange(2);coin=0;a,b,axis=oe[eid];pc=qc=0;outcome='skip';events=[]
            def inner_step(eid,coin):
                nonlocal ipos,inner,calls
                a,b,_=ie[eid];j=int(ipos==a)-coin if ipos in (a,b) else 0
                if j:ipos=b if j==1 else a
                events.append([eid,coin,j]);inner+=1;calls+=1;return j
            if axis==-1 and reader in (a,b):
                target=oc[reader//12][:3];failed=False
                for ax in range(3):
                    while ic[ipos//12][ax]!=target[ax]:
                        z=list(ic[ipos//12]);z[ax]+=1 if target[ax]>z[ax] else -1
                        dest=12*ii[tuple(z)]+ipos%12;pe=ip[frozenset((ipos,dest))]
                        for attempt in range(128):
                            j=inner_step(pe,rng.randrange(2));pc+=1
                            if j:break
                        else:failed=True
                        if failed:break
                    if failed:break
                if failed:outcome='preparation_timeout'
                else:
                    for attempt in range(257):
                        z=ic[ipos//12]
                        if z==(0,0,0):outcome='hit';break
                        if max(map(abs,z))==2:outcome='escape';break
                        if attempt==256:outcome='probe_timeout';break
                        pe=rng.choice(ia[ipos]);inner_step(pe,rng.randrange(2));qc+=1
            previous=reader;j=int(reader==a) if reader in (a,b) else 0
            if j:reader=b
            outer+=1;calls+=1;occ[reader]+=1
            if j and reader%12!=previous%12:
                internal_moves+=1;downward+=int(reader%12<previous%12)
            tick=None
            if oc[reader//12][3]==(origin+ticks+1)%3:ticks+=1;tick=[step,ticks,reader];clock.append(tick)
            expected=dict(proposal=eid,coin=0,current=j,auxiliary_hold=False,target_coupled=False,
                inner_events=events,preparation_calls=pc,probe_calls=qc,outcome=outcome,tick=tick)
            assert entry==expected
            if 'timeout' in outcome:censors.append([step,outcome])
        assert ch['occupation_counts']==occ and ch['clock']==clock and ch['censors']==censors
        assert calls==ch['native_calls'] and n-(clock[-1][0] if clock else 0)==ch['final_open_interval']
        counts.append(occ);stats.append({'calls':calls,'ticks':ticks,'open_interval':ch['final_open_interval'],'censors':censors})
    assert outer==131072 and outer+inner==receipt['biased_tape_native_calls']
    assert internal_moves==receipt['biased_observed_internal_moves'] and downward==receipt['biased_observed_downward_port_moves']==0
    assert receipt['total_native_calls']==347533+2*len(oe)+outer+inner
    ys=np.array([math.sqrt(NV/(4*W*n))*(np.array(c)-n/NV) for c in counts]);av=(ys.reshape(4,81,12)@np.ones(12)/math.sqrt(12)).reshape(4,3,3,3,3)
    dirs=(3,0,1,2);curl=[np.roll(av[nu],-1,dirs[mu])-av[nu]-np.roll(av[mu],-1,dirs[nu])+av[mu] for mu in range(4) for nu in range(mu+1,4)]
    assert np.max(abs(np.array(curl)-r['field_strength_readouts']))<1e-14
    # Derive all directed branches and exact cell lump, without producer imports.
    L=np.zeros((NV,NV),dtype=np.int64)
    for a,b,axis in oe:L[a,a]+=1;L[a,b]-=1
    cellbasis=np.kron(np.eye(81,dtype=np.int64),np.ones((12,1),dtype=np.int64));C=L@cellbasis;cell=C[::12]
    assert np.array_equal(C,np.repeat(cell,12,axis=0)) and np.all(cell.sum(axis=0)==0)
    for pos,z in enumerate(oc):
        expected={pos:4}
        for ax in range(4):
            y=list(z);y[ax]+=1
            if ax==3:y[ax]%=3
            elif y[ax]>1:y[ax]=-1
            expected[oi[tuple(y)]]=-1
        assert {j:int(v) for j,v in enumerate(cell[pos]) if v}==expected
    sinks=[i for i in range(12) if not any(a==i for a,b in internal)]
    assert sinks==receipt['biased_absorbing_internal_ports']==[9,11]
    exact=[]
    for bsq,saved in zip((3,0),receipt['changed_law_same_F12_population']):
        covariance=F(3,2*(9+bsq))-F(1,4*W);fair=F(1,6)-F(1,4*W)
        assert abs(float(covariance)-saved['biased_covariance'])<1e-12
        assert abs(float(6*covariance)-saved['biased_F12_covariance'])<1e-12
        assert abs(float(fair)-saved['fair_covariance'])<1e-12
        exact.append(covariance)
    assert 6*(exact[1]-exact[0])==F(1,4)
    # Independently check all 144 graph-label entries against frozen Lean text.
    iso=json.loads((P/'PORT_ISOMORPHISM.json').read_text());lean=ROOT/'reverse-engineering-reality/Lean/Screen/DiscreteCoulombGreen.lean'
    assert hashlib.sha256(lean.read_bytes()).hexdigest()==iso['Lean_source_sha256']
    block=lean.read_text().split('def laplacianEntriesZ')[1].split('def greenNumZ')[0]
    rows=re.findall(r'!\[([^]\[]+)\]',block);mat=np.array([[int(x.strip()) for x in row.split(',')] for row in rows[:12]],dtype=np.int64)
    lint=np.zeros((12,12),dtype=np.int64)
    for a,b in internal:lint[a,a]+=1;lint[b,b]+=1;lint[a,b]-=1;lint[b,a]-=1
    perm=[iso['source_port_to_Lean_port'][str(i)] for i in range(12)]
    assert np.array_equal(lint,mat[np.ix_(perm,perm)])
    o={'pass':True,'same_group':True,'G6_pass':False,'source_ancestor_records_hash':oldhash,
       'control_native_events_independently_replayed':outer+inner,'control_macro_records':outer,'control_inner_records':inner,
       'all_rng_routes_current_cost_clock_and_censor_fields_checked':True,'full_directed_cell_lump_verified':True,
       'equal_norm_F12_population_contrast_exact':'1/4','graph_label_entries_checked':144,'stats':stats,
       'population_is_short_tape_estimate':False,'moving_target_covariance_proved':False,'minimal_electric_coupling_proved':False,
       'receipt_sha256':hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest(),'control_archive_sha256':rawhash}
    (HERE/'VERIFICATION.json').write_text(json.dumps(o,indent=2,sort_keys=True)+'\n');print(json.dumps(o))

if __name__=='__main__':main()
