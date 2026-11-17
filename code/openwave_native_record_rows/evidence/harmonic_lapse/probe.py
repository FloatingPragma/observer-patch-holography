"""Native plus-hit metric gates and nested reciprocal hitting experiment."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,math,random,struct,sys,time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
engine=module('gr7_frozen_probe_engine',HERE/'inputs/native_probe_engine.py');engine.HERE=HERE
STARTS=((1,0,0),(1,1,0),(1,1,1));K=512
METRIC=(('internal_A',0,'internal'),('internal_B',1,'internal'),('internal_C',2,'internal'),
 ('temporal_A',0,'temporal'),('temporal_B',1,'temporal'),('temporal_C',2,'temporal'),
 ('spatial',0,'spatial'),('ignore_hits_A',0,'spatial'),('hold',0,'hold'))
OUTER=(('weighted_A',0,1024,True),('weighted_B',1,512,True),
       ('weighted_C',2,512,True),('ungated_A',0,512,False))
def dump(name,obj): (HERE/name).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def solve(A,b):
    q=[list(row)+[b[i]] for i,row in enumerate(A)]
    for j in range(len(q)):
        k=next(k for k in range(j,len(q)) if q[k][j]);q[j],q[k]=q[k],q[j]
        d=q[j][j];q[j]=[x/d for x in q[j]]
        for i in range(len(q)):
            if i!=j:
                f=q[i][j];q[i]=[x-f*y for x,y in zip(q[i],q[j])]
    return [row[-1] for row in q]

def population(cap,tail):
    hit=[1-cap[i]-tail[i] for i in range(3)];b=[(1+x)/2 for x in hit]
    neighbors=(((1,4),(3,1),(4,1)),((0,2),(2,2),(4,2)),((1,3),(4,3)))
    Q=[[F(i==j) for j in range(3)] for i in range(3)];target=[F(0)]*3
    for i in range(3):
      for j,count in neighbors[i]:
        bj=b[j] if j<3 else (F(1) if j==3 else F(1,2))
        move=F(count,12)*b[i]*bj;Q[i][i]-=move
        if j<3:Q[i][j]+=move
        elif j==3:target[i]+=move
    H=solve([[F(i==j)-Q[i][j] for j in range(3)] for i in range(3)],target)
    uncapped=(F(10,27),F(1,6),F(2,23));hit_exact=(F(5,22),F(1,11),F(1,22))
    inner_error=max(hit_exact[i]-hit[i] for i in range(3))
    assert all(abs(H[i]-uncapped[i])<=24*inner_error for i in range(3))
    outer_tail=F(47,48)**K
    assert outer_tail<F(1,40000)
    dump('POPULATION.json',{'hit_cap':[str(x) for x in hit],'plus_cap':[str(x) for x in b],
      'weighted_Q':[[str(x) for x in row] for row in Q],'weighted_b':[str(x) for x in target],
      'weighted_hit_infinite':[str(x) for x in H],'uncapped_reciprocal_hit':[str(x) for x in uncapped],
      'inner_cap_error':str(inner_error),'reciprocal_cap_error_bound':str(24*inner_error),
      'outer_cap_tail_bound':str(outer_tail)})
    return b,H,outer_tail

def plus(native,rng,pos,start,raw):
    pos=engine.event(native,rng,pos,6,raw);coin=(raw[-1]>>2)&1
    pos,rec=engine.probe_word(native,rng,pos,start,raw)
    return pos,{'coin':coin,'probe':rec,'plus':bool(coin==0 or rec['outcome']=='hit')}

def reset(native,rng,pos,target,raw):
    failed=False;moves=0
    for _ in range(9):
        if failed or pos==target:
            for _ in range(16):pos=engine.event(native,rng,pos,6,raw)
        else:
            axis=next(i for i in range(3) if pos[i]!=target[i]);direction=2*axis+int(target[axis]<pos[axis])
            origin=pos
            for _ in range(16):pos=engine.event(native,rng,pos,6 if pos!=origin else direction,raw)
            if pos==origin:failed=True
            else:moves+=1
    assert failed or pos==target
    return pos,failed,moves

def metric(native,b):
    inner=STARTS[0];main=0;records=[];results=[];checks=[]
    predictions=[b[0]**2/4,b[1]**2/4,b[2]**2/4,b[0]**4,b[1]**4,b[2]**4,F(1,16),F(1,16),F(0)]
    tolerance=F(47,1000);reset_error=F(36,65536)
    with gzip.open(HERE/'METRIC_EVENTS.bin.gz','wb') as tape:
      for ci,(name,where,rule) in enumerate(METRIC):
        rng=random.Random(202611150727+ci);accepted=0;moved=0
        for trial in range(2048):
            raw=bytearray();initial_inner=inner;initial_main=main;slots=[]
            for _ in range(4):inner,s=plus(native,rng,inner,STARTS[where],raw);slots.append(s)
            if rule=='internal':gate=slots[0]['plus'] and slots[1]['plus'] and slots[2]['coin']==0 and slots[3]['coin']==0
            elif rule=='temporal':gate=all(s['plus'] for s in slots)
            elif rule=='spatial':gate=all(s['coin']==0 for s in slots)
            else:gate=False
            coin=rng.randrange(2);values=tuple(int(i==main) for i in range(4))
            out,(_,j)=engine.call(native,values,(0,1) if gate else (2,3),coin)
            main=out.index(1);assert j*(1 if initial_main==0 else -1)==abs(j)
            raw.append(((7 if gate else 6)<<3)|(coin<<2)|(j+1));assert len(raw)==837
            tape.write(raw);accepted+=gate;moved+=abs(j)
            records.append({'case':name,'index':trial,'inner_start':initial_inner,'inner_end':inner,
                'main_start':initial_main,'main_end':main,'slots':slots,'gate':bool(gate),'current':j})
        p=predictions[ci]
        ok={'acceptance':abs(F(accepted,2048)-p)<=tolerance+reset_error,
            'main_current':abs(F(moved,2048)-p/2)<=tolerance+reset_error/2}
        checks.extend(ok.values());results.append({'case':name,'accepted':accepted,'moved':moved,
             'target':str(p),'checks':ok});print(json.dumps(results[-1]),flush=True)
    contrast=F(results[3]['accepted']-results[7]['accepted'],2048)>F(1,50)
    checks.append(contrast)
    dump('METRIC_TRIALS.json',records);dump('METRIC_RESULTS.json',{'cases':results,'changed_ignore_hit_contrast':contrast,
         'martingale_tolerance':str(tolerance),'reset_offset':str(reset_error),'all_pass':all(checks)})
    return inner,all(checks)

def outer_experiment(native,inner,H,tail):
    outer=STARTS[0];records=[];results=[];checks=[];calls=0
    deterministic=F(225,65536)+tail
    with gzip.open(HERE/'OUTER_EVENTS.bin.gz','wb') as tape:
      for ci,(name,where,n,weighted) in enumerate(OUTER):
        rng=random.Random(202611150827+ci);hits=0
        for trial in range(n):
            raw=bytearray();initial_outer=outer;initial_inner=inner
            outer,failed,reset_moves=reset(native,rng,outer,STARTS[where],raw);calls+=144
            status='reset_failed' if failed else 'active';steps=[]
            for step in range(K):
                if status!='active':break
                direction=rng.randrange(6);neighbor=engine.advance(outer,direction)
                raw.append(direction);slots=[]
                for start in (outer,neighbor):inner,s=plus(native,rng,inner,start,raw);slots.append(s)
                gate=all(s['plus'] for s in slots) if weighted else True
                coin=rng.randrange(2)
                out,(_,j)=engine.call(native,(1,0) if gate else (0,0),(0,1),coin)
                assert out==((1-j,j) if gate else (0,0))
                raw.append(((direction if gate else 6)<<3)|(coin<<2)|(j+1));calls+=419
                if j:outer=neighbor
                steps.append({'direction':direction,'slots':slots,'gate':bool(gate),'coin':coin,'current':j,'end':outer})
                if outer==(0,0,0):status='hit'
                elif max(map(abs,outer))==2:status='escape'
            if status=='active':status='timeout'
            hits+=status=='hit';assert len(raw)==144+420*len(steps)
            tape.write(struct.pack('<I',len(raw)));tape.write(raw)
            records.append({'case':name,'index':trial,'inner_start':initial_inner,'inner_end':inner,
                'outer_start':initial_outer,'outer_end':outer,'reset_failed':failed,'reset_moves':reset_moves,
                'outcome':status,'steps':steps})
        target=H[where] if weighted else F(5,22);eps=F(65,1000) if n==1024 else F(9,100)
        ok=abs(F(hits,n)-target)<=eps+deterministic
        checks.append(ok);results.append({'case':name,'trials':n,'hits':hits,'target':str(target),
             'martingale_tolerance':str(eps),'deterministic_error':str(deterministic),'pass':ok})
        print(json.dumps(results[-1]),flush=True)
    changed=F(results[0]['hits'],1024)-F(results[3]['hits'],512)>F(1,25)
    null=F(results[0]['hits'],1024)>F(5,22)+F(1,10)
    checks.extend([changed,null])
    dump('OUTER_TRIALS.json',records);dump('OUTER_RESULTS.json',{'cases':results,'changed_ungated_contrast':changed,
        'ungated_null_rejected':null,'all_pass':all(checks),'native_calls':calls,
        'reset_failures':sum(r['outcome']=='reset_failed' for r in records),
        'timeouts':sum(r['outcome']=='timeout' for r in records),'steps':sum(len(r['steps']) for r in records)})
    return all(checks),calls

def run():
    started=time.time();native_module=module('gr7_native_pin',HERE/'inputs/native.py')
    assert hashlib.sha256((HERE/'inputs/native.py').read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    cap,tail=engine.source_census(native_module.native);engine.embedding(native_module.native);engine.internal_spectrum(native_module.native)
    b,H,outer_tail=population(cap,tail)
    inner,ok1=metric(native_module.native,b);ok2,outer_calls=outer_experiment(native_module.native,inner,H,outer_tail)
    assert engine.CALLS==1824+9*2048*837+outer_calls
    r={'all_pass':ok1 and ok2,'literal_native_calls':engine.CALLS,'metric_calls':9*2048*837,
       'outer_calls':outer_calls,'pretrial_calls':1824,'metric_proposals':9*2048,'outer_trials':2560,
       'initial_units':3,'fresh_preparations':0,'class':'N+','curved_Maxwell_dynamics':False,
       'physical_G':False,'G6':False,'wall_seconds':time.time()-started}
    dump('RECEIPT.json',r);print(json.dumps(r,sort_keys=True),flush=True)
    return r
if __name__=='__main__':raise SystemExit(0 if run()['all_pass'] else 1)
