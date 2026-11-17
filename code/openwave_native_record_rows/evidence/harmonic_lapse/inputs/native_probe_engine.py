"""Actual native reset/probe/gate controller; no potential lookup in policy."""
from __future__ import annotations
import gzip,hashlib,importlib.util,itertools,json,math,random,sys,time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
T=64; RETRIES=16; RESET_STEPS=9; TRIALS=2048
SLOTS=4*(RESET_STEPS*RETRIES+T)+1
DIRS=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
STARTS=((1,0,0),(1,1,0),(1,1,1))
CASES=(('internal_A',(0,0,0,0),'two'),('internal_B',(1,1,1,1),'two'),
 ('internal_C',(2,2,2,2),'two'),('spatial_AB',(0,0,1,1),'four'),
 ('spatial_BC',(1,1,2,2),'four'),('temporal',(0,0,0,0),'all'),
 ('one_miss_A',(0,0,0,0),'one'),('hold',(0,0,0,0),'hold'))
CALLS=0

def dump(name,obj): (HERE/name).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def call(native,loads,seam,coin):
    global CALLS
    CALLS+=1
    out,record=native(loads,seam,coin)
    assert record[0]==coin
    return out,record
def advance(pos,direction): return tuple(a+b for a,b in zip(pos,DIRS[direction]))
def kind(pos): return sum(v!=0 for v in pos)-1
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def source_census(native):
    rows=[];class_counts=[None]*3
    for pos in itertools.product((-1,0,1),repeat=3):
        if pos==(0,0,0):continue
        count=[0]*5 # A,B,C,target,escape
        for direction in range(6):
          for coin in (0,1):
            out,(_,j)=call(native,(1,0),(0,1),coin)
            assert out==(1-j,j)
            end=advance(pos,direction) if j else pos
            terminal=3 if end==(0,0,0) else (4 if max(map(abs,end))==2 else kind(end))
            count[terminal]+=1;rows.append([pos,direction,coin,j,end,terminal])
        k=kind(pos)
        if class_counts[k] is None:class_counts[k]=count
        assert class_counts[k]==count
    assert class_counts==[[6,4,0,1,1],[2,6,2,0,2],[0,3,6,0,3]]
    Q=[[F(v,12) for v in row[:3]] for row in class_counts]
    r=[F(17,22),F(10,11),F(21,22)]
    for i in range(3): assert r[i]==F(class_counts[i][4],12)+sum(Q[i][j]*r[j] for j in range(3))
    QT=[[F(i==j) for j in range(3)] for i in range(3)]
    for _ in range(T): QT=mm(QT,Q)
    tail=[sum(row) for row in QT]
    cap=[r[i]-sum(QT[i][j]*r[j] for j in range(3)) for i in range(3)]
    assert max(tail)<F(1,100000)
    dump('HIT_CENSUS.json',{'rows':rows,'class_counts':class_counts,'Q':[[str(x) for x in row] for row in Q],
         'r':[str(x) for x in r],'r_cap':[str(x) for x in cap],'tail':[str(x) for x in tail]})
    return cap,tail

def embedding(native):
    rows=[]
    for a in range(12):
      for b in range(a+1,12):
        outside=next(x for x in range(12) if x not in (a,b))
        for site in (a,b,outside):
          loads=tuple(int(i==site) for i in range(12))
          for coin in (0,1):
            full,rec=call(native,loads,(a,b),coin)
            pair,rp=call(native,(loads[a],loads[b]),(0,1),coin)
            expected=list(loads);expected[a],expected[b]=pair
            assert tuple(expected)==full and rec==rp
            rows.append([a,b,site,coin,rec[1],list(pair)])
    dump('EMBEDDING.json',rows)

def internal_spectrum(native):
    faces=json.loads((HERE/'inputs/faces.json').read_text())
    edges=sorted({tuple(sorted((face[j],face[(j+1)%3]))) for face in faces for j in range(3)})
    assert len(edges)==30
    C=[[0]*12 for _ in range(12)];rows=[]
    for x in range(12):
      initial=tuple(int(i==x) for i in range(12))
      for ei,e in enumerate(edges):
        for coin in (0,1):
            out,(_,j)=call(native,initial,e,coin);y=out.index(1)
            C[x][y]+=1;rows.append([x,ei,coin,j,y])
    L=[[60*int(i==j)-C[i][j] for j in range(12)] for i in range(12)]
    assert all(sum(row)==0 for row in L) and all(L[i][i]==5 for i in range(12))
    L2=mm(L,L)
    # Source polynomial z(z-6)(z^2-10z+20) vanishes exactly.
    factor=[[L2[i][j]-10*L[i][j]+20*int(i==j) for j in range(12)] for i in range(12)]
    product=mm(mm(L,[[L[i][j]-6*int(i==j) for j in range(12)] for i in range(12)]),factor)
    assert all(x==0 for row in product for x in row)
    dump('INTERNAL_SPECTRUM.json',{'edges':edges,'branch_rows':rows,'counts':C,'L':L,
         'annihilating_polynomial':[0,-120,80,-16,1],
         'positive_values':['5-sqrt(5)','6','5+sqrt(5)'],'polynomial_verified':True})

def event(native,rng,pos,direction,raw):
    coin=rng.randrange(2)
    if direction==6:
        out,(_,j)=call(native,(0,0),(0,1),coin);assert out==(0,0) and j==0
    else:
        out,(_,j)=call(native,(1,0),(0,1),coin);assert out==(1-j,j)
        if j:pos=advance(pos,direction)
    raw.append((direction<<3)|(coin<<2)|(j+1))
    return pos

def probe_word(native,rng,pos,target,raw):
    reset_failed=False;reset_moves=0
    for slot in range(RESET_STEPS):
        if reset_failed or pos==target:
            for _ in range(RETRIES):pos=event(native,rng,pos,6,raw)
            continue
        axis=next(i for i in range(3) if pos[i]!=target[i])
        direction=2*axis+int(target[axis]<pos[axis]);old=pos;moved=False
        for _ in range(RETRIES):
            pos=event(native,rng,pos,6 if moved else direction,raw)
            if pos!=old:moved=True
        if not moved:reset_failed=True
        else:reset_moves+=1
    assert reset_failed or pos==target
    outcome=('reset_failed' if reset_failed else ('hit' if pos==(0,0,0) else ('escape' if max(map(abs,pos))==2 else 'active')));active_calls=0
    for _ in range(T):
        direction=rng.randrange(6) if outcome=='active' else 6
        if outcome=='active':active_calls+=1
        pos=event(native,rng,pos,direction,raw)
        if outcome=='active':
            if pos==(0,0,0):outcome='hit'
            elif max(map(abs,pos))==2:outcome='escape'
    if outcome=='active':outcome='timeout'
    return pos,{'outcome':outcome,'miss':int(outcome=='escape'),'reset_moves':reset_moves,'active_probe_calls':active_calls}

def experiment(native,cap):
    probe_pos=STARTS[0];main_pos=0;results=[];trials=[]
    reset_error=F(36,1<<RETRIES);epsilon=F(9,200)
    all_checks=[];reset_failures=0;timeouts=0
    predictions=(cap[0]**2,cap[1]**2,cap[2]**2,cap[0]**2*cap[1]**2,
                 cap[1]**2*cap[2]**2,F(1),cap[0],F(0))
    with gzip.open(HERE/'EVENTS.bin.gz','wb') as stream:
      for case_index,(name,starts,rule) in enumerate(CASES):
        rng=random.Random(202611150626+case_index);accepts=0;moves=0
        for trial in range(TRIALS):
            before_probe=probe_pos;before_main=main_pos;raw=bytearray();probes=[]
            for which in starts:
                probe_pos,rec=probe_word(native,rng,probe_pos,STARTS[which],raw);probes.append(rec)
                reset_failures+=int(rec['outcome']=='reset_failed');timeouts+=int(rec['outcome']=='timeout')
            bits=[q['miss'] for q in probes]
            accepted=(bool(bits[0] and bits[1]) if rule=='two' else
                     (all(bits) if rule=='four' else (True if rule=='all' else (bool(bits[0]) if rule=='one' else False))))
            coin=rng.randrange(2);initial=tuple(int(i==main_pos) for i in range(4))
            seam=(0,1) if accepted else (2,3)
            out,(_,current)=call(native,initial,seam,coin)
            # Continuity on the active pair, including either current sign.
            main_pos=out.index(1)
            assert main_pos in (0,1) and abs(current)<=1 and (accepted or current==0)
            aligned=current*(1 if before_main==0 else -1)
            assert aligned==abs(current)
            raw.append(((7 if accepted else 6)<<3)|(coin<<2)|(current+1))
            assert len(raw)==SLOTS
            stream.write(raw);accepts+=accepted;moves+=abs(current)
            trials.append({'case':name,'index':trial,'probe_start':before_probe,'probe_end':probe_pos,
                'main_start':before_main,'main_end':main_pos,'probes':probes,'accepted':accepted,'current':current})
        p=predictions[case_index]
        checks={'acceptance':abs(F(accepts,TRIALS)-p)<=epsilon+reset_error,
                'main_current':abs(F(moves,TRIALS)-p/2)<=epsilon+reset_error/2}
        all_checks.extend(checks.values())
        results.append({'case':name,'trials':TRIALS,'accepts':accepts,'moves':moves,
                        'target_acceptance':str(p),'target_move':str(p/2),'checks':checks})
        print(json.dumps(results[-1]),flush=True)
    changed=F(results[6]['accepts']-results[0]['accepts'],TRIALS)>F(2,25)
    all_checks.append(changed)
    dump('TRIALS.json',trials)
    dump('CONTROLLER_RESULTS.json',{'cases':results,'changed_one_miss_separated':changed,
        'reset_failures':reset_failures,'probe_timeouts':timeouts,'conditional_reset_error':str(reset_error),
        'martingale_epsilon':str(epsilon),'simultaneous_failure_bound':'32*exp(-2*2048*(9/200)^2)<0.01',
        'record_bytes':len(trials)*SLOTS,'initial_probe':STARTS[0],'initial_main':0})
    return all(all_checks),reset_failures,timeouts

def clock_calculation(cap):
    output=[]
    for cls,r in enumerate(cap):
      for a in (F(1),F(1,2),F(1,4),F(1,8),F(1,16)):
        row=[]
        for lam in (5-math.sqrt(5),6,5+math.sqrt(5)):
            ratio=math.asinh(float(a*r)*math.sqrt(lam)/2)/math.asinh(float(a)*math.sqrt(lam)/2)
            row.append(ratio)
        output.append({'class':cls,'a':str(a),'r_cap':str(r),'clock_ratios_display_only':row})
    bound=F(10,24*256-10)
    assert bound<F(1,500)
    dump('CLOCK_POPULATION.json',{'values':output,'final_uniform_rational_error_bound':str(bound),
         'source':'same gated native scalar kernel and previously proved one-plane first-hit recurrence',
         'empirical_clock_trajectories':False,'finite_spacing_error_retained':True})

def run():
    start=time.time();pin=hashlib.sha256((HERE/'inputs/native.py').read_bytes()).hexdigest()
    assert pin=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('native_controller_pin',HERE/'inputs/native.py')
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    cap,tail=source_census(mod.native);embedding(mod.native);internal_spectrum(mod.native)
    passed,failures,timeouts=experiment(mod.native,cap);clock_calculation(cap)
    assert CALLS==13649696
    receipt={'all_pass':passed,'literal_native_calls':CALLS,'literal_proposals':len(CASES)*TRIALS,
       'calls_per_proposal':SLOTS,'native_reset_retry_cap':RETRIES,'native_probe_cap':T,
       'initial_preparations':{'main_units':1,'probe_units':1},'fresh_units_per_proposal':0,
       'reset_failures_retained':failures,'probe_timeouts_retained':timeouts,
       'source_pin':pin,'class':'N+','empirical_large_field':False,'physical_Newton_constant':False,
       'physical_Maxwell_transport':False,'G6':False,'wall_seconds':time.time()-start}
    dump('RECEIPT.json',receipt);print(json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=='__main__': raise SystemExit(0 if run()['all_pass'] else 1)
