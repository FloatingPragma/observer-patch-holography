"""Independent GR7 raw-record reader; imports only frozen native endpoint."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,random,struct,sys,time
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
DIRECTIONS=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
PREPARATIONS=((1,0,0),(1,1,0),(1,1,1))
NCALLS=0
def load(name):return json.loads((ROOT/name).read_text())
def emit(name,obj):(ROOT/name).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def verify():
    global NCALLS
    began=time.time()
    assert hashlib.sha256((ROOT/'inputs/native.py').read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('independent_gr7_native_reader',ROOT/'inputs/native.py')
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    def apply(values,seam,coin):
        global NCALLS
        NCALLS+=1
        return mod.native(values,seam,coin)
    fine=load('HIT_CENSUS.json');per_start={}
    for pos,dirn,coin,j,end,terminal in fine['rows']:
        actual,(rc,rj)=apply((1,0),(0,1),coin)
        assert (rc,rj)==(coin,j) and actual==(1-j,j)
        dest=tuple(pos[k]+j*DIRECTIONS[dirn][k] for k in range(3))
        status=3 if dest==(0,0,0) else (4 if max(map(abs,dest))==2 else sum(bool(x) for x in dest)-1)
        assert list(dest)==end and status==terminal
        key=tuple(pos);per_start.setdefault(key,[0]*5)[status]+=1
    assert len(per_start)==26
    by_class={}
    for pos,count in per_start.items():
        k=sum(bool(x) for x in pos)-1
        if k in by_class:assert count==by_class[k]
        by_class[k]=count
    Q=[[Fraction(by_class[i][j],12) for j in range(3)] for i in range(3)]
    escape=[Fraction(by_class[i][4],12) for i in range(3)]
    # Forward finite-absorption recursion, independent of producer's Q^T r.
    cap=[Fraction(0)]*3;survive=[Fraction(1)]*3
    for _ in range(64):
        cap=[escape[i]+sum(Q[i][j]*cap[j] for j in range(3)) for i in range(3)]
        survive=[sum(Q[i][j]*survive[j] for j in range(3)) for i in range(3)]
    assert [str(x) for x in cap]==fine['r_cap'] and [str(x) for x in survive]==fine['tail']
    for a,b,site,coin,j,pair in load('EMBEDDING.json'):
        initial=tuple(int(i==site) for i in range(12))
        out,r=apply(initial,(a,b),coin);two,r2=apply((initial[a],initial[b]),(0,1),coin)
        expected=list(initial);expected[a],expected[b]=two
        assert out==tuple(expected) and r==r2==(coin,j) and list(two)==pair
    spectral=load('INTERNAL_SPECTRUM.json');counts=[[0]*12 for _ in range(12)]
    for x,ei,coin,j,y in spectral['branch_rows']:
        out,r=apply(tuple(int(i==x) for i in range(12)),spectral['edges'][ei],coin)
        assert r==(coin,j) and out.index(1)==y;counts[x][y]+=1
    assert counts==spectral['counts']
    # Separate exact elimination uses source-derived transition counts.
    hit=[1-cap[i]-survive[i] for i in range(3)];b=[(1+h)/2 for h in hit]
    population=load('POPULATION.json')
    assert population['hit_cap']==list(map(str,hit)) and population['plus_cap']==list(map(str,b))
    C=[[Fraction(0)]*3 for _ in range(3)];rhs=[Fraction(0)]*3
    for i in range(3):
        for j,count in enumerate(by_class[i]):
            if i==j:continue
            bj=b[j] if j<3 else Fraction(1) if j==3 else Fraction(1,2)
            weight=Fraction(count,12)*b[i]*bj
            C[i][i]+=weight
            if j<3:C[i][j]-=weight
            elif j==3:rhs[i]+=weight
    augmented=[C[i]+[rhs[i]] for i in range(3)]
    for j in range(3):
        pivot=augmented[j][j]
        augmented[j]=[x/pivot for x in augmented[j]]
        for i in range(3):
            if i!=j:
                factor=augmented[i][j]
                augmented[i]=[augmented[i][k]-factor*augmented[j][k] for k in range(4)]
    H=[row[3] for row in augmented]
    assert population['weighted_hit_infinite']==list(map(str,H))
    for i in range(3):
        for j in range(3):assert Fraction(population['weighted_Q'][i][j])==int(i==j)-C[i][j]
        assert Fraction(population['weighted_b'][i])==rhs[i]
    exact=[Fraction(10,27),Fraction(1,6),Fraction(2,23)]
    inner_bound=24*max(Fraction(x)-y for x,y in zip(('5/22','1/11','1/22'),hit))
    assert all(abs(x-y)<=inner_bound for x,y in zip(H,exact))
    tail=Fraction(47,48)**512
    assert Fraction(population['outer_cap_tail_bound'])==tail<Fraction(1,40000)
    class Reader:
        def __init__(self,data,rng):self.data=data;self.rng=rng;self.at=0
        def byte(self):
            value=self.data[self.at];self.at+=1;return value
        def event(self,point,direction):
            value=self.byte();coin=(value>>2)&1;current=(value&3)-1
            assert value>>3==direction and coin==self.rng.randrange(2)
            active=direction<6
            out,record=apply((1,0) if active else (0,0),(0,1),coin)
            assert record==(coin,current) and out==((1-current,current) if active else (0,0))
            if active:point=tuple(point[k]+current*DIRECTIONS[direction][k] for k in range(3))
            return point,coin
        def reset(self,point,wanted):
            failed=False;moves=0
            for _ in range(9):
                if failed or point==wanted:
                    for _ in range(16):point,_=self.event(point,6)
                else:
                    axis=next(k for k in range(3) if point[k]!=wanted[k])
                    direction=2*axis+int(point[axis]>wanted[axis]);origin=point
                    for _ in range(16):point,_=self.event(point,6 if point!=origin else direction)
                    if point==origin:failed=True
                    else:moves+=1
            assert failed or point==wanted
            return point,failed,moves
        def plus(self,point,wanted):
            point,flagcoin=self.event(point,6)
            point,failed,moves=self.reset(point,wanted)
            outcome='reset_failed' if failed else 'hit' if point==(0,0,0) else 'escape' if max(map(abs,point))==2 else None
            active=0
            for _ in range(64):
                direction=self.rng.randrange(6) if outcome is None else 6
                active+=outcome is None
                point,_=self.event(point,direction)
                if outcome is None:
                    if point==(0,0,0):outcome='hit'
                    elif max(map(abs,point))==2:outcome='escape'
            if outcome is None:outcome='timeout'
            return point,{'coin':flagcoin,'plus':bool(flagcoin==0 or outcome=='hit'),
                'probe':{'outcome':outcome,'miss':int(outcome=='escape'),'reset_moves':moves,'active_probe_calls':active}}
    metric_programs=(('internal_A',0,2),('internal_B',1,2),('internal_C',2,2),
      ('temporal_A',0,4),('temporal_B',1,4),('temporal_C',2,4),('spatial',0,0),('ignore_hits_A',0,0),('hold',0,-1))
    targets=[b[0]**2/4,b[1]**2/4,b[2]**2/4,b[0]**4,b[1]**4,b[2]**4,Fraction(1,16),Fraction(1,16),Fraction(0)]
    records=load('METRIC_TRIALS.json');cursor=0;point=(1,0,0);main=0;summary=[];inner_fail=0;inner_timeout=0
    with gzip.open(ROOT/'METRIC_EVENTS.bin.gz','rb') as tape:
      for ci,(label,location,power) in enumerate(metric_programs):
        rng=random.Random(202611150727+ci);yes=0;movement=0
        for trial in range(2048):
            data=tape.read(837);assert len(data)==837;reader=Reader(data,rng)
            initial=point;oldmain=main;slots=[]
            for _ in range(4):
                point,slot=reader.plus(point,PREPARATIONS[location]);slots.append(slot)
                inner_fail+=slot['probe']['outcome']=='reset_failed';inner_timeout+=slot['probe']['outcome']=='timeout'
            gate=False if power<0 else all(s['plus'] if k<power else s['coin']==0 for k,s in enumerate(slots))
            value=reader.byte();coin=(value>>2)&1;j=(value&3)-1
            assert value>>3==(7 if gate else 6) and coin==rng.randrange(2)
            out,record=apply(tuple(int(k==main) for k in range(4)),(0,1) if gate else (2,3),coin)
            assert record==(coin,j);main=out.index(1)
            assert j*(1 if oldmain==0 else -1)==abs(j) and reader.at==837
            expected={'case':label,'index':trial,'inner_start':list(initial),'inner_end':list(point),
              'main_start':oldmain,'main_end':main,'slots':slots,'gate':gate,'current':j}
            assert records[cursor]==expected;cursor+=1;yes+=gate;movement+=abs(j)
        p=targets[ci];tol=Fraction(47,1000);reset=Fraction(36,65536)
        assert abs(Fraction(yes,2048)-p)<=tol+reset
        assert abs(Fraction(movement,2048)-p/2)<=tol+reset/2
        summary.append({'case':label,'accepted':yes,'moved':movement,'target':str(p)})
        print(json.dumps(summary[-1]),flush=True)
      assert tape.read(1)==b''
    assert cursor==len(records)==18432
    assert Fraction(summary[3]['accepted']-summary[7]['accepted'],2048)>Fraction(1,50)
    for row,produced in zip(summary,load('METRIC_RESULTS.json')['cases']):
        assert all(produced[k]==v for k,v in row.items())
    outer_programs=(('weighted_A',0,1024,True),('weighted_B',1,512,True),('weighted_C',2,512,True),('ungated_A',0,512,False))
    records=load('OUTER_TRIALS.json');cursor=0;outer=(1,0,0);outer_summary=[];nsteps=0;failures=0;timeouts=0;outer_raw=0
    with gzip.open(ROOT/'OUTER_EVENTS.bin.gz','rb') as tape:
      for ci,(label,location,n,weighted) in enumerate(outer_programs):
        rng=random.Random(202611150827+ci);hits=0
        for trial in range(n):
            header=tape.read(4);assert len(header)==4
            size=struct.unpack('<I',header)[0];data=tape.read(size);assert len(data)==size
            reader=Reader(data,rng);initial_inner=point;initial_outer=outer
            outer,failed,moves=reader.reset(outer,PREPARATIONS[location])
            outcome='reset_failed' if failed else None;steps=[]
            for _ in range(512):
                if outcome is not None:break
                direction=reader.byte();assert direction==rng.randrange(6)
                neighbor=tuple(outer[k]+DIRECTIONS[direction][k] for k in range(3));slots=[]
                for desired in (outer,neighbor):
                    point,slot=reader.plus(point,desired);slots.append(slot)
                    inner_fail+=slot['probe']['outcome']=='reset_failed';inner_timeout+=slot['probe']['outcome']=='timeout'
                gate=all(s['plus'] for s in slots) if weighted else True
                value=reader.byte();coin=(value>>2)&1;j=(value&3)-1
                assert value>>3==(direction if gate else 6) and coin==rng.randrange(2)
                out,record=apply((1,0) if gate else (0,0),(0,1),coin)
                assert record==(coin,j) and out==((1-j,j) if gate else (0,0))
                if j:outer=neighbor
                steps.append({'direction':direction,'slots':slots,'gate':gate,'coin':coin,'current':j,'end':list(outer)})
                if outer==(0,0,0):outcome='hit'
                elif max(map(abs,outer))==2:outcome='escape'
            if outcome is None:outcome='timeout'
            expected={'case':label,'index':trial,'inner_start':list(initial_inner),'inner_end':list(point),
              'outer_start':list(initial_outer),'outer_end':list(outer),'reset_failed':failed,'reset_moves':moves,
              'outcome':outcome,'steps':steps}
            assert records[cursor]==expected;cursor+=1
            assert reader.at==size==144+420*len(steps)
            hits+=outcome=='hit';failures+=outcome=='reset_failed';timeouts+=outcome=='timeout'
            nsteps+=len(steps);outer_raw+=size
        target=H[location] if weighted else Fraction(5,22)
        epsilon=Fraction(65,1000) if n==1024 else Fraction(9,100)
        assert abs(Fraction(hits,n)-target)<=epsilon+Fraction(225,65536)+tail
        outer_summary.append({'case':label,'trials':n,'hits':hits,'target':str(target)})
        print(json.dumps(outer_summary[-1]),flush=True)
      assert tape.read(1)==b''
    assert cursor==len(records)==2560
    assert Fraction(outer_summary[0]['hits'],1024)-Fraction(outer_summary[3]['hits'],512)>Fraction(1,25)
    assert Fraction(outer_summary[0]['hits'],1024)>Fraction(5,22)+Fraction(1,10)
    produced=load('OUTER_RESULTS.json')
    for row,result in zip(outer_summary,produced['cases']):assert all(result[k]==v for k,v in row.items())
    assert produced['steps']==nsteps and produced['reset_failures']==failures and produced['timeouts']==timeouts
    assert NCALLS==1824+18432*837+2560*144+419*nsteps==load('RECEIPT.json')['literal_native_calls']
    result={'all_pass':True,'imports_producer':False,'native_calls':NCALLS,'metric_trials':18432,'outer_trials':2560,
       'outer_steps':nsteps,'raw_native_event_bytes':18432*837+2560*144+419*nsteps,
       'outer_bytes_including_address_headers':outer_raw,'outer_reset_failures':failures,'outer_timeouts':timeouts,
       'inner_reset_failures':inner_fail,'inner_timeouts':inner_timeout,'exact_forward_cap_recursion':True,
       'exact_reciprocal_elimination':True,'every_rng_and_current_checked':True,'metric_cases':summary,
       'outer_cases':outer_summary,'wall_seconds':time.time()-began}
    emit('VERIFICATION.json',result);print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':verify()
