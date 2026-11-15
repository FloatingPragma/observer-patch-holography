#!/usr/bin/env python3
"""Independent full adaptive scheduler/coin/current replay."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import gzip,hashlib,json,math,random
import numpy as np
HERE=Path(__file__).resolve().parent
P=HERE.parent.parent/'work/2026-11-14_codex_probe_binding'

def solve(a,b):
    m=[list(row)+[v] for row,v in zip(a,b)]
    for i in range(len(m)):
        j=next(j for j in range(i,len(m)) if m[j][i]);m[i],m[j]=m[j],m[i]
        pivot=m[i][i];m[i]=[v/pivot for v in m[i]]
        for j in range(len(m)):
            if j!=i:
                c=m[j][i];m[j]=[v-c*w for v,w in zip(m[j],m[i])]
    return [row[-1] for row in m]

def main():
    raw=gzip.decompress((P/'RECORDS.json.gz').read_bytes());r=json.loads(raw);result=json.loads((P/'RECEIPT.json').read_text())
    assert hashlib.sha256(raw).hexdigest()==result['raw_uncompressed_sha256']
    cells=list(product(range(-2,3),repeat=3));ci={v:i for i,v in enumerate(cells)}
    faces=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
           (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
    intr=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in faces for j in range(3)})
    edges=[]
    for v in cells:
        c=ci[v];edges += [(12*c+a,12*c+b) for a,b in intr]
        for axis in range(3):
            if v[axis]<2:
                w=list(v);w[axis]+=1;edges += [(12*c+p,12*ci[tuple(w)]+p) for p in range(12)]
    assert r['graph']==dict(cells=[list(v) for v in cells],edges=[list(e) for e in edges])
    adj=[[] for _ in range(1500)];lookup={}
    for ei,(a,b) in enumerate(edges):adj[a].append(ei);adj[b].append(ei);lookup[frozenset((a,b))]=ei
    branches=iter(r['branches'])
    for ei in range(len(edges)):
        for slot in range(3):
            for c in (0,1):assert next(branches)==[ei,slot,c,int(slot==0)-c if slot<2 else 0]
    assert next(branches,None) is None
    Q=[[F(8,11),F(2,11),F(0)],[F(1,11),F(8,11),F(1,11)],[F(0),F(3,22),F(8,11)]]
    hit=[F(1,22),F(0),F(0)];escape=[F(1,22),F(2,22),F(3,22)]
    # Re-derive every orbit row at every interior cell/port from the endpoint identity.
    for xyz in product(range(-1,2),repeat=3):
        if xyz==(0,0,0):continue
        k=sum(v!=0 for v in xyz)-1
        for p in range(12):
            x=12*ci[xyz]+p;row=[F(0)]*3;h=F(0);e=F(0);assert len(adj[x])==11
            for ei in adj[x]:
                a,b=edges[ei]
                for coin in (0,1):
                    j=int(x==a)-coin;y=b if j==1 else a if j==-1 else x;v=cells[y//12]
                    if v==(0,0,0):h+=F(1,22)
                    elif max(map(abs,v))==2:e+=F(1,22)
                    else:row[sum(n!=0 for n in v)-1]+=F(1,22)
            assert row==Q[k] and h==hit[k] and e==escape[k]
    imq=[[F(i==j)-Q[i][j] for j in range(3)] for i in range(3)]
    assert solve(imq,hit)==[F(5,22),F(1,11),F(1,22)]
    assert solve(imq,[F(1)]*3)==[F(31,3),F(10),F(26,3)]
    # Exact finite cap probabilities, without floating matrix powers.
    hc=[F(0)]*3;ec=[F(0)]*3;tail=[F(1)]*3
    for _ in range(256):
        hc=[hit[i]+sum(Q[i][j]*hc[j] for j in range(3)) for i in range(3)]
        ec=[escape[i]+sum(Q[i][j]*ec[j] for j in range(3)) for i in range(3)]
        tail=[sum(Q[i][j]*tail[j] for j in range(3)) for i in range(3)]
    assert all(h+e+t==1 for h,e,t in zip(hc,ec,tail))
    assert max(tail)<=F(10,11)**256
    assert max(abs(float(a)-b) for a,b in zip(hc,result['capped_hit']))<1e-14
    assert max(abs(float(a)-b) for a,b in zip(ec,result['capped_escape']))<1e-14
    pos={'i':12*ci[(0,0,0)],'o':12*ci[(0,0,0)]};rng=random.Random(2026111431)
    trials=iter(r['trials']);checked=0;stats=[]
    for target,feedback in [((0,0,0),True),((1,0,0),True),((1,1,0),True),((1,1,1),True),((1,0,0),False)]:
        hits=moves=censors=cost=0
        for _ in range(4096):
            trial=next(trials);records=iter(trial['records']);used=0
            def step(which,eid):
                nonlocal checked,used
                coin=rng.randrange(2)
                if which=='z':j=0
                else:
                    a,b=edges[eid];s=pos[which];j=int(s==a)-coin if s in (a,b) else 0
                    if j:pos[which]=b if j==1 else a
                assert next(records)==[which,eid,coin,j];checked+=1;used+=1;return j
            def route(which,zero=False):
                if zero and pos[which]%12:
                    assert pos[which]%12==1;c=pos[which]//12;eid=lookup[frozenset((12*c,12*c+1))]
                    for _ in range(128):
                        step(which,eid)
                        if pos[which]%12==0:break
                    else:return False
                for axis in range(3):
                    while cells[pos[which]//12][axis]!=target[axis]:
                        s=pos[which];v=list(cells[s//12]);v[axis]+=1 if v[axis]<target[axis] else -1
                        dest=12*ci[tuple(v)]+s%12;eid=lookup[frozenset((s,dest))]
                        for _ in range(128):
                            step(which,eid)
                            if pos[which]==dest:break
                        else:return False
                return True
            ok=route('o',True) and route('i');prep=used;probe=0;outcome='preparation_timeout'
            if ok:
                for k in range(257):
                    v=cells[pos['i']//12]
                    if v==(0,0,0):outcome='hit';break
                    if max(map(abs,v))==2:outcome='escape';break
                    if k==256:outcome='probe_timeout';break
                    opts=adj[pos['i']];step('i',opts[rng.randrange(len(opts))]);probe+=1
            held=not ok or feedback and outcome!='escape'
            if held:move=step('z',-1)
            else:
                base=12*ci[target];assert pos['o']==base;move=step('o',lookup[frozenset((base,base+1))])
            assert next(records,None) is None
            assert trial['x']==list(target) and trial['feedback']==feedback
            assert trial['preparation_events']==prep and trial['probe_events']==probe
            assert trial['outcome']==outcome and trial['held']==held
            hits+=outcome=='hit';moves+=move;censors+=outcome.endswith('timeout');cost+=used
        stats.append(dict(hit_count=hits,move_count=moves,censors=censors,cost=cost))
    assert next(trials,None) is None and checked+len(r['branches'])==result['actual_native_calls']==305157
    for s,t in zip(stats,result['summaries']):
        assert s['hit_count']/4096==t['hit_mean'] and s['move_count']/4096==t['move_mean']
        assert s['censors']==t['censors']==0 and s['cost']==t['actual_native_calls']
        assert abs(t['hit_mean']-t['hit_target'])<=math.sqrt(math.log(4000)/8192)+25/2**128
        assert abs(t['move_mean']-t['move_target'])<=math.sqrt(math.log(4000)/8192)+25/2**128
    # The exact uniform trial bound uses 1 +6h100+12h110+8h111=42/11.
    assert (1+6*F(5,22)+12*F(1,11)+8*F(1,22))/27==F(14,99)
    for w in result['finite_wells']:
        v=np.array(w['ground_profile']);assert abs(v@v-1)<1e-12
        assert 0<w['lowest_squared_energy']<F(85,99)*w['native_lambda']
    out={'pass':True,'same_group':True,'G6_pass':False,'branch_records_checked':len(r['branches']),
         'adaptive_records_checked':checked,'all_source_calls_accounted':305157,'all_adaptive_rng_choices_replayed':True,
         'all_current_only_routes_and_feedback_replayed':True,'stats':stats,'exact_orbit_probabilities_verified':True,
         'exact_probe_cap_probabilities_verified':True,'exact_uniform_well_bound':'85 lambda / 99',
         'preparation_failure_scope':'25 * 2^-128 per attempt; exact ideal capped-probe kernel is not the history-dependent fully capped marginal',
         'packet_receipt_sha256':hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest(),
         'packet_records_sha256':hashlib.sha256((P/'RECORDS.json.gz').read_bytes()).hexdigest()}
    (HERE/'VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
