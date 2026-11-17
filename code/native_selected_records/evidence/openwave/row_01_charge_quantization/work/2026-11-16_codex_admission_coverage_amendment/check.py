"""Coverage-amended read-only checker; no producer imports or evidence writes.

Domain identity is derived from census currents and observed public coin tapes.
Seed/total/filename/order audit labels are never used to identify a case.
"""
import argparse,gzip,hashlib,json,math,time
import random
from collections import Counter
from pathlib import Path
import numpy as np
OLD=Path(__file__).parent.parent/'2026-11-16_codex_record_only_admission_v2'

def require(x,why):
    if not x:raise ValueError(why)

def array(x):return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)
def near(a,b,why):require(np.max(abs(np.asarray(a)-np.asarray(b)))<1e-10,why)

def census(probe,n):
    relative=[None]*n;relative[0]=0;change=[0]*n;at=0;edges=0
    while at<len(probe):
        a,b,c,j=probe[at];correction=change[a]-change[b];first=(c,j)
        require(relative[a] is not None and relative[b] is None,'census tree ordering')
        change[a]-=j;change[b]+=j;at+=1
        while True:
            require(at<len(probe),'truncated census');aa,bb,c,j=probe[at];at+=1
            require((aa,bb)==(a,b),'census seam changed')
            change[a]-=j;change[b]+=j
            if c!=first[0]:break
        relative[b]=relative[a]-(2*first[1]+j-correction);edges+=1
    require(edges==n-1 and all(x is not None for x in relative),'incomplete census')
    initial=[x-relative[12] for x in relative]
    return initial,[x+d for x,d in zip(initial,change)]

def check(path):
    data=json.loads(gzip.decompress(path.read_bytes()));require(data['passed'] is True,'saved scientific verdict failed')
    public=json.loads((Path(__file__).parent/'PUBLIC_GRAPHS.json').read_text())
    require(data['rows_count']==1752 and len(data['rows'])==1752,'missing input cases')
    calls=0;native_success=0;zeros=0;controls=0;global_columns=0;source_by_prefix={}
    # Audit-only provenance check. Both public streams are used on EVERY input;
    # identifying a stream from its observed coins therefore reveals no total.
    length=max(len(r['audit']['events']) for r in data['rows'])
    streams=[]
    for public_seed in (72361000,72362000):
        rng=random.Random(public_seed);streams.append(tuple(rng.randrange(2) for _ in range(length)))
    domain=Counter()
    expected_domain=Counter({(b,m,q,s,c):1 for b in (2,3,4) for m in (4,8)
                            for q in range(-12*b,12*b+1) for s in (0,1) for c in (False,True)})
    expected_columns={(p,e,c) for p in range(3) for e in range(5) for c in (0,1)}
    for r in data['rows']:
        b,m=r['B'],r['M'];require(b in (2,3,4) and m in (4,8) and type(r['control']) is bool,'undeclared public setting')
        g=public[f'{b}:{m}'];n=12*b;obs=r['controller']
        initial,current=census(obs['probe'],n)
        require(initial==obs['initial_core']==r['initial_audit_only'][:n],'initial current reconstruction')
        require(current==obs['post_probe_core'],'probe disturbance')
        total=sum(initial);require(-12*b<=total<=12*b,'decoded input outside registered domain')
        sign=1 if total>=0 else -1
        declared=[sign*min(b,max(0,abs(total)-b*i)) for i in range(12)]
        empty=[i for i,x in enumerate(declared) if not x]
        for a,c in zip(empty[::2],empty[1::2]):declared[a],declared[c]=b,-b
        require(initial==declared+[0]*(n-12),'decoded input differs from registered preparation')
        observed_coins=tuple(rec[2] for _,rec in r['audit']['events'])
        matches=[i for i,stream in enumerate(streams) if observed_coins==stream[:len(observed_coins)]]
        require(len(matches)==1,'observed coin tape does not identify one declared public stream')
        domain[b,m,total,matches[0],r['control']]+=1
        declared_global=initial+[0]*(g['N']-n)
        for ring in g['rings']:declared_global[ring[2*m]]=1
        require(r['initial_audit_only']==declared_global,'declared ancilla preparation differs')
        q=list(r['initial_audit_only']);allowed={tuple(sorted(x)) for x in g['edges']}
        stage_records={x:[] for x in ('probe','balance','route')}
        for stage,rec in r['audit']['events']:
            a,c,coin,j=rec;require(tuple(sorted((a,c))) in allowed,'forbidden graph edge')
            d=q[a]-q[c];wanted=coin*d if r['control'] and stage=='balance' else (d-(2*coin-1)*(abs(d)%2))//2
            require(coin in (0,1) and j==wanted,'native current mismatch');q[a]-=j;q[c]+=j;calls+=1
            stage_records[stage].append(rec)
        require(q==r['audit']['final']==obs['observed_final'],'final state mismatch')
        require(all(stage_records[s]==obs[s] for s in stage_records),'missing stage outcomes')
        if r['control']:
            nonunit=max(current)-min(current)>1
            require(not nonunit or obs['outcome']=='repair_failed','changed-law decrease credited');controls+=int(nonunit)
        elif total==0:
            require(obs['outcome']=='no_unit' and not any(q[:n]),'zero source became a unit');zeros+=1
        else:
            sign=1 if total>0 else -1
            require(obs['outcome']=='signed_target_unit' and obs['sign']==sign,'nonzero source not extracted')
            require([q[v] for v in g['target']]==[0,sign,0],'target count or sign wrong')
            require(sum(q[:n])==total-sign and set(q[:n])<={0,sign},'excess reservoir changed')
            require(all(after<=before-2 for before,after in obs['strict_decreases']),'potential did not decrease')
            native_success+=1;source_by_prefix.setdefault(f'{b}:{m}:{sign}',[]).append(q)
            columns=r['consumer_global_columns']
            require(len(columns)==30 and all(len(x)==6 for x in columns),'incomplete global consumer catalog')
            require({tuple(x[:3]) for x in columns}==expected_columns,'missing or duplicated global consumer key')
            for p,e,c,j,y,cc in r['consumer_global_columns']:
                edge=((0,1),(0,1),(0,2),(0,2),(1,2))[e];x=[0,0,0];x[p]=sign
                d=x[edge[0]]-x[edge[1]];wanted=(d-(2*c-1)*(abs(d)%2))//2
                x[edge[0]]-=wanted;x[edge[1]]+=wanted
                require(j==wanted and x[y]==sign and cc==c^bool(j),'global consumer column');global_columns+=1
    require(domain==expected_domain,'decoded input/public-stream multiplicity differs from registered domain')
    require(global_columns==25920,'incomplete total global consumer catalog')
    require(native_success==864 and zeros==12 and controls>=10,'case/control counts')
    for key,pref in data['complete_prefixes'].items():
        b,m,sign=map(int,key.split(':'));g=public[f'{b}:{m}']
        require(pref['actual_source_before'] in source_by_prefix[key],'consumer not on actual first source')
        R=np.eye(33);fine=[]
        for e,(a,b_) in enumerate(((0,1),(0,1),(0,2),(0,2),(1,2))):
            for c in (0,1):
                maps=np.zeros((3,3,3))
                for p in range(3):
                    x=[0,0,0];x[p]=sign;d=x[a]-x[b_];j=(d-(2*c-1)*(abs(d)%2))//2;x[a]-=j;x[b_]+=j;y=x.index(sign)
                    maps[j+1,y,p]=1;R[:,11*p+2*e+c]=0;R[11*y+2*e+(c^bool(j)),11*p+2*e+c]=1
                fine.extend(maps)
        fine=np.asarray(fine)/math.sqrt(10)
        def eraser(x):
            count=len(x);H=np.zeros((count,count));H[0,:]=1/math.sqrt(count)
            for j in range(1,count):H[j,:j]=1/math.sqrt(j*(j+1));H[j,j]=-j/math.sqrt(j*(j+1))
            return np.einsum('ij,jab->iab',H,x)
        continuation=eraser(fine);blank=np.zeros_like(fine)
        for j in range(10):blank[3*j+1]=np.eye(3)/math.sqrt(10)
        reference=eraser(blank);contrast=np.array([(x+s*y)/2 for y,x in zip(continuation,reference) for s in (1,-1)])
        a=np.r_[np.full(10,1/math.sqrt(10)),0.];A=np.kron(np.eye(3),a[:,None]);W=R@(2*A@A.T-np.eye(33));near(W,pref['W'],'native word')
        waves=[];live=A.astype(complex)
        for _ in range(256):waves.append(live.copy());live=W@live
        counter=np.fft.fft(np.asarray(waves),axis=0)/256;rho=np.diag([0.,1.,0.]);reach=1.
        ops_list=[contrast]+[continuation]*8+[counter]
        for rec,ops in zip(pref['records'],ops_list):
            near(rho,array(rec['input_density']),'actual prefix input')
            outs=ops@rho@ops.conj().transpose(0,2,1);probs=np.trace(outs,axis1=1,axis2=2).real
            near(probs,rec['all_probabilities'],'complete prefix outcomes');require(min(probs)>-1e-12 and abs(sum(probs)-1)<1e-10,'prefix probability law')
            keep=sum(outs[k] for k in rec['accepted']);p=float(np.trace(keep).real);reach*=p;rho=keep/p
        near(rho,array(pref['clock_density']),'actual finite counter density');near(reach,pref['complete_prefix_probability'],'prefix cost')
    require(set(data['complete_prefixes'])=={f'{b}:{m}:{s}' for b in (2,3,4) for m in (4,8) for s in (-1,1)},'missing signed consumer prefixes')
    require(data['timeout_control']['controller']['outcome']=='probe_timeout','timeout became success')
    return {'passed':True,'native_extractions':native_success,'zero_outcomes':zeros,'control_challenges':controls,
            'record_events_checked':calls,'global_columns':global_columns,'complete_prefixes':12,'no_producer_imports':True,
            'decoded_domain_cases':sum(domain.values()),'exact_input_multiplicities':True,'metadata_labels_used':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--evidence',type=Path,default=OLD/'records/EVIDENCE.json.gz');args=parser.parse_args()
    try:
        if args.evidence==OLD/'records/EVIDENCE.json.gz':
            receipt=json.loads((args.evidence.parent/'RECEIPT.json').read_text())
            require(receipt['passed'] is True,'saved receipt reports scientific failure')
            require(hashlib.sha256(args.evidence.read_bytes()).hexdigest()==receipt['evidence_sha256'],'saved evidence custody mismatch')
        start=time.perf_counter();result=check(args.evidence);result['seconds']=time.perf_counter()-start
        print(json.dumps(result));return 0
    except Exception as e:print(json.dumps({'passed':False,'error':repr(e)}));return 1

if __name__=='__main__':raise SystemExit(main())
