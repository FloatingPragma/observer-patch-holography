"""Independent read-only checker. No producer-module imports and no evidence writes."""
import argparse,gzip,hashlib,json,math,time
from pathlib import Path
import numpy as np

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
    for r in data['rows']:
        b,m=r['B'],r['M'];g=public[f'{b}:{m}'];n=12*b;obs=r['controller']
        initial,current=census(obs['probe'],n)
        require(initial==obs['initial_core']==r['initial_audit_only'][:n],'initial current reconstruction')
        require(current==obs['post_probe_core'],'probe disturbance')
        q=list(r['initial_audit_only']);allowed={tuple(sorted(x)) for x in g['edges']}
        stage_records={x:[] for x in ('probe','balance','route')}
        for stage,rec in r['audit']['events']:
            a,c,coin,j=rec;require(tuple(sorted((a,c))) in allowed,'forbidden graph edge')
            d=q[a]-q[c];wanted=coin*d if r['control'] and stage=='balance' else (d-(2*coin-1)*(abs(d)%2))//2
            require(coin in (0,1) and j==wanted,'native current mismatch');q[a]-=j;q[c]+=j;calls+=1
            stage_records[stage].append(rec)
        require(q==r['audit']['final']==obs['observed_final'],'final state mismatch')
        require(all(stage_records[s]==obs[s] for s in stage_records),'missing stage outcomes')
        total=sum(initial)
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
            native_success+=1;source_by_prefix.setdefault(f'{b}:{m}:{sign}',q)
            for p,e,c,j,y,cc in r['consumer_global_columns']:
                edge=((0,1),(0,1),(0,2),(0,2),(1,2))[e];x=[0,0,0];x[p]=sign
                d=x[edge[0]]-x[edge[1]];wanted=(d-(2*c-1)*(abs(d)%2))//2
                x[edge[0]]-=wanted;x[edge[1]]+=wanted
                require(j==wanted and x[y]==sign and cc==c^bool(j),'global consumer column');global_columns+=1
    require(native_success==864 and zeros==12 and controls>=10,'case/control counts')
    for key,pref in data['complete_prefixes'].items():
        b,m,sign=map(int,key.split(':'));g=public[f'{b}:{m}']
        require(pref['actual_source_before']==source_by_prefix[key],'consumer not on actual first source')
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
    require(len(data['complete_prefixes'])==12,'missing signed consumer prefixes')
    require(data['timeout_control']['controller']['outcome']=='probe_timeout','timeout became success')
    return {'passed':True,'native_extractions':native_success,'zero_outcomes':zeros,'control_challenges':controls,
            'record_events_checked':calls,'global_columns':global_columns,'complete_prefixes':12,'no_producer_imports':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--evidence',type=Path,default=Path(__file__).parent/'records/EVIDENCE.json.gz');args=parser.parse_args()
    try:
        if args.evidence==Path(__file__).parent/'records/EVIDENCE.json.gz':
            receipt=json.loads((args.evidence.parent/'RECEIPT.json').read_text())
            require(receipt['passed'] is True,'saved receipt reports scientific failure')
            require(hashlib.sha256(args.evidence.read_bytes()).hexdigest()==receipt['evidence_sha256'],'saved evidence custody mismatch')
        start=time.perf_counter();result=check(args.evidence);result['seconds']=time.perf_counter()-start
        print(json.dumps(result));return 0
    except Exception as e:print(json.dumps({'passed':False,'error':repr(e)}));return 1

if __name__=='__main__':raise SystemExit(main())
