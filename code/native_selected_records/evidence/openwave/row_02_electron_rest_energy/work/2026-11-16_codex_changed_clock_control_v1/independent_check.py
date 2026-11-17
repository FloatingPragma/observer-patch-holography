"""Verification-only independent current/full-W audit; no producer imports."""
import argparse,gzip,json,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def require(ok,why):
    if not ok:raise ValueError(why)

def array(x):return np.asarray(x['real'])+1j*np.asarray(x['imag'])

def step(q,r,law):
    a,b,c,j=r;require(type(c) is int and c in (0,1),'actual coin domain');new=q[b] if law=='swap' else (q[a]+q[b]+c)//2
    require(q[a]-new==j,'literal source current');q[a]=new;q[b]+=j

def verify(data):
    start=time.monotonic();require(data['passed'],'current verdict');g=data['source_join']['graph'];n=g['N'];counts={'events':0,'joint_columns':0,'counter_bins':0,'returned_states':0,'ramsey_laws':0};worst=0.
    def close(x,y,label,tol=2e-9):
        nonlocal worst
        error=float(np.max(abs(np.asarray(x)-np.asarray(y))));worst=max(worst,error);require(error<tol,label)
    rows=data['source_join']['pairs'];require(len(rows)==32,'pair domain')
    require(len({(r['clock_law'],r['pattern'],r['public_stream_pair'],tuple(r['signs'])) for r in rows})==32,'unique pair grid')
    for row in rows:
        for side,s in enumerate(row['sources']):
            q=list(s['initial_audit_only']);require(s['outcome']=='ready','continuing outcome')
            for stage,r in s['events']:step(q,r,row['clock_law'] if stage in ('clock_probe','clock_blank') else 'native');counts['events']+=1
            require(q==s['final']==s['observer']['admission']['observed_final'],'continuing final/background')
            sign=s['observer']['admission']['sign'];require(sign==row['signs'][side],'observed sign')
            for label in ('prep','clock'):
                observed=s['observer'][label];position=1;found={};tape=[]
                for e in observed['events']:
                    r=e['record'];require(e['before_position']==position,'record-only position');local=[sign*int(i==position) for i in range(3)]
                    i,j=map(g['target'].index,r[:2]);local[i]-=r[3];local[j]+=r[3]
                    require(local.count(sign)==1 and local.count(0)==2,'signed unit invariant');position=local.index(sign);require(position==e['after_position'],'decoded position');tape.append(r)
                    if e['stage']=='calibration_probe':
                        key=(e['before_position'],i,j,r[2]);value={'input_position':e['before_position'],'edge':[i,j],'coin':r[2],'current':r[3],'decoded_output':position}
                        require(key not in found or found[key]==value,'conditional stationarity');found[key]=value
                require([r for st,r in s['events'] if st in (label+'_route',label+'_probe')]==tape and position==1 and len(found)==18,'same-particle table records/return')
                require(sorted(found.values(),key=str)==sorted(observed['columns'],key=str),'actual18 column table')
                blank=[r for st,r in s['events'] if st==label+'_blank'];require(blank==s['observer'][label+'_blank']['records'] and {r[2] for r in blank}=={0,1} and all(r[3]==0 for r in blank),'actual blank records')
        base=row['sources'][0]['final']+row['sources'][1]['final'];require(len(row['joint_columns'])==2178,'joint column count');seen=set()
        for rec in row['joint_columns']:
            x,y=rec['input'];key=(x,y,tuple(rec['order']));require(key not in seen,'duplicate joint column');seen.add(key);q=list(base)
            for side,col in enumerate((x,y)):q[side*n+g['target'][1]]-=row['signs'][side];q[side*n+g['target'][col//11]]+=row['signs'][side]
            for event in rec['records']:step(q,event,row['clock_law'])
            for side,col in enumerate(rec['output']):q[side*n+g['target'][col//11]]-=row['signs'][side];q[side*n+g['target'][1]]+=row['signs'][side]
            require(q==base,'full joint spectator state');counts['joint_columns']+=1
        require(seen=={(x,y,o) for x in range(33) for y in range(33) for o in ((0,1),(1,0))},'complete joint domain')
    clocks={}
    for key,record in data['clocks'].items():
        law,sign,k=key.split(':');sign=int(sign);k=int(k);c=record['clock'];cols=data['source_join']['representatives'][f'{law}:{sign}']['clock_columns'];tab={(r['input_position'],tuple(r['edge']),r['coin']):r for r in cols};edges=((0,1),(0,1),(0,2),(0,2),(1,2));R=np.zeros((33,33))
        for p in range(3):
            for h in range(11):
                if h==10:out=11*p+10
                else:
                    r=tab[p,edges[h//2],h%2];out=11*r['decoded_output']+(h^int(r['current']!=0))
                R[out,11*p+h]=1
        eps=2.**-k;a=np.r_[np.full(10,eps/np.sqrt(10)),np.sqrt(1-eps*eps)];A=np.kron(np.eye(3),a[:,None]);W=R@(2*A@A.T-np.eye(33))
        close(R,c['R'],'observed R');close(A,c['A'],'declared A');close(W,c['W'],'actual full W');close(A.T@R@A,c['P'],'source transfer');q=256*2**k;live=A.astype(complex);words=[]
        for _ in range(q):words.append(live.copy());live=W@live
        ops=np.fft.fft(np.asarray(words),axis=0)/q;native=array(data['native_preparation'][str(sign)]['final_factor']);waves=ops@native;prob=np.sum(abs(waves)**2,axis=(1,2));entry=record['input'];close(prob,entry['counter_probabilities'],'full initial counter bins');counts['counter_bins']+=q
        positive=waves[1:q//2];success=sum(prob[1:q//2]);v=a.copy();v[0]-=1;eraser=np.eye(11)-2*np.outer(v,v)/(v@v);decode=np.asarray([np.kron(np.eye(3),eraser[j:j+1]) for j in range(11)]);selected=positive.transpose(1,0,2).reshape(33,-1)/np.sqrt(success);decoded=decode@selected
        close(np.sum(abs(decoded)**2,axis=(1,2)),entry['decode_probabilities'],'all11 decoder flags');clocks[key]=(W,ops,q)
    require(len(data['cases'])==16*len(data['k_values']),'quantum case count');seen=set()
    for row in data['cases']:
        key=(row['clock_law'],row['charge_sign'],row['k'],row['turns'],row['direction']);require(key not in seen,'unique quantum cases');seen.add(key);require(row['preparation_law']==row['frame_law']=='native','fixed preparation/frame')
        W,ops,q=clocks[f'{key[0]}:{key[1]}:{key[2]}'];factor=array(row['rotation_output_factor']);E=np.asarray([[0,-2/np.sqrt(6)],[1/np.sqrt(2),1/np.sqrt(6)],[-1/np.sqrt(2),1/np.sqrt(6)]])
        left=ops@(E@factor[:2]);right=ops@(E@factor[2:]);allp=np.sum(abs(left)**2+abs(right)**2,axis=(1,2));close(allp,row['returned']['positive_reencode_probabilities'],'full returned counter law');counts['counter_bins']+=q
        for index,branch in enumerate(row['returned']['conditional_clocks']):
            v=(left[1:q//2]+(1 if index==0 else -1)*right[1:q//2])/np.sqrt(2);rho=np.einsum('kxa,kya->xy',v,v.conj());rho/=np.trace(rho);close(rho,array(branch['arc_density']),'actual returned full state');counts['returned_states']+=1
            for stride in (1,8):
                U=np.linalg.matrix_power(W,stride*2**key[2]);power=np.eye(33,dtype=complex)
                for obs in branch['ramsey_stride'+str(stride)]:
                    z=np.trace(power@rho);close([(1+z.real)/2,(1+z.imag)/2],[obs['x_plus'],obs['y_plus']],'all Ramsey binary laws');counts['ramsey_laws']+=1;power=power@U
                if index==(1 if row['turns']%2 else 0):
                    z=np.trace(U@rho)**2;rec=row['independent_pair']['records'][str(stride)];close([(1+z.real)/2,(1-z.real)/2,(1+z.imag)/2,(1-z.imag)/2],[rec[n] for n in ('x_plus','x_minus','y_plus','y_minus')],'actual product pair record')
    return {'passed':True,'verification_only':True,'independent_producer_imports':False,'counts':counts,'maximum_float_difference':worst,'seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path,nargs='?',default=HERE/'EVIDENCE.json.gz');args=p.parse_args();print(json.dumps(verify(json.loads(gzip.decompress(args.source.read_bytes()))),indent=2))
