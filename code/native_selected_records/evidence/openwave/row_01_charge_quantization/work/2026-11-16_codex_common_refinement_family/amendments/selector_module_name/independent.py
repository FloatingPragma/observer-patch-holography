"""Separate standard-library/NumPy verifier. Imports no producer code."""
import argparse
from collections import deque
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np


def arr(x):
    return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)


def close(x,y,tol=2e-10):
    error=float(np.max(abs(np.asarray(x)-np.asarray(y))))
    assert error<tol,error
    return error


def pair(a,b,c,law):
    if law=='native':
        low=(a+b)//2;return (low,a+b-low) if c==0 else (a+b-low,low)
    if law=='swap':return b,a
    return (b,a) if c else (a,b)


def trace(q,tape,law):
    q=list(q)
    for a,b,c,j in tape:
        aa,bb=pair(q[a],q[b],c,law)
        assert j==q[a]-aa
        q[a],q[b]=aa,bb
    return q


def current_reader(record,n,edges):
    body,probe=record['body'],record['probe'];delta=np.zeros(n,dtype=int)
    for a,b,c,j in body:delta[a]-=j;delta[b]+=j
    # Simultaneous edge equations + an independent reference-current anchor.
    equations=[];rhs=[];change=np.zeros(n,dtype=int);start=0
    while start<len(probe):
        a,b,c,j=probe[start];correction=change[a]-change[b];end=start
        while end<len(probe):
            aa,bb,cc,jj=probe[end];assert (a,b)==(aa,bb)
            change[a]-=jj;change[b]+=jj;end+=1
            if cc!=c:
                line=np.zeros(n);line[a]=1;line[b]=-1
                equations.append(line);rhs.append(2*j+jj-correction);break
        else:raise AssertionError('probe incomplete')
        start=end
    anchor=np.zeros(n);anchor[12]=1;equations.append(anchor);rhs.append(delta[12])
    before=np.linalg.solve(np.asarray(equations),rhs)
    close(before,np.rint(before),1e-9);before=np.rint(before).astype(int)
    initial=before-delta;post=before+change
    close(initial,record['reading']['initial']);close(before,record['reading']['before_probe'])
    close(post,record['reading']['after_probe'])
    # Hidden initial audit is compared only after the observer has finished.
    close(initial,record['initial_audit_only'][:n])
    return initial.tolist(),post.tolist()


def native_ops(edges,law,vacuum=False):
    fine=[]
    for a,b in edges:
        for c in (0,1):
            outcomes=np.zeros((3,3,3))
            for x in range(3):
                q=[int(i==x) for i in range(3)]
                if vacuum:j,y=0,x
                else:
                    aa,bb=pair(q[a],q[b],c,law);j=q[a]-aa;q[a],q[b]=aa,bb;y=q.index(1)
                outcomes[j+1,y,x]=1
            fine.extend(outcomes)
    fine=np.array(fine)/math.sqrt(2*len(edges));n=len(fine)
    H=np.zeros((n,n));H[0]=1/math.sqrt(n)
    for i in range(1,n):H[i,:i]=1/math.sqrt(i*(i+1));H[i,i]=-i/math.sqrt(i*(i+1))
    close(H@H.T,np.eye(n))
    return np.einsum('ij,jab->iab',H,fine)


def contrast(edges,law):
    a=native_ops(edges,law);b=native_ops(edges,law,True)
    return np.array([(y+s*x)/2 for x,y in zip(a,b) for s in (1,-1)])


def stage(rho,ops,stored,index):
    ids=set(index if isinstance(index,list) else [index]);selected=np.zeros((ops.shape[1],ops.shape[1]),complex)
    probabilities=[]
    for offset in range(0,len(ops),256):
        chunk=ops[offset:offset+256];outputs=(chunk@rho)@np.swapaxes(chunk.conj(),1,2)
        probabilities.extend(np.trace(outputs,axis1=1,axis2=2).real.tolist())
        for j,out in enumerate(outputs):
            if offset+j in ids:selected+=out
    close(probabilities,stored);close(sum(probabilities),1)
    p=float(np.trace(selected).real)
    return selected/p,p


def run(folder,k):
    dyadic_k=k
    begin=time.perf_counter();blob=(folder/f'k{k}/EVIDENCE.json.gz').read_bytes()
    envelope=json.loads(gzip.decompress(blob));d=envelope['base_execution'];receipt=json.loads((folder/f'k{k}/RECEIPT.json').read_text())
    epsilon=2.**(-k);Q=256*2**k;stride=2**k
    assert hashlib.sha256(blob).hexdigest()==receipt['evidence_sha256']
    prereg=json.loads((folder/'PREREG.json').read_text())
    for name,digest in prereg['files'].items():
        if name=='family.py':digest=json.loads((folder/'AMENDMENT.json').read_text())['after_sha256']
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest
    graphs={len(g['rings'][0])//3:g for g in d['graphs']};readings=0;route_calls=0
    for r in d['charge_records']:
        g=graphs[r['M']];n=g['core_N'];initial,post=current_reader(r,n,g['core']);readings+=1
        full=r['initial_audit_only'];after=trace(full,r['body'],r['law']);after=trace(after,r['probe'],'native')
        close(after[:n],post);live=initial.copy();first=None
        for t,record in enumerate(r['body']+[None]):
            if first is None and (set(live)<={0,1} or set(live)<={-1,0}):first=t
            if record:
                a,b,c,j=record;live[a]-=j;live[b]+=j
        assert first==r['reading']['first_unit_event']
        if r['outcome']=='target_basis_unit':
            assert post.count(1)==1 and set(post)<={0,1}
            address=post.index(1)
            for a,b,c,j in r['route']:
                assert address in (a,b)
                if j:address=b if j==1 else a
            final=trace(after,r['route'],'native');assert address==g['target'][1]
            assert final[address]==1 and sum(final)==3
            assert sum(final[v] for v in set(range(g['N']))-set(sum(g['rings'],[])))==1
            route_calls+=len(r['route'])
        else:assert r['outcome']=='ineligible_target_sector'
        all_edges={tuple(sorted(e)) for e in g['edges']}
        assert all(tuple(sorted((a,b))) in all_edges for a,b,c,j in r['body']+r['probe']+r['route'])
        seen={0};queue=[0]
        for v in queue:
            for a,b in all_edges:
                w=b if v==a else a if v==b else None
                if w is not None and w not in seen:seen.add(w);queue.append(w)
        assert len(seen)==g['N']
        for vertices in (g['target'],g['copy']):
            for a,b in ((0,1),(0,2),(1,2)):
                assert tuple(sorted((vertices[a],vertices[b]))) in all_edges
        for a,b in zip(g['target'],g['copy']):assert tuple(sorted((a,b))) in all_edges
        assert tuple(sorted(g['blank'])) in all_edges
        for ring in g['rings']:
            assert all(tuple(sorted((a,b))) in all_edges for a,b in zip(ring,ring[1:]+ring[:1]))
    catalog={k:arr(v) for k,v in d['catalog'].items()};tri=((0,1),(1,2),(0,2));weighted=((0,1),(0,1),(0,2),(0,2),(1,2))
    max_catalog=0.
    for name,ops in catalog.items():
        close(np.einsum('kxa,kxb->ab',ops.conj(),ops),np.eye(ops.shape[2]))
        if ':target:' in name:
            m,law,_,slot=name.split(':');m=int(m);v,k=divmod(int(slot),m)
            edges=[e for e in tri if v in e]*(m-k)+[e for e in tri if (v+1)%3 in e]*k
            max_catalog=max(max_catalog,close(ops,native_ops(edges,law)))
        elif name.endswith(':plane'):max_catalog=max(max_catalog,close(ops,contrast(tri,name.split(':')[0])))
    close(catalog['cooling'],native_ops(weighted,'native'))
    close(catalog['cooling_contrast'],contrast(weighted,'native'))
    R=np.eye(33)
    for p in range(3):
        for e,(a,b) in enumerate(weighted):
            for c in (0,1):
                q=[int(x==p) for x in range(3)];aa,bb=pair(q[a],q[b],c,'native');j=q[a]-aa
                q[a],q[b]=aa,bb;y=q.index(1);cc=c^int(j!=0)
                R[:,11*p+2*e+c]=0;R[11*y+2*e+cc,11*p+2*e+c]=1
    a=np.r_[np.full(10,epsilon/math.sqrt(10)),math.sqrt(1-epsilon**2)];A=np.kron(np.eye(3),a[:,None])
    W=R@np.kron(np.eye(3),2*np.outer(a,a)-np.eye(11));close(W,arr(d['clock']['W']))
    close(R@R,np.eye(33));close(W.T@W,np.eye(33));close(catalog['decode'][0],A.T)
    waves=[];live=A.copy()
    for _ in range(Q):waves.append(live);live=W@live
    close(np.fft.fft(waves,axis=0)/Q,catalog['qpe'])
    rho=np.diag([0.,1.,0.]);prep_log=0.
    for record in d['preparation_records']:
        close(rho,arr(record['input_density']));close(prep_log,record['log10_reach_before'])
        rho,p=stage(rho,catalog[record['catalog']],record['probabilities'],record['accepted']);prep_log+=math.log10(p)
    close(rho,arr(d['decoded_position_density']))
    perm=np.eye(3)[[1,2,0]];spin_input=perm@rho@perm.T;event_count=0;marker_histories=0;max_frequency=0.
    for case in d['cases']:
        m=case['M'];law=case['law'];g=graphs[m];rho=np.kron(np.ones((2,2))/2,spin_input);logp=0.;markers={'rot':2*m,'ref':2*m}
        for rec in case['records']:
            close(logp,rec['log10_reach_before'],1e-9)
            if rec['kind']=='marker':
                arm=rec['arm'];old=markers[arm];edge=rec['edge'];ring=g['rings'][0 if arm=='rot' else 1]
                assert [ring[x] for x in edge]==rec['physical_edge']
                assert (edge[1]-edge[0]) in (1,3*m-1)
                probability=0.;success=0.;decoded=set()
                for h in rec['histories']:
                    probability+=h['probability'];close(h['probability'],2.**(-len(h['coins'])))
                    q=[int(x==old) for x in range(3*m)]
                    tape=[[edge[0],edge[1],c,j] for c,j in zip(h['coins'],h['currents'])]
                    out=trace(q,tape,'native');assert out.index(1)==h['decoded']
                    if h['outcome']=='hop':success+=h['probability'];decoded.add(h['decoded'])
                    else:assert len(h['coins'])==8 and h['decoded']==old
                    marker_histories+=1
                close(probability,1);assert len(decoded)==1;markers[arm]=decoded.pop()
                step=rec['step'];k=step%m
                desired=(2*m+case['direction']*(step if arm=='rot' else min(k,m-k)))%(3*m)
                assert markers[arm]==desired;logp+=math.log10(success);continue
            close(rho,arr(rec['input_density']))
            if rec['kind']=='target':
                assert rec['rot']==markers['rot'] and rec['ref']==markers['ref']
                rops=catalog[f'{m}:{law}:target:{rec["rot"]}'];bops=catalog[f'{m}:{law}:target:{rec["ref"]}']
            else:rops=bops=catalog[f'{law}:plane']
            ops=np.zeros((len(rops),6,6),complex);ops[:,:3,:3]=rops;ops[:,3:,3:]=bops
            rho,p=stage(rho,ops,rec['probabilities'],rec['accepted']);logp+=math.log10(p);event_count+=1
        close(rho,arr(case['spin_final_density']))
        change=np.kron(np.eye(2),perm.T);rho=change@rho@change.T;close(rho,arr(case['reencode_input']))
        ops=np.zeros((Q,66,6),complex);ops[:,:33,:3]=catalog['qpe'];ops[:,33:,3:]=catalog['qpe']
        rho,p=stage(rho,ops,case['all_reencode_bin_probabilities'],list(range(1,Q//2)))
        close(rho,arr(case['final_density']));logp+=math.log10(p)
        close(logp+prep_log,case['log10_total_with_preparation'],1e-9)
        blocks=rho.reshape(2,33,2,33);path=[];states=[]
        for sign in (1,-1):
            state=(blocks[0,:,0,:]+blocks[1,:,1,:]+sign*(blocks[0,:,1,:]+blocks[1,:,0,:]))/2
            probability=np.trace(state).real;path.append(probability);states.append(state/probability if probability>1e-12 else None)
        close(path,case['path_probabilities']);close(2*path[0]-1,case['visibility'])
        arc=states[case['clock_path_flag']];close(arc,arr(case['clock_density']))
        power=np.eye(33,dtype=complex);step=np.linalg.matrix_power(W,stride)
        z=[]
        for record in case['ramsey_records']:
            moment=np.trace(arc@power);close(moment,2*record['x_plus']-1+1j*(2*record['y_plus']-1))
            z.append(moment);power=step@power
        # Stable interval uses only the first binary record plus fixed eta.
        eta=.002;real=2*case['ramsey_records'][1]['x_plus']-1
        low=math.acos(min(1.,max(-1.,(real+eta)/(1-eta))))
        high=math.acos(min(1.,max(-1.,(real-eta)/(1-eta))))
        close([low,high],[case['stable_reader']['lower'],case['stable_reader']['upper']])
        assert low<=d['expected_frequency']<=high and high-low<.01
        f=np.array([0.,1.,-1.])/math.sqrt(2);lam=1-.4*epsilon**2;af=A@f
        jf=(R@af-lam*af)/math.sqrt(1-lam*lam);plus=(af-1j*jf)/math.sqrt(2)
        positive=float(np.vdot(plus,arc@plus).real);negative=float(np.vdot(plus.conj(),arc@plus.conj()).real)
        off=1-positive-negative;assert -1e-10<=off<=.002+1e-10
        close(off,case['audited_off_compact_band_weight'],1e-9)
        pair_record=next(p for p in envelope['actual_returned_pair_records'] if all(p[key]==case[key] for key in ('M','turns','direction','law')))
        word=np.linalg.matrix_power(W,stride);joint=np.kron(arc,arc)
        pair_z=complex(np.einsum('ab,ba->',joint,np.kron(word,word)))
        close(pair_z,complex(*pair_record['pair_ramsey']))
        eta_pair=.01;pair_low=math.acos(min(1.,max(-1.,(pair_z.real+eta_pair)/(1-eta_pair))))
        pair_high=math.acos(min(1.,max(-1.,(pair_z.real-eta_pair)/(1-eta_pair))))
        assert pair_low<=2*d['expected_frequency']<=pair_high and pair_high-pair_low<.03
        pair_off=1-positive**2-negative**2;assert -1e-10<=pair_off<=eta_pair+1e-10
        close(pair_off,pair_record['audited_pair_off_band_weight'],1e-9)
        max_frequency=max(max_frequency,high-low)
    # Native work word carries the whole choice/reference tensor, not a scalar label.
    work=d['transport'];mapping=[]
    for p0,tr in enumerate(work['traces']):
        p,c=p0,0
        assert [row[0] for row in tr]==work['word']
        for k,pp,cc,j,y,newc in tr:
            assert (pp,cc)==(p,c)
            column=work['columns'][12*k+2*p+c]
            assert column[1:]==[p,c,j,y,newc]
            a,b=column[0];q=[int(x==p) for x in range(6)];aa,bb=pair(q[a],q[b],c,'native')
            assert j==q[a]-aa;q[a],q[b]=aa,bb;assert q.index(1)==y
            p,c=y,newc
        assert p==(p0+3)%6 and c==0;mapping.append([p,c])
    close(mapping,work['final_mapping'])
    # Reader sensitivity controls are synthetic; never scientific law alternatives.
    plants=[-.6,0.,.6];detected=[2*((1+v)/2)-1 for v in plants];close(plants,detected)
    # A representable higher optical mode must survive the same re-encoder as a
    # DIFFERENT frequency. This detects a reader/re-encoder that always emits f.
    higher=np.array([-2.,1.,1.])/math.sqrt(6);source=np.outer(higher,higher)
    wave=catalog['qpe'][1:Q//2]@higher;higher_rho=wave.T@wave.conj();higher_rho/=np.trace(higher_rho)
    theta=2*math.asin(epsilon*math.sqrt(.3));moment=np.trace(higher_rho@np.linalg.matrix_power(W,stride))
    # Recover two signed conjugate poles using three observed complex moments.
    moments=np.array([np.trace(higher_rho@np.linalg.matrix_power(W,stride*t)) for t in range(4)])
    coefficient=np.linalg.solve(np.array([[moments[0],moments[1]],[moments[1],moments[2]]]),-moments[2:4])
    poles=np.roots(np.r_[1.,coefficient[::-1]]);higher_frequency=max(np.angle(poles))
    close(higher_frequency,stride*theta,1e-8)
    assert abs(higher_frequency-d['expected_frequency'])>.1
    first=arr(d['cases'][0]['final_density']);first[:33,33:]=0;first[33:,:33]=0
    assert abs(2*np.trace(first[:33,33:]).real)<1e-15
    return {'k':dyadic_k,'same_group':True,'external_G6':False,'evidence_sha256':hashlib.sha256(blob).hexdigest(),
        'current_only_reconstructions':readings,'route_native_calls':route_calls,
        'complete_catalog_count':len(catalog),'catalog_rebuild_error':max_catalog,
        'rotation_instrument_stages':event_count,'retained_marker_histories':marker_histories,
        'same_object_cases':len(d['cases']),'maximum_stable_interval_width':max_frequency,
        'higher_mode_reader_frequency':float(higher_frequency),'synthetic_reader_plants':detected,'committed_path_control_visibility':0.,
        'passed':True,'seconds':time.perf_counter()-begin}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--k',type=int,required=True)
    parser.add_argument('folder',type=Path,nargs='?',default=Path(__file__).resolve().parent)
    args=parser.parse_args();print(json.dumps(run(args.folder,args.k),indent=2,sort_keys=True))
