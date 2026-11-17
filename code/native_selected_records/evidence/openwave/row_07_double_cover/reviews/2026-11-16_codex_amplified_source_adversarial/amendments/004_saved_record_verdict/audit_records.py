"""Independent full-state record audit. Imports reviewer helpers only.

Complete source catalogs have their own independently reconstructed audit.
This reader propagates density matrices, including the outer-arm crossblocks,
from those bound full catalogs; it never imports producer Python modules.
"""
import argparse,gzip,hashlib,json,math,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from circuits import array,need,near,reversal,EDGES

HERE=Path(__file__).resolve().parent
OPT=np.array([[0,-2/math.sqrt(6)],[1/math.sqrt(2),1/math.sqrt(6)],[-1/math.sqrt(2),1/math.sqrt(6)]])
EMBED=np.kron(np.eye(2),OPT)
COUNTS={};DEFECT=0.

def check(a,b,label,tol=2e-9):
    global DEFECT
    e=near(a,b,label,tol);DEFECT=max(DEFECT,e)

def rel(a,b,label,tol=2e-8):
    need(math.isfinite(a) and math.isfinite(b) and (a==0)==(b==0),label+' lost exact zero/positive')
    need(a==b or abs(a-b)<=tol*max(abs(a),abs(b)),label+' relative error')

def count(key,n=1):COUNTS[key]=COUNTS.get(key,0)+n

def load(path):return json.loads(gzip.decompress(path.read_bytes())) if path.suffix=='.gz' else json.loads(path.read_text())

def density(f):return f@f.conj().T

def trace(x):return np.trace(x,axis1=-2,axis2=-1).real

def scalar(compiled):
    p=[[F(v)for v in row]for row in compiled['counts']]
    const=sum(p[0]);tr=sum(p[i][i]for i in range(3))-const
    need(tr==0,'Swap source optical trace must vanish')
    u=(sum(p[i][j]*p[j][i]for i in range(3)for j in range(3))-const*const)/2
    if compiled['preboost']:u=u*(4*u-3)**2
    alpha=2*(1-F(compiled['phase_record']['real']))
    return 1-(1-u)*(1-alpha*u)**2

def prepare(rec):
    rho=np.diag([0.,1.,0.]).astype(complex);reach=1.;failed=0.;rcount=0
    need(len(rec['stages'])==50,'Preparation stages omitted')
    for n,row in enumerate(rec['stages']):
        name='contrast'if n<2 else'cool';orientation=('forward','inverse')[n%2]
        c=rec[name];ops=array(c[orientation]);check(density(array(row['input_factor'])),rho,'Preparation actual input')
        outputs=ops@rho@ops.conj().transpose(0,2,1);probs=trace(outputs)
        check(probs,row['probabilities'],'Preparation all flags');check(sum(probs),1,'Preparation complete')
        check(row['reach_before'],reach,'Preparation reaching cost')
        fail=reach*sum(probs[1:]);check(fail,row['first_failure_mass'],'Preparation first failure');failed+=fail
        nr=c['native_reversals_per_filter']//2;rcount+=nr
        need(row['source_forward_calls']==nr*(2 if n%2==0 else 1)and row['source_inverse_calls']==nr*(1 if n%2==0 else 2),'Preparation directed query counts')
        need(row['kind']==name+':'+orientation and row['accepted_flags']==[0],'Preparation program/flag')
        rho=outputs[0]/probs[0];reach*=probs[0];count('preparation_flags',len(probs))
    check(density(array(rec['final_factor'])),rho,'Prepared actual density')
    check(rec['success'],reach,'Preparation success');check(rec['first_failure_mass'],failed,'Preparation full failure')
    check(failed+reach,1,'Preparation terminal mass');need(rec['R_count']==rcount,'Preparation R count')
    optical=OPT.T@rho@OPT;check(OPT@optical@OPT.T,rho,'Preparation invariant coordinates')
    check(np.diag(optical).real,rec['optical_weights'],'Preparation mode weights')
    return optical,reach

def clock(k,sign,rec,prep):
    width=11;perm,_,cols=reversal(EDGES,width,'native',sign);r=np.eye(33)[perm]
    eps=2.**-k;a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps*eps)]
    aa=np.kron(np.eye(3),a[:,None]);w=r@np.kron(np.eye(3),2*np.outer(a,a)-np.eye(11))
    v=a.copy();v[0]-=1;hh=np.eye(11)-2*np.outer(v,v)/(v@v)
    dec=np.array([np.kron(np.eye(3),hh[j:j+1])for j in range(11)])
    c=rec['clock'];q=256*2**k
    for name,actual in [('R',r),('A',aa),('P',aa.T@r@aa),('W',w),('a',a),('decode',dec)]:check(actual,array(c[name]),'Clock '+name)
    need(c['Q']==q and c['charge_sign']==sign and c['blank_predicate']=={'bits':2*k,'words':4**k,'active':1},'Clock source/refinement declaration')
    need(len(c['columns'])==30,'Clock current columns')
    for row,(pos,ticket,edge,coin,j,out,cc)in zip(c['columns'],cols):
        need(row['input']==pos and row['ticket']==ticket and row['forward']==[edge[0],edge[1],coin,j]and row['output']==out and row['feedback_coin']==cc,'Clock literal signed reversal')
    # Independent full source orbit, restricted only by the exact optical
    # coordinate identity; both components and all counter bins remain.
    live=aa@OPT;orbit=np.empty((q,33,2),complex)
    for i in range(q):orbit[i]=live;live=w@live
    ops=np.fft.fft(orbit,axis=0)/q
    gram=ops.conj().transpose(0,2,1)@ops
    probabilities=np.einsum('kij,ji->k',gram,prep).real
    s=rec['source'];check(probabilities,s['counter_probabilities'],'Initial complete counter bins')
    need(s['counter_accepted']==list(range(1,q//2)),'Counter admitted bins')
    p=sum(probabilities[1:q//2]);selected=ops[1:q//2]
    superop=np.einsum('kia,kjb->ijab',selected,selected.conj(),optimize=True)
    arc=np.einsum('ijab,ab->ij',superop,prep)/p
    decouts=dec@arc@dec.conj().transpose(0,2,1);dprobs=trace(decouts)
    check(dprobs,s['decode_probabilities'],'Complete current-choice eraser')
    decoded=decouts[0]/dprobs[0];rho=OPT.T@decoded@OPT
    check(density(array(s['output_optical_factor'])),rho,'Actual mixed counter/decode state')
    check(density(array(s['input_factor'])),OPT@prep@OPT.T,'Counter actual prepared input')
    check(s['counter_success'],p,'Initial counter success');check(s['decode_success'],dprobs[0],'Choice success')
    count('counter_bins',q);count('decoder_flags',11)
    return {'W':w,'Q':q,'gram':gram,'T':superop,'rho':rho,'success':p*dprobs[0]}

def marker(row,step,turns,direction,reach,markerreach):
    old=(64+direction*(step-1))%96;wanted=(64+direction*step)%96;edge=sorted((old,wanted));e0,e1=edge
    need(row['step']==step and row['old']==old and row['decoded_slot']==wanted,'Marker decoded schedule')
    need(row['edge']==edge and row['physical_edge']==[31+e0,31+e1]and row['reference_hold_port']==191,'Actual marker edge/reference hold')
    need(len(row['histories'])==17,'Marker complete stopped-coin partition')
    seen=set();total=0.;good=0.
    for h in row['histories']:
        coins=h['coins'];need(1<=len(coins)<=16 and len(coins)==len(h['currents']),'Marker record shape')
        need(tuple(coins)not in seen,'Duplicate marker history');seen.add(tuple(coins))
        pos=old
        for i,(c,j)in enumerate(zip(coins,h['currents'])):
            need(c in(0,1),'Marker coin alphabet');d=int(pos==e0)-int(pos==e1);observed=(d-(2*c-1)*(abs(d)%2))//2
            need(j==observed,'Marker literal current');pos=e1 if j==1 else e0 if j==-1 else pos
            need(i==len(coins)-1 or j==0,'History continues after observed hop')
            count('marker_literal_events')
        need(h['decoded']==pos and h['outcome']==('hop'if pos==wanted else'timeout'),'Marker observer outcome')
        need(pos==wanted or len(coins)==16,'Premature timeout');need(h['probability']==2.**-len(coins),'Marker probability')
        total+=h['probability'];good+=h['probability']*(pos==wanted)
    need(total==1 and good==1-2.**-16,'Marker exhaustive probability partition')
    check(row['success'],good,'Marker success');check(row['reach_before'],reach,'Joint marker reach')
    check(row['marker_only_reach_before'],markerreach,'Marker-only reach')
    check(row['first_failure_mass'],reach*(1-good),'Marker first failure')
    return wanted,good

def ramsey(rec,rho,w,k):
    for name,stride in [('ramsey_stride1',2**k),('ramsey_stride8',8*2**k)]:
        need(len(rec[name])==17,'Ramsey observation count');unit=np.linalg.matrix_power(w,stride);power=np.eye(33,dtype=complex)
        for n,row in enumerate(rec[name]):
            z=np.trace(rho@power);need(row['time']==n,'Ramsey time');check(row['x_plus'],(1+z.real)/2,'Actual X binary law');check(row['y_plus'],(1+z.imag)/2,'Actual Y binary law');power=unit@power;count('ramsey_binary_laws',2)
    need(rec['actual_native_W_steps_per_readout_tick']==[2**k,8*2**k],'Native clock stride scaling')

def audit_case(row,cl,catalog,prefix):
    k,turns,sign,law,stride,h,direction=(row[x]for x in('k','turns','charge_sign','law','filter_stride','h','direction'))
    need(row['M']==32 and h==2*stride and stride in(1,2),'Applied frame program')
    rho=np.kron(np.ones((2,2))/2,cl['rho']);check(density(array(row['rotation_input_factor'])),cl['rho'],'No ideal ray substitution')
    reach=1.;failed=0.;markerreach=1.;query=0;index=0;records=row['target_first_failures']
    need(len(row['marker_records'])==96*turns,'Marker frame coverage')
    a=F(1);b=F(1);scalars=[scalar(c)for c in catalog]if law!='native'else None
    for step,mark in enumerate(row['marker_records'],1):
        slot,p=marker(mark,step,turns,direction,reach,markerreach);failed+=reach*(1-p);reach*=p;markerreach*=p
        if step%stride:
            record=records[index];index+=1
            need(record['kind']=='declared_filter_hold'and record['step']==step and record['identity_on_both_arms']is True and record['target_source_calls']==0,'Missing actual filter hold')
            check(record['reach_before'],reach,'Hold reach');continue
        rot,ref=catalog[slot],catalog[64]
        if scalars is not None:a*=scalars[slot]**h;b*=scalars[64]**h
        for repetition in range(h):
            for orientation in('forward','inverse'):
                record=records[index];index+=1
                need((record['step'],record['repetition'],record['orientation'],record['rot_slot'],record['ref_slot'])==(step,repetition,orientation,slot,64),'Target program addressing')
                full=EMBED@rho@EMBED.T;check(density(array(record['input_factor'])),full,'Full actual arm density/crossblock')
                ro,rf=array(rot[orientation]),array(ref[orientation]);probs=trace(ro@full[:3,:3]@ro.conj().transpose(0,2,1))+trace(rf@full[3:,3:]@rf.conj().transpose(0,2,1))
                need(len(probs)==128 and len(record['probabilities'])==128,'Complete target flags')
                check(probs,record['probabilities'],'All actual target outcomes');check(sum(probs),1,'Complete target law')
                need(record['accepted_flags']==[0]and record['failure_output_rule']=='diag(K_rot_flag,K_ref_flag) @ input_factor','Complete failure-state rule')
                check(record['reach_before'],reach,'Target prefix cost');mass=reach*sum(probs[1:]);failed+=mass;check(record['first_failure_mass'],mass,'Target first-failure mass')
                aa=OPT.T@ro[0]@OPT;bb=OPT.T@rf[0]@OPT;op=np.zeros((4,4),complex);op[:2,:2]=aa;op[2:,2:]=bb
                rho=op@rho@op.conj().T/probs[0];reach*=probs[0]
                nr=rot['native_reversals_per_filter']//2;nf=ref['native_reversals_per_filter']//2;maximum=max(nr,nf)
                need((record['rot_R'],record['ref_R'],record['joint_R_schedule'],record['rot_controlled_idle_R_slots'],record['ref_controlled_idle_R_slots'])==(nr,nf,maximum,maximum-nr,maximum-nf),'Which-arm padding budget')
                fwd=maximum*(2 if orientation=='forward'else 1);inv=maximum*(1 if orientation=='forward'else 2)
                need(record['forward_source_calls']==fwd and record['inverse_source_calls']==inv,'Directed source/inverse schedule');query+=fwd+inv;count('target_flags',len(probs));count('full_outer_crossblocks')
    need(index==len(records),'Extra or omitted target stages')
    outfactor=array(row['rotation_output_factor']);check(density(outfactor),rho,'Final mixed state and arm crossblocks')
    check(row['marker_target_first_failure_mass'],failed,'Full target/marker failure sum');check(row['raw_float_rotation_marker_success'],reach,'Raw target/marker survival');check(failed+reach,1,'Full first-failure mass',1e-8)
    returned=row['returned'];rawposition=rho[:2,:2]+rho[2:,2:]
    rawprobs=np.einsum('kij,ji->k',cl['gram'],rawposition).real
    check(rawprobs,returned['raw_float_propagated_reencode_probabilities'],'All raw returned counter bins')
    actualposition=cl['rho']if law!='native'else rawposition
    probs=np.einsum('kij,ji->k',cl['gram'],actualposition).real;p=sum(probs[1:cl['Q']//2])
    check(probs,returned['positive_reencode_probabilities'],'Authoritative actual counter bins');check(p,returned['reencode_success'],'Returned clock counter success')
    pathrho=[(rho[:2,:2]+rho[2:,2:]+s*(rho[:2,2:]+rho[2:,:2]))/2 for s in(1,-1)]
    arcs=[np.einsum('ijab,ab->ij',cl['T'],r)for r in pathrho];rawbinary=[np.trace(r).real/p for r in arcs]
    check(rawbinary,returned['raw_float_path_probabilities'],'Raw path laws')
    if law!='native':
        den=2*(a*a+b*b);exact=[(a+b)**2/den,(a-b)**2/den];binary=[float(x)for x in exact]
        ec=row['exact_scalar_control'];need(ec['scalars']==[str(s)for s in scalars]and ec['minus_exact_zero']==(exact[1]==0),'Exact actual finite scalar path')
        for i in range(2):rel(returned['path_probabilities_plus_minus'][i],binary[i],'Exact positive branch probability',1e-10)
        reach=markerreach*float((a*a+b*b)/2)
        actualarc=np.einsum('ijab,ab->ij',cl['T'],cl['rho']);arcs=[actualarc,actualarc]
    else:binary=rawbinary;check(binary,returned['path_probabilities_plus_minus'],'Native binary outcomes')
    for i,(arc,rec)in enumerate(zip(arcs,returned['conditional_clocks'])):
        actual=arc/np.trace(arc).real;check(array(rec['arc_density']),actual,'Actual returned mixed arc state')
        need(rec['available']==(binary[i]>0)and rec['path']==('plus','minus')[i],'Physical conditional clock branch')
        ramsey(rec,actual,cl['W'],k)
    check(row['rotation_marker_success'],reach,'Physical rotation reach')
    selected=turns%2;need(row['selected_path']==('plus','minus')[selected]and row['selected_clock_available']==(binary[selected]>0),'Same fixed parity reader across laws')
    selector=1-F(3,8)**8;need(F(row['selector_success_each'])==selector,'Selector timeout omitted')
    cost=prefix*reach*p*binary[selected]*float(selector)**2
    rel(row['selected_clock_probability_given_entry'],cost,'Full selected clock cost')
    lower=1-(1 if law=='native'else 2)*F(597,2**64)
    need(F(row['source_calibration_success_lower'])==lower and row['calibration_inclusive_cost_is_lower_bound']is True,'Calibration success bound')
    rel(row['selected_clock_probability_including_calibration_given_entry'],cost*float(lower),'Calibration-inclusive selected clock cost')
    need(row['source_accesses_rotation_including_padding']==query and row['source_accesses_counter']==6*(cl['Q']-1),'Complete source-query accounting')
    check(row['visibility'],binary[0]-binary[1],'Fixed interference visibility')
    count('returned_counter_bins',cl['Q']);count('cases')
    return {'sign':sign,'k':k,'stride':stride,'turns':turns,'direction':direction,'law':law,'visibility':binary[0]-binary[1],'selected_clock_cost':cost,'minus_probability':binary[1]}

def calculate(folder,ks):
    index=load(folder/'INDEX.json');names=[x['path']for x in index['files']];need(len(names)==len(set(names)),'Duplicate custody artifacts')
    for x in index['files']:need(hashlib.sha256((folder/x['path']).read_bytes()).hexdigest()==x['sha256'],'Producer custody changed')
    catalog_receipt=load(HERE/'CATALOG_AUDIT.json')
    for name,expected in catalog_receipt['input_hashes'].items():need(hashlib.sha256((folder/name).read_bytes()).hexdigest()==expected,'Previously independent catalog changed')
    prep=load(folder/'PREPARATION.json.gz');catalog=load(folder/'CATALOG.json.gz');prepared={}
    for sign in(-1,1):
        for law in('native','swap'):prepared[sign,law]=prepare(prep[str(sign)][law])
    rows=[];keys=[]
    for sign in(-1,1):
        for k in ks:
            rho,success=prepared[sign,'native'];cl=clock(k,sign,load(folder/f'INPUT_s{sign}_k{k}.json.gz'),rho)
            for stride in(1,2):
                for turns in(1,2):
                    for direction in(-1,1):
                        for law in(('native','swap','swap_fixed_native_calibration')if k in(0,6)else('native','swap')):
                            name=f'CASE_s{sign}_k{k}_r{stride}_t{turns}_d{direction}_{law}.json.gz'
                            row=load(folder/name);key=(sign,k,stride,turns,direction,law)
                            need(tuple(row[x]for x in('charge_sign','k','filter_stride','turns','direction','law'))==key,'Case source index mismatch')
                            rows.append(audit_case(row,cl,catalog[str(sign)][law],success*cl['success']));keys.append(key)
            print(json.dumps({'progress_sign':sign,'k':k,'verified_cases':len(rows)}),flush=True)
    return {'passed':True,'same_group':True,'external_G6':False,'no_producer_imports':True,'index_sha256':hashlib.sha256((folder/'INDEX.json').read_bytes()).hexdigest(),'k_grid':ks,'counts':COUNTS,'maximum_absolute_defect':DEFECT,'rows':rows}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('folder',type=Path);ap.add_argument('--k',type=int,nargs='+',default=list(range(7)));ap.add_argument('--record-to',type=Path);ap.add_argument('--verify',type=Path,default=HERE/'RECORD_AUDIT.json');args=ap.parse_args()
    if args.record_to:need(not args.record_to.exists(),'Exclusive first-record path exists')
    start=time.perf_counter();result=calculate(args.folder,args.k)
    if args.record_to:
        with args.record_to.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:
        saved=load(args.verify);need(result['counts']==saved['counts']and result['index_sha256']==saved['index_sha256']and result['k_grid']==saved['k_grid'],'Saved verdict custody/coverage')
        for a,b in zip(result['rows'],saved['rows']):
            for k in a:
                if isinstance(a[k],float):rel(a[k],b[k],'Saved result '+k)
                else:need(a[k]==b[k],'Saved result identity')
        need(len(result['rows'])==len(saved['rows']),'Saved row count')
    print(json.dumps({'passed':True,'counts':COUNTS,'maximum_defect':DEFECT,'seconds':time.perf_counter()-start}))
