"""Complete finite stationary-reference source with first-failure ledgers.

No reader receives the source-law label. All continuation settings depend only
on the declared/record-checked number of frame circuits.
"""
import math
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import native_source as old
from amplified_native import require,OPTICAL
from finite_source import weighted,instrument,star,scalar_control,current_reversal
from selected_source import contrast_dilation

PERM=np.eye(3)[[1,2,0]]
COOL_H=24
MARKER_CAP=16


def packed(array):return old.encode(array)


def check_probabilities(p):
    require(np.min(p)>=0 and abs(float(sum(p))-1)<2e-10,'Complete probability law failed')


def apply_single(factor,ops,name,reach):
    outs=np.einsum('kxa,ab->kxb',ops,factor)
    probabilities=np.sum(abs(outs)**2,axis=(1,2))
    check_probabilities(probabilities)
    p=float(probabilities[0]);require(p>0,'Required successful source flag absent')
    record={'kind':name,'input_factor':packed(factor),'probabilities':probabilities.tolist(),
        'reach_before':reach,'accepted_flags':[0],'failure_output_rule':'K_flag @ input_factor'}
    return outs[0]/math.sqrt(p),reach*p,record


def prepare(law,columns=None,sign=1,blank_evidence=None):
    raw=weighted(law,columns,sign,blank_evidence);c=instrument(contrast_dilation(raw));p=instrument(raw)
    factor=np.array([[0.],[1.],[0.]],complex);reach=1.;records=[];failed=0.
    for name,compiled in [('contrast',c)]+[('cool',p)]*COOL_H:
        for orientation in ('forward','inverse'):
            factor,reach,record=apply_single(factor,compiled[orientation],name+':'+orientation,reach)
            count=compiled['native_reversals_per_filter']//2
            record['source_forward_calls']=count*(2 if orientation=='forward'else 1)
            record['source_inverse_calls']=count*(1 if orientation=='forward'else 2)
            record['first_failure_mass']=record['reach_before']*sum(record['probabilities'][1:])
            failed+=record['first_failure_mass']
            records.append(record)
    optical=OPTICAL.T@factor
    residual=float(np.linalg.norm(factor-OPTICAL@optical))
    require(residual<1e-11,'Exact constant-sector removal identity failed')
    weights=np.sum(abs(optical)**2,axis=1)
    require(abs(failed+reach-1)<1e-9,'Preparation first-failure ledger incomplete')
    return optical,reach,{'law':law,'charge_sign':sign,'stages':records,'success':reach,
        'first_failure_mass':failed,'complete_outcome_mass':failed+reach,
        'optical_weights':weights.tolist(),'constant_identity_residual':residual,
        'final_factor':packed(factor),'contrast':c,'cool':p,
        'R_count':c['native_reversals_per_filter']+COOL_H*p['native_reversals_per_filter']}


def clock(k,columns,sign):
    edges=((0,1),(0,1),(0,2),(0,2),(1,2))
    permutation,_,records=current_reversal(edges,'native',11,columns,sign)
    reversal=np.eye(33)[permutation]
    # The declared coherent blank-bit predicate has exactly one active word.
    word_count=4**k;active=sum(word==0 for word in range(word_count))
    require(active==1,'Clock blank-word predicate')
    epsilon=math.sqrt(active/word_count)
    choice=np.r_[np.full(10,epsilon/math.sqrt(10)),math.sqrt(1-epsilon**2)]
    a=np.kron(np.eye(3),choice[:,None])
    reflection=np.kron(np.eye(3),2*np.outer(choice,choice)-np.eye(11))
    w=reversal@reflection
    v=choice.copy();v[0]-=1;eraser=np.eye(11)-2*np.outer(v,v)/(v@v)
    decode=np.asarray([np.kron(np.eye(3),eraser[j:j+1])for j in range(11)])
    q=256*2**k;live=a.astype(complex);waves=[]
    for _ in range(q):waves.append(live.copy());live=w@live
    ops=np.fft.fft(np.asarray(waves),axis=0)/q
    return {'R':reversal,'A':a,'P':a.T@reversal@a,'W':w,'a':choice,'decode':decode,
        'qpe':ops,'Q':q,'epsilon':epsilon,'columns':records,'edges':edges,'charge_sign':sign,
        'blank_predicate':{'bits':2*k,'words':word_count,'active':active}}


def counter_input(prepared,cl):
    factor=OPTICAL@prepared
    waves=np.einsum('kxa,ab->kxb',cl['qpe'],factor)
    probabilities=np.sum(abs(waves)**2,axis=(1,2));check_probabilities(probabilities)
    positive=waves[1:cl['Q']//2];p=float(sum(probabilities[1:cl['Q']//2]))
    selected=positive.transpose(1,0,2).reshape(33,-1)/math.sqrt(p)
    decoded=np.einsum('kxa,ab->kxb',cl['decode'],selected)
    decode_prob=np.sum(abs(decoded)**2,axis=(1,2));check_probabilities(decode_prob)
    q=float(decode_prob[0]);out=decoded[0]/math.sqrt(q)
    optical=OPTICAL.T@out
    residual=float(np.linalg.norm(out-OPTICAL@optical))
    require(residual<1e-11,'Counter violated exact optical invariant')
    # QR compresses a purification without selecting an eigenvector or deleting
    # an eigenvalue. Both optical components of the actual mixed input remain.
    factor=np.linalg.qr(optical.T,mode='reduced')[1].T
    return factor,p*q,{'counter_probabilities':probabilities.tolist(),
        'counter_accepted':list(range(1,cl['Q']//2)),'decode_probabilities':decode_prob.tolist(),
        'counter_success':p,'decode_success':q,'input_factor':packed(OPTICAL@prepared),
        'output_optical_factor':packed(factor),'invariant_residual':residual}


def coordinate(slot,m):
    vertices=((1.,0.),(-.5,math.sqrt(3)/2),(-.5,-math.sqrt(3)/2))
    vertex,step=divmod(slot%(3*m),m)
    return tuple((1-step/m)*vertices[vertex][i]+step/m*vertices[(vertex+1)%3][i]for i in (0,1))


def marker(m,turns,direction):
    n=31+6*m;ring=list(range(31,31+3*m));reference=list(range(31+3*m,n))
    old_slot=2*m;records=[];angles=[];reach=1.
    for step in range(1,3*m*turns+1):
        wanted=(2*m+direction*step)%(3*m)
        p,record=old.hop(old_slot,wanted,ring,n,MARKER_CAP)
        require(all(r['decoded']==wanted for r in record['histories']if r['outcome']=='hop'),'Marker current decoder failed')
        a,b=coordinate(old_slot,m),coordinate(wanted,m)
        angle=math.atan2(a[0]*b[1]-a[1]*b[0],sum(x*y for x,y in zip(a,b)))
        records.append({'step':step,'old':old_slot,'decoded_slot':wanted,'marker_only_reach_before':reach,
            'success':p,'frame_angle_increment':angle,'decoded_coordinate':list(b),
            'reference_hold_port':reference[2*m],**record})
        reach*=p;angles.append(angle);old_slot=wanted
    winding=sum(angles)/(2*math.pi)
    require(abs(winding-direction*turns)<1e-12,'Observed frame winding')
    require(max(abs(x)for x in angles)<math.pi/2,'Frame sampling aliases half turns')
    return records,reach,{'N':n,'target_energy_order':[0,24,25],'target_spin_order':[24,25,0],
        'rotating_ring':ring,'stationary_ring':reference,'winding':winding,
        'maximum_increment_degrees':math.degrees(max(abs(x)for x in angles)),
        'reference_winding':0.,'reference_native_hops':0}


def catalog(m,law,native_catalog=None,columns=None,sign=1):
    endpoint='swap'if law!='native'else'native'
    return [star(j,m,endpoint,native_catalog[j]if law=='swap_fixed_native_calibration'else None,columns,sign)for j in range(3*m)]


def paired_stage(factor,rot,ref,orientation,reach,step,repetition):
    require(rot['coordinate_system']==ref['coordinate_system']=='energy_target_order','Frame source coordinate mismatch')
    rops=rot[orientation];fops=ref[orientation]
    full=np.vstack((OPTICAL@factor[:2],OPTICAL@factor[2:]))
    rout=np.einsum('kxa,ab->kxb',rops,full[:3])
    fout=np.einsum('kxa,ab->kxb',fops,full[3:])
    probabilities=np.sum(abs(rout)**2+abs(fout)**2,axis=(1,2));check_probabilities(probabilities)
    p=float(probabilities[0]);require(p>0,'No surviving rotation prefix')
    out=np.vstack((OPTICAL.T@rout[0],OPTICAL.T@fout[0]))/math.sqrt(p)
    residual=max(np.linalg.norm(rout[0]-OPTICAL@(OPTICAL.T@rout[0])),np.linalg.norm(fout[0]-OPTICAL@(OPTICAL.T@fout[0])))
    require(residual<1e-11,'Accepted instrument violated exact optical invariant')
    nr=rot['native_reversals_per_filter']//2;nf=ref['native_reversals_per_filter']//2
    return out,reach*p,{'kind':'target','step':step,'repetition':repetition,'orientation':orientation,
        'rot_slot':rot['slot'],'ref_slot':ref['slot'],'input_factor':packed(full),
        'probabilities':probabilities.tolist(),'reach_before':reach,'accepted_flags':[0],
        'failure_output_rule':'diag(K_rot_flag,K_ref_flag) @ input_factor',
        'invariant_residual':float(residual),'joint_R_schedule':max(nr,nf),
        'rot_R':nr,'ref_R':nf,'rot_controlled_idle_R_slots':max(nr,nf)-nr,
        'ref_controlled_idle_R_slots':max(nr,nf)-nf,
        'padding_source':'Each unused arm R slot addresses admitted vacuum edge[29,30], with V,V,V† and zero current workspaces.',
        'forward_source_calls':max(nr,nf)*(2 if orientation=='forward'else 1),
        'inverse_source_calls':max(nr,nf)*(1 if orientation=='forward'else 2)}


def exact_control(catalog,m,h,turns,direction,filter_stride=1):
    scalars=tuple(scalar_control(c)for c in catalog)
    return scalar_path(scalars,m,h,turns,direction,filter_stride)


@lru_cache(maxsize=32)
def scalar_path(scalars,m,h,turns,direction,filter_stride):
    a=F(1);b=F(1)
    for step in range(1,3*m*turns+1):
        if step%filter_stride:continue
        a*=scalars[(2*m+direction*step)%(3*m)]**h
        b*=scalars[2*m]**h
    denominator=2*(a*a+b*b)
    plus=(a+b)**2/denominator;minus=(a-b)**2/denominator
    require(plus+minus==1 and plus>=minus>=0,'Exact finite scalar path law')
    return a,b,plus,minus


def returned_records(factor,cl,exact=None,control_input=None):
    left=OPTICAL@factor[:2];right=OPTICAL@factor[2:]
    x=np.einsum('kxa,ab->kxb',cl['qpe'],left)
    y=np.einsum('kxa,ab->kxb',cl['qpe'],right)
    qprobs=np.sum(abs(x)**2+abs(y)**2,axis=(1,2));check_probabilities(qprobs)
    raw_qprobs=qprobs.copy()
    if exact is not None:
        unrotated=np.einsum('kxa,ab->kxb',cl['qpe'],OPTICAL@control_input)
        qprobs=np.sum(abs(unrotated)**2,axis=(1,2));check_probabilities(qprobs)
    pos=slice(1,cl['Q']//2);p=float(sum(qprobs[pos]))
    plus=(x[pos]+y[pos])/math.sqrt(2);minus=(x[pos]-y[pos])/math.sqrt(2)
    raw=[float(np.sum(abs(z)**2)/p)for z in (plus,minus)]
    require(abs(sum(raw)-1)<1e-10,'Joint path law')
    binary=[float(exact[2]),float(exact[3])]if exact is not None else raw
    clock_records=[]
    for index,branch in enumerate((plus,minus)):
        if exact is not None:
            # Both scalar-control path branches carry this SAME actual mixed
            # input. The tiny rational difference is never subtracted in float.
            branch=np.einsum('kxa,ab->kxb',cl['qpe'][pos],OPTICAL@control_input)
        density=np.einsum('kxa,kya->xy',branch,branch.conj(),optimize=True)
        probability=float(np.trace(density).real)
        require(probability>0,'Clock branch representation vanished')
        density/=probability
        clock_records.append({'path':'plus'if index==0 else'minus','available':binary[index]>0,
            'arc_density':packed(density),'ramsey_stride1':old.ramsey(density,cl['W'],round(1/cl['epsilon'])),
            'ramsey_stride8':old.ramsey(density,cl['W'],round(8/cl['epsilon'])),
            'actual_native_W_steps_per_readout_tick':[round(1/cl['epsilon']),round(8/cl['epsilon'])],
            'scope':'Conditional nominal finite circuit; a rare control branch has no guaranteed implementation acceptance.'})
    return p,binary,raw,{'positive_reencode_probabilities':qprobs.tolist(),
        'raw_float_propagated_reencode_probabilities':raw_qprobs.tolist(),
        'reencode_success':p,'path_probabilities_plus_minus':binary,'raw_float_path_probabilities':raw,
        'conditional_clocks':clock_records,
        'all_path_arc_factors_rule':'(QPE_q @ F_rot ± QPE_q @ F_ref)/sqrt(2), retain every q and both signs'}


def run_case(input_factor,cl,m,h,turns,direction,law,compiled,prefix,filter_stride=1):
    marks,marker_success,geometry=marker(m,turns,direction)
    applied=[sum(r['frame_angle_increment']for r in marks[i:i+filter_stride])for i in range(0,len(marks),filter_stride)]
    geometry['applied_target_orientations_per_turn']=3*m//filter_stride
    geometry['maximum_applied_increment_degrees']=math.degrees(max(abs(x)for x in applied))
    require(geometry['maximum_applied_increment_degrees']<90,'Applied rotation sampling too coarse')
    factor=np.vstack((input_factor,input_factor))/math.sqrt(2)
    reach=1.;records=[];query_count=0;failed=0.
    for step,mark in enumerate(marks,1):
        # Marker timeout aborts before this target stage. The stationary marker
        # remains fixed; successful current records supply the moving address.
        mark['reach_before']=reach
        mark['first_failure_mass']=reach*sum(r['probability']for r in mark['histories']if r['outcome']=='timeout')
        failed+=mark['first_failure_mass']
        reach*=mark['success']
        slot=mark['decoded_slot']
        if step%filter_stride:
            records.append({'kind':'declared_filter_hold','step':step,'identity_on_both_arms':True,
                'target_source_calls':0,'reach_before':reach})
            continue
        for repetition in range(h):
            for orientation in ('forward','inverse'):
                factor,reach,record=paired_stage(factor,compiled[slot],compiled[2*m],orientation,reach,step,repetition)
                record['first_failure_mass']=record['reach_before']*sum(record['probabilities'][1:])
                failed+=record['first_failure_mass']
                records.append(record);query_count+=record['forward_source_calls']+record['inverse_source_calls']
    require(abs(failed+reach-1)<1e-8,'Joint marker/target first-failure ledger incomplete')
    exact=exact_control(compiled,m,h,turns,direction,filter_stride)if law!='native'else None
    p,binary,raw,returned=returned_records(factor,cl,exact,input_factor)
    visibility=binary[0]-binary[1]
    selected=1 if turns%2 else 0  # Same fixed analyzer on every source law.
    probability=binary[selected]
    raw_reach=reach
    if exact is not None:
        a,b,_,_=exact
        reach=marker_success*float((a*a+b*b)/2)
    total=prefix*reach*p*probability
    # The old clock choice selectors are retained, including their timeout.
    selector=F(1)-F(3,8)**8
    total*=float(selector)**2
    exact_record=None
    if exact is not None:
        a,b,plus,minus=exact
        # Huge rational products remain checkable without printing a million
        # decimal digits; individual rational stage scalars define them exactly.
        exact_record={'definition':'a=product moving scalar^h on filter frames; b=stationary scalar^(3Mturns h/filter_stride)',
            'scalars':[str(scalar_control(c))for c in compiled],
            'plus':float(plus),'minus':float(minus),'minus_exact_zero':minus==0,
            'a_log10':(math.log(a.numerator)-math.log(a.denominator))/math.log(10),
            'b_log10':(math.log(b.numerator)-math.log(b.denominator))/math.log(10)}
    return {'M':m,'h':h,'filter_stride':filter_stride,'charge_sign':cl['charge_sign'],
        'k':int(round(-math.log2(cl['epsilon']))),'turns':turns,'direction':direction,
        'law':law,'geometry':geometry,'marker_records':marks,'target_first_failures':records,
        'rotation_marker_success':reach,'raw_float_rotation_marker_success':raw_reach,
        'marker_target_first_failure_mass':failed,'raw_complete_marker_target_mass':failed+raw_reach,
        'rotation_input_factor':packed(input_factor),
        'rotation_output_factor':packed(factor),'returned':returned,'visibility':visibility,
        'selected_path':'minus'if selected else'plus','selected_clock_probability_given_entry':total,
        'selected_clock_available':probability>0,'selector_success_each':str(selector),
        'exact_scalar_control':exact_record,'source_accesses_rotation_including_padding':query_count,
        'source_accesses_counter':6*(cl['Q']-1),
        'ideal_compiler_versus_finite_note':'Finite rational matching phases; nonzero rare control ports are retained.'}
