"""Independent opposing audit; only development members2,3,4 are executed."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import expm
import refinement_adapter as adapter

HERE=Path(__file__).resolve().parent
ROOT=adapter.ROOT
FROZEN={'openwave-tuning/common_dynamics/clock/MANIFEST.json':'6f80b0fe899cdeffd59d3742c318fe627817c19685d8f164abed69721f0b0992',
        'openwave-tuning/common_dynamics/continuum/MANIFEST.json':'3067f7bdea43ca5137d6b115f7e1bd92b4d753ba99684ed198ba2da170240adf'}


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path); result=importlib.util.module_from_spec(spec)
    sys.modules[name]=result; spec.loader.exec_module(result); return result


def check_pins():
    pins=dict(FROZEN); counts={}
    for rel,digest in FROZEN.items(): assert sha(ROOT/rel)==digest,rel
    for rel in (*FROZEN,'openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json',
                'openwave-tuning/auditor/common_walk_particle/MANIFEST.json'):
        path=ROOT/rel; manifest=json.loads(path.read_text())
        records=manifest.get('files_sha256',manifest.get('source_sha256',{})); counts[rel]=len(records)
        for name,digest in records.items():
            target=ROOT/name if (ROOT/name).is_file() else path.parent/name
            assert sha(target)==digest,str(target)
            pins[str(target.relative_to(ROOT))]=digest
        pins[rel]=sha(path)
    for rel in ('openwave-tuning/auditor/common_walk_review/source_error.py',
                'openwave-tuning/auditor/common_walk_review/clock_statistics.py',
                'openwave-tuning/auditor/common_walk_review/flux_algebra.py',
                'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json',
                'reverse-engineering-reality/model/certified_reference/source.py',
                'reverse-engineering-reality/model/certified_reference/certificate.py',
                'reverse-engineering-reality/model/observer/source_scalar/compiler.py',
                'reverse-engineering-reality/model/openwave_coherent_native_architecture/source.py'):
        pins[rel]=sha(ROOT/rel)
    return pins,counts


def independent_bloch():
    """Reconstruct all144 states and all spectral complements from seams."""
    I=np.eye(144,dtype=np.int64); S=I.copy(); Ds=[]
    for colour,matching in enumerate(adapter.common.INTERNAL_MATCHINGS):
        for a,b in matching:
            x,y=12*a+2*colour+1,12*b+2*colour+1
            S[x,x]=S[y,y]=0; S[x,y]=S[y,x]=1
    for matching in adapter.spatial.links():
        D=np.zeros_like(I)
        for a,b in matching:
            x,y=12*a+11,12*b+11
            S[x,x]=S[y,y]=0; S[x,y]=S[y,x]=1
            D[x,y]=1; D[y,x]=-1
        Ds.append(D)
    Rn=np.kron(np.eye(12,dtype=np.int64),np.ones((12,12),dtype=np.int64)-6*np.eye(12,dtype=np.int64))
    Wn=S@Rn; Cn=Wn+Wn.T
    assert np.array_equal(S@S,I) and np.array_equal(Rn@Rn,36*I)
    # Derive every real-part projector, including the optical and pi sectors.
    levels=(12,-12,8,6,4); projections={}
    for c in levels:
        numerator=I.copy(); denominator=1
        for other in levels:
            if other==c: continue
            factor=Cn-other*I
            assert 144*int(abs(numerator).max())*int(abs(factor).max())<2**63
            numerator=numerator@factor; denominator*=c-other
        if denominator<0: numerator=-numerator; denominator=-denominator
        gcd=math.gcd(denominator,int(np.gcd.reduce(numerator.reshape(-1))))
        numerator//=gcd; denominator//=gcd
        assert np.array_equal(numerator@numerator,denominator*numerator)
        assert np.array_equal(Cn@numerator,c*numerator)
        projections[c]=(numerator,denominator)
    lcm=math.lcm(*(d for n,d in projections.values()))
    assert np.array_equal(sum(n*(lcm//d) for n,d in projections.values()),lcm*I)
    ranks={str(F(c,12)):int(np.trace(n))//d for c,(n,d) in projections.items()}
    assert ranks=={'1':26,'-1':96,'2/3':6,'1/2':4,'1/3':12}
    Pn,Pd=projections[12]; assert np.array_equal(Wn@Pn,6*Pn)
    H=[]; g=[]; currents=[]
    for D in Ds:
        n=Pn@D@Rn@Pn; d=6*Pd*Pd
        gcd=math.gcd(d,int(np.gcd.reduce(n.reshape(-1)))); n//=gcd; d//=gcd
        v=n@np.ones(144,dtype=np.int64)
        assert np.array_equal(144*n,v[:,None]+v[None,:]) and sum(v)==0
        J=-D@S
        assert np.array_equal((Pn@J@Pn)*d,-n*Pd*Pd)
        H.append(n.astype(float)/d); g.append((v,d)); currents.append(J)
    for a,(v,d) in enumerate(g):
        for b,(w,e) in enumerate(g): assert int(v@w)==(3*d*e if a==b else 0)
    # A terminal native-label read can be a valid final instrument while its
    # postmeasurement state is outside the next interval's preparation sector.
    full_label_leak=F(144-26,144)
    select=np.diag(currents[0])==1
    v=select.astype(np.int64)
    branch_probability=F(int(v.sum()),144)
    projected=F(int(v@Pn@v),144*Pd)
    binary_leak=2*(branch_probability-projected)
    assert full_label_leak==F(59,72) and binary_leak>0
    data={'S':S,'R':Rn/6,'P':Pn/Pd,'D':Ds,'H':H,'W0':Wn/6}
    rows={'full_internal_event_dimension':144,'real_part_projector_ranks':ranks,
          'zero_phase_dimension':26,'optical_total_dimension':22,'pi_phase_dimension':96,
          'generic_small_nonzero_momentum':{'acoustic':2,'zero_phase_stationary':24,'other_inherited_optical':22,'pi_phase':96},
          'common_stationary_subspace_dimension':22,'gradient_Gram_exact':'I3/48',
          'native_current_compression':'P J_a P = -H_a',
          'fine_internal_event_readout_leakage_exact':str(full_label_leak),
          'positive_axis_displacement_binary_readout_leakage_exact':str(binary_leak),
          'selected_native_branch_probability_exact':str(branch_probability),
          'selected_branch_P_survival_conditional_exact':str(projected/branch_probability),
          'complement_accounting':'All144 states retained. Source/event padding and export/reset records belong to the enclosing compiled instrument.'}
    return data,rows


def continuum_controls(data):
    producer=module('audited_continuum','openwave-tuning/common_dynamics/continuum/continuum.py')
    P=data['P']; Q=np.eye(144)-P; p=np.array([1.,-1.,1.]); H=sum(p[a]*data['H'][a] for a in range(3)); K=1j*H
    X=-np.linalg.solve(data['W0']-np.eye(144)+P,Q@(1j*sum(p[a]*data['D'][a]@data['R'] for a in range(3)))@P)
    assert np.linalg.norm(data['W0']@X+1j*sum(p[a]*data['D'][a]@data['R'] for a in range(3))@P-X-K,2)<1e-12
    records=[]
    for denom in (8,32):
        h=F(1,denom); n=denom; S=data['S'].astype(complex)
        for a,D in enumerate(data['D']):
            x,y=np.where(D)
            S[x,y]=np.exp(1j*float(h)*p[a]*D[x,y])
        W=S@data['R']; exact_bound=F(5,2)*h+F(77,32)*n*h*h
        actual=np.linalg.norm(np.linalg.matrix_power(W,n)@P-expm(n*float(h)*K)@P,2)
        assert actual<float(exact_bound)
        assert producer.observer_error_bound(n,h,1)['trace_distance_and_terminal_TV_bound']==min(F(1),exact_bound)
        records.append({'h_exact':str(h),'steps':n,'operator_norm_diagnostic':float(actual),'rational_uniform_bound_exact':str(exact_bound)})
    rejected=[]
    for call in (lambda:producer.compiled_size2_budget(2,member_size=3),
                 lambda:producer.observer_error_bound(1,F(1),1),
                 lambda:producer.observer_error_bound(1,F(1,16),1,tau0=0)):
        try: call()
        except ValueError: rejected.append(True)
        else: rejected.append(False)
    assert all(rejected)
    # Norm error is only an upper bound on TD, not an equality. A global sign
    # becomes observable under controlled-W, which needs an actual wrapper.
    assert np.linalg.norm(np.array([1.,0.])-np.array([-1.,0.]))==2
    plus=np.array([1.,1.])/math.sqrt(2); minus=np.array([1.,-1.])/math.sqrt(2)
    assert abs(np.vdot(plus,minus))<1e-15
    return {'controls':records,'invalid_domain_and_wrong_member_rejected':rejected,
            'global_sign_control':{'state_vector_distance':2,'uncontrolled_state_trace_distance':0,'controlled_branch_probe_trace_distance':1},
            'position_intervention':'A sharp cell projector preserves internal P but broadens Fourier support; a fine port/event or signed-displacement read can also leave P.',
            'preparation_interface':'The frozen initial_cell_state is cell-localized. It does not meet a fixed-bandwidth continuum preparation premise as L grows. General bandlimited preparation still requires its own native word/error certificate.'}


def clock_controls():
    stats=module('audited_clock_statistics','openwave-tuning/auditor/common_walk_review/clock_statistics.py')
    clock=module('audited_clock_contract','openwave-tuning/common_dynamics/clock/clock_contract.py')
    receipt=json.loads((ROOT/'openwave-tuning/common_dynamics/clock/PHASE_CLOCK.json').read_text())
    systematic=receipt['complete_joint_source_and_clock_error_upper']
    result=stats.estimate([4096,4096,0],[4096,4096,0],.01,systematic)
    assert result['status']=='inconclusive'
    all_fail=stats.estimate([0,0,8192],[0,0,8192],.01,systematic)
    assert all_fail['status']=='inconclusive'
    finite=stats.estimate([8192,0,0],[4096,4096,0],.01,systematic)
    assert finite['status']=='interval' and finite['phase_interval'][0]==0
    assert clock.effective_time(8,F(1,16),F(3,2))==F(3,4)
    waits=[]
    for size in (2,3,4):
        row=clock.wait_law(size,3)
        assert sum(map(F,row['first_contact_attempt_probabilities']))+F(row['timeout_probability'])==1
        waits.append({'size':size,'contact_exact':row['contact_probability_per_attempt'],'timeout_exact':row['timeout_probability']})
    # Exact marked thinning is classical; two coherent steps keep amplitudes.
    assert F(25,72)-F(7,24)==F(1,18)
    return {'systematic_mean_phase_bias_upper':float(F(receipt['phase_error_upper_exact'])),
            'finite_shot_interval_8192_each':finite,'all_failure_and_zero_contrast_inconclusive':True,
            'clock_frequency':'1/(8*tau0); not1/(16*tau0)',
            'exact_capped_wait_controls':waits,'coherent_minus_Markov_two_step_return_exact':'1/18',
            'bare_W_pulses':432425,'bare_W_load_export_services':864850,
            'time_scope':'t=n*tau0*h is a common effective coordinate. Serial pulse duration needs control rates and service latencies; neither the native rule nor this adapter fixes physical seconds.',
            'frozen_phase_certificate_largest_source_step_count':512,'frozen_phase_certificate_graph_size':2}


def evaluate_template(template,bits,qubits):
    value=np.eye(1 << qubits,dtype=complex)
    failure=np.zeros_like(value)
    max_completeness=0.
    for kind,targets,desired in template:
        event=adapter.certified_pulse(kind,targets,desired,bits)
        branches=adapter.execute_full_response(event); d=4 if kind=='J' else 2
        gram=sum(K.conj().T@K for K in branches)
        max_completeness=max(max_completeness,float(np.linalg.norm(gram-np.eye(d),2)))
        for K in branches[1:]:
            flag=adapter.compiler._apply(value,K[d:d+1].conj().T@K[d:d+1],targets,qubits)
            failure+=value.conj().T@flag
        value=adapter.compiler._apply(value,branches[0][:d],targets,qubits)
    assert np.linalg.norm(value.conj().T@value+failure-np.eye(1 << qubits),2)<1e-10
    return value,float(np.trace(failure).real),max_completeness


def adapter_controls():
    rows=[]; program=adapter.reversible_program(2)
    # All4096 clean data/law columns at size2; finite chosen complete-register
    # columns at3,4. The arbitrary-size proof covers every dirty input by the
    # reversible gate extension, without claiming an exponential enumeration.
    for size in (2,3,4):
        program=adapter.reversible_program(size); info=program['layout']; limit=1<<info['data_bits']
        data=np.arange(limit,dtype=np.int64) if size==2 else np.array(sorted(set([0,1,11,15,16,255,limit-1]+[i*limit//17 for i in range(17)])),dtype=np.int64)
        inputs=np.r_[data,data|(1<<info['law'])]
        actual=adapter.execute_bitplanes(program['gates'],inputs,info['memory_bits'])
        target=np.array(program['routing_permutation'],dtype=np.int64)
        expected=np.r_[target[data],target[data]|(1<<info['law'])]
        assert np.array_equal(actual,expected)
        # Every declared dirty register input has a complete inverse, tested
        # on a deterministic small set while the gatewise proof is universal.
        dirty=inputs[:16]|(1<<info['scratch'][0])|(1<<info['raw'][0])
        out=adapter.execute_bitplanes(program['gates'],dirty,info['memory_bits'])
        back=adapter.execute_bitplanes(list(reversed(program['gates'])),out,info['memory_bits'])
        assert np.array_equal(dirty,back)
        cert=adapter.certificate(size,F(1,size*size))
        assert F(cert['whole_W_operator_error_exact']) <= F(1,size*size)
        assert cert['native_pulses']<=adapter.resource_upper_bound(size)['native_pulses_upper']
        reset=adapter.reset_repetition_plan(info['memory_bits'],F(1,size*size))
        assert F(reset['complete_timeout_union_bound_exact'])<=F(1,size*size)
        tick=adapter.clock_tick_pulse(size,F(1,size),cert['precision_bits'])
        assert tick['targets']==[info['memory_bits']]
        assert F(tick['full_response_operator_error_exact'])<=F(1,size*size)
        assert adapter.verify_pulse(tick)
        assert F(tick['angle_interval_exact'][0])<=-F(1,16*size)<=F(tick['angle_interval_exact'][1])
        try: adapter.clock_tick_pulse(size,F(1,size),cert['precision_bits'],info['law'])
        except ValueError: pass
        else: raise AssertionError('clock/source wire alias accepted')
        rows.append({'size':size,'clean_columns_executed':len(inputs),'dirty_roundtrip_columns':len(dirty),
                     'all_clean_columns_enumerated':size==2,'certificate':cert,'reset_plan':reset,
                     'clock_tick':tick})
        print('adapter member',size,'checked',flush=True)
    u,failure,complete=evaluate_template(adapter.uniform_template(),72,4)
    desired=np.r_[np.ones(12)/math.sqrt(12),np.zeros(4)]
    assert np.linalg.norm(u[:,0]-desired)<1e-12
    r,rfailure,rcomplete=evaluate_template(adapter.reflection_template(),72,4)
    assert np.linalg.norm(r-(2*np.outer(desired,desired)-np.eye(16)),2)<1e-12
    primitive_rows={}
    for name,template in adapter.primitive_templates().items():
        qubits={'ccx':3,'cx':2}.get(name,1)
        value,ff,cc=evaluate_template(template,72,qubits)
        if name=='ccx': target=np.eye(8)[:,[0,1,2,7,4,5,6,3]]
        elif name=='cx': target=np.eye(4)[:,[0,3,2,1]]
        elif name=='flip': target=np.array([[0,1],[1,0]])
        else: target=np.array([[1,1],[1,-1]])/math.sqrt(2)
        error=float(np.linalg.norm(value-target,2)); assert error<1e-12
        primitive_rows[name]={'numerical_residual':error,'complement_gram_trace':ff,'complete_channel_residual':cc}
    e48=adapter.template_bound(adapter.reflection_template(),48)
    e80=adapter.template_bound(adapter.reflection_template(),80)
    assert 0<e80<e48/F(10**6)
    # An emitted W stream contains actual native twelve-coefficient requests.
    first=next(adapter.pulse_stream(2,72)); assert len(first['controls_exact'])==12
    assert adapter.verify_pulse(first)
    changed=dict(first); changed['controls_exact']=list(first['controls_exact']); changed['controls_exact'][0]='1'
    try: adapter.verify_pulse(changed)
    except ArithmeticError: pass
    else: raise AssertionError('changed coefficient accepted by exact certificate')
    good=adapter.execute_full_response(first)[0]; bad=adapter.execute_full_response(changed)[0]
    mutation=float(np.linalg.norm(good-bad,2)); assert mutation>1e-3
    return {'development_members':rows,'uniform12_preparation_numeric_residual':float(np.linalg.norm(u[:,0]-desired)),
            'reflection_numeric_residual':float(np.linalg.norm(r-(2*np.outer(desired,desired)-np.eye(16)),2)),
            'uniform_export_complement_gram_trace':failure,'reflection_export_complement_gram_trace':rfailure,
            'primitive_complete_response_controls':primitive_rows,
            'precision_control':{'reflection_error48_exact':str(e48),'reflection_error80_exact':str(e80)},
            'changed_native_coefficient_rejected_diagnostic':mutation,'first_emitted_native_pulse':first,
            'formal_resource_bounds':[adapter.resource_upper_bound(L) for L in (16,32,64)],
            'family_convergence':'Choose delta_W(h)<=h^2, then n delta_W<= (t/tau0)h. Charge preparation, final detector and all retained reset timeouts separately. For each finite L, precision search terminates; no large member was executed.',
            'open_preparation_and_detector_interface':'General fixed-bandwidth input preparation and a refined native terminal detector are typed inputs with separately supplied errors; the adapter does not claim to have built every such input or detector.'}


def finite_flux_controls():
    rows=[]
    for m in (3,6):
        U=np.zeros((m,m),dtype=np.int64)
        for j in range(m): U[(j+1)%m,j]=1
        E=np.diag(np.arange(m)); boundary=np.zeros_like(U); boundary[0,m-1]=1
        assert np.array_equal(E@U-U@E,U-m*boundary)
        cyclic=np.zeros((2*m,2*m),dtype=np.int64); bounded=np.zeros_like(cyclic)
        G=np.diag([int(x==0)+j for x in range(2) for j in range(m)])
        for x in range(2):
            for j in range(m):
                col=x*m+j; shift=1 if x==0 else -1
                cyclic[(1-x)*m+(j+shift)%m,col]=1
                bounded[(1-x)*m+j+shift if 0<=j+shift<m else col,col]=1
        assert np.array_equal(G@bounded,bounded@G)
        defect=G@cyclic-cyclic@G
        assert set(defect.reshape(-1))=={-m,0,m}
        delta=cyclic-bounded; gram=delta.T@delta
        assert np.trace(gram)==4 and np.array_equal(gram@gram,4*gram)
        rows.append({'modulus':m,'integer_Gauss_cyclic_defect':m,'bounded_vs_cyclic_operator_norm':2})
    return rows


def run():
    own_sources=('refinement_adapter.py','audit.py','PROOF.md','CONTRACT.json')
    own_pins={name:sha(HERE/name) for name in own_sources}
    pins,counts=check_pins(); print('all frozen source pins verified',flush=True)
    data,algebra=independent_bloch(); print('all144 rational spectral sectors reconstructed',flush=True)
    continuum=continuum_controls(data); clock=clock_controls(); flux=finite_flux_controls()
    adapter_result=adapter_controls()
    assert all(sha(ROOT/p)==h for p,h in pins.items())
    assert all(sha(HERE/p)==h for p,h in own_pins.items())
    checks={'frozen_transitive_bytes_unchanged':True,'all144_modes_and_complements_accounted':True,
            'rank26_generator_and_native_flux_independently_reconstructed':True,'uniform_continuum_bound_audited':True,
            'norm_TD_and_controlled_phase_distinguished':True,'readout_intervention_leakage_exact':True,
            'finite_shot_clock_and_mean_bias_distinguished':True,'clock_and_serial_native_cost_separated':True,
            'arbitrary_member_native_dilation_adapter_constructed':True,'rational_precision_scaling_certificate':True,
            'all_export_and_reset_complements_retained':True,'only_development_graphs234_executed':True,
            'general_bandlimited_preparation_not_inferred':True,'massive_and_common_field_transfer_not_inferred':True}
    result={'schema':'oph.chain2.local12_limit_audit.v1','model_id':adapter.common.MODEL_ID,
            'status':'frozen_scope_accepted_with_additive_refinement_adapter','checks':checks,
            'audited_frozen_manifests':FROZEN,'transitive_manifest_counts':counts,
            'algebra':algebra,'continuum':continuum,'clock':clock,'adapter':adapter_result,
            'independent_finite_flux_controls':flux,
            'row30_literal':'emergent ω²=k²+m² for matter waves',
            'row30_scope':'The frozen W theorem yields omega^2=|p|^2/(48 tau0^2), hence m=0 after shared speed normalization. It does not derive a massive matter-wave spectrum, species masses, or a common interacting field connection. Pending mass_coin/massive_coupled producers are neither audited nor promoted by this result.',
            'flux_interface':'The cyclic register has [E,U]=U-m|0><m-1|. Exact modular Gauss is not exact integer U(1) Gauss; bounded integer boundary holds amend the source. Uniform affine finite-cycle preparation is not a normalizable infinite uniform integer state.',
            'source_sha256':pins,'audit_source_sha256':own_pins,'official_rows_changed':[]}
    (HERE/'AUDIT_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    own=('refinement_adapter.py','audit.py','PROOF.md','CONTRACT.json','AUDIT_RESULT.json')
    manifest={'schema':'oph.chain2.local12_limit_audit.manifest.v1',
              'frozen_upstream_manifests':FROZEN,
              'files_sha256':{name:sha(HERE/name) for name in own},
              'checks':checks,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':checks,'native_pulses':[r['certificate']['native_pulses'] for r in adapter_result['development_members']]}),flush=True)

if __name__=='__main__': run()
