#!/usr/bin/env python3
"""One native history, two local clocks, and the existing product-field code.

The finite executed domain and the conditional continuum Lorentzian theorem
are typed separately. No physical metric, radial apparatus or particle law
is inferred merely from the finite clock's successful execution.
"""
from __future__ import annotations
from evidence_contract import public_receipt
from collections import Counter
from fractions import Fraction as F
import hashlib,json,math,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import sympy as s
from scipy.linalg import expm,logm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PUBLIC=ROOT
sys.path.insert(0,str(PUBLIC))
from model.observer import native
from model.observer.source_scalar import compiler
from model.observer_time_fields import continuum_dynamics as field
from model.physical_space_time import source_green,observation_process
from model.relative_observer_geometry import intrinsic_geometry

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,s.Basic):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def finite_apparatus(native_eigenvalue=6):
    dictionary=observation_process.product_mode_dictionary(native_eigenvalue,radial_mode=1,radial_points=3,length=F(2))
    q=630;theta=s.acos(1-s.Rational(native_eigenvalue)/(2*q));a=s.sqrt(q)*theta;c=s.sqrt(dictionary['radial_eigenvalue'])
    omega=s.sqrt(a*a+c*c);beta=F(1);angle=s.atan2(c,a)
    r=field._source_dyadic_scalar(s.tanh(omega/2));rf=float(r)
    rotate=sum(([compiler.pulse('Y',-float(field._source_dyadic_scalar(angle/2)),[j])] for j in range(2)),[])
    step=[compiler.pulse('Z',float(field._source_dyadic_scalar(theta/2)),[0]),
          compiler.pulse('X',float(field._source_dyadic_scalar(c/(2*s.sqrt(q)))),[0])]
    prep=field._thermal_code_preparation(r,2)+rotate
    blank=np.zeros(32,complex);blank[0]=1
    prepared=compiler.execute(prep,4,blank[:16])
    x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]],complex);z=np.diag([1.,-1.])
    af,cf,of=map(lambda v:float(s.N(v,70)),(a,c,omega))
    normal=(cf*x+af*z)/of;transverse=(af*x-cf*z)/of
    rho=(np.eye(2)+rf*normal)/2
    reduced=prepared.reshape(4,4).T@prepared.reshape(4,4).conj()
    assert np.linalg.norm(reduced-np.kron(rho,rho))<3e-11
    def weyl(strength):
        normalized=float(s.N(s.Rational(strength)/s.sqrt(s.Rational(r)*omega),70))
        word=field.spin_weyl_word([[normalized,0]],2,controlled=True)
        return compiler.inverse(rotate)+compiler.remap(word,[0,1,4])+rotate
    first_strength,second_strength=F(1,2),F(1,3)
    outputs=[];branches=[];records=[]
    step_matrix=compiler.execute(step,1)
    for count in (0,1,2):
        evolve=sum((compiler.remap(step*count,[j]) for j in range(2)),[])
        by_axis=[];by_branch=[]
        for axis in ('X','Y'):
            word=prep+[compiler.pulse('Y',-math.pi/4,[4])]
            word+=weyl(second_strength)+evolve+weyl(first_strength)+compiler.inverse(evolve)
            word+=[compiler.pulse('Y' if axis=='X' else 'X',math.pi/4 if axis=='X' else -math.pi/4,[4])]
            vector=compiler.execute(word,5,blank)
            law=abs(vector)**2;assert abs(law.sum()-1)<3e-11
            by_axis.append(vector);by_branch.append(vector.reshape(2,16))
            records.append({'incident_count':count,'setting':axis,'source_pulses':len(word),
                'word_sha256':hashlib.sha256(json.dumps(word,sort_keys=True).encode()).hexdigest(),
                'complete32_fine_refinement_probabilities':law.tolist(),
                'actual_pointer_only_probabilities':law.reshape(2,16).sum(axis=1).tolist(),
                'pointer_mean':float(law[:16].sum()-law[16:].sum())})
        outputs.append(by_axis);branches.append(by_branch)
    # Owning public producer independently computes ordered matrix laws,
    # modular matching, signed kicks and split-error contracts on the same
    # state, band and controls. Include a held-out25-trigger experiment.
    public={n:field.product_mode_source_control(native_eigenvalue,radial_mode=1,radial_points=3,length=F(2),
              inverse_temperature=beta,steps=n) for n in (1,2,25)}
    errors=[]
    for n in (1,2):
        for k,axis in enumerate(('X','Y')):
            oracle=public[n]['complete_ordered_two_time_pointer_records'][axis]['all_system_reference_pointer_probabilities']
            errors.append(float(np.max(np.abs(abs(np.array(outputs[n][k]))**2-oracle))))
    assert max(errors)<3e-11
    # An independent finite matrix formula also checks the zero-trigger law.
    def ordered_for_u(u):
        def w(strength):
            t=float(strength)/math.sqrt(4*rf*of)
            return math.cos(t)*np.eye(2)+1j*math.sin(t)*transverse
        return complex(np.trace(rho@u.conj().T@w(first_strength)@u@w(second_strength))**2)
    initial_ordered=ordered_for_u(np.eye(2))
    assert abs(records[0]['pointer_mean']-initial_ordered.real)<3e-11
    assert abs(records[1]['pointer_mean']-initial_ordered.imag)<3e-11
    h=(of*np.eye(2)-af*z-cf*x)/2
    return {'dictionary':dictionary,'a':a,'c':c,'omega':omega,'beta':beta,'polarization':r,
        'rho_one':rho,'h_one':h,'prepared_state':prepared,'step':step_matrix,
        'records':records,'outputs':np.array(outputs),'branches':np.array(branches),
        'public_controls':public,'crosscheck_error':max(errors),'ordered_for_u':ordered_for_u}

def native_joint_history(apparatus):
    graph,lap,_,_,_=source_green.matrices(0);l=np.asarray(lap.toarray(),dtype=np.int64)
    assert graph.ports==240 and len(graph.seams)==630
    fixture=field.product_native_mode_control()
    profile=np.zeros(graph.ports,dtype=np.int64)
    profile[:12]=[int(s.sympify(v)) for v in fixture['exact_first_carrier_profile']]
    assert np.array_equal(l@profile,6*profile) and sum(profile)==0
    assert np.all(profile.reshape(20,12).sum(axis=1)==0)
    clock_a=next(p for p,nbs in enumerate(graph.neighbours) if len(nbs)==5)
    clock_b=next(p for p in graph.neighbours[clock_a] if len(graph.neighbours[p])==6)
    da,db=len(graph.neighbours[clock_a]),len(graph.neighbours[clock_b])
    initial_ports=(next(i for i,v in enumerate(profile) if v==1),12)
    initial=tuple(int(i in initial_ports) for i in range(graph.ports));f0=int(profile@initial)
    catalogue=[initial];lookup={initial:0};tables={}
    def transitions(source_id):
        if source_id in tables:return tables[source_id]
        before=catalogue[source_id];rows=[]
        for seam_id,seam in enumerate(graph.seams):
            for coin in (0,1):
                after,(_,current)=native.native(before,seam,coin)
                assert sum(after)==2 and max(after)==1
                if after not in lookup:lookup[after]=len(catalogue);catalogue.append(after)
                target=lookup[after]
                decoded=int(profile@before)+current*int(profile[seam[1]]-profile[seam[0]])
                assert decoded==int(profile@after)
                rows.append((seam_id,coin,current,target,int(clock_a in seam),int(clock_b in seam),decoded))
        tables[source_id]=np.array(rows,dtype=np.int32);return tables[source_id]
    first=transitions(0)
    first_parents=sorted(set(int(r[3]) for r in first))
    for source_id in first_parents:transitions(source_id)
    first_summary=Counter((int(r[3]),int(r[4]),int(r[5]),int(r[6])) for r in first)
    joint=Counter()
    for (state,na,nb,fa),multiplicity in first_summary.items():
        for row in tables[state]:
            dest,ma,mb,fb=map(int,row[3:])
            joint[(dest,na+ma,nb+mb,fb)]+=multiplicity
    total=1260**2;assert sum(joint.values())==total
    # Full fine record law is represented by exact first/conditional tables;
    # every one of1260^2 histories retains its native state and record pair.
    count_law=Counter();field_mean=F();clock_mismatch=0
    for (state,na,nb,f),n in joint.items():
        count_law[na,nb]+=n;field_mean+=F(n*f,total)
        if db*na!=da*nb:clock_mismatch+=n
        assert f==int(profile@catalogue[state])
    assert field_mean==F(f0)*F(1-F(6,1260))**2
    mean_a=sum((F(n*na,total) for (na,nb),n in count_law.items()),F())
    mean_b=sum((F(n*nb,total) for (na,nb),n in count_law.items()),F())
    assert mean_a==F(2*da,630) and mean_b==F(2*db,630)
    assert mean_a/mean_b==F(da,db)
    # Native fine records control the same prepared quantum code, including
    # zero-current incident events. No successful-current conditioning.
    combined=[];conditional_pointer=np.zeros((3,2,2))
    for na in range(3):
        for axis in range(2):conditional_pointer[na,axis]=np.sum(abs(apparatus['branches'][na,axis])**2,axis=1)
    for (state,na,nb,f),n in sorted(joint.items()):
        combined.append({'native_terminal_state_id':state,'clock_A':na,'clock_B':nb,'decoded_field':f,
            'native_numerator':n,'native_denominator':total,
            'X_pointer_joint_probabilities':(float(F(n,total))*conditional_pointer[na,0]).tolist(),
            'Y_pointer_joint_probabilities':(float(F(n,total))*conditional_pointer[na,1]).tolist()})
    for axis in ('X','Y'):
        assert abs(sum(sum(row[f'{axis}_pointer_joint_probabilities']) for row in combined)-1)<3e-11
    exact_mean=sum(float(F(n,total))*(conditional_pointer[na,0,0]-conditional_pointer[na,0,1]) for (na,nb),n in count_law.items())
    mean_clock_u=expm(float(mean_a)*logm(apparatus['step']))
    mean_clock_prediction=apparatus['ordered_for_u'](mean_clock_u).real
    mean_clock_defect=abs(exact_mean-mean_clock_prediction)
    assert mean_clock_defect>1e-9
    # Each observer's classical view is a deterministic restriction of this
    # same joint instrument; overlap-field and shared pointer laws coincide.
    localA=Counter();localB=Counter()
    for (_,na,nb,f),n in joint.items():localA[na,f]+=n;localB[nb,f]+=n
    assert {f:sum(n for (na,g),n in localA.items() if g==f) for f in set(g for _,g in localA)}=={
        f:sum(n for (nb,g),n in localB.items() if g==f) for f in set(g for _,g in localB)}
    arrays={'native_seams':np.array(graph.seams,dtype=np.int16),'native_laplacian':l,
        'native_exact_angular_profile':profile,'native_first_records':first,
        'native_second_parent_ids':np.array(first_parents,dtype=np.int16),
        'native_second_records':np.array([tables[i] for i in first_parents]),
        'native_state_catalogue_occupied_pairs':np.array([[i for i,x in enumerate(row) if x] for row in catalogue],dtype=np.int16),
        'quantum_prepared_work_reference_state':apparatus['prepared_state'],
        'quantum_preterminal32_amplitudes':apparatus['outputs'],
        'quantum_pointer_branch16_amplitudes':apparatus['branches'],
        'work_single_code_hamiltonian':apparatus['h_one'],'work_single_code_reference':apparatus['rho_one'],
        'actual_source_split_step':apparatus['step']}
    return {'source_ports':240,'source_seams':630,'native_attempts':2,'all_fine_native_histories':total,
        'conditional_native_parent_states':len(first_parents),'initial_occupied_ports':initial_ports,
        'initial_preparation_disclosed':True,'native_record_row_format':['seam_index','coin','signed_current','destination_state_id','clock_A_incident','clock_B_incident','decoded_profile_value'],
        'clock_ports':[clock_a,clock_b],'clock_degrees':[da,db],
        'zero_current_incidents_retained':True,'mean_clock_counts':[str(mean_a),str(mean_b)],
        'B_to_A_mean_display_factor':str(F(da,db)),'finite_pathwise_clock_display_mismatch_probability':str(F(clock_mismatch,total)),
        'decoded_field_mean_after_two_attempts':str(field_mean),'same_native_density_formula':'(1-6/1260)^2 times initial profile',
        'shared_field_marginals_agree':True,'joint_record_rows':combined,
        'native_history_joint_law_complete':True,'replace_actual_incident_count_by_mean_pointer_error':mean_clock_defect,
        'timing_contract':'Actual N_A controls protected two-axis source pulses; t=N_A/sqrt(630). Repairs pause during atomic coherent services.',
        'field_setting_contract':'X orY is a classical requested setting; incompatible pointer quadratures are not asserted jointly measured.',
        'retained_quantum_output_contract':'Only pointerbit is committed. The unnormalized16-dimensional work/reference vector remains for downstream measurements. Full32 computational probabilities are an optional refinement, not an already executed destructive work read.'},arrays

def low_angular_control(joint,arrays):
    """Actual low angular eigenmode, with explicitly numerical membership.

    The source code uses a dyadic approximation of the native eigenvalue.
    This does not compile its 240-component spatial preparation/encoder.
    All original clock-prefix weights are reused, without a fresh clock law.
    """
    lap=arrays['native_laplacian'];values,vectors=np.linalg.eigh(lap.astype(float))
    lam=float(values[1]);v=vectors[:,1]
    if v[np.argmax(abs(v))]<0:v=-v
    rounded=F(round(lam*2**48),2**48)
    residual=float(np.linalg.norm(lap@v-float(rounded)*v))
    compressed=v.reshape(20,12).mean(axis=1)
    cell_weight=float(12*compressed@compressed)
    assert 0<lam<.05 and residual<5e-13 and cell_weight>.98
    apparatus=finite_apparatus(rounded)
    weights=[F(0)]*3
    for row in joint['joint_record_rows']:
        weights[row['clock_A']]+=F(row['native_numerator'],row['native_denominator'])
    laws=np.array([abs(x)**2 for x in apparatus['outputs']])
    mixed=np.einsum('n,nak->ak',np.array(list(map(float,weights))),laws)
    assert np.max(abs(mixed.sum(axis=1)-1))<3e-11
    arrays['low_angular_profile_numerical']=v
    arrays['low_angular_quantum_preterminal32_amplitudes']=apparatus['outputs']
    return {'actual_source_operator':'full240 native Laplacian from source_green.matrices(0)',
        'lowest_nonzero_eigenvalue_numerical':lam,'used_dyadic_native_eigenvalue':str(rounded),
        'eigenpair_residual_numerical':residual,'membership_is_a_numerical_diagnostic':True,
        'cell_average_weight_numerical':cell_weight,
        'internal_weight_bound_at_numeric_eigenvalue':lam/(5-math.sqrt(5)),
        'first_nonzero_spectrum_numerical':values[1:8].tolist(),
        'radial_dictionary':apparatus['dictionary'],'source_angular_frequency':str(apparatus['a']),
        'source_joint_frequency':str(apparatus['omega']),
        'same_native_clock_prefix_weights':[str(x) for x in weights],
        'complete_source_records':apparatus['records'],
        'same_band_state_clock_response_error_controls':apparatus['public_controls'],
        'owning_public_control_crosscheck_error':apparatus['crosscheck_error'],
        'complete_clock_mixed32_laws_by_setting':mixed.tolist(),
        'new_spatial_eigenvector_encoder_executed':False,
        'spatial_eigenmode_and_mode_code_join':'The eigenvalue is computed from the actual source L and passed through the existing declared mode-code dictionary; this numeric control compiles the mode code, not a new full240 spatial encoder.',
        'single_finite_low_mode_declared_a_continuum_eigenfunction':False,
        'scope':'The low finite source band is cell-visible, unlike the exact lambda6 fixture. Continuum membership still uses the stated common-subsequence and bounded-band construction, not this one finite diagonalization.'}

def metric_and_limits():
    mu,a,delta=s.symbols('mu a delta',positive=True);b=s.symbols('b',real=True)
    screen=s.Matrix([[a,b],[b,(b*b+delta)/a]])
    spatial=s.diag(1,screen);coefficient=mu*s.diag(-1,spatial)
    volume=mu**2*s.sqrt(delta);metric=s.simplify(volume*coefficient.inv())
    assert s.simplify(metric.det()-coefficient.det())==0
    assert s.simplify(volume*metric.inv()-coefficient)==s.zeros(4)
    k=s.symbols('k',positive=True);jac=s.diag(k,1,1,1)
    transformed=s.simplify(jac.T*metric*jac)
    clock_coefficient=k*mu*s.diag(-1/k**2,spatial)
    assert s.simplify(s.sqrt(-clock_coefficient.det())*clock_coefficient.inv()-transformed)==s.zeros(4)
    bounds=[field.product_clock_bound(l,degree=5,tolerance=F(1,2**(l//2+3)),step_error=F(1,2**(44+4*l))) for l in (0,4,8,12)]
    assert bounds[-1]['combined_single_code_unitary_error_upper']<bounds[0]['combined_single_code_unitary_error_upper']
    assert bounds[-1]['clock_failure_probability_upper']<bounds[0]['clock_failure_probability_upper']
    return {'same_limit_operator':'B=-d_r^2+(A direct_sum0), A is the inverse of the native angular Green limit on its selected common spatial subsequence',
        'same_field_state':'Z^-1 exp(-beta dGamma(sqrt(B)))','field_beta':1,
        'same_ordered_covariance':'[coth(beta sqrt(B)/2) cos(t sqrt(B))-i sin(t sqrt(B))]/(2 sqrt(B))',
        'same_retarded_response':'B^(-1/2) sin(t sqrt(B))','same_modular_clock':'sigma_s=alpha_(-beta s)',
        'incident_clock_limit':'t=(2d/sqrt(630))*tau; local normalizedcounts sharetau in the stated concentration/refinement limit',
        'action_density':str(coefficient),'metric':str(metric),'lapse':str(s.sqrt(mu)*delta**s.Rational(1,4)),
        'Lorentz_signature':[-1,1,1,1],'density_to_metric_and_clock_covariance_exact':True,
        'coefficient_scope':'Symbolic measurable limiting form with mu>0 and positive2D screen tensor. No finite measured continuum tensor is substituted.',
        'comparison_support_cone':'|t-u|^2 >= |r-rprime|^2+2 d_screen^2',
        'full_radial_constant_restored_and_Dirichlet_gap_positive':True,
        'smooth_metric_or_curvature_selected':False,
        'smooth_effective_description':{
            'quantifiers':'For each fixed finite admissible menu, bounded wave-time interval and positive observation tolerance, one sufficiently fine member of the fixed smooth coefficient-approximation sequence describes the entire menu.',
            'admissible_menu':'Finite smeared thermal Weyl words and finite adaptive Ramsey instruments formed from finite Weyl sums/feedback, with every classical record and any fixed finite detector/reference output. Momentum smearings lie in the fixed H1 form domain with zero radial-end traces.',
            'common_support_preserving_dictionary':'U_epsilon f=sqrt(rho_epsilon/rho) f; smooth-theory smearing U_epsilon^-1 f',
            'same_native_clock_tapes_and_controller_retained':True,
            'errors':'Source synthesis, copy/radial/angular refinement, clock concentration/timeouts and uniform smooth-description errors add within the proved contracts.',
            'one_exact_smooth_metric_for_all_observations_claimed':False,
            'geometry_selected_separately_per_readout':False,
            'computable_smoothing_rate_for_unknown_limit_claimed':False,
            'arbitrary_sharp_projectors_or_full_surviving_field_trace_norm_claimed':False,
            'exact_native_density_smooth':False,
            'source':'model/relative_observer_geometry/README.md:844-936'},
        'finite_240_port_band_itself_is_a_local3D_continuum':False,
        'growing_encoder_execution_claimed':False,
        'public_clock_error_contracts':bounds,
        'passive_density_vs_product_field':{'passive':'native density heat semigroup exp(-tau A), Gaussian path OS generator dGamma(A)',
            'protected':'separately encoded product field wavegenerator sqrt(B), faithful thermal state and active incidentclock',
            'same_native_angular_operator_family':True,'same_instrument_or_same_state_asserted':False},
        'remaining_declared_inputs':['radial interval/grid/boundary and itsX control','mode access and independently prepared copies',
            'common beta and protected code transport/storage','complete native current and zero-inclusive clock records',
            'uniform IID scheduler with paused coherent services','controlled synthesis precision in refinement'],
        'common_source_ProductBand_morphism':{'finite_common_source_oscillator_histories_identified_with_this_product_field':False,
            'existing_compiler_shared':True,'intervention_preserving_refinement_embedding_emitted_here':False}}

def main():
    apparatus=finite_apparatus();joint,arrays=native_joint_history(apparatus);limits=metric_and_limits()
    low=low_angular_control(joint,arrays)
    artifact=HERE/'results/spacetime_correspondence_arrays.npz';np.savez_compressed(artifact,**arrays)
    paths=[]
    for module in tuple(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f and Path(f).is_file() and Path(f).is_relative_to(PUBLIC):paths.append(Path(f))
    paths += [PUBLIC/p for p in compiler.source_summary()['source_sha256']]
    paths += [PUBLIC/p for p in ('model/CONSTRUCTION.md','model/observer_time_fields/REFINEMENT.md',
        'model/relative_observer_geometry/README.md','model/physical_space_time/OCCUPATION_FIELD.md',
        'model/OBSERVER_RELATIVE_RECOVERY.md')]
    result={'status':'PASS_SAME_HISTORY_NATIVE_FIELD_CLOCK_AND_CONDITIONAL_LORENTZIAN_PRODUCT_CORRESPONDENCE',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))},
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'finite_product_domain':{'native_angular_eigenvalue':6,'radial_mode':1,'radial_points':3,'radial_length':2,
            'beta':1,'Q':630,'angular_frequency':str(apparatus['a']),'radial_frequency':str(apparatus['c']),
            'joint_frequency':str(apparatus['omega']),'polarization':str(apparatus['polarization']),
            'quantum_copy_count':2,'work_qubits':[0,1],'reference_qubits':[2,3],'pointer_qubit':4,
            'finite_spin_fluctuations_are_exact_CCR_field':False,
            'actual_lambda6_profile_is_internal_with_zero_cell_compression':True,
            'profile_survives_bounded_nonconstant_continuum_band':False},
        'complete_source_records':apparatus['records'],'owning_public_control_crosscheck_error':apparatus['crosscheck_error'],
        'same_band_state_clock_response_error_controls':apparatus['public_controls'],
        'joint_native_history':joint,'continuum_correspondence':limits,'low_angular_mode_source_control':low,
        'downstream_state_interface':{'artifact':'spacetime_correspondence_arrays.npz',
            'branch_array':'quantum_pointer_branch16_amplitudes[count0to2,settingXorY,pointerbit,work_reference_address]',
            'native_weights':'joint_native_history.joint_record_rows; exactnumerator/denominator',
            'matrix_bit_order':'littleendian; workbits0,1, references2,3',
            'branch_states_unnormalized_and_no_postselection':True,
            'downstream_effects_may_be_appended_without_reset':True},
        'integrated_physical_claim':{'conditional_free_product_field_with_Lorentz_signature_supported':True,
            'all_same_history_measurements_already_identified_as_SM_fields_and_GR':False,
            'particle_or_Einstein_dynamics_selected':False,'fundamental_spacetime_required':False}}
    (HERE/'results/spacetime_correspondence.json').write_text(json.dumps(serial(public_receipt(result)),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'fine_native_histories':joint['all_fine_native_histories'],
        'field_clock_joint_rows':len(joint['joint_record_rows']),'mean_clock_substitution_error':joint['replace_actual_incident_count_by_mean_pointer_error'],
        'source_law_crosscheck_error':apparatus['crosscheck_error']},indent=2))

if __name__=='__main__':main()
