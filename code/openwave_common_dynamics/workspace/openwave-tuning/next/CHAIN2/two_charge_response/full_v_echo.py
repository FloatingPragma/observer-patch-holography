"""Promote the two-charge echo to one unchanged full one-face V_m.

Inverse local coins are preparation; inverse transport and local event
analysis are declared terminal apparatus with native words and errors.
"""
from fractions import Fraction as F
import hashlib
import json
import math
import numpy as np
import two_charge as t
import native_contract as nc


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    edge=json.loads((t.HERE/'EDGE_ECHO.json').read_text())
    mb=t.af.mb;compiler=mb.compiler
    face=t.af.data()[2][:,0];geometry=t.af.current_walk.geometry()
    candidates=[]
    for left,right,e in mb.mode_program()['transport_pairs']:
        a=mb.common.decode_position(left//16);b=mb.common.decode_position(right//16)
        particular=geometry['P'][:,a]+geometry['P'][:,b]
        if not face[e] and a==0 and np.count_nonzero(particular)==1 and particular[e]==1:
            candidates.append((left,right,e,a,b,particular))
    left,right,e,a,b,particular=candidates[0]
    assert not np.any(particular[np.flatnonzero(face)])
    assert np.count_nonzero(face)==3 and face[e]==0
    permutation=mb.common.expected_swap_permutation()
    assert int(permutation[left-1])==left-1 and int(permutation[right-1])==right-1
    assert int(permutation[left])==right and int(permutation[right])==left
    field=bytes((particular%16).astype(np.uint8))
    # Exact common field factor of the actual nonzero-beta source. On the
    # face orbit all three magnitudes equal the encoded centered label.
    program=json.loads((nc.BASE/'vector_field/FIELD_PROGRAM.json').read_text())
    word=program['plaquette_core_word']+3*program['electric_word']
    seed=np.zeros(16,complex);seed[0]=1
    actual_field=compiler.execute(word,4,seed)
    _,_,_,kick=t.face_data(.41)
    centered=(np.arange(16)+8)%16-8
    ideal_field=np.exp(-1j*.17*3*centered**2/2)*kick[:,0]
    field_error=float(np.linalg.norm(actual_field-ideal_field));assert field_error<2e-11
    # Executed controls retain their field factors and actual native phase
    # arrays, so control classification does not depend on prose alone.
    beta_only=compiler.execute(program['plaquette_core_word'],4,seed)
    alpha_off_seed=np.zeros(16,complex);alpha_off_seed[(1<<0)|(1<<2)]=1
    event_prepare=mb.givens_word(0,1,math.pi/4)+mb.givens_word(2,3,math.pi/4)
    alpha_off_events=compiler.execute(event_prepare+compiler.inverse(event_prepare),4,alpha_off_seed)
    alpha_off_event_probs=[float(abs(alpha_off_events[(1<<lm)|(1<<(2+rm))])**2)
                           for lm,rm in ((0,0),(0,1),(1,0),(1,1))]
    alpha_off_joint=np.outer(abs(beta_only)**2,alpha_off_event_probs)
    assert abs(alpha_off_joint.sum()-1)<2e-11 and abs(alpha_off_event_probs[0]-1)<2e-11
    assert np.linalg.norm(beta_only-kick[:,0])<2e-11
    electric_actual=compiler.execute(program['electric_word'],4)
    separate=[electric_actual[lm,lm]*electric_actual[1-rm,1-rm]
              for lm,rm in ((0,0),(0,1),(1,0),(1,1))]
    separate_connected=separate[0]*separate[3]/(separate[1]*separate[2])
    assert abs(separate_connected-1)<2e-11
    omitted_label=1;lost=float(abs(ideal_field[omitted_label])**2)
    retained=float(sum(abs(ideal_field[q])**2 for q in range(16)if q!=omitted_label))
    assert lost>1e-6 and abs(1-retained-lost)<2e-12
    executed_controls={
        'alpha_off':{'beta':.41,'all16_native_common_field_probabilities':(abs(beta_only)**2).tolist(),
                     'all4_native_event_probabilities':alpha_off_event_probs,
                     'all64_native_joint_probabilities':alpha_off_joint.tolist(),
                     'common_field_native_residual':float(np.linalg.norm(beta_only-kick[:,0])),
                     'joint_norm_squared':float(alpha_off_joint.sum())},
        'independent_electric_registers':{'all4_native_phase_values':[[float(z.real),float(z.imag)]for z in separate],
                     'initial_edge_fields':[0,1],
                     'preparation':'Separate P0 and P7 fields; their sum equals the shared initial P0+P7 field. Each individual charge has Gauss background delta0.',
                     'connected_phase':[float(separate_connected.real),float(separate_connected.imag)],
                     'connected_phase_residual_from_one':float(abs(separate_connected-1)),
                     'two_distinct_field_registers':True},
        'dropped_face_label':{'omitted_label':omitted_label,'omitted_probability':lost,
                              'retained_unconditional_norm':retained,'renormalized':False,'rejected':True},
        'geometry':{'face0_support_edges':np.flatnonzero(face).tolist(),'selected_edge_outside_face':True,
                    'held_modes_fixed_by_actual_S':True,'moving_modes_swapped_by_actual_S':True}}
    ideal_events=[(10+6*math.cos(.17))/16]+[(2-2*math.cos(.17))/16]*3
    law=[]
    for q in range(16):
        fieldq=t.af.add_field(field,q*face)
        for k,(lm,rm) in enumerate(( (0,0),(0,1),(1,0),(1,1) )):
            modes=[left-1+lm,right-1+rm]
            n=t.af.mode_field.port_loads(modes)
            assert not np.any((n+geometry['B']@np.frombuffer(fieldq,dtype=np.uint8)-2*np.eye(96,dtype=int)[:,0])%16)
            law.append({'occupied_modes':modes,'all288_field_hex':fieldq.hex(),
                        'probability':float(ideal_events[k]*abs(ideal_field[q])**2)})
    assert abs(sum(r['probability']for r in law)-1)<2e-12
    base=json.loads((nc.BASE/'vector_field/COMPOSITION.json').read_text())
    massdata=json.loads((nc.BASE/'particles/MASSIVE_COUPLED.json').read_text())
    root=json.loads((t.INPUT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    compiler_programs=[];bounds=[]
    for mass in (0,1):
        local_coin=nc.massive.fock_kick_word(F(1,8),F(mass))+mb.fock_reflection_word()
        inverse_words=[]
        for vertex in (a,b):
            group=list(range(16*mb.common.encode_position(vertex),16*mb.common.encode_position(vertex)+16))
            inverse_words+=compiler.remap(compiler.inverse(local_coin),group)
        # Local event qubits are (held,moving) at each selected endpoint.
        support=[left-1,left,right-1,right]
        prep=[]
        for mode in (left-1,right-1,2048+4*e):
            prep+=compiler.remap(mb.architecture.primitive_routines()['flip'],[mode])
        prep+=mb.givens_word(left-1,left,math.pi/4)+mb.givens_word(right-1,right,math.pi/4)
        prep+=inverse_words
        rephase=(mb.occupation_phase_word(left,.17*1.5)+mb.occupation_phase_word(right,-.17*.5))
        analysis=mb.givens_word(left-1,left,-math.pi/4)+mb.givens_word(right-1,right,-math.pi/4)
        massrow=next(r for r in massdata['experiments']if r['mass_parameter']==str(mass))
        onsite=(F(base['onsite_GammaR_error_exact'])+
                F(massrow['native_complete_error']['onsite_Gamma_mass_kick_operator_error_exact']))
        # The exact literal inverse has the same operator error as forward.
        helper=nc.massive.mass_coin.source_error.exact_template_bound
        prep_error=3*F(root['primitive_operator_errors_exact']['flip'])+2*onsite
        analysis_error=F()
        for x in (left,right):
            prep_error+=helper(mb.givens_word(x-1,x,math.pi/4),nc.composition.givens_template(x-1,x,nc.composition.cert._pi(1,4)))
            analysis_error+=helper(mb.givens_word(x-1,x,-math.pi/4),nc.composition.givens_template(x-1,x,nc.composition.cert._pi(-1,4)))
        for mode,angle in ((left,F(.17)*F(3,2)),(right,-F(.17)/2)):
            analysis_error+=helper(mb.occupation_phase_word(mode,float(angle)),nc.composition.phase_template(mode,angle))
        apparatus=nc.errors(1,True,mass,0)
        total=(F(apparatus['complete_native_TV_upper_exact'])-F(apparatus['preparation_operator_error_exact'])+
               prep_error+F(massrow['native_complete_error']['whole_step_operator_error_exact'])+
               F(base['complete_matter_current_layer_error_exact'])+analysis_error)
        inverse_S_pulses=base['source_pulses_per_matter_current_layer']-128*len(mb.fock_reflection_word())
        source_step_pulses=nc.errors(1,True,mass,1)['evolution_native_pulses_per_step']
        total_pulses=len(prep)+source_step_pulses+inverse_S_pulses+len(rephase)+len(analysis)+19200
        bounds.append({'mass':mass,'h':'1/8','full_native_TV_upper_exact':str(total),
                       'full_native_TV_upper':math.nextafter(float(total),math.inf),
                       'complete_native_pulses_excluding_resets':total_pulses,
                       'preparation_native_pulses':len(prep),'full_V_native_pulses':source_step_pulses,
                       'inverse_full_S_native_pulses':inverse_S_pulses,'analysis_native_pulses':len(rephase)+len(analysis),
                       'readout_native_pulses':19200,
                       'blank_resets':3203,'ready_resets':3200,
                       'all_reset_native_attempts_max':apparatus['all_reset_native_attempts_max'],
                       'source_bits':3203,'pointer_records':3200,'field_bits':1152,'matter_bits':2048})
        compiler_programs.append({'mass':mass,'source_preparation_word':prep,
             'full_V_word':'Prepend128remapped local mass-kick words to frozen native_field.complete_word(.17,[.41]+[0]*207); no clock spacing appended.',
             'terminal_inverse_S_word_generator':'Literal inverse of all288field-current-block-plus-full-CAR-swap words in frozen transport order; no onsite reflection.',
             'terminal_rephase_word':rephase,'terminal_analysis_word':analysis})
    result={'schema':'oph.two_charge.unchanged_full_V_echo.v1','candidate_id':t.MODEL,
            'status':'executed_exact_full_source_factorization_and_native_factors',
            'alpha':.17,'beta_face0':.41,'other207_beta':0,'edge':e,'vertices':[a,b],
            'global_qubit_modes':support,'initial_G':'2delta0','initial_field_hex':field.hex(),
            'source_update':'Exactly the frozen V_m=UE UB Gamma(S)_F Gamma(R exp(-ihmR)); h1/8,m0or1. Neither beta nor the source transport is removed.',
            'preparation':'Two equal held/moving local orbitals followed by inverse local Gamma(C_m) at each endpoint. The full source coin restores these two event qubits before transport. This is explicitly different from a position-delta uniform-event preparation.',
            'apparatus':'Full native inverse S/current word, explicit one-body phase calibration, local native Givens analysis, then complete3200-bit native occupation/field readout.',
            'factorization':'The selected transported edge is outside face0. All source current is on that edge. The remaining field evolves by UE_rest UB_face0 as one common factor, including every16cycle label. Inverse S returns the transported field to its initial value. Local inverse coin preparation cancels the unchanged full source coin.',
            'two_charge_connected_phase':'exp(i*alpha), independent of supplied mass under the specified mass-dependent preparation.',
            'event_TV_vs_additive_control':float(3*(1-math.cos(.17))/8),
            'all64_joint_event_field_outcomes':law,'ideal_norm_squared':sum(r['probability']for r in law),
            'native_common_field_factor_residual':field_error,
            'local_native_echo_evidence':{'file':'EDGE_ECHO.json','sha256':sha(t.HERE/'EDGE_ECHO.json')},
            'native_complete_contracts':bounds,'emitted_programs':compiler_programs,
            'controls':{'coupling_alpha_off':'connected phase1 and event00 probability1; UB common field factor still retained.',
                        'independent_field_substitution':'Two separate electric registers yield connected phase1 after their own one-body calibration.',
                        'CAR_phase':'Every transported pair has canonical signs in native source and inverse; same-Gauss coherent mutation witness is in RESULT.',
                        'dropped_field_outcomes':'All16plaquette labels appear in every event branch; none is omitted or conditioned.'},
            'executed_control_arrays':executed_controls,
            'clock':'No optional phase-clock word or clock/probe readout is included. Macrostep semantics only; clock relocation and errors require their separate actual composition.',
            'claim':'Concrete native two-charge connected interaction under one unchanged full source step and a charged, declared inverse/analysis instrument.',
            'nonclaim':'This interference record is not a freely observed first-particle displacement, unperturbed multistep force trajectory, Coulomb law, reciprocal acceleration, physical mass or continuum result.',
            'input_snapshot_sha256':sha(t.HERE/'INPUT_SNAPSHOT.json')}
    (t.HERE/'FULL_V_ECHO.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'edge':e,'vertices':[a,b],'event_TV':result['event_TV_vs_additive_control'],
                      'native_TV':[r['full_native_TV_upper']for r in bounds],
                      'common_field_residual':field_error}),flush=True)


if __name__=='__main__':run()
