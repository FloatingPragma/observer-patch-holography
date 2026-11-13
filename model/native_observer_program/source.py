"""Executable complete observation policies on one padded two-load carrier."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import pi
import json
import numpy as np
from .contract import (manifest,validate_manifest,validate_state,source_addresses,
                       decoded_basis,matching_vector,memory_character,descriptor_sha256)
from .reflection import plan,run_reflections,complete_bound
from model.observer.source_scalar import compiler
from model.common_source import record_code,parameter_selection
from model.openwave_native_process.source import raw_local_branches
from model.openwave_active_native_charge import source as active
from model.openwave_native_coherent_clock import source as clock
from model.common_source import repair_continuation as continuation
from model.native_charged_propagation import source as magnetic


@lru_cache(maxsize=1)
def compiled_event_matrices():
    p=plan(60)
    return (compiler.execute(p['preparation'],p['bits']),
            compiler.execute(p['reflection'],p['bits']))


@lru_cache(maxsize=1)
def compiled_inverse_reflection():
    """Execute requested inverse native pulses, not an adjoint of a result."""
    from model.native_effect_synthesis.synthesis import certify_primitive
    p=plan(60);word=compiler.inverse(p['reflection'])
    for kind,angle in {(x['kind'],x['angle'])for x in word}:certify_primitive(kind,angle)
    return compiler.execute(word,p['bits'])


def event_action(state,matrix):
    value=validate_state(state)
    if np.asarray(matrix).shape!=(64,64):raise ValueError('complete64-address event operation required')
    return np.einsum('ef,smfr->smer',matrix,value,optimize=True)


def prepare(program=None):
    validate_manifest(manifest()if program is None else program)
    psi=(matching_vector(0)+matching_vector(1))/np.sqrt(2)
    state=np.zeros((512,8,64,1),complex)
    state[:,:6,0,0]=psi[:,None]*memory_character(1)[None,:]
    return event_action(state,compiled_event_matrices()[0])


def prepare_idle_mode(program=None):
    """Explicit alternative apparatus-mode preparation; not the clock input."""
    validate_manifest(manifest()if program is None else program)
    state=np.zeros((512,8,64,1),complex)
    state[:,:6,0,0]=matching_vector(0)[:,None]*memory_character(1)[None,:]
    word=[compiler.pulse('H',pi/2,[1]),compiler.pulse('Z',pi/2,[1])]
    from model.native_effect_synthesis.synthesis import certify_primitive
    for row in word:certify_primitive(row['kind'],row['angle'])
    return event_action(state,compiler.execute(word,6))


def idle_mode_control(program=None):
    state=prepare_idle_mode(program);moved=spatial_step(state,stored=True)
    overlap=np.einsum('smer,e->smr',state,np.r_[np.ones(60)/np.sqrt(60),np.zeros(4)])
    error=float(np.linalg.norm(moved+state))
    assert error<2e-10 and np.linalg.norm(overlap)<2e-11
    return {'preparation':'matching word0, memory character1, event(|0>-|2>)/sqrt2 up to the recorded native word global phase',
      'stored_native_preparation_pulses':2,'source_memory_untouched_by_event_preparation':True,
      'all64_event_columns_executed':True,'full_W_eigenvalue_target':'-1',
      'stored_full_step_residual_diagnostic':error,'uniform_A_overlap_diagnostic':float(np.linalg.norm(overlap)),
      'complete_terminal_records':terminal_records(moved),
      'preparation_error_bound':complete_bound(2,F(2,10**15)),
      'physical_particle_or_rest_energy_identification':False}


def event_swap(state):
    """Complete reversible source-memory feedback, including all padding."""
    value=validate_state(state);out=value.copy();g=magnetic.carrier();codes=source_addresses()
    for seam_index,edges in enumerate(g.seam_edges):
        event=2*seam_index+1
        for i,j,chi in edges:
            out[codes[j],:6,event]=np.roll(value[codes[i],:6,event],chi,axis=0)
            out[codes[i],:6,event]=np.roll(value[codes[j],:6,event],-chi,axis=0)
    return out


def spatial_step(state,*,inverse=False,stored=False):
    if type(inverse)is not bool or type(stored)is not bool:raise ValueError('explicit direction and implementation flags required')
    if stored:R=compiled_inverse_reflection()if inverse else compiled_event_matrices()[1]
    else:
        u=np.zeros(64);u[:60]=1/np.sqrt(60);R=2*np.outer(u,u)-np.eye(64)
    return event_action(event_swap(state),R)if inverse else event_swap(event_action(state,R))


def terminal_records(state):
    value=validate_state(state);codes=source_addresses()
    admitted=value[np.array(codes),:6,:60,:]
    char=np.stack([memory_character(k)for k in range(6)],axis=1)
    amps=np.einsum('km,nmer->nker',char.conj().T,admitted,optimize=True)
    probabilities=np.sum(abs(amps)**2,axis=(2,3))
    # Compute the orthogonal failure amplitude directly, never by clamping.
    outside=value.copy();outside[np.array(codes),:6,:60,:]=0
    failure=float(np.vdot(outside,outside).real)
    total=float(np.vdot(value,value).real)
    assert abs(float(probabilities.sum())+failure-total)<3e-11
    return {'source_pair_character_probabilities':probabilities.tolist(),
            'complement_probability':failure,'input_norm_squared':total,
            'outcomes':397,'normalization_performed':False}


@lru_cache(maxsize=1)
def energy_diagonal():
    value=np.zeros(512)
    for state,energy in zip(parameter_selection.code_geometry()['states'],parameter_selection.diagonal_cost(2)['values'],strict=True):
        from model.observer.source_discrimination import circuit
        value[circuit.encode(state)]=float(energy)
    return value


def source_energy(state):
    a=validate_state(state)
    return float(np.sum(energy_diagonal()*np.sum(abs(a)**2,axis=(1,2,3))))


def observe_spatial(state,steps,program=None):
    validate_manifest(manifest()if program is None else program)
    if type(steps)is not int or not 0<=steps<=8:raise ValueError('declared zero through8 coherent steps required')
    initial=validate_state(state);actual=initial.copy();ideal=initial.copy();work=[]
    for _ in range(steps):
        before=source_energy(actual);actual=spatial_step(actual,stored=True);ideal=spatial_step(ideal)
        work.append(source_energy(actual)-before)
    records=terminal_records(actual)
    records.update({'policy':'coherent_spatial','steps':steps,'source_work_per_step_diagnostic':work,
        'actual_stored_reflection_vs_exact_service_vector_error':float(np.linalg.norm(actual-ideal)),
        'all_full_event_source_memory_coordinates_retained':True})
    return records,actual


def _matching_rows():
    from model.observer.source_discrimination import circuit
    return tuple(circuit.encode(state)for support in record_code.SUPPORTS for state in support)


@lru_cache(maxsize=1)
def fixed_clock_error():
    return F(continuation.error_bound('C',F(64),4096)['trace_norm_upper_exact'])


def observe_fixed_clock(state,axis,program=None):
    """Two complete records on the actual36-row invariant source, with Q retained.

    The source density is obtained by tracing retained memory/event/reference
    axes only for these local probabilities. No source outcome is selected or
    divided by its success probability. The inherited channel proof retains
    arbitrary untouched references.
    """
    validate_manifest(manifest()if program is None else program)
    if axis not in ('X','Y'):raise ValueError('declared complete X or Y instrument required')
    value=validate_state(state);rows=np.array(_matching_rows())
    amplitudes=value[rows,:6,:,:].reshape(36,-1);rho=amplitudes@amplitudes.conj().T
    outside=value.copy();outside[rows,:6,:,:]=0
    failure=float(np.vdot(outside,outside).real)
    law='C';time=F(64);steps=4096;P=clock.effects(0,axis)
    first=clock.evolve(rho,law,time,steps);joint=np.zeros((3,3))
    for i,p in enumerate(P):
        second=clock.evolve(p@first@p,law,time,steps)
        for j,q in enumerate(P):joint[i,j]=np.trace(q@second).real
    after1=sum(p@first@p for p in P);second=clock.evolve(after1,law,time,steps);after2=sum(p@second@p for p in P)
    E=np.diag(energy_diagonal()[rows]);work=[float(np.trace(E@(x-y)).real)for x,y in ((first,rho),(after1,first),(second,after1),(after2,second))]
    total=float(np.vdot(value,value).real)
    assert joint.min()>-2e-11 and abs(float(joint.sum())+failure-total)<2e-10
    error=fixed_clock_error()
    return {'policy':'fixed_compensated_clock','axis':axis,'joint_probabilities':joint.tolist(),
            'initial_complement_probability':failure,'input_norm_squared':total,'outcomes':10,
            'source_repreparation_between_blocks':False,'normalization_performed':False,
            'source_work_components_diagnostic':work,'inherited_complete_two_time_L1_bound_exact':str(2*error),
            'stored_feedback_and_readout_service_error_budget_additional':True,
            'inherited_bound_hypothesis':'unit-trace initial state in the six matching-code words, with arbitrary untouched reference; not arbitrary36-row input',
            'failure_preinstrument':'P36 source support times valid memory, and its exact complement',
            'event_register_is_retained_spectator':True}


def observe_committed_event(state,program=None):
    """Execute every original raw record, keeping the coherent event as spectator."""
    validate_manifest(manifest()if program is None else program)
    value=validate_state(state);g=magnetic.carrier();codes=np.array(source_addresses())
    amps=value[codes,:6,:,:];outside=value.copy();outside[codes,:6,:,:]=0
    failure=float(np.vdot(outside,outside).real)
    characters=np.stack([memory_character(k)for k in range(6)],axis=1)
    joint=[];work=0.;energies=energy_diagonal()[codes]
    for seam in g.seams:
        for branch in raw_local_branches('mixture'):
            out=np.zeros_like(amps)
            for i,pair in enumerate(g.pairs):
                target=active.target(pair,seam,branch)
                if target is not None:
                    j=g.index[target];chi=active.voltage(pair,target,seam)
                    out[j]=np.roll(amps[i],chi,axis=0)
                    work+=float(branch.weight)/30*(energies[j]-energies[i])*float(np.vdot(amps[i],amps[i]).real)
            char=np.einsum('km,nmer->nker',characters.conj().T,out,optimize=True)
            joint.append(float(branch.weight)/30*np.sum(abs(char)**2,axis=(2,3)))
    joint=np.asarray(joint);total=float(np.vdot(value,value).real)
    assert abs(joint.sum()+failure-total)<2e-10
    return {'policy':'committed_event','raw_records':750,'terminal_source_character_outcomes_per_record':396,
      'complete_outcomes_including_initial_failure':297001,'raw_record_probabilities':joint.sum(axis=(1,2)).tolist(),
      'complete_joint_probability_sha256':sha256(json.dumps(joint.tolist(),separators=(',',':')).encode()).hexdigest(),
      'initial_complement_probability':failure,'input_norm_squared':total,
      'source_work_diagnostic':work,'normalization_performed':False,
      'event_register_role':'retained spectator; original raw coin is not the coherent swap selector'}


def electric_extension_protocol(program=None,descriptor=None):
    """Declared carrier extension; execute actual raw current, not a label join."""
    from model.native_electric_gauss import source as electric
    validate_manifest(manifest()if program is None else program)
    expected=electric.extension_descriptor()
    if descriptor is not None and descriptor!=expected:raise ValueError('undeclared integer extension substitution')
    state=electric.pair_string();rows=[]
    for seam in electric.LOOP:
        events=electric.complete_raw_events(state,seam,limit=4)
        assert len(events)==50 and sum(e['probability']for e in events)==1
        for e in events:
            if e['state'] is not None:assert not any(electric.gauss(e['state']))
        hop=electric.apply_native_current(state,seam,event_bit=1,limit=4)
        rows.append({'seam':list(seam),'all_fifty_raw_probabilities':[str(e['probability'])for e in events],
            'selected_coherent_event_bit':1,'overflow':hop.overflow,
            'source_pair_before':list(state.pair),'source_pair_after':list(hop.state.pair),
            'source_address_before':source_addresses()[magnetic.carrier().index[state.pair]],
            'source_address_after':source_addresses()[magnetic.carrier().index[hop.state.pair]],
            'memory_before':state.memory,'memory_after':hop.state.memory,
            'realized_signed_current':hop.realized_signed_current,
            'electric_work_exact':str(electric.electric_energy(hop.state)-electric.electric_energy(state))})
        state=hop.state
    return {'extension_descriptor':expected,'steps':rows,'final_full_integer_flux':list(state.flux),
      'final_gauss_residual':list(electric.gauss(state)),
      'preparation':'declared neutral pair/string service on added integer registers; separate from matching-clock preparation',
      'raw_probabilities_are_controls':'all50 records are evaluated at each input; the coherent bit1 word is a distinct declared intervention',
      'integer_carrier_is_inside_base_512x8':False,'physical_electric_unit_derived':False}


def compensated_branch_identity():
    g=magnetic.carrier();checked=0;slots=0
    for seam in g.seams:
        for branch in raw_local_branches('mixture'):
            slots+=1
            for pair in g.pairs:
                target=active.target(pair,seam,branch)
                if target is None:continue
                chi=active.voltage(pair,target,seam)
                for memory in range(6):
                    after=(memory+chi)%6;restored=(after-chi)%6
                    assert restored==memory;checked+=1
    return {'all_seam_raw_slots':slots,'nonzero_structural_source_memory_columns':checked,
      'identity':'(R_r tensorI)V_r*V_r(K_r tensorI)=(R_r K_r)tensorI',
      'source_and_memory_complement':'retained initial failure, no repaired discarded branch',
      'scope':'Exact service contract; complete compiled primitive precision is additional.'}


def run():
    p=manifest();initial=prepare(p);spatial,end=observe_spatial(initial,2,p)
    X=observe_fixed_clock(initial,'X',p);Y=observe_fixed_clock(initial,'Y',p)
    committed=observe_committed_event(initial,p)
    back=spatial_step(spatial_step(initial),inverse=True)
    assert np.linalg.norm(back-initial)<2e-12
    matrices=compiled_event_matrices();prep_norm=float(np.vdot(initial,initial).real)
    return {'schema':'oph.single_native_observer_program.results.v1',
      'verdict':'PASS_SHARED_TWO_LOAD_PROGRAM_COMPLETE_RECORDS_AND_NATIVE_REFLECTION',
      'program':p,'program_sha256':descriptor_sha256(),'decoded_H66_basis':[[a,list(pair)]for a,pair in decoded_basis()],
      'shared_input':{'same_full_state_used_in_each_protocol':True,'norm_squared_diagnostic':prep_norm,
          'shape':list(initial.shape),'preparation_error_and_export_bound':complete_bound(len(plan(60)['preparation']),F(6,1<<41)),
          'event_padding_mass_diagnostic':float(np.sum(abs(initial[:,:,60:,:])**2))},
      'complete_protocols':{'coherent_spatial':spatial,'fixed_clock_X':X,'fixed_clock_Y':Y,'committed_event':committed},
      'integer_extension_protocol':electric_extension_protocol(p),
      'explicit_idle_apparatus_mode':idle_mode_control(p),
      'compensated_branch_identity':compensated_branch_identity(),
      'inverse_full_state_error_diagnostic':float(np.linalg.norm(back-initial)),
      'stored_inverse_execution':{'implementation':'actual emitted chronological inverse native word, with negative requested angles',
          'all64_event_columns_executed':True,'additional_forward_inverse_operator_error_diagnostic':float(np.linalg.norm(compiled_inverse_reflection()@compiled_event_matrices()[1]-np.eye(64))),
          'full_inverse_word_bound':complete_bound(len(plan(60)['reflection']),F(13,1<<41))},
      'event_reflection_compilation':run_reflections(),
      'physical_boundaries':{'one_local_two_load_construction_executed':True,
          'moving_clock_bundle_implemented_by_this_package':False,'local_autonomous_clock_proved':False,
          'physical_rest_masses_or_electric_units_derived':False,'continuum_or_gravity_limit':False,
          'official_OpenWave_awards':0,'all_observables_reproduced':False}}
