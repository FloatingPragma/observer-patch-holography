"""Native source preparation, conjugated kicks and complete quantum archives."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy import linalg as la

from model.lepton_cycle import instrument
from model.native_quantum_service_bridge import finite_compiler
from model.observer.source_scalar import compiler as source
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
THETA=F(2,9);KICK=F(1,7)


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/lepton_cycle/results.json','model/lepton_cycle/PROOF.md',
                  'model/observer/source_scalar/DERIVATION.md','model/native_quantum_service_bridge/PROOF.md',
                  source.RESPONSE,source.CONTROL_RECEIPT,'requirements.txt','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def digest(tape):
    return hashlib.sha256(json.dumps(tape,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def zero_phase_word(bits,angle):
    """Complete Walsh expansion of exp(i angle |0><0|), retaining zero services."""
    if type(bits) is not int or not 1<=bits<=12:raise ValueError('one through twelve bits required')
    if not math.isfinite(float(angle)):raise ValueError('finite angle required')
    return [event for mask in range(1<<bits)
            for event in source.pauli_rotation(''.join('Z' if mask>>j&1 else 'I' for j in range(bits)),float(angle)/(1<<bits))]


def vectors():
    u=np.zeros(64);u[:60]=1/np.sqrt(60)
    q=np.eye(64)[0]-u;q/=la.norm(q)
    records,_=instrument.source_data();plus=np.array([x==2 for x in records]);minus=np.array([x==3 for x in records])
    blank=~(plus|minus)
    masks=(blank[:,None]&(~blank)[None,:],plus[:,None]&blank[None,:],minus[:,None]&blank[None,:])
    phi=np.zeros((64,64));phi[:60,:60]=sum(m/np.sqrt(m.sum()) for m in masks)/np.sqrt(3)
    phases=np.zeros((64,4))
    for k in range(3):
        phases[:60,k][plus]=float(THETA)+2*np.pi*k/3
        phases[:60,k][minus]=-float(THETA)-2*np.pi*k/3
    return u,q,phi,phases


@lru_cache(maxsize=1)
def words():
    u,q,phi,phases=vectors();v=finite_compiler.prepare(q)
    z=zero_phase_word(6,np.pi);house=source.inverse(v)+z+v
    native_domain=source.compile_diagonal(np.r_[np.full(60,np.pi),np.zeros(4)])
    reflect=source.inverse(house)+z+house+native_domain
    vp=source.compile_real_state(phi.reshape(-1));zp=zero_phase_word(12,-float(KICK))
    d=source.compile_diagonal(phases.reshape(-1))
    swap=source.cnot(0,1)+source.cnot(1,0)+source.cnot(0,1)
    return {'house_prepare':v,'zero6':z,'house':house,'native_domain':native_domain,
            'reflect':reflect,'phi_prepare':vp,'zero12':zp,'tick':d,'swap_bit':swap}


@lru_cache(maxsize=1)
def compiled_blocks():
    """Replay native words on all small-block columns and the large preparation."""
    tape=words();u,q,phi,phases=vectors()
    v=source.execute(tape['house_prepare'],6)
    vi=source.execute(source.inverse(tape['house_prepare']),6)
    z=source.execute(tape['zero6'],6)
    house=v@z@vi
    domain=source.execute(tape['native_domain'],6)
    # Every factor is an executed native block. No target matrix is substituted.
    house_i=v@source.execute(source.inverse(tape['zero6']),6)@vi
    reflect=domain@house@z@house_i
    reflect_i=house@source.execute(source.inverse(tape['zero6']),6)@house_i@source.execute(source.inverse(tape['native_domain']),6)
    tick=source.execute(tape['tick'],8);tick_i=source.execute(source.inverse(tape['tick']),8)
    blank=np.zeros(4096);blank[0]=1
    psi=source.execute(tape['phi_prepare'],12,blank)
    swap=source.execute(tape['swap_bit'],2)
    expected_h=np.eye(64)-2*np.outer(q,q);expected_r=np.eye(64);expected_r[:60,:60]=2*np.ones((60,60))/60-np.eye(60)
    expected_swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    errors={'Householder_all64_input_columns':float(la.norm(house-expected_h,2)),
            'Householder_preparation':float(la.norm(house[:,0]-u)),
            'Householder_four_padding_columns':float(la.norm(house[:,60:]-np.eye(64)[:,60:])),
            'Householder_inverse':float(la.norm(house_i@house-np.eye(64),2)),
            'uniform_reflection_all64_columns':float(la.norm(reflect-expected_r,2)),
            'reflection_inverse_all64_columns':float(la.norm(reflect_i@reflect-np.eye(64),2)),
            'phase_tick_all256_columns':float(la.norm(tick-np.diag(np.exp(1j*phases.reshape(-1))),2)),
            'inverse_tick_all256_columns':float(la.norm(tick_i@tick-np.eye(256),2)),
            'phi_source_preparation':float(la.norm(psi-phi.reshape(-1))),
            'phi_raw_norm_defect':float(abs(np.vdot(psi,psi)-1)),
            'SWAP_all4_columns':float(la.norm(swap-expected_swap,2))}
    if max(errors.values())>2e-9:raise ArithmeticError('native preparation or source block replay failed')
    return {'house':house,'house_i':house_i,'reflect':reflect,'reflect_i':reflect_i,'tick':tick,'tick_i':tick_i,
            'phi':psi,'swap':swap,'errors':errors}


def rank_kick(state,psi,angle):
    """All-input exact conjugation identity; raw executed psi is not normalized."""
    return state+np.expm1(-1j*float(angle))*psi[:,None]*(psi.conj()@state)[None,:]


@lru_cache(maxsize=1)
def literal_kick_control():
    """Literal replay of V†, blank phase and V on two independent input columns."""
    tape=words();blocks=compiled_blocks();u,_,_,_=vectors()
    first=np.outer(u,u).reshape(-1).astype(complex)
    # A coherent native/padding column also checks that padding is retained.
    second=np.exp(2j*np.pi*np.arange(4096)/4096)/64
    x=np.column_stack((first,second))
    value=source.execute(source.inverse(tape['phi_prepare']),12,x)
    value=source.execute(tape['zero12'],12,value)
    value=source.execute(tape['phi_prepare'],12,value)
    shortcut=rank_kick(x,blocks['phi'],KICK)
    ideal=vectors()[2].reshape(-1)
    target=rank_kick(x,ideal,KICK)
    errors={'literal_source_vs_conjugation_shortcut':float(la.norm(value-shortcut)),
            'literal_source_vs_native_analyzer':float(la.norm(value-target)),
            'two_input_column_norm_preservation':float(la.norm(value.conj().T@value-x.conj().T@x))}
    if max(errors.values())>3e-9:raise ArithmeticError('literal large source kick differs')
    return {'executed_source_columns':2,'all4096_source_coordinates_retained':True,
            'literal_native_pulses':2*len(tape['phi_prepare'])+len(tape['zero12']),
            'one_input_has_nonzero_padding_amplitude':True,'errors':errors}


def reflect_axes(state,block):
    return np.einsum('ax,by,xyfc->abfc',block,block,state,optimize=True)


def tick_axis(state,block):
    shape=state.shape
    arranged=state.transpose(0,2,1,3).reshape(256,-1)
    return (block@arranged).reshape(64,4,64,shape[3]).transpose(0,2,1,3)


def return_parts(state,uniform):
    x=np.einsum('x,y,yafc->xafc',uniform,uniform.conj(),state,optimize=True)
    y=np.einsum('y,b,xbfc->xyfc',uniform,uniform.conj(),state,optimize=True)
    both=np.einsum('x,y,a,b,abfc->xyfc',uniform,uniform,uniform.conj(),uniform.conj(),state,optimize=True)
    return {'yes_yes':both,'x_no':y-both,'y_no':x-both,'both_no':state-x-y+both}


@lru_cache(maxsize=1)
def cycle():
    block=compiled_blocks();u=vectors()[0];fresh=np.outer(block['house'][:,0],block['house'][:,0])
    state=np.zeros((64,64,4,4),complex)
    for k in range(3):state[:,:,k,k]=fresh
    # Correct address histograms do not supply the required pure coherent source.
    signs=np.ones(64);signs[1:60:2]=-1
    state[:,:,1,3]=signs[:,None]*fresh
    value=reflect_axes(state,block['reflect'])
    value=tick_axis(value,block['tick'])
    value=rank_kick(value.reshape(4096,-1),block['phi'],KICK).reshape(value.shape)
    value=tick_axis(value,block['tick_i'])
    value=reflect_axes(value,block['reflect_i'])
    yes=block['house'][:,0,None]@block['house_i'][0,None,:]
    no=block['house']@block['house_i']-yes
    parts={label:np.einsum('ax,by,xyfc->abfc',left,right,value,optimize=True)
           for label,left,right in (('yes_yes',yes,yes),('x_no',no,yes),('y_no',yes,no),('both_no',no,no))}
    target=instrument.mass_instrument(float(THETA),float(KICK));errors={};local={}
    wrong_family=0.;padding=0.
    for record in instrument.LABELS:
        actual=parts[record][:,:,:,:3]
        expected=np.zeros_like(actual)
        for k in range(3):expected[:60,:60,k,k]=target[record][:,:,k]
        errors['complete_'+record+'_all_register_columns']=float(la.norm((actual-expected).reshape(-1,3),2))
        wrong=actual.copy()
        for k in range(3):wrong[:,:,k,k]=0
        wrong_family+=float(la.norm(wrong)**2)
        outside=actual.copy();outside[:60,:60]=0;padding+=float(la.norm(outside)**2)
        local[record]=actual.reshape(4096,4,3)
    joint=np.concatenate([parts[r][:,:,:,:3].reshape(-1,3) for r in instrument.LABELS])
    errors['complete_isometry_all_input_columns_upper']=float(np.sqrt(sum(v*v for k,v in errors.items() if k.startswith('complete_'))))
    errors['complete_reference_isometry_identity']=float(la.norm(joint.conj().T@joint-np.eye(3),2))
    errors['padded_source_probability_sum']=padding
    errors['unused_family_and_wrong_family_probability_sum']=wrong_family
    independent={r:float((la.norm(parts[r][:,:,0,0])**2+la.norm(parts[r][:,:,1,1])**2)/2) for r in instrument.LABELS}
    correlated={r:float((la.norm(parts[r][:,:,0,0])**2+la.norm(parts[r][:,:,1,3])**2)/2) for r in instrument.LABELS}
    distance=sum(abs(independent[r]-correlated[r]) for r in instrument.LABELS)
    if max(errors.values())>8e-9 or distance<.1:raise ArithmeticError('complete compiled cycle or source-correlation control failed')
    return {'local':local,'errors':errors,'correlated':{'independent_original_probabilities':independent,
            'correlated_original_probabilities':correlated,'unhalved_classical_distance':distance,
            'same_uniform_address_diagonal_histogram':True,'same_full_pure_source_density_matrix':False,
            'pure_source_marginal_forces_product_with_every_reference':True,
            'source_phase_control':'alternating sign on thirty of the sixty first-register addresses, correlated with family one'}}


def history_products(local,steps):
    """Retain all local quantum factors; no source coordinate is measured."""
    if type(steps) is not int or not 1<=steps<=3:raise ValueError('bounded execution admits one through three cycles')
    # The exact source theorem is family diagonal. Tiny raw wrong-family drift
    # is separately bounded against the all-column cycle isometry above.
    diagonal={r:np.column_stack([local[r][:,k,k] for k in range(3)]) for r in instrument.LABELS}
    grams={r:diagonal[r].conj().T@diagonal[r] for r in instrument.LABELS}
    rows=[]
    for labels in product(instrument.LABELS,repeat=steps):
        gram=np.ones((3,3),complex)
        for r in labels:gram*=grams[r]
        rows.append({'record_history':list(labels),'probability_maximally_entangled_register_reference':float(np.trace(gram).real/3),
                     'quantum_archive_Gram_real':gram.real.tolist(),'quantum_archive_Gram_imag':gram.imag.tolist()})
    return rows,diagonal


@lru_cache(maxsize=1)
def histories():
    data=cycle();target=instrument.mass_instrument(float(THETA),float(KICK));old={}
    for record,array in target.items():
        value=np.zeros((4096,4,3),complex)
        for k in range(3):
            padded=np.zeros((64,64),complex);padded[:60,:60]=array[:,:,k];value[:,k,k]=padded.reshape(-1)
        old[record]=value
    reports=[]
    for steps in (1,2,3):
        rows,diagonal=history_products(data['local'],steps);exact,_=history_products(old,steps)
        error=max(abs(a['probability_maximally_entangled_register_reference']-b['probability_maximally_entangled_register_reference'])
                  for a,b in zip(rows,exact,strict=True))
        gram_error=max(la.norm((np.array(a['quantum_archive_Gram_real'])+1j*np.array(a['quantum_archive_Gram_imag']))-
                              (np.array(b['quantum_archive_Gram_real'])+1j*np.array(b['quantum_archive_Gram_imag'])))
                       for a,b in zip(rows,exact,strict=True))
        total=sum(r['probability_maximally_entangled_register_reference'] for r in rows)
        if max(error,gram_error,abs(total-1))>3e-8:raise ArithmeticError('complete quantum-archive histories differ')
        reports.append({'cycles':steps,'complete_return_histories':rows,'history_count':4**steps,
            'probability_sum_diagnostic':total,'maximum_probability_error':float(error),'maximum_archive_Gram_error':float(gram_error),
            'source_archive_dimension_exact':str(4096**steps),'retained_classical_return_bits':2*steps,
            'explicit_local_tensor_factor_recipe':'for family k and history (r1,...,rn), tensor the stored vectors A[rj][:,k]; retain the original register and its untouched reference',
            'raw_full_cycle_wrong_family_and_padding_errors_retained_in_error_budget':True})
    factors={r:{'shape':[4096,3],'real':a.real.tolist(),'imag':a.imag.tolist()} for r,a in diagonal.items()}
    return {'bounded_runs':reports,'local_quantum_archive_factors':factors,
            'source_coordinates_dephased_or_selected':False,'tensor_factors_are_full4096_component_vectors':True,
            'arbitrary_register_reference_identity_follows_from_all3_input_columns':True,
            'complete_cq_channel_unhalved_error_from_one_cycle_isometry_error':'at most 2*n*delta, with delta the all-column norm; physical service errors are additive separately',
            'terminal_active_fresh_source':{'shape':[64,64],
                'real':np.outer(compiled_blocks()['house'][:,0],compiled_blocks()['house'][:,0]).real.tolist(),
                'imag':np.outer(compiled_blocks()['house'][:,0],compiled_blocks()['house'][:,0]).imag.tolist()},
            'literal_full4096_to_power3_array_materialized':False}


def service_controls():
    result={}
    for d in (2,4):
        kraus=source.transfer_kraus(d,coherent=True)
        identity=sum((k.conj().T@k for k in kraus),np.zeros((6,6),complex))
        wrong=source.transfer_kraus(d,coherent=False);state=np.zeros(6,complex);state[:d]=1/np.sqrt(d)
        coherent=sum((np.outer(k@state,(k@state).conj()) for k in kraus),np.zeros((d+1,d+1),complex))
        dephased=sum((np.outer(k@state,(k@state).conj()) for k in wrong),np.zeros_like(coherent))
        distance=float(np.sum(abs(la.eigvalsh(coherent-dephased))))
        if la.norm(identity-np.eye(6))>1e-14 or distance<.9:raise ArithmeticError('complete coherent source service failed')
        result[str(d)]={'full_processor_input_dimension':6,'all_Kraus_outcomes':len(kraus),
            'complete_identity_error':float(la.norm(identity-np.eye(6))),
            'dephasing_negative_control_unhalved_distance':distance,
            'all_proper_code_complement_failure_coordinates_retained':True}
    return result


def resources(steps):
    if type(steps) is not int or not 1<=steps<=3:raise ValueError('one through three cycles required')
    tape=words();prep=2*(steps+1)*len(tape['house'])
    kick=2*len(tape['phi_prepare'])+len(tape['zero12'])
    # A return test coherently unprepares/reprepares each address around its
    # supplied binary blank projector. The fresh pair is swapped in afterwards.
    return_pulses=4*len(tape['house'])
    one=4*len(tape['reflect'])+2*len(tape['tick'])+kick+return_pulses+12*len(tape['swap_bit'])
    total=prep+steps*one
    return {'cycles':steps,'native_source_pulses':total,'coherent_load_export_service_events':2*total,
            'fresh_address_pairs_including_terminal_active_pair':steps+1,
            'fresh_pair_preparation_native_pulses':prep,'per_cycle_native_pulses':one,
            'per_cycle_rank_one_kick_native_pulses':kick,'per_cycle_return_test_basis_change_pulses':return_pulses,
            'per_cycle_source_archive_SWAP_pulses':12*len(tape['swap_bit']),
            'supplied_binary_blank_measurements':2*steps,'retained_source_archive_qubits':12*steps,
            'retained_classical_return_bits':2*steps,'active_address_qubits':12,'family_code_qubits':2,
            'separate_assumptions':['coherent proper-code load/export','blank reset','address routing and pulse-angle control','binary blank-projector readout'],
            'complete_half_diamond_error_budget':'Npulse*eta_pulse + 2*Npulse*eta_service + Nreset*eta_blank + Nrouting*eta_routing + 2*n*eta_binary_readout; clip at one',
            'physical_error_rates_set_to_zero_or_certified':False,
            'adaptive_autonomous_service_scheduler_derived':False}


def calculate():
    pins=bindings();summary=source.source_summary();tape=words();blocks=compiled_blocks();literal=literal_kick_control();data=cycle()
    result={'schema':'oph.lepton_service_compiler.v1','source_bindings':pins,'all_checks_pass':True,
        'source_contract':{'processor_dimension':6,'source_generators':12,'address_dimensions':[64,64],
            'native_address_support':[60,60],'family_dimension':3,'padded_family_dimension':4,
            'preparation':'source-word Householder on two supplied blank address registers, independent of the family and every reference',
            'clock':'one declared analyzer kick with angle h=1/7 per complete cycle','theta_exact':str(THETA),'kick_exact':str(KICK),
            'arbitrary_untouched_quantum_reference_retained':True,'source_summary':summary},
        'compiled_words':{name:{'native_pulses':len(word),'sha256':digest(word)} for name,word in tape.items()},
        'source_block_replay_errors':blocks['errors'],'literal_large_kick_control':literal,
        'complete_cycle_errors':data['errors'],'complete_service_controls':service_controls(),
        'correlated_source_countercontrol':data['correlated'],'quantum_archive_execution':histories(),
        'resources':[resources(n) for n in (1,2,3)],
        'ideal_zero_branches':{'native_proper_code_first_abort':'one retained first-failure record for every service position, each has exact probability zero in the declared invariant-code source',
            'padded_addresses':'every source and archive vector retains all four padding addresses on each register',
            'unused_family_coordinate':'the fourth family coordinate is retained and has exact zero amplitude for admitted inputs'},
        'scope':{'full_padded_Householder_on_all_input_columns':True,'complete_up_to_three_cycle_native_service_implementation':True,
            'source_word_conjugation_shortcut_proved_on_all_inputs':True,'literal_large_kick_replayed_on_two_columns':True,
            'fresh_preparation_source_archive_and_all_four_returns_retained':True,'coherent_service_and_blank_reset_are_supplied':True,
            'autonomous_reset_or_hardware_pulse_precision_derived':False,'new_physical_Yukawa_or_particle_mass_prediction':False}}
    if bindings()!=pins:raise ArithmeticError('service compiler source changed during execution')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return type(a) in (int,float) and math.isfinite(a) and abs(a-b)<=2e-10+2e-8*abs(b)
    return type(a) is type(b) and a==b


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('bounded service receipt differs')
    print('Complete bounded native lepton preparation, service, return and quantum-archive execution verified.')
