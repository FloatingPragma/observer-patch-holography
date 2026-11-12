"""Native source spans, certified ranks, streaming unitaries and full histories."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy import linalg as la
import sympy as sy

from model.lepton_cycle import instrument,composition
from model.lepton_mass_compiler import analyzer
from model.native_flavor_orbit import intervals as iv
from model.observer.source_scalar import compiler as native
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
THETA=F(2,9);TIME=F(1,7);LABELS=instrument.LABELS
ACTIVE=('x_no','both_no')


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/lepton_cycle/results.json','model/lepton_cycle/PROOF.md',
        'model/lepton_mass_compiler/results.json','model/lepton_mass_compiler/PROOF.md',
        'model/observer/source_scalar/DERIVATION.md',native.RESPONSE,native.CONTROL_RECEIPT,
        'requirements.txt','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def array_json(value):
    value=np.asarray(value);return {'shape':list(value.shape),'real':value.real.tolist(),'imag':value.imag.tolist()}


@lru_cache(maxsize=1)
def local_frames():
    records,_=instrument.source_data();b=np.array([r not in (2,3) for r in records],float)/np.sqrt(50)
    p=np.array([r==2 for r in records],float)/np.sqrt(5);m=np.array([r==3 for r in records],float)/np.sqrt(5)
    d=(p+m)/np.sqrt(2);c=(p-m)/np.sqrt(2);s=np.sqrt(5/6);t=np.sqrt(1/6)
    u=s*b+t*d;v=t*b-s*d
    frame=np.column_stack((u,v,c))
    if la.norm(frame.T@frame-np.eye(3))>2e-14:raise ArithmeticError('native category frame is not orthonormal')
    bases={'yes_yes':np.kron(u,u)[:,None],
        'x_no':np.column_stack((np.kron(v,u),1j*np.kron(c,u))),
        'y_no':np.kron(u,v)[:,None],
        'both_no':np.column_stack((np.kron(v,v),1j*np.kron(c,v)))}
    angles=float(THETA)+2*np.pi*np.arange(3)/3;cs=np.cos(angles);sn=np.sin(angles)
    lam=np.sqrt(5/108)+2*np.sqrt(5/216)*cs
    a=np.vstack(((np.sqrt(2)*s*s*cs-t*t)/np.sqrt(3),np.sqrt(2/3)*s*sn))
    y=(s*s-np.sqrt(2)*t*t*cs)/np.sqrt(3)
    z=np.vstack((-s*t*(1+np.sqrt(2)*cs)/np.sqrt(3),-np.sqrt(2/3)*t*sn))
    geometry={};lengths={}
    for label,value in [('x_no',a),('both_no',z)]:
        lengths[label]=la.norm(value,axis=0);geometry[label]=value/lengths[label]
    return dict(frame=frame,bases=bases,lambda_values=lam,rank_two_geometry=geometry,
        rank_two_lengths=lengths,y_geometry=y,unscaled_rank_two={'x_no':a,'both_no':z})


@lru_cache(maxsize=1)
def exact_rank_certificate():
    c,s=sy.symbols('c s',real=True);root=sy.sqrt
    cs=[c,-c/2-root(3)*s/2,-c/2+root(3)*s/2]
    sn=[s,-s/2+root(3)*c/2,-s/2-root(3)*c/2]
    geometry={'x_no':sy.Matrix([[root(6)*sy.Rational(5,18)*x-root(3)/18 for x in cs],
                               [root(5)*x/3 for x in sn]]),
        'both_no':sy.Matrix([[-root(15)*(1+root(2)*x)/18 for x in cs],[-x/3 for x in sn]])}
    sine=iv.sine(iv.I(THETA));cosine=1-2*iv.sine(iv.I(THETA)/2)**2
    substitutions={c:cosine,s:sine};witnesses=[]
    for label,matrix in geometry.items():
        for i,j in combinations(range(3),2):
            determinant=sy.expand(matrix[:,[i,j]].det())
            determinant=sy.Poly(determinant,s).rem(sy.Poly(s*s+c*c-1,s)).as_expr().expand()
            enclosure=iv.rounded_out(iv.expression(determinant,substitutions),36)
            if enclosure.lo<=0<=enclosure.hi:raise ArithmeticError('projective family separation is not certified')
            witnesses.append({'outcome':label,'family_pair':[i,j],'exact_minor':str(determinant),
                              'rational_minor_interval':enclosure.strings()})
    lam=analyzer.rational_effect_intervals()['lambda_intervals_k0_k1_k2']
    if any(F(lo)<=0 for lo,hi in lam):raise ArithmeticError('native family amplitude vanished')
    return {'cos_theta_interval':iv.rounded_out(cosine,36).strings(),'sin_theta_interval':iv.rounded_out(sine,36).strings(),
        'formal_cosine_symbol':'c','formal_sine_symbol':'s','polynomial_relation':'c²+s²=1',
        'unscaled_real_rank_two_columns':{k:[[str(x) for x in row] for row in v.tolist()] for k,v in geometry.items()},
        'all_six_nonzero_pair_minors':witnesses,'positive_lambda_intervals':lam,
        'y_no_nonzero_scalar_lower':'(5-sqrt(2))/(6*sqrt(3))>0',
        'yes_yes_nonzero_modulus_lower':'1-2*max(e_k)>4/9',
        'rank_rule':'1 if zero rank-two outcomes; 2 if exactly one; 3 if at least two',
        'rank_is_selected_by_floating_threshold':False}


def rank_for(history):
    count=sum(r in ACTIVE for r in history)
    return min(3,count+1)


def local_columns(steps):
    if type(steps) is not int or steps not in (1,2,4):raise ValueError('one, two or four complete rounds required')
    d=local_frames();z=np.expm1(-1j*float(TIME)/steps);lam=d['lambda_values']
    scalar={'yes_yes':1+z*lam*lam,'y_no':z*lam*d['y_geometry']}
    scalar.update({r:z*lam*d['rank_two_lengths'][r] for r in ACTIVE})
    geometry={r:d['rank_two_geometry'][r] if r in ACTIVE else np.ones((1,3)) for r in LABELS}
    columns={r:geometry[r]*scalar[r] for r in LABELS}
    original=instrument.mass_instrument(float(THETA),float(TIME)/steps)
    errors={r:float(la.norm(d['bases'][r]@columns[r]-original[r].reshape(3600,3))) for r in LABELS}
    if max(errors.values())>2e-12:raise ArithmeticError('actual native archive does not equal its explicit local span')
    return dict(geometry=geometry,scalar=scalar,columns=columns,original=original,errors=errors)


def complete_frame(columns,rank):
    """A full orthogonal basis; rank comes from the exact minor theorem."""
    columns=np.asarray(columns,float);dimension=columns.shape[0];basis=[]
    candidates=[columns[:,k]/la.norm(columns[:,k]) for k in range(rank)]+list(np.eye(dimension))
    for value in candidates:
        vector=value.copy()
        for _ in range(2):
            for previous in basis:vector-=previous*np.dot(previous,vector)
        norm=la.norm(vector)
        if norm>1e-10:basis.append(vector/norm)
        if len(basis)==dimension:break
    if len(basis)!=dimension:raise ArithmeticError('unitary complement construction failed')
    frame=np.column_stack(basis)
    if la.norm(frame.T@frame-np.eye(dimension))>3e-13 or la.norm(frame[:,rank:].T@columns)>3e-12:
        raise ArithmeticError('certified active columns were lost')
    return frame


@lru_cache(maxsize=None)
def streaming_sequence(sequence):
    if any(r not in ACTIVE for r in sequence) or not 1<=len(sequence)<=4:raise ValueError('one through four rank-two outcomes required')
    first=local_frames()['rank_two_geometry'][sequence[0]]
    if len(sequence)==1:return {'columns':first.copy(),'original':first.copy(),'embedding':np.eye(2),'unitary':None}
    previous=streaming_sequence(sequence[:-1]);current=local_frames()['rank_two_geometry'][sequence[-1]]
    before=np.column_stack([np.kron(previous['columns'][:,k],current[:,k]) for k in range(3)])
    frame=complete_frame(before,3);dimension=len(frame)
    destinations=[0,1,2] if dimension==4 else [0,2,4]
    others=[j for j in range(dimension) if j not in destinations]
    target=np.eye(dimension)[:,destinations+others];unitary=target@frame.T
    output=unitary@before
    columns=output if dimension==4 else output.reshape(4,2,3)[:,0,:]
    blank_error=la.norm(output[3:]) if dimension==4 else la.norm(output.reshape(4,2,3)[:,1,:])
    original=np.column_stack([np.kron(previous['original'][:,k],current[:,k]) for k in range(3)])
    embedding=original@la.inv(columns[:3,:])
    errors={'full_orthogonal_unitary':float(la.norm(unitary.T@unitary-np.eye(dimension))),
            'work_register_blank':float(blank_error),'full_archive_embedding':float(la.norm(embedding.T@embedding-np.eye(3))),
            'all3_family_columns_reconstructed':float(la.norm(embedding@columns[:3]-original))}
    if max(errors.values())>2e-11:raise ArithmeticError('streaming archive unitary or reconstruction failed')
    return {'columns':columns,'original':original,'embedding':embedding,'unitary':unitary,'errors':errors,
            'input_dimension':dimension,'output_active_rank':3,'fixed_blank_work_qubits':int(dimension==8)}


def native_small_services():
    rows=[]
    for length in (2,3,4):
        for sequence in product(ACTIVE,repeat=length):
            data=streaming_sequence(sequence);matrix=data['unitary'];bits=len(matrix).bit_length()-1
            tape=native.compile_real_orthogonal(matrix);actual=native.execute(tape,bits)
            error=float(la.norm(actual-matrix,2))
            if error>3e-10:raise ArithmeticError('native small compression service failed')
            rows.append({'rank_two_history':list(sequence),'input_dimension':len(matrix),'rank':3,
                'native_pulses':len(tape),'coherent_load_export_events':2*len(tape),
                'all_input_column_error_diagnostic':error,'unitary':array_json(matrix),
                'complete_native_word':tape,'reconstruction_errors':data['errors']})
    return rows


def native_address_service():
    """Actual native preparation columns with exact Householder conjugation."""
    target=np.zeros((64,3));target[:60]=local_frames()['frame'];working=target.copy();exact=np.eye(64);actual=np.eye(64)
    zword=native.compile_diagonal(np.r_[np.pi,np.zeros(63)]);blank=np.eye(64)[:,0];rows=[]
    for j in range(3):
        q=working[:,j]-np.eye(64)[:,j];q/=la.norm(q)
        tape=native.compile_real_state(q);raw=native.execute(tape,6,blank)
        house=np.eye(64)-2*np.outer(q,q);rawhouse=np.eye(64)+np.expm1(1j*np.pi)*np.outer(raw,raw.conj())
        exact=house@exact;actual=rawhouse@actual;working=house@working
        row={'column':j,'preparation_pulses':len(tape),'Householder_pulses':2*len(tape)+len(zword),
            'native_prepared_column_error':float(la.norm(raw-q)),'raw_prepared_norm_defect':float(abs(np.vdot(raw,raw)-1)),
            'preparation_word_sha256':hashlib.sha256(json.dumps(tape,separators=(',',':')).encode()).hexdigest()}
        if j==0:
            literal=native.execute(native.inverse(tape)+zword+tape,6)
            row['literal_all64_Householder_input_columns_error']=float(la.norm(literal-house,2))
            row['literal_vs_conjugation_all64_error']=float(la.norm(literal-rawhouse,2))
        rows.append(row)
    phase=native.compile_diagonal(np.r_[0.,0.,-np.pi/2,np.zeros(61)])
    ph=native.execute(phase,6);permutation=np.eye(4)[[1,2,0,3],:]
    route_word=native.compile_real_orthogonal(permutation);route=native.execute(route_word,2)
    flip_word=[native.pulse('X',np.pi/2,[0]),native.pulse('I',-np.pi/2,[0])]
    flip=native.execute(flip_word,1)
    errors={'source_category_all3_columns':float(la.norm(actual@target-np.eye(64)[:,:3])),
        'full64_unitary_diagnostic':float(la.norm(actual.conj().T@actual-np.eye(64),2)),
        'full64_target_unitary_error':float(la.norm(actual-exact,2)),
        'all4_padding_columns':float(la.norm(actual[:,60:]-np.eye(64)[:,60:])),
        'native_phase_all64_columns':float(la.norm(ph-np.diag(np.r_[1.,1.,-1j,np.ones(61)]),2)),
        'conditional_rank2_routing_all4_columns':float(la.norm(route-permutation,2)),
        'conditional_second_return_flip_all2_columns':float(la.norm(flip-np.array([[0,1],[1,0]]),2))}
    if max(errors.values())>3e-10:raise ArithmeticError('native address compression service differs')
    return {'Householder_columns':rows,'left_and_right_address_basis_pulses':2*sum(r['Householder_pulses'] for r in rows)+len(phase),
        'conditional_rank2_routing_pulses':len(route_word),'conditional_second_return_bit_flip_pulses':2,
        'all_original_source_padding_coordinates_retained':True,'errors':errors,
        'ideal_real_angle_recipe_proves_exact_target':True,'binary64_tape_or_roundoff_is_certified_exact_target':False,
        'blank_recycling_on_admitted_span_is_unitary':True,'autonomous_thermal_reset_or_work_budget_derived':False}


def history_run(steps):
    d=local_columns(steps);raw_grams={r:v.reshape(-1,3).conj().T@v.reshape(-1,3) for r,v in d['original'].items()}
    rows=[];ranks={1:0,2:0,3:0};complete=np.zeros((3,3),complex);diagnostic_misclassified=0
    wrong_collapsed=0.;wrong_orthogonal=0.;maximum=0.
    for history in product(LABELS,repeat=steps):
        active=tuple(r for r in history if r in ACTIVE);rank=rank_for(history);ranks[rank]+=1
        scalar=np.prod(np.vstack([d['scalar'][r] for r in history]),axis=0)
        stream=streaming_sequence(active) if active else {'columns':np.ones((1,3)),'original':np.ones((1,3)),'embedding':np.eye(1)}
        columns=stream['columns'][:rank]*scalar;original=stream['original']*scalar
        reconstructed=stream['embedding']@columns
        error=float(la.norm(reconstructed-original));maximum=max(maximum,error)
        relative=max(float(la.norm(reconstructed[:,k]-original[:,k])/la.norm(original[:,k])) for k in range(3))
        if relative>3e-11:raise ArithmeticError('a rare family archive was lost in compression')
        gram=columns.conj().T@columns;expected=np.ones((3,3),complex)
        for record in history:expected*=raw_grams[record]
        if la.norm(gram-expected)>3e-12 or error>3e-12:raise ArithmeticError('an original history or reference coherence was changed')
        complete+=gram;probability=float(np.trace(gram).real/3)
        eigen=np.linalg.eigvalsh(gram);naive=int(np.sum(eigen>max(eigen)*1e-10));diagnostic_misclassified+=naive!=rank
        lengths=np.sqrt(np.maximum(np.diag(gram).real,0))
        if rank>1:
            phases=np.exp(1j*np.angle(columns[0]));false_columns=lengths*phases
            wrong_collapsed+=float(np.sum(abs(la.eigvalsh((np.outer(false_columns.conj(),false_columns)-gram)/3))))
        wrong_orthogonal+=float(np.sum(abs(la.eigvalsh((np.diag(lengths*lengths)-gram)/3))))
        rows.append({'original_history':list(history),'rank_two_subsequence':list(active),'minimal_quantum_rank':rank,
            'minimal_binary_quantum_qubits':math.ceil(math.log2(rank)),'original_probability_maximally_entangled_input':probability,
            'compressed_environment_columns':array_json(columns),'quantum_archive_decoder':array_json(stream['embedding']),
            'raw_floating_relative_eigenvalue_rank_control':naive,'reference_extended_column_reconstruction_error':error,'maximum_relative_family_column_error':relative})
    expected=sum(raw_grams.values())**steps
    errors={'complete_original_history_sum':float(la.norm(complete-expected)),
        'trace_preservation_all_family_inputs':float(np.max(abs(np.diag(complete)-1))),
        'all_original_archive_columns_reconstructed':maximum}
    if max(errors.values())>5e-12 or wrong_collapsed<1e-7 or wrong_orthogonal<.1:
        raise ArithmeticError('complete history or lossy countercontrol failed')
    return {'rounds':steps,'total_kick_time_exact':str(TIME),'per_round_angle_exact':str(TIME/steps),
        'all_original_histories':rows,'rank_counts':{str(k):v for k,v in ranks.items()},
        'original_raw_source_environment_dimension_exact':str(3600**steps),'padded_archive_qubits_without_compression':12*steps,
        'retained_classical_original_record_bits':2*steps,'worst_case_quantum_archive_qubits':2 if steps>1 else 1,
        'complete_probability_sum':sum(r['original_probability_maximally_entangled_input'] for r in rows),
        'wrong_single_environment_vector_unhalved_reference_output_distance':wrong_collapsed,
        'wrong_orthogonalization_unhalved_reference_output_distance':wrong_orthogonal,
        'relative_eigenvalue_threshold_misclassified_histories':diagnostic_misclassified,'errors':errors}


def calculate():
    pins=bindings();rank=exact_rank_certificate();runs=[history_run(n) for n in (1,2,4)]
    local=native_address_service();services=native_small_services();by_sequence={tuple(r['rank_two_history']):r['native_pulses'] for r in services}
    resources=[]
    for run in runs:
        n=run['rounds'];costs=[]
        for row in run['all_original_histories']:
            history=row['original_history'];active=tuple(r for r in history if r in ACTIVE)
            cost=n*local['left_and_right_address_basis_pulses']+len(active)*local['conditional_rank2_routing_pulses']
            cost+=sum(r in ('y_no','both_no') for r in history)*local['conditional_second_return_bit_flip_pulses']
            cost+=sum(by_sequence[active[:j]] for j in range(2,len(active)+1))
            row['additional_compression_native_pulses']=cost;costs.append(cost)
        resources.append({'rounds':n,'additional_native_compression_pulses_worst_case':max(costs),
            'additional_native_compression_pulses_minimum':min(costs),
            'additional_coherent_load_export_service_events_worst_case':2*max(costs),
            'retained_classical_history_bits':2*n,'retained_quantum_archive_qubits':run['worst_case_quantum_archive_qubits'],
            'distinct_small_streaming_unitaries_needed':sum(2**j for j in range(2,n+1)),
            'source_preparation_cycle_measurement_and_archive_transport_costs_are_additional':True,
            'program_storage_classical_addressing_and_analog_precision_are_supplied':True,
            'information_dimension_equal_to_physical_work_claimed':False})
    result={'schema':'oph.lepton_archive_compression.v1','source_bindings':pins,'all_checks_pass':True,
        'source_contract':{'source':'actual lepton_cycle mass instrument','theta_exact':str(THETA),
            'independent_fresh_source_pairs':True,'arbitrary_family_and_untouched_reference_inputs':True,
            'original_classical_outcomes':list(LABELS),'original_source_dimension_per_round':3600,
            'source_padding_dimension_per_round':4096,'committed_classical_history_domain':True},
        'exact_rank_certificate':rank,'local_quantum_archive_basis_vectors':{r:array_json(v) for r,v in local_frames()['bases'].items()},
        'complete_bounded_runs':runs,'native_address_compression_service':local,'native_streaming_compression_services':services,
        'resources':resources,'scope':{'exact_ideal_instrument_and_full_archive_isometry':True,
            'original_record_histories_erased_or_coarsened':False,'any_rejected_outcome_removed':False,
            'lossless_on_all_admitted_family_reference_states':True,'arbitrary_unprepared_raw_source_states_compressed':False,
            'physical_work_or_autonomous_reset_derived':False,'binary64_native_pulses_have_certified_precision':False,
            'same_microphysical_controls_supplied_as_native_service_compiler':True}}
    if bindings()!=pins:raise ArithmeticError('archive source changed during calculation')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return type(a) in (int,float) and math.isfinite(a) and abs(a-b)<2e-10+2e-8*abs(b)
    return type(a) is type(b) and a==b


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('archive compression receipt differs')
    print('PASS: exact native archive ranks, complete histories, reference isometries and compression services.')
