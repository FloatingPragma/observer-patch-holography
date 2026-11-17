"""Native analyzer span, full-occupation archive isometry and leakage bound."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy import linalg as la
import sympy as s

from model.lepton_mass_compiler import analyzer
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/lepton_mass_compiler/results.json','model/lepton_mass_compiler/PROOF.md',
                  'model/lepton_mass_compiler/compiler.py','model/charged_lepton_response/CAR_KOIDE.md',
                  'requirements.txt','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def authenticate_parent():
    p=json.loads((ROOT/'model/lepton_mass_compiler/results.json').read_text())
    for name,digest in p['source_bindings'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('changed full-CAR compiler source: '+name)
    return p


@lru_cache(maxsize=1)
def native_spans():
    """Pull the actual Q columns through U; no fitted mass enters the span."""
    parent=authenticate_parent()
    blank=np.zeros((3600,3,3),complex);blank[0]=np.eye(3)
    vectors=analyzer.apply(blank,adjoint=True)
    forward=analyzer.apply(blank)
    rows=[]
    for k in range(3):
        vector=vectors[:,k,k]
        lam=vector[0];effect=float(abs(lam)**2)
        physical=np.zeros(3600,complex);physical[0]=1
        ancilla=(vector-lam*physical)/np.sqrt(1-effect)
        frame=np.column_stack((physical,ancilla))
        reduced=np.array([[effect,lam*np.sqrt(1-effect)],
                          [lam.conjugate()*np.sqrt(1-effect),1-effect]],complex)
        # Matrix-unit attachment of the source projector to the two-column span.
        injected=np.zeros((3600,3,2),complex);injected[:,k,:]=frame
        projected=analyzer.apply(injected);projected[1:]=0
        direct=analyzer.apply(projected,adjoint=True)[:,k,:]
        maximum=0.
        for start in range(0,3600,180):
            stop=min(3600,start+180)
            original=vector[:,None]*vector[start:stop].conj()[None,:]
            reconstructed=frame@reduced@frame[start:stop].conj().T
            maximum=max(maximum,float(np.max(abs(original-reconstructed))))
        errors={'orthonormal_frame':float(la.norm(frame.conj().T@frame-np.eye(2))),
                'ancilla_has_zero_physical_component':float(abs(ancilla[0])),
                'literal_U_projector_span':float(la.norm(direct-frame@reduced)),
                'all3600_projector_columns':maximum,
                'real_positive_blank_amplitude':float(abs(lam.imag)),
                'parent_effect':abs(effect-parent['actual_native_analyzer']['effect_values_k0_k1_k2'][k]),
                'forward_blank_amplitude':float(abs(forward[0,k,k]-lam))}
        if max(errors.values())>3e-12 or not 0<effect<1 or lam.real<=0:
            raise ArithmeticError('literal native active span failed')
        rows.append({'character':k,'effect':effect,'blank_amplitude':float(lam.real),
                     'raw_address_dimension':3600,'retained_address_dimension':2,
                     'all_projector_columns_checked':3600,'errors':errors})
    return rows


def car(modes,index):
    a=np.zeros((1<<modes,1<<modes),complex)
    for word in range(1<<modes):
        if word>>index&1:
            a[word^(1<<index),word]=(-1)**((word&((1<<index)-1)).bit_count())
    return a


@lru_cache(maxsize=1)
def exact_leakage():
    """Polynomial identity on every physical occupation of one directed pair."""
    e,r=s.symbols('e r',real=True)
    a=[s.Matrix(car(4,i)).applyfunc(s.nsimplify) for i in range(4)]
    w=s.Matrix([[e,r],[r,1-e]])
    pair=sum((w[i,j]*a[i].H*a[2+j].H for i in range(2) for j in range(2)),s.zeros(16))
    v=pair+pair.H;physical=[0,1,4,5]
    compressed=v.extract(physical,physical)
    leakage=(v*v).extract(physical,physical)-compressed*compressed
    leakage=leakage.applyfunc(lambda x:s.expand(x).subs(r*r,e*(1-e)).expand())
    expected=s.diag(1-e*e,1-e,1-e,(1-e)**2)
    if any(s.simplify(x)!=0 for x in leakage-expected):
        raise ArithmeticError('exact physical-vacuum leakage operator differs')
    return {'four_mode_Fock_dimension':16,'all_physical_occupations':4,
            'physical_basis':['blank','left','right','left_and_right'],
            'leakage_eigenvalues':['1-e^2','1-e','1-e','(1-e)^2'],
            'identity':'J†V²J−(J†VJ)²=(1−e²)I−e(1−e)(n_left+n_right)',
            'validity':'0<e<1; r²=e(1−e)',
            'maximum_exact':'1-e^2','maximum_attained_by':'physical pair vacuum',
            'source_pairs_per_character_per_site':2,
            'all_site_leakage_norm_squared':'V_sites*(6-2*sum_k e_k^2)',
            'all_site_maximum_attained':True}


def wedge(frame):
    """Full exterior lift of a number-preserving mode isometry."""
    raw,small=frame.shape;out=np.zeros((1<<raw,1<<small),complex)
    for column in range(1<<small):
        chosen=[j for j in range(small) if column>>j&1]
        for occupied in combinations(range(raw),len(chosen)):
            row=sum(1<<j for j in occupied)
            out[row,column]=la.det(frame[np.ix_(occupied,chosen)]) if chosen else 1
    return out


def pairing_hamiltonian(w,cost=F(1,3)):
    half=w.shape[0];modes=2*half;ann=[car(modes,i) for i in range(modes)]
    pair=sum((w[i,j]*ann[i].conj().T@ann[half+j].conj().T
              for i in range(half) for j in range(half)),np.zeros((1<<modes,1<<modes),complex))
    number=sum((a.conj().T@a for a in ann),np.zeros_like(pair))
    return pair+pair.conj().T+float(cost)*number


def kraus(u,physical,ancillary):
    pwords=[sum(1<<position for j,position in enumerate(physical) if word>>j&1)
            for word in range(1<<len(physical))]
    awords=[sum(1<<position for j,position in enumerate(ancillary) if word>>j&1)
            for word in range(1<<len(ancillary))]
    operators=[]
    for mask in awords:
        signs=np.array([(-1)**sum((mask&((1<<position)-1)).bit_count()
                                  for position in physical if word>>position&1) for word in pwords])
        operators.append(signs[:,None]*u[np.ix_([word|mask for word in pwords],pwords)])
    return np.asarray(operators)


def histories(operators,steps):
    result=np.eye(operators.shape[1],dtype=complex)[None,:,:]
    for _ in range(steps):
        result=np.einsum('aij,bjk->baik',operators,result,optimize=True).reshape(-1,*operators.shape[1:])
    return result


@lru_cache(maxsize=8)
def complete_dilation_control(effect,steps=3,time=F(1,8)):
    """Independent six-mode raw embedding versus the four-mode active source."""
    if type(steps) is not int or not 1<=steps<=3:raise ValueError('control cap is three rounds')
    e=float(effect)
    if not 0<e<1:raise ValueError('strict effect required')
    w=np.array([[e,np.sqrt(e*(1-e))],[np.sqrt(e*(1-e)),1-e]])
    # Complex ancillary rotation prevents an accidental real-frame comparison.
    b=np.array([[1,0],[0,3/5],[0,4j/5]],complex)
    raw_w=b@w@b.conj().T
    one=la.block_diag(b,b.conj());embedding=wedge(one)
    hc=pairing_hamiltonian(w);hr=pairing_hamiltonian(raw_w)
    h=float(time)/steps;uc=la.expm(-1j*h*hc);ur=la.expm(-1j*h*hr)
    small=kraus(uc,(0,2),(1,3));raw=kraus(ur,(0,3),(1,2,4,5))
    env_one=np.array([[3/5,0],[4j/5,0],[0,3/5],[0,-4j/5]],complex)
    env=wedge(env_one)
    small_hist=histories(small,steps);raw_hist=histories(raw,steps)
    history_embedding=env
    for _ in range(steps-1):history_embedding=np.kron(history_embedding,env)
    predicted=np.einsum('ab,bij->aij',history_embedding,small_hist,optimize=True)
    errors={'full16_input_occupation_intertwiner':float(la.norm(embedding.conj().T@embedding-np.eye(16))),
            'Hamiltonian_all16_columns':float(la.norm(hr@embedding-embedding@hc)),
            'unitary_all16_columns':float(la.norm(ur@embedding-embedding@uc)),
            'quantum_environment_one_round_all_input_operators':float(la.norm(raw-np.einsum('ab,bij->aij',env,small))),
            'quantum_archive_all_histories_all_input_operators':float(la.norm(raw_hist-predicted)),
            'complete_small_Kraus_identity':float(la.norm(np.einsum('aij,aik->jk',small.conj(),small)-np.eye(4))),
            'complete_raw_Kraus_identity':float(la.norm(np.einsum('aij,aik->jk',raw.conj(),raw)-np.eye(4)))}
    if max(errors.values())>3e-12:raise ArithmeticError('full-CAR archive isometry failed')
    wrong=pairing_hamiltonian(b@w@b.T)
    wrong_frame_error=float(la.norm(wrong@embedding-embedding@hc))
    if wrong_frame_error<.1:raise ArithmeticError('conjugate-leg negative control lost sensitivity')
    coarse={''.join(map(str,bits)):0. for bits in product((0,1),repeat=steps)}
    fine=[];max_coherence=0.
    for index,labels in enumerate(product(range(4),repeat=steps)):
        key=''.join(str(int(j!=0)) for j in labels)
        probability=float(la.norm(small_hist[index],'fro')**2/4)
        coarse[key]+=probability;fine.append({'ancillary_history':list(labels),'probability':probability})
    for bits in product((0,1),repeat=steps):
        selected=[row['probability'] for row in fine if tuple(int(j!=0) for j in row['ancillary_history'])==bits]
        max_coherence=max(max_coherence,float(np.sqrt(max(0.,sum(selected)**2-sum(x*x for x in selected)))))
    physical=[0,1,4,5];hp=hc[np.ix_(physical,physical)];ideal=la.expm(-1j*float(time)*hp)
    aa=small_hist[0].reshape(-1)/2;bb=ideal.reshape(-1)/2;weight=float(np.vdot(aa,aa).real)
    complete_error=float(np.sqrt(max(0.,(1+weight)**2-4*abs(np.vdot(bb,aa))**2))+1-weight)
    sharper=2*(1-e*e)*float(time)**2/steps
    if complete_error>sharper+2e-12 or abs(sum(coarse.values())-1)>2e-12:
        raise ArithmeticError('sharper complete-history bound failed its independent control')
    return {'effect':e,'raw_Fock_dimension':64,'compressed_Fock_dimension':16,
            'physical_input_dimension':4,'rounds':steps,'clock_exact':str(time),
            'raw_ancillary_outcomes_per_round':16,'compressed_ancillary_outcomes_per_round':4,
            'raw_quantum_archive_coordinates':16**steps,'compressed_quantum_archive_coordinates':4**steps,
            'committed_classical_histories':2**steps,'complete_fine_diagonal_records':fine,
            'complete_coarse_history_probabilities':coarse,
            'maximum_retained_offdiagonal_archive_norm':max_coherence,
            'entangled_input':'maximally entangled physical/reference input; equality also checked on every input operator',
            'complete_output_trace_error':complete_error,'sharper_unhalved_bound':sharper,
            'wrong_equal_leg_frame_Hamiltonian_error':wrong_frame_error,
            'all_input_quantum_archive_isometry':True,'fine_environment_coordinates_measured':False,
            'spectator_archive_coordinates_are_fixed_vacuum_under_inverse_isometry':True,'errors':errors}


def resource_bound(steps,time=F(1),sites=1):
    if type(steps) is not int or steps<1 or type(sites) is not int or sites<1:
        raise ValueError('positive integer rounds and sites required')
    if type(time) not in (int,F) or time<0:raise ValueError('nonnegative exact time required')
    intervals=[tuple(map(F,x)) for x in analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    lower=F(sites)*(6-2*sum(hi*hi for lo,hi in intervals))
    upper=F(sites)*(6-2*sum(lo*lo for lo,hi in intervals))
    epsilon=min(F(2),2*upper*F(time)**2/steps)
    original=F(108*sites*sites)*F(time)**2/steps
    return {'sites':sites,'rounds':steps,'time_exact':str(time),
            'leakage_norm_squared_interval':[str(lower),str(upper)],
            'complete_unhalved_trace_error_upper_exact':str(epsilon),
            'complete_unhalved_trace_error_upper_numeric':float(epsilon),
            'old_bound_divided_by_new_bound_before_clipping':float(original/(2*upper*F(time)**2/steps)) if time else None,
            'physical_CAR_modes':90*sites,'active_ancillary_CAR_modes':12*sites,
            'total_CAR_modes':102*sites,'original_total_CAR_modes':324000*sites,
            'original_ancillary_CAR_modes':323910*sites,
            'ancillary_mode_compression_factor_exact':'53985/2',
            'quantum_archive_modes':12*sites*steps,
            'classical_flag_bits':steps,'original_fine_addresses_reconstructible_by_isometry':True,
            'precision_round_count_executed':False}


def calculate():
    pins=bindings();spans=native_spans();leakage=exact_leakage()
    base=resource_bound(1)
    coefficient=2*F(base['leakage_norm_squared_interval'][1])
    examples=[]
    for sites in (1,64):
        rounds=math.ceil(coefficient*sites*10**12)
        row=resource_bound(rounds,sites=sites)
        row['target_unhalved_error_exact']='1/1000000000000'
        row['old_certified_rounds']=108*sites**2*10**12
        row['round_reduction_factor']=row['old_certified_rounds']/rounds
        row['total_exported_ancillary_reduction_factor']=323910*sites*row['old_certified_rounds']/row['quantum_archive_modes']
        examples.append(row)
    result={'schema':'oph.lepton_compiler_compression.v1','source_bindings':pins,'all_checks_pass':True,
            'native_active_spans':spans,'exact_leakage_operator':leakage,
            'independent_full_CAR_controls':[complete_dilation_control(row['effect']) for row in spans],
            'precision_resource_examples':examples,
            'exact_reduction':{'physical_modes_per_site':90,'ancillary_modes_per_site':12,
                'total_modes_per_site':102,'physical_spectator_modes_per_site':78,
                'active_pair_blocks_per_site':6,
                'full_original_source_intertwiner':'Gamma(T) with identity on physical modes and native active vectors on ancillary modes; unused directions stay in vacuum',
                'full_kinetic_intertwiner':'the same spin/site Weyl hopping is identity on address and family, hence preserves the selected spans',
                'record_intertwiner':'ancillary vacuum is mapped to vacuum and its orthogonal complement to the nonvacuum sector; each committed flag is unchanged',
                'arbitrary_physical_occupation_and_untouched_reference':True,
                'minimality_domain':'number-preserving CAR mode embeddings fixing the90 physical modes and reproducing the first-order source coupling',
                'minimal_ancillary_modes_in_that_domain':12,
                'universal_minimal_environment_Hilbert_dimension_claimed':False,
                'neutral_background_fixed':True,'all_Higgs_orientations_simultaneously_compressed':False},
            'archive_bound':{'offdiagonal_block_norm_squared':'B²=V_sites*(6−2 sum_k e_k²)',
                'single_round_vacuum_operator_error':'||P exp(−ihH)P−exp(−ihPHP)||≤B²h²/2',
                'complete_committed_history_unhalved_error':'min(2,2B²t²/N)',
                'uniform_simplification':'min(2,12 V_sites t²/N)',
                'any_nonvacuum_probability_upper':'min(1,B²t²/N)',
                'additional_complete_dilation_error':'2N eta when eta is per-round half-trace/half-diamond error',
                'kinetic_norm_enters_bound':False,'coherent_unmeasured_flag_error_has_same_rate':False},
            'scope':{'same_declared_neutral_compiler_complete_instrument':True,
                'source_phase_or_Higgs_background_changed':False,'any_outcomes_or_archives_discarded':False,
                'positive_compression_and_resource_theorem':True,'physical_Standard_Model_attachment_closed':False,
                'passive_native_clock_or_hardware_implementation_derived':False,
                'compression_is_a_declared_exact_coordinate_implementation':True}}
    if bindings()!=pins:raise ArithmeticError('compression source changed during execution')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return type(a) in (int,float) and math.isfinite(a) and abs(a-b)<=2e-12+2e-9*abs(b)
    return type(a) is type(b) and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('compression receipt differs')
    print('Native active spans, complete quantum archives and linear-volume compiler bounds verified.')
