"""Complete current records for the actual elementary/Higgs shared excitation."""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import hashlib,json
from . import local
from model.native_charge_resolution import bounds as resolution
from model.native_higgs_filtered_state import model as filtered
from model.native_weak_current_response.bounds import divide_positive, expand, probability, encode
from model.native_weak_current_response.certificate import certificate as replay_weak
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENTS={
    'model/native_higgs_filtered_state/results.json':'19527188c4ffa855df73330c5fc1d9efead08308e0c7719e6f67fd845e901029',
    'model/native_weak_current_response/results.json':'8415fe269dfa3c18b2136c3505e1fd2c701456120b1dd395dc6bd1499a7dd3d1'}


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def parents():
    out={}
    for name,expected in PARENTS.items():
        if digest(ROOT/name)!=expected: raise ValueError('frozen source differs: '+name)
        value=json.loads((ROOT/name).read_text())
        filtered.higgs.validate(value['input_sha256'])
        for owned,sha in value.get('owned_source_sha256',{}).items():
            if digest((ROOT/name).parent/owned)!=sha: raise ValueError('owned source differs: '+owned)
        out[name]=value
    assert filtered.certificate()==out['model/native_higgs_filtered_state/results.json']
    assert replay_weak()==out['model/native_weak_current_response/results.json']
    return out


@lru_cache(None)
def calculate():
    parents(); atom=local.calculate(); data=resolution.data(); window=resolution.window(128)
    prep=filtered.preparation(3); w0=data['weights'][0]
    w=tuple(map(F,window['unnormalized_bin_weight_interval']))
    delay=F(2)
    error=filtered.rounded(
        filtered.root(w0[1],True)*data['v']/F(window['selected_atom_distance_to_bin_complement_lower'])
        +F(window['complete_other_atom_leakage_norm_upper'])+data['eta']
        +delay*data['v']*filtered.root(w0[1],True),True)
    epsilon=F(prep['retained_purification_error_upper']); Q=F(prep['complete_preparation_success_lower'])
    Hnorm=8*(F(600)+F(360,2**34))
    normsum=filtered.root(w[1],True)+filtered.root(w0[1],True)
    atomic_insert=tuple(map(F,atom['atomic_unnormalized_W_insertion_interval']))
    exact_insert=expand(atomic_insert,2*error*normsum)
    elementary_mean=expand(divide_positive(exact_insert,w),3*epsilon/w[0])
    joint=filtered.filtered_bounds(window['source_excitation_interval_relative_actual_native_ground'],w,epsilon)
    distance=F(joint['actual_same_history_total_variation_upper'])
    composite_mean=expand(elementary_mean,4*distance)
    assert elementary_mean[1]<0 and composite_mean[1]<0
    pA=F(joint['actual_elementary_in_bin_probability_lower'])
    pB=F(joint['actual_composite_B_over_2_in_bin_probability_lower'])
    histogram=[]
    for eig,row in atom['atomic_addition_branch_W'].items():
        if not eig.startswith('probability_'): continue
        p=tuple(map(F,row['outward_interval'])); atomic_joint=(w0[0]*p[0],w0[1]*p[1])
        source_joint=probability(expand(atomic_joint,error*normsum))
        actual_joint=probability(expand(source_joint,epsilon))
        conditional=probability(expand(divide_positive(source_joint,w),epsilon/w[0]))
        composite=probability(expand(conditional,distance))
        histogram.append({'W_value':eig.removeprefix('probability_'),
            'atomic_addition_bin_conditional_probability_interval':encode(p),
            'actual_elementary_joint_bin_current_probability_interval':encode(actual_joint),
            'actual_elementary_bin_conditional_probability_interval':encode(conditional),
            'actual_composite_bin_conditional_probability_interval':encode(composite),
            'original_elementary_attempt_joint_probability_lower':str(Q*actual_joint[0]),
            'original_composite_attempt_joint_probability_lower':str(Q*pB*composite[0]),
            'fair_setting_original_elementary_joint_probability_lower':str(Q*actual_joint[0]/2),
            'fair_setting_original_composite_joint_probability_lower':str(Q*pB*composite[0]/2)})
    negative=next(x for x in histogram if x['W_value']=='-1/2')
    assert F(negative['original_composite_attempt_joint_probability_lower'])>0
    return {'schema':'oph.native_higgs_current_response.v1',
        'source':{'Hamiltonian':'unchanged full eight-cell h=H0+2^-34 T and complete old/dual/port references',
            'reference_energy_on_selected_even_ports':'0',
            'full_source_reference_complements_retained':True,
            'elementary_probe':'A_alpha=c_(alpha,ec)^dagger',
            'Higgs_probe':'B_alpha=-epsilon_(alpha,beta) sum_i phi_i^dagger c_(beta,10+i)',
            'current':'original quantum dressed weak density W=sum_a tr(sigma_a phi phi^dagger) j_a',
            'original_W_unchanged':True,'current_is_not_bare_Gauss_generator':True},
        'conjugation':{'full_native_H_preserved':True,'source_time_reversed':True,
            'formula':'Theta W Theta^-1=-sum_a(sum_ij conjugate(sigma_a[j,i]) phi_i^dagger phi_j) j_a',
            'W_invariant':False,'same_untransformed_current_on_opposite_charge_branches_equal':False,
            'full_equal_cost_opposite_charge_reference_required':True,
            'old_reference_4N_particle_hole_symmetry_used':False,
            'charge_conjugation_CP_violation_or_particle_identity_inferred':False},
        'atomic_calculation':atom,'preparation':prep,
        'energy_admission':{'radius_in_units_of_360lambda':128,
            'absolute_native_source_bin':[str(F(x)-4) for x in window['absolute_source_plus_reference_bin']],
            'excitation_energy_enclosure':window['source_excitation_interval_relative_actual_native_ground'],
            'full_source_bin_weight_interval':encode(w),
            'full_atomic_packet_amplitude_error_after_two_source_units':str(error),
            'common_elapsed_time_before_current_service':'3 source units including preparation',
            'exact_finite_spectral_service_declared':True,'bin_QPE_implemented':False},
        'prepared_current':{'elementary_W_mean_interval':encode(elementary_mean),
            'Higgs_composite_W_mean_interval':encode(composite_mean),
            'same_history_TV_upper':str(distance),'five_complete_current_records':histogram,
            'elementary_original_attempt_bin_probability_lower':str(Q*pA),
            'composite_original_attempt_bin_probability_lower':str(Q*pB),
            'each_setting_probability':'1/2','actual_preparation_purification_error_upper':str(epsilon)},
        'apparatus':{'complete_setting_pump_energy_current_alphabet_size':80,
            'current_pointer_dimension':5,'full_fine_preparation_records_retained':True,
            'additional_zero_energy_current_pointer_capacity_nats':'ln5',
            'complete_register_capacity_nats_formula':'8*(240 ln2+8 ln5+24 ln119+ln(2^31+1)+16 ln2)+7lnM+ln48+4ln2+ln5',
            'controller_coefficient_and_synthesis_workspace_included':False,
            'failures_and_outside_energy_bin_retained':True,
            'complete_law':['1-Q','Q/2 times each pump failure and outside-bin record for C=A,B/2',
                'Q/2 tr(P_w U(2) P_I C rho C^dagger P_I U(2)^dagger P_w), every current label'],
            'current_service_duration':'1 source unit; current refers to its interval start',
            'common_clock_service':'V(s)=exp(-i Hext s) L exp(i Hext s), exp(-iL)=U_read',
            'constant_drive_alternative':'L=i Log(exp(-i Hext) U_read), V=L-Hext',
            'complete_elapsed_source_time':'4 source units',
            'full_raw_Hext_norm_upper':str(Hnorm),
            'constant_control_norm_upper':str(Hnorm+F(22,7)),
            'net_bare_source_work_absolute_upper':str(2*Hnorm),
            'gross_four_constant_pulse_switching_work_upper':str(8*(Hnorm+F(22,7))),
            'full_readout':'exp(-iHext) U_read',
            'current_measurement_can_change_admitted_energy':True,
            'all_pulse_and_failure_branches_retain_Hext':True,
            'programmable_matrix_service_is_autonomous_local_native_compiler':False},
        'scope':{'same_H_same_current_response_of_actual_elementary_and_Higgs_preparations':True,
            'nonzero_negative_current_outcome_for_both_preparations':True,
            'equal_unconditional_A_and_B_probabilities':False,
            'quantum_scalar_replaced_by_classical_background':False,
            'weak_boson_chiral_Lorentz_interaction_or_SM_identification':False,
            'GeV_rest_mass_or_measured_lepton_ratio':False}}


def certificate():
    value=calculate(); seeds=set(PARENTS)
    for parent in parents().values(): seeds.update(parent['input_sha256'])
    seeds.update(str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md'))
    seeds.add('model/particle_common_contract/provenance.py')
    value=dict(value)
    value['input_sha256']={name:digest(ROOT/name) for name in sorted(import_closure(seeds))}
    value['owned_source_sha256']={path.name:digest(path) for path in sorted(HERE.iterdir()) if path.suffix in ('.py','.md')}
    return value
