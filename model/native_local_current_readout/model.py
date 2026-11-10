"""Complete native current records from a local pulse with source drift retained."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json
import sympy as s
from . import apparatus
from model.native_higgs_current_response import model as response
from model.coherent_defect_wave.chiral_joint.atomic_reference_control import lattice
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENT='model/native_higgs_current_response/results.json'
PARENT_SHA='55fea7f49b3683fb34f3a72aa686cab58646608a3267db02f7f6d3a943b9352a'
TAU=F(1,2**20)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def require(test,message):
    if not test:raise ValueError(message)

def encode(pair):return list(map(str,pair))

def widen_probability(pair,error):
    lo,hi=map(F,pair);return(max(F(),lo-error),min(F(1),hi+error))


@lru_cache(None)
def parent():
    require(digest(ROOT/PARENT)==PARENT_SHA,'frozen same-source current parent changed')
    value=json.loads((ROOT/PARENT).read_text())
    for path,sha in value['input_sha256'].items():require(digest(ROOT/path)==sha,'source parent input differs: '+path)
    require(response.certificate()==value,'same-source current parent replay differs')
    return value


@lru_cache(None)
def local_source_bound():
    config=Configuration(side=2)
    assert (config.y_up,config.y_down,config.y_lepton,config.scalar_penalty,config.link_penalty)==(1,1,1,1,1)
    assert config.spacing==s.pi
    rows=kinetic_terms(config)
    touched=[t for t in rows if 0 in(t.left_cell,t.right_cell)]
    assert len(rows)==24 and len(touched)==6
    weights=sum(abs(t.coefficient)for t in touched)
    assert weights==s.Rational(3,2)
    assert lattice.local_norm()==57
    # The actual scalar penalty commutes with W; all up/down/lepton Yukawa
    # terms, including their shared scalar, enter the remaining norm56.
    W,f,_=apparatus.current_algebra()
    scalar_penalty=s.kronecker_product(s.diag(0,1,1,1,1),s.eye(16))
    assert scalar_penalty*f==f*scalar_penalty
    from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
    native_W,*_=apparatus.dressed_weak_polynomial()
    parts,ys=car.quantum_yukawa();values={y:1 for y in ys.values()}
    witnesses={}
    for name,part in parts.items():
        commutator={key:s.expand(v.subs(values))for key,v in car.commutator(part,native_W).items()}
        commutator={key:v for key,v in commutator.items()if v}
        assert commutator
        witnesses[name]=len(commutator)
    normY=F(56);normT=F(60)*F(weights);lam=F(1,2**34)
    comm=2*F(2)*(normY+lam*normT)
    assert comm==224+F(360,2**34)
    return {'measured_source_copy':0,'measured_native_cell':0,
            'native_source_hopping_coefficient':str(lam),
            'all_periodic_directed_terms':len(rows),'touched_directed_terms':len(touched),
            'actual_touched_paths':[{'left_cell':t.left_cell,'right_cell':t.right_cell,'path':list(t.path),
                                      'coefficient':str(t.coefficient)}for t in touched],
            'touched_coefficient_absolute_sum':str(weights),
            'complete_local_Yukawa_norm_upper':str(normY),
            'all_up_down_and_lepton_Yukawa_terms_included':True,
            'nonzero_actual_unit_Yukawa_W_commutator_terms':witnesses,
            'scalar_penalty_and_link_energy_commutators_zero':True,
            'all_other_source_copies_and_reference_commutators_zero':True,
            'touched_full_Fock_hopping_norm_upper':str(normT),
            'full_source_commutator_H_fW_norm_upper':str(comm),
            'source_operator_or_reference_clock_replaced':False}


@lru_cache(None)
def calculate():
    require(type(TAU)is F and TAU>0,'positive exact rational source duration required')
    inherited=parent();local=apparatus.local_execution();bound=local_source_bound()
    Pnorm=F(local['pointer_generator_norm_upper']);Gnorm=2*Pnorm
    quantization=F(local['ideal_readout_error_from_pointer_quantization_upper'])
    comm=F(bound['full_source_commutator_H_fW_norm_upper'])*Pnorm
    drift=TAU*comm/2;distance=drift+quantization
    assert 0<distance<F(27,100000)
    means={}
    old=inherited['prepared_current']
    for key in ('elementary_W_mean_interval','Higgs_composite_W_mean_interval'):
        lo,hi=map(F,old[key]);means[key]=encode((lo-4*distance,hi+4*distance))
        assert F(means[key][1])<0
    rows=[]
    for oldrow in old['five_complete_current_records']:
        row={'W_value':oldrow['W_value']}
        for short,prefix in (('elementary','actual_elementary'),('composite','actual_composite')):
            interval=widen_probability(oldrow[prefix+'_bin_conditional_probability_interval'],distance)
            weight=F(old[short+'_original_attempt_bin_probability_lower'])
            row[prefix+'_bin_conditional_probability_interval']=encode(interval)
            row['original_'+short+'_attempt_joint_probability_lower']=str(weight*interval[0])
            row['fair_setting_original_'+short+'_joint_probability_lower']=str(weight*interval[0]/2)
        rows.append(row)
    minus=next(row for row in rows if row['W_value']=='-1/2')
    assert all(F(minus['original_'+pump+'_attempt_joint_probability_lower'])>0 for pump in('elementary','composite'))
    Hnorm=F(inherited['apparatus']['full_raw_Hext_norm_upper'])
    pi_upper=F(apparatus.iv.pi().hi,apparatus.iv.SCALE)
    drive=Gnorm/TAU
    return {'schema':'oph.native_local_current_readout.v1',
        'source':{'Hamiltonian':'the unchanged complete native Hext of the actual8-cell preparation and all retained copies/references',
            'source_hamiltonian_continues_through_final_pulse':True,'full_Hext_cancellation_used_for_final_readout':False,
            'full_Hext_matrix_exponentiated_in_this_local_evaluator':False,
            'first_three_preparation_pump_bin_services_unchanged':True,
            'exact_interacting_ground_projection_added':False,
            'parent_current_response':PARENT,'parent_current_response_sha256':PARENT_SHA},
        'local_apparatus':local,'source_locality_bound':bound,
        'complete_clock_and_error':{
            'earlier_service_source_durations':['1','1','1'],'final_pulse_source_duration':str(TAU),
            'total_elapsed_source_time':str(3+TAU),'old_time4_result_or_receipt_modified':False,
            'centered_current_polynomial':'f(W)=(4W^4-4W^3-13W^2+25W)/15',
            'local_drive':'V=success_flag tensor f(W) tensor P40/tau',
            'actual_evolution':'exp[-i tau (Hext+V)]',
            'ideal_comparison':'exp(-i tau Hext) times the complete controlled current-pointer shift',
            'full_source_drive_commutator_norm_upper':str(comm),
            'uninterrupted_source_product_error_upper':str(drift),
            'pointer_quantization_error_upper':str(quantization),
            'complete_output_operator_and_half_diamond_distance_upper':str(distance),
            'all_retained_history_total_variation_error_upper':str(distance),
            'conditional_current_mean_absolute_error_upper':str(4*distance),
            'arbitrary_entangled_partners_and_references_in_bound':True,
            'post_readout_success_selection_or_normalization_used':False},
        'complete_prepared_records':{
            'current_means':means,
            'intervals_describe_pointer_decoded_current_means':True,
            'instantaneous_post_measurement_W_expectation_asserted':False,'five_complete_current_records':rows,
            'all_fine_preparation_flags_and_source_outputs_retained':True,
            'complete_setting_pump_energy_current_alphabet_size':80,
            'same_apparatus_A_B_conditional_history_distance_upper':old['same_history_TV_upper'],
            'A_B_history_comparison_conditions_only_on_retained_bin_success':True,
            'unconditional_A_and_B_rates_equated':False,
            'disabled_drive_branches_freely_evolve_for_tau':True,
            'disabled_pointer_zero_is_not_a_W_measurement':True,
            'every_current_label_retained_on_admitted_branch':True,
            'ideal_branch_law':'Q/2 tr(P_w U(2) P_I C rho C^dagger P_I U(2)^dagger P_w); actual complete law differs by the stated bound',
            'current_refers_to_the_interval_start_at_source_time3':True},
        'control_resources':{
            'explicit_additional_primitive':'bounded-support gauge-invariant polynomial-current coupling to a neutral five-state pointer, controlled by retained apparatus flags',
            'coupling_is_declared_not_compiled_from_native_C_E':True,
            'represented_pointer_coefficients_are_exact_rationals':True,
            'analog_realization_accuracy_or_efficient_synthesis_proved':False,
            'full_Gauss_constraints_preserved_throughout_drive':True,
            'source_raw_Hext_norm_upper':str(Hnorm),'local_final_drive_norm_upper':str(drive),
            'work_is_unconditioned_mean_bare_source_energy_accounting':True,
            'a_separate_work_measurement_instrument_added':False,
            'final_pulse_net_bare_source_work_absolute_upper':str(min(2*Hnorm,comm)),
            'final_pulse_gross_switching_work_upper':str(2*drive),
            'complete_four_interval_net_bare_source_work_absolute_upper':str(2*Hnorm),
            'complete_four_interval_gross_switching_work_upper':str(6*(Hnorm+pi_upper)+2*drive),
            'retained_register_capacity_nats':inherited['apparatus']['complete_register_capacity_nats_formula'],
            'additional_pointer_beyond_parent_current_package':False,
            'controller_coefficient_and_synthesis_workspace_in_capacity':False,
            'full_matrix_preparation_and_exact_spectral_bin_premises_removed':False},
        'scope':{'local_same_source_original_W_readout_with_nonzero_response':True,
            'whole_lattice_wavefunction_numerically_evolved':False,
            'autonomous_native_repair_compiler':False,'weak_boson_or_chiral_Lorentz_interaction':False,
            'Standard_Model_species_or_GeV_mass_identified':False}}


def certificate():
    result=calculate();saved=parent()
    paths=set(saved['input_sha256'])|{PARENT,'model/particle_common_contract/provenance.py'}
    paths.update(str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md'))
    result=dict(result)
    result['input_sha256']={p:digest(ROOT/p)for p in sorted(import_closure(paths))}
    result['owned_source_sha256']={p.name:digest(p)for p in sorted(HERE.iterdir())if p.suffix in('.py','.md')}
    return result
