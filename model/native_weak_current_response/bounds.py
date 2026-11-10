"""Full native weak-current records from exact local contractions."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from . import local,charged
from model.native_charge_resolution.bounds import root,round_up
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as charged_data
from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.model import certificate as preparation

def interval(record):return tuple(map(F,record['outward_interval']))
def down(value):return -round_up(-F(value))
def box(pair):return down(pair[0]),round_up(F(pair[1]))
def encode(pair):return list(map(str,pair))
def expand(pair,error):return pair[0]-error,pair[1]+error
def probability(pair):return max(F(0),pair[0]),min(F(1),pair[1])
def square(pair):
    lo,hi=pair
    return (F(0) if lo<=0<=hi else min(lo*lo,hi*hi),max(lo*lo,hi*hi))
def divide_positive(pair,denominator):
    assert 0<denominator[0]<=denominator[1]
    values=[x/y for x in pair for y in denominator]
    return min(values),max(values)

@lru_cache(None)
def calculate():
    atomic=local.calculate();target=charged.calculate();source=charged_data('ec');prep=preparation()
    eta=round_up(F(source['eta']));v=F(source['v']);filter=prep['finite_energy_agreement_filter']
    delta=round_up(F(filter['conditional_one_source_trace_distance_upper']))
    pair_fidelity=down(filter['conditional_ground_pair_fidelity_lower'])
    assert delta*delta>=1-pair_fidelity
    q=down(filter['success_probability_lower']);length=int(filter['explicit_integer_M'])
    leak=round_up(F(filter['mixed_ground_excited_amplitude_upper']));assert leak<=F(1,1024)
    original=atomic['actual_local_ground_moments'];m=interval(original['W']);second=interval(original['W2'])
    ground_m=expand(m,3*eta);ground_second=expand(second,4*eta)
    ground_var=(ground_second[0]-square(ground_m)[1],ground_second[1]-square(ground_m)[0])
    assert 0<ground_m[0]<ground_m[1] and ground_var[0]>0
    prepared_m=expand(ground_m,3*delta);prepared_second=expand(ground_second,4*delta)
    histogram=[]
    for eigenvalue in local.EIGENVALUES:
        atom=interval(original['probability_'+str(eigenvalue)])
        exact=probability(expand(atom,eta));actual=probability(expand(exact,delta))
        assert actual[0]>0
        histogram.append({'W_value':str(eigenvalue),'atomic_probability_interval':encode(atom),
          'actual_ground_probability_interval':encode(exact),'prepared_probability_interval':encode(actual),
          'original_attempt_probability_lower':str(q*actual[0])})
    elastic=tuple(x/4 for x in square(ground_m));inelastic=tuple(x/4 for x in ground_var)
    pump=tuple(x/4 for x in ground_second)
    phase_nonzero=probability((inelastic[0]*(1-leak*leak)-delta,inelastic[1]+delta))
    phase_zero=probability((elastic[0]-delta,elastic[1]+leak*leak*inelastic[1]+delta))
    fail=probability((1-pump[1]-delta,1-pump[0]+delta))
    assert min(phase_nonzero[0],phase_zero[0],fail[0])>0
    # The charged instrument records W after one unit for the pump and one
    # for bin admission. Both source Hamiltonians continue throughout.
    delay=F(2);weight=box(source['local']['rows'][0]['weight_interval'])
    norm_error=round_up(F(source['window_error']))+v*delay*root(weight[1],True)
    bin_weight=box(source['probabilities']['line'])
    normsum=root(bin_weight[1],True)+root(weight[1],True)
    target_m=interval(target['charged_target_ground_current']['W'])
    atom_insert=tuple(map(F,target['atomic_unnormalized_charged_current_insertion_interval']))
    exact_insert=expand(atom_insert,2*norm_error*normsum)
    exact_conditional=divide_positive(exact_insert,bin_weight)
    prepared_conditional=expand(exact_conditional,3*delta/bin_weight[0])
    actual_bin=(bin_weight[0]-delta,bin_weight[1]+delta)
    prepared_insert=expand(exact_insert,4*delta)
    assert prepared_conditional[0]>0 and prepared_insert[0]>0
    charged_histogram=[]
    for eigenvalue in local.EIGENVALUES:
        atom=interval(target['charged_target_ground_current']['probability_'+str(eigenvalue)])
        atom_joint=(weight[0]*atom[0],weight[1]*atom[1])
        joint=probability(expand(atom_joint,norm_error*normsum))
        conditional=probability(divide_positive(joint,bin_weight))
        actual=probability(expand(conditional,delta/bin_weight[0]))
        joint_actual=probability(expand(joint,delta))
        charged_histogram.append({'W_value':str(eigenvalue),'atomic_conditional_probability_interval':encode(atom),
          'exact_ground_input_joint_bin_current_probability_interval':encode(joint),
          'actual_prepared_joint_bin_current_probability_interval':encode(joint_actual),
          'actual_prepared_bin_conditional_probability_interval':encode(actual),
          'original_attempt_joint_probability_lower':str(q*joint_actual[0])})
    return {'schema':'oph.actual_native_weak_current_response.v1',
      'atomic_local_calculation':atomic,'atomic_charged_calculation':target,
      'source':{'Hamiltonian':'h_active+h_partner+4N_reference, both h are the full eight-cell H0+2^-34 T',
        'source_dimension':'D^2, D=2^240 5^8 119^24','charge_reference_dimension':'2^30',
        'position_reference_dimension':8,'position_reference_H':'0',
        'source_only_ground_unique':True,'joint_ground_unique_with_position_reference':False,
        'charged_setting_position_reference':'origin basis state; static charge reference is at the recorded position',
        'neutral_setting_position_reference':'uniform superposition of eight positions',
        'full_source_complements_retained':True,'same_source_clock_on_all_branches':True},
      'actual_preparation':{'pair_success_lower':str(q),'pair_ground_fidelity_lower':str(pair_fidelity),
        'pair_trace_distance_upper':str(delta),'single_source_trace_distance_upper':str(delta),
        'exact_P0_added':False,'actual_prepared_state_asserted_stationary':False,
        'all_source_reference_and_preparation_filter_records_retained':True},
      'ground_current':{'ground_trial_vector_error_upper':str(eta),
        'actual_ground_W_mean_interval':encode(ground_m),'actual_ground_W2_mean_interval':encode(ground_second),
        'actual_ground_variance_interval':encode(ground_var),
        'actual_prepared_W_mean_interval':encode(prepared_m),'actual_prepared_W2_mean_interval':encode(prepared_second),
        'five_complete_current_records':histogram,
        'source_energy_spectral_inelastic_weight_for_W_over_two_interval':encode(inelastic),
        'source_energy_spectral_elastic_weight_for_W_over_two_interval':encode(elastic)},
      'finite_inelastic_experiment':{'current_pump':'W_tilde/2 on uniform position reference',
        'Cayley':'C=(I-i h/B)(I+i h/B)^(-1), Z=C_active tensor C_partner^dagger',
        'B':filter['full_absolute_H_norm_upper'],'source_gap_lower':filter['spectral_gap_lower'],
        'additional_pointer_dimension':length,'all_additional_pointer_labels_retained':True,
        'nonzero_phase_amplitude_on_equal_ground_pair':'0',
        'mixed_ground_excited_zero_phase_amplitude_upper':str(leak),
        'conditional_pump_failure_interval':encode(fail),
        'conditional_zero_phase_interval':encode(phase_zero),
        'conditional_nonzero_phase_interval':encode(phase_nonzero),
        'original_attempt_nonzero_phase_probability_lower':str(q*phase_nonzero[0]),
        'complete_law':['1-q','q tr(F rho_pair F^dagger)',
          'q tr(K_k (W_tilde/2) rho_pair (W_tilde/2) K_k^dagger), every k=0,...,M-1'],
        'F':'sqrt(I-W_tilde^2/4)','K_k':'M^-1 sum_j exp(2 pi i j k/M) Z^j',
        'unknown_ground_energy_as_free_clock_phase':False,
        'exact_finite_algebraic_Cayley_service_declared':True,
        'large_source_pointer_matrix_executed':False},
      'charged_current':{'family':'ec','reference_cost':'4',
        'absolute_total_source_plus_reference_bin':list(map(str,source['window'])),
        'source_excitation_bin_after_declared_reference_subtraction':list(map(str,source['source_excitation'])),
        'bin_selected_on':'A=h_active+4N_reference, exact full-source spectral service',
        'source_delay_before_current_record':str(delay),
        'norm_error_after_delay_upper':str(norm_error),
        'unnormalized_exact_ground_bin_weight_interval':encode(bin_weight),
        'actual_prepared_bin_weight_interval':encode(actual_bin),
        'original_attempt_bin_probability_lower':str(q*actual_bin[0]),
        'exact_ground_input_unnormalized_P_I_W_P_I_insertion_interval':encode(exact_insert),
        'actual_prepared_unnormalized_insertion_interval':encode(prepared_insert),
        'exact_ground_input_bin_conditional_W_mean_interval':encode(exact_conditional),
        'actual_prepared_bin_conditional_W_mean_interval':encode(prepared_conditional),
        'actual_prepared_bin_conditional_E_plus_probability_interval':encode(tuple(F(1,2)+x/4 for x in prepared_conditional)),
        'five_complete_current_records':charged_histogram,
        'conditioned_preparation_trace_bound':str(delta/bin_weight[0]),
        'complete_law':['1-q','q tr(K_fail rho K_fail^dagger)',
          'q tr((I-P_I) B rho B^dagger (I-P_I))',
          'q tr(P_w U(2) P_I B rho B^dagger P_I U(2)^dagger P_w), every W eigenvalue w'],
        'current_measurement_changes_the_admitted_energy_distribution':True,
        'outside_bin_charged_spectrum_or_ground_identified':False},
      'apparatus':{'current_pointer_dimension':5,'binary_E_plus_minus_is_alternative_setting':True,
        'current_readout_duration':'1 source unit, position/current refers to interval start',
        'interaction_picture_completion':'V(s)=exp(-i Htotal s) L exp(i Htotal s), exp(-iL)=U_read',
        'pump_and_bin_duration':'one source unit each; all no-op branches retain the same autonomous Htotal',
        'source_work_upper':'2 ||Htotal|| per completed unitary service',
        'retained_storage_upper':'D^2 2^30 2^32 Mprep^2 8 5 2^3; Mprep is the reported Cayley length',
        'controller_and_gate_synthesis_workspace_included':False,
        'efficient_native_gate_synthesis_claimed':False},
      'scope':{'nonzero_same_H_dressed_weak_density_and_inelastic_record':True,
        'nonzero_same_H_energy_selected_charged_density_insertion':True,
        'bare_Gauss_generator_measured':False,'W_boson_or_Lorentz_vector_identified':False,
        'weak_scattering_cross_section_or_SM_matching':False,'continuum_or_GeV_calibration':False}}
