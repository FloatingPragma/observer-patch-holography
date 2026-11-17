"""Kernel-controlled charged arrivals including physical preparation error."""
import json
from pathlib import Path
import sympy as s
from .algebra import coefficient,routes
from ..collective_spin_response.transfer import data as neutral_data
from ..interacting_charged_response.transfer import data as energy_data
from ..higgs_composites.certificate import validate_pins
from ..atomic_reference_control.lattice import local_norm
from ..atomic_momentum_response.model import reference_data
from ..interacting_charged_response.model import FAMILIES
from ..family_cell.quotient import sqrt_interval,outward

ROOT=Path(__file__).resolve().parents[4]
PREPARATION='model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'
TIME=s.Rational(1,2**24)


def arrival_bounds(moment,coupling,time,edge_norm,trial_error,tau,centered_atomic_norm,preparation_vector_error=0):
    lam,t,k,eta,norm,b,prep=map(s.sympify,(coupling,time,edge_norm,trial_error,tau,centered_atomic_norm,preparation_vector_error))
    if any(value.is_Rational is not True for value in (lam,t,k,eta,norm,b,prep)) or not(
        lam>0 and t>=0 and k>0 and eta>=0 and norm>=k and b>0 and prep>=0):
        raise ValueError('exact positive source norm controls and nonnegative time/vector errors required')
    lo,hi=map(s.Rational,moment)
    if not 0<lo<=hi<=k*k:raise ValueError('a positive complete-route moment interval is required')
    error=k*(eta+prep)+norm*(2*b+lam*norm)*t/2
    root_low,root_high=sqrt_interval((lo,hi),100)
    lower=max(0,root_low-error);upper=root_high+error
    return {'first_order_scaled_amplitude_error':s.cancel(error),
        'probability_interval':(s.cancel((lam*t*lower)**2),s.cancel((lam*t*upper)**2)),
        'strict_scaled_amplitude_margin':s.cancel(root_low-error)}


def certificate(family):
    neutral=neutral_data();local=coefficient(family);energy=energy_data(family)
    reference,_=reference_data()
    assert reference['native_irrep_ground_occupations'][family]['exact_ground_mean_occupation']==str(len(FAMILIES[family]))
    assert reference['mean_occupations_also_hold_on_each_pure_ground_spin_state']
    receipt=json.loads((ROOT/PREPARATION).read_text());validate_pins(receipt['input_sha256'])
    preparation=receipt['preparation']['finite_energy_agreement_filter']
    ground_fidelity=s.Rational(preparation['conditional_ground_pair_fidelity_lower'])
    sqrt_error=s.Rational(preparation['conditional_one_source_trace_distance_upper'])
    assert sqrt_error**2>=1-ground_fidelity
    prep_error=2*sqrt_error
    graph=neutral['errors']['invariant_graph_norm']
    eta=neutral['errors']['full_ground_vector_error_after_phase_alignment']
    centered_atomic=8*(local_norm()+12)+24
    assert centered_atomic==576
    controls=arrival_bounds(local['moment_interval'],neutral['coupling'],TIME,15,eta,
        neutral['tau'],centered_atomic,prep_error)
    exact=arrival_bounds(local['moment_interval'],neutral['coupling'],TIME,15,eta,
        neutral['tau'],centered_atomic,0)
    assert controls['strict_scaled_amplitude_margin']>0
    lower,upper=controls['probability_interval']
    qlower=s.Rational(preparation['success_probability_lower'])
    atomic_failure_upper=(graph+prep_error)**2
    atomic_success_lower=1-atomic_failure_upper
    total_arrival=(6*lower,6*upper)
    nonarrival_lower=atomic_success_lower/2-total_arrival[1]
    assert nonarrival_lower>0
    record_rows=[{'record':'physical_ground_preparation_failure',
        'exact_probability':'1-q','probability_interval':['0',str(1-qlower)]},
        {'record':'physical_preparation_success_atomic_localization_failure',
        'exact_probability':'q(1-a)','probability_interval':['0',str(atomic_failure_upper)]},
        {'record':'localized_pump_failure','exact_probability':'q a/2',
        'probability_interval':outward((qlower*atomic_success_lower/2,s.Rational(1,2)),160)}]
    for route in routes():
        record_rows.append({'record':'charged_arrival_with_selected_flux',**route,
            'exact_probability':f'q p_edge_{route["edge"]}',
            'probability_interval':outward((qlower*lower,upper),192),
            'static_reference_cell':0,'source_origin_and_other_nonarrival_cells_in_true_local_ground':True})
    record_rows.append({'record':'localized_pump_success_complete_nonarrival',
        'exact_probability':'q(a/2-sum_edge p_edge)',
        'probability_interval':outward((qlower*nonarrival_lower,s.Rational(1,2)),160)})
    return {'family':family,'source':{
            'Hamiltonian':'same full native N2 H0+lambda*T, plus the fixed static reference mu N_ref',
            'hopping_coefficient':str(neutral['coupling']),'reference_particle_cost':'4',
            'source_time':str(TIME),'source_clock':'exp(-it H_lambda)',
            'reference_clock':'exp(-4it N_ref)','initial_stationary_target':'actual unique full interacting neutral ground',
            'actual_input':'the inherited physical algebraic preparation with its complete success/failure law',
            'Hamiltonian_extended_by_a_hand_selected_charge_hop':False},
        'localization_and_probe':{
            'P0':'tensor of all eight actual local ground-doublet projectors, with every link blank',
            'complete_localization_effects':['P0','I-P0'],'failure_retained_and_subsequent_controls_off':True,
            'successful_probe':'the same native B_r applied after P0',
            'exact_P0_probe_effect':'P0 B_r^dagger B_r P0=P0/2',
            'exact_pump_success_and_failure_conditioned_on_localization':['1/2','1/2'],
            'atomic_localization_success_lower':str(atomic_success_lower),
            'atomic_localization_failure_upper':str(atomic_failure_upper),
            'interacting_ground_silently_replaced_by_atomic_trial':False},
        'arrival_detector':{
            'six_pairwise_orthogonal_effects':'selected incident link excited, all other links blank; actual local ground at every nonarrival cell; identity on every arrival-cell CAR/scalar state',
            'complete_complement':'I-sum_of_the_six_arrival_effects on pump success',
            'reference_cell_fixed':0,'arrival_cells':[1,2,4],
            'explicit_quantum_flux_and_reference_retained':True,
            'all_uncompensated_source_Gauss_and_joint_origin_constraints_preserved':True,
            'exact_no_operation_arrival_probability':'0 for every preparation input after P0',
            'exact_no_hopping_arrival_probability':'0 at every source time',
            'arrival_is_odd_source_charge_motion':True},
        'full_H_time_and_preparation_control':{
            'centered_full_atomic_H_norm_upper':str(centered_atomic),
            'full_hopping_norm_upper':str(neutral['tau']),'selected_edge_norm_upper':'15',
            'exact_kernel':'P_arrival exp(-it H0) B_r P0=0 on the complete input carrier',
            'ground_trial_vector_error':str(eta),'inherited_ground_fidelity_lower':str(ground_fidelity),
            'retained_environment_purification_vector_error_upper':str(prep_error),
            'preparation_error_subtracted_as_global_probability_error':False,
            'scaled_arrival_amplitude_error':str(controls['first_order_scaled_amplitude_error']),
            'strict_scaled_amplitude_margin':str(controls['strict_scaled_amplitude_margin']),
            'exact_ground_per_route_arrival_interval':outward(exact['probability_interval'],192),
            'prepared_per_route_arrival_interval_conditioned_on_physical_preparation_success':outward((lower,upper),192),
            'per_route_arrival_interval_from_original_physical_vacuum_input':outward((qlower*lower,upper),192),
            'total_six_route_arrival_interval_from_original_physical_vacuum_input':outward((qlower*total_arrival[0],total_arrival[1]),192)},
        'complete_original_physical_records':{'records':record_rows,'coarse_record_count':10,
            'definitions':'q=physical preparation success, a=P0 success conditional on q, p_edge=route probability conditional on q before normalization by localization/pump',
            'exact_completeness':'(1-q)+q(1-a)+qa/2+q sum_edge p_edge+q(a/2-sum_edge p_edge)=1',
            'all_preparation_and_detector_failure_outputs_retained':True,
            'all_source_copies_and_quantum_references_retained':True},
        'same_preparation_energy_setting':{
            'preparation_operator':'B_r P0 on the same prepared interacting source',
            'same_total_energy_bin':list(map(str,energy['window'])),
            'source_energy_difference_bin_after_explicit_reference_cost_subtraction':outward(energy['source_excitation'],64),
            'conditional_line_interval':outward((max(0,energy['probabilities']['line'][0]-sqrt_error),
                min(1,energy['probabilities']['line'][1]+sqrt_error)),64),
            'conditional_outside_interval':outward((max(0,energy['probabilities']['outside'][0]-sqrt_error),
                min(1,energy['probabilities']['outside'][1]+sqrt_error)),64),
            'energy_and_arrival_are_distinct_observer_settings':True,
            'exact_completeness':'(1-q)+q(1-a)+qa/2+q(p_line+p_outside)=1, with p_line+p_outside=a/2',
            'joint_commuting_energy_location_measurement_or_energy_resolved_arrival_inferred':False},
        'scope':{'strict_unconditional_same_H_charge_arrival_lower_bound':True,
            'source_clock_and_static_compensator_preserved':True,
            'complete_flux_route_and_nonarrival_outputs':True,
            'free_charged_dispersion_momentum_band_or_continuum_limit':False,
            'rest_mass_or_relativistic_velocity_identification':False}}
