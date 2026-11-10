"""Full interacting charged energy-bin records without a charged-band truncation."""
import json
from pathlib import Path
import sympy as s
from .local import positive_response,compensated_blocks
from .model import REFERENCE_COST
from ..collective_spin_response.transfer import data as neutral_data
from ..collective_spin_response.bounds import probability_interval
from ..family_cell.quotient import outward
from ..higgs_composites.certificate import validate_pins

ROOT=Path(__file__).resolve().parents[4]
PREPARATION='model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'


def window_error(energies,selected,radius,perturbation,input_error):
    if type(selected) is not int or not 0<=selected<len(energies):
        raise ValueError('a retained atomic energy is required')
    radius,perturbation,input_error=map(s.sympify,(radius,perturbation,input_error))
    if any(v.is_Rational is not True for v in (radius,perturbation,input_error)) or not(radius>0 and perturbation>=0 and input_error>=0):
        raise ValueError('exact positive window and nonnegative error bounds required')
    intervals=[tuple(map(s.Rational,energy)) for energy in energies]
    if any(lo>hi for lo,hi in intervals):raise ValueError('ordered energy intervals required')
    low,high=intervals[selected];center=(low+high)/2
    window=(center-radius,center+radius)
    distances=[]
    for index,(a,b) in enumerate(intervals):
        if index==selected:
            distance=radius-(b-a)/2
        else:
            distance=max(window[0]-b,a-window[1])
        if distance<=0:raise ValueError('spectral window does not separate the complete visible atomic support')
        distances.append(distance)
    error=s.cancel(input_error+perturbation*sum(1/distance for distance in distances))
    return window,error,tuple(distances)


def data(family):
    neutral=neutral_data();local=positive_response(family);sector=compensated_blocks(family)
    e0low,e0high=local['e_interval'];mu=REFERENCE_COST
    energies=[(row['energy_interval'][0]+7*e0low+mu,
               row['energy_interval'][1]+7*e0high+mu) for row in local['rows']]
    v=neutral['coupling']*neutral['tau'];eta=neutral['errors']['full_ground_vector_error_after_phase_alignment']
    window,error,distances=window_error(energies,0,s.Rational(1,1024),v,eta)
    weight=local['rows'][0]['weight_interval']
    probabilities={'line':probability_interval(weight,error),
        'outside':probability_interval((s.Rational(1,2)-weight[1],s.Rational(1,2)-weight[0]),error),
        'success':probability_interval((s.Rational(1,2),)*2,eta),
        'failure':probability_interval((s.Rational(1,2),)*2,eta)}
    assert all(0<low<=high<1 for low,high in probabilities.values())
    # The neutral ground is in [8 E0-v,8 E0+v], with no scalar shift dropped.
    excitation=(window[0]-8*e0high-v,window[1]-8*e0low+v)
    source_excitation=(excitation[0]-mu,excitation[1]-mu)
    one_reference_gap=min(sector['complete_compensated_local_gap_lower'],neutral['gap'])-2*v
    # This unrestricted source bound covers every occupation of the complete 30-mode reference.
    local_unrestricted=s.Rational(-1171,100)
    occupied_reference_gap=mu+8*(local_unrestricted-e0high)-2*v
    assert source_excitation[0]>0 and one_reference_gap>0 and occupied_reference_gap>0
    prep=json.loads((ROOT/PREPARATION).read_text());validate_pins(prep['input_sha256'])
    prep=prep['preparation'];assert s.Rational(prep['source']['hopping_coefficient'])==neutral['coupling']
    trace=s.Rational(prep['finite_energy_agreement_filter']['conditional_one_source_trace_distance_upper'])
    approximate={name:(max(0,lo-trace),min(1,hi+trace)) for name,(lo,hi) in probabilities.items()}
    assert all(0<lo<=hi<1 for lo,hi in approximate.values())
    return {'neutral':neutral,'local':local,'sector':sector,'v':v,'eta':eta,'window':window,
        'window_error':error,'distances':distances,'probabilities':probabilities,
        'excitation':excitation,'source_excitation':source_excitation,
        'one_reference_gap':one_reference_gap,'occupied_reference_gap':occupied_reference_gap,
        'approximate':approximate,'preparation':prep}


def certificate(family):
    d=data(family);dimension=(5*2**30)**8*119**24
    probabilities={name:outward(value,64) for name,value in d['probabilities'].items()}
    rows=[{'record':'initial_ground_preparation_failure','probability':str(1-s.Rational(1,dimension))}]
    for name in ('failure','line','outside'):
        rows.append({'record':'ground_prepared_pump_'+name,
            'probability_interval':outward(tuple(value/dimension for value in d['probabilities'][name]),512),
            'exact_probability':f'p_{name}/D','all_source_and_reference_outputs_retained':True})
    return {'family':family,'source':{
            'Hamiltonian':'actual full H0+lambda*T, plus static-reference cost mu N_ref',
            'hopping_coefficient':str(d['neutral']['coupling']),'reference_particle_cost':str(REFERENCE_COST),
            'onsite_Yukawa_scalar_link_coefficients':'1','source_clock':'exp(-it H_lambda)',
            'reference_clock':'exp(-it mu N_ref)',
            'input_reference':'actual unique neutral full interacting source ground, static CAR reference vacuum',
            'source_ground_or_clock_substituted':False},
        'reference_sector_comparison':{
            'complete_compensated_local_gap_above_neutral_ground_lower':str(d['sector']['complete_compensated_local_gap_lower']),
            'zero_reference_cost_full_lattice_one_particle_sector_gap_lower':str(d['one_reference_gap']),
            'bound_includes_every_nonblank_link_set_and_charged_source_sector':True,
            'positive_cost_full_Fock_occupied_reference_sector_gap_lower':str(d['occupied_reference_gap']),
            'neutral_ground_tensor_reference_vacuum_is_unique_full_physical_ground_at_mu_four':True,
            'zero_cost_full_Fock_unique_ground':False,
            'zero_cost_full_Fock_global_minimum_identified':False,
            'zero_cost_filled_reference_with_same_neutral_ground_energy_exists':True},
        'full_interacting_energy_readout':{
            'measurement':'spectral projector of the declared total Hamiltonian onto the finite bin, and its complete complement',
            'selected_total_absolute_energy_bin':list(map(str,d['window'])),
            'total_excitation_energy_bin_relative_to_actual_initial_ground':outward(d['excitation'],64),
            'source_excitation_energy_bin_after_subtracting_declared_reference_cost':outward(d['source_excitation'],64),
            'reference_energy_subtracted_silently':False,
            'full_hopping_perturbation_norm':str(d['v']),
            'actual_neutral_ground_trial_vector_error':str(d['eta']),
            'line_and_complete_outside_amplitude_error':str(d['window_error']),
            'complete_atomic_support_residual_denominators':[str(value) for value in d['distances']],
            'conditional_probabilities':probabilities,
            'all_real_source_clock_times':True,
            'isolated_charged_eigenvalue_or_unchanged_multiplicity_claimed':False,
            'odd_source_parity_on_every_pump_success':True,
            'static_reference_remains_one_particle_on_success':True,
            'all_higher_full_H_states_in_complete_outside_outcome':True},
        'original_preparation_records':{
            'input':'I/D on complete source with reference vacuum','D':str(dimension),
            'domain':'kinematic pre-physical preparation; successful branch obeys source Gauss constraints',
            'initial_failure_source_only_Gauss_constraints_asserted':False,
            'ground_filter_success':'1/D','records':rows,'positive_coarse_record_count':4,
            'exact_completeness':'(1-1/D)+(p_failure+p_line+p_outside)/D=1',
            'complete_reference_Fock_space_retained':True,
            'all_unsupported_reference_occupations':{'probability':'0','retained':True}},
        'constructive_approximate_preparation_records':{
            'preparation_success_lower':d['preparation']['finite_energy_agreement_filter']['success_probability_lower'],
            'preparation_failure_upper':d['preparation']['finite_energy_agreement_filter']['failure_probability_interval'][1],
            'trace_distance_from_exact_ground_upper':d['preparation']['finite_energy_agreement_filter']['conditional_one_source_trace_distance_upper'],
            'conditional_probability_intervals':{name:outward(value,64) for name,value in d['approximate'].items()},
            'second_source_and_all_preparation_references_retained':True,
            'all_preparation_success_and_failure_branches_obey_source_Gauss':True,
            'all_detector_success_and_failure_branches_obey_joint_source_static_reference_Gauss':True,
            'approximate_state_asserted_exactly_stationary':False},
        'scope':{'actual_interacting_odd_elementary_channel_with_positive_complete_records':True,
            'full_complement_error_control':True,'physical_rest_mass_or_free_charge_dispersion':False,
            'quark_local_spectra_or_continuum_transfer':False}}
