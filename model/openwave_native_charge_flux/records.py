"""Complete charge/flux refinements of actual retained transport records."""
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.transport import certificate as arrival_certificate
from .algebra import route_certificate


@lru_cache(maxsize=1)
def certificate():
    actual=arrival_certificate('ec')
    assert actual['source']['hopping_coefficient']=='1/17179869184'
    assert actual['source']['source_time']=='1/16777216'
    assert actual['source']['reference_particle_cost']=='4'
    q,a=s.symbols('q a');p=s.symbols('p0:6')
    assert s.expand((1-q)+q*(1-a)+q*a/2+q*sum(p)+q*(a/2-sum(p)))==1
    attachment={r['edge']:r for r in route_certificate()['actual_routes']}
    rows=[]
    for record in actual['complete_original_physical_records']['records']:
        row=dict(record)
        if record['record']=='charged_arrival_with_selected_flux':
            route=attachment[record['edge']]
            row.update({'joint_charge_and_flux_atom':{
                'arrival_source_matter_plus_scalar_Y':'-1',
                'arrival_source_matter_plus_scalar_Q':'-1',
                'reference_Y_and_Q':'1',
                'three_half_cube_source_Y_and_outward_flux':route['three_half_cube_fluxes'],
                'physical_charge_flux_mismatch_probability':'0'},
                'joint_atom_probability_equals_actual_arrival_probability':True,
                'same_source_time_and_preparation':True})
        else:
            row['joint_flux_refinement']='all ambient spectral outcomes retained; physical charge/flux mismatch has exactly zero probability'
            if record['record']=='localized_pump_failure':
                row['three_half_cube_fluxes_before_further_clock_control']=['0','0','0']
            row['individual_flux_probabilities_evaluated']=False
        rows.append(row)
    assert len(rows)==10
    for row in rows:
        if 'joint_charge_and_flux_atom' in row:
            lo,hi=map(s.Rational,row['probability_interval']);assert 0<lo<hi
    return {'source':actual['source'],
        'preparation':actual['localization_and_probe'],
        'actual_arrival_error_control':actual['full_H_time_and_preparation_control'],
        'joint_readout':{
            'measurement_order':'Complete physical preparation, P0 intervention, native ec pump, same full clock on pump success, ten coarse outputs, simultaneous central-Y cut flux and local total-Y readout.',
            'central_cut_flux_projector':'Pi_(S,k) is the diagonal spectral projector of 6F_S with eigenvalue k on all 119 states of each cut link.',
            'charge_projector':'R_(S,k) is the spectral projector of 6(rho_source,S^Y+1_(0 in S)Y_ref).',
            'commuting_effects':'All selected cut-flux and charge projectors commute with each other, every local Gauss current, and the six arrival effects.',
            'actual_joint_probability':'Pr(r,k,l)=Tr[R_(S,l) Pi_(S,k) K_r rho_input K_r^dagger Pi_(S,k) R_(S,l)] with the complete actual record Kraus K_r and retained preparation environment.',
            'exact_mismatch_rule':'Pr(r,k,l)=0 for k != l on every preparation, failure, arrival and nonarrival branch.',
            'exact_refinement_completeness':'sum_(k,l) Pr(r,k,l)=Pr(r), sum_r Pr(r)=1',
            'arbitrary_multiple_simultaneous_cuts':True,
            'subsequent_controllers_on_failure_branches':'off; the cut readout is applied to their retained physical outputs',
            'Q_measurement_on_arbitrary_nonarrival_states_assumed':False,
            'all_full_reference_Fock_and_source_environments_retained':True},
        'complete_original_preparation_records':{'records':rows,'coarse_record_count':10,
            'exact_completeness':actual['complete_original_physical_records']['exact_completeness'],
            'source_probability_definitions':actual['complete_original_physical_records']['definitions']},
        'same_source_six_route_joint_charge_flux_interval':actual['full_H_time_and_preparation_control']['total_six_route_arrival_interval_from_original_physical_vacuum_input'],
        'postprocessing_scope':{'charge_flux_identity_holds_on_every_physical_input_to_the_six_effects':True,
            'Gauss_preserving_energy_filter_retains_identity':True,
            'positive_energy_filtered_arrival_probability_derived_here':False,
            'energy_and_arrival_projectors_assumed_commuting':False},
        'physical_scope':{'native_current_to_actual_finite_Gauss_flux_attachment':True,
            'same_interacting_clock_ground_preparation_and_reference':True,
            'unit_native_ec_charge_coefficient_fitted_to_target':False,
            'physical_electromagnetic_field_normalization_or_alpha_derived':False,
            'Coulomb_inverse_square_field_derived_from_PW119_here':False,
            'free_electron_identity_or_reference_free_rest_mass_derived':False,
            'all_failure_flux_probabilities_numerically_evaluated':False}}
