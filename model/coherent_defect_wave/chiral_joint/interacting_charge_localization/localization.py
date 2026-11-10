"""Rational full-H bounds on single-unscreened-charge records."""
from functools import lru_cache
import sympy as s
from .geometry import graph
from ..interacting_charged_response.transfer import data as charged_data
from ..interacting_charged_transport.transport import certificate as transport_certificate
from ..interacting_charged_response.model import REFERENCE_COST
from ..family_cell.quotient import outward


def ordered_projection_bound(atomic_lower,energy_upper,perturbation):
    a,b,v=map(s.sympify,(atomic_lower,energy_upper,perturbation))
    if any(x.is_Rational is not True for x in (a,b,v)) or not(a>b and v>=0):
        raise ValueError('ordered exact spectral bounds and a nonnegative exact perturbation required')
    return min(s.Integer(1),(v/(a-b))**2)


@lru_cache(maxsize=2)
def data(family):
    d=charged_data(family);ground=d['local']['e_interval'];v=d['v'];mu=REFERENCE_COST
    unrestricted=s.Rational(-1171,100);delta=ground[1]-unrestricted
    charged_gap=d['sector']['complete_compensated_local_gap_lower']
    assert 0<delta<1 and 0<charged_gap<1-delta
    gamma=1-2*delta
    assert gamma>0
    _,distance=graph();rows=[]
    for site,ell in sorted(distance.items()):
        if not site:continue
        lower=8*ground[0]+ell+charged_gap+mu
        bound=ordered_projection_bound(lower,d['window'][1],v)
        rows.append({'cell':site,'distance':ell,'absolute_atomic_energy_lower':lower,
            'total_energy_cost_above_actual_neutral_ground_lower':mu+ell+charged_gap-v,
            'source_energy_cost_after_explicit_reference_subtraction_lower':ell+charged_gap-v,
            'spectral_separation_lower':lower-d['window'][1],
            'conditional_probability_upper':bound})
    return {'charged':d,'local_unrestricted_lower':unrestricted,'delta':delta,'gamma':gamma,
        'charged_gap':charged_gap,'rows':rows,'remote_upper':sum(row['conditional_probability_upper'] for row in rows)}


def certificate(family):
    d=data(family);source=d['charged'];prep=transport_certificate(family)
    energy=prep['same_preparation_energy_setting'];lower,upper=map(s.Rational,energy['conditional_line_interval'])
    q=s.Rational(source['preparation']['finite_energy_agreement_filter']['success_probability_lower'])
    a=s.Rational(prep['localization_and_probe']['atomic_localization_success_lower'])
    rows=[{'record':'physical_ground_preparation_failure','probability':'1-q'},
        {'record':'P0_localization_failure','probability':'q(1-a)'},
        {'record':'localized_pump_failure','probability':'qa/2'},
        {'record':'energy_bin_complement','probability':'q(a/2-b)'},
        {'record':'single_unscreened_charge_at_origin','probability':'qb w_0'}]
    rows.extend({'record':'single_unscreened_charge_at_remote_cell','cell':row['cell'],
        'probability':f'qb w_{row["cell"]}',
        'unconditional_probability_upper':outward((0,upper*row['conditional_probability_upper']),128)[1]}
        for row in d['rows'])
    rows.append({'record':'complete_multiple_charged_cell_complement','probability':'qb w_cloud'})
    return {'family':family,'source':{'Hamiltonian':'same full native N2 H0+lambda*T+mu N_ref',
        'lambda':str(source['neutral']['coupling']),'mu':str(REFERENCE_COST),'full_hopping_perturbation_norm_upper':str(source['v']),
        'reference_family_one_particle_sector':family,'reference_cell':0,
        'clock':'exp(-it(H_lambda+mu N_ref))','clock_times':'every real t'},
        'single_unscreened_charge_projectors':{
            'definition':'Pi_y requires full gauge-singlet local matter at every x != y, retaining every link and the entire local CAR/scalar space at y',
            'domain':'joint physical source plus fixed one-reference-family sector',
            'all_eight_Pi_y_pairwise_orthogonal':True,'complete_complement_retained':True,
            'commutes_with_H0':True,'exact_compressed_hopping':'Pi_y T Pi_y=0',
            'matter_at_y_hypercharge':'minus the fixed reference hypercharge',
            'origin_source_only_Gauss_zero_asserted':False},
        'full_source_energy_bounds':{
            'local_unrestricted_lower':str(d['local_unrestricted_lower']),
            'neutral_minus_unrestricted_upper':str(d['delta']),
            'complete_compensated_local_gap_lower':str(d['charged_gap']),
            'extra_link_margin':str(1-d['delta']-d['charged_gap']),
            'single_unscreened_charge_compression':'Pi_y H0 Pi_y >= (8E0+distance(0,y)+delta_r) Pi_y, y != 0',
            'all_state_link_count_slope_lower':str(d['gamma']),
            'all_state_atomic_bound':'H0 >= 8E0+gamma N_nonblank+delta_r I_origin_has_no_incident_nonblank_link',
            'all_state_interacting_relative_bound':'H_lambda+mu-E_ground >= mu+gamma N_nonblank-2v',
            'all_state_distance_from_compensator_bound_asserted':False},
        'energy_filtered_localization':{
            'same_absolute_total_energy_bin':list(map(str,source['window'])),
            'projection_estimate':'||Pi_y P_I(H)||^2 <= (v/(atomic_lower_y-upper(I)))^2',
            'finite_bin_complement_or_higher_states_truncated':False,
            'every_normalized_state_in_bin_and_every_time':True,
            'per_cell':[dict(cell=row['cell'],distance=row['distance'],
                **{name:str(value) for name,value in row.items() if name not in ('cell','distance')},
                conditional_probability_upper_outward=outward((0,row['conditional_probability_upper']),128)[1]) for row in d['rows']],
            'union_remote_single_unscreened_charge_probability_upper':outward((0,d['remote_upper']),128)[1],
            'origin_probability_lower_inferred_by_discarding_screened_complement':False},
        'sequential_complete_original_preparation_records':{
            'measurement_order':'physical preparation; P0; native pump; total-energy bin; same clock; complete single-charge-location measurement',
            'joint_commuting_energy_and_position_measurement_claimed':False,
            'preparation_success_lower':str(q),'P0_success_conditional_lower':str(a),
            'energy_bin_success_conditional_on_preparation_interval':list(map(str,(lower,upper))),
            'unconditional_energy_bin_success_lower':str(q*lower),
            'records':rows,'coarse_record_count':len(rows),
            'definitions':'b is joint P0/pump/bin probability conditional on physical preparation; w_0+...+w_7+w_cloud=1 after bin selection',
            'exact_completeness':'(1-q)+q(1-a)+qa/2+q(a/2-b)+qb(sum_y w_y+w_cloud)=1',
            'all_preparation_failures_and_reference_systems_retained':True,
            'localization_bound_depends_on_purity_or_preparation_fidelity':False},
        'scope':{'linear_unscreened_string_energy_cost_on_actual_finite_source':True,
            'uniform_time_low_bin_charge_localization':True,'dynamical_screening_excluded':False,
            'thermodynamic_confinement_or_asymptotically_free_charged_particle_proved':False,
            'finite_regulator_mobile_charged_band_proved':False,
            'sharp_continuum_charged_mass_pole_required_for_Coulomb_physics':False}}
