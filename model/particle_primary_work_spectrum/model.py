"""Actual phase work, modular fluctuation bounds and finite distinguishability."""
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_critical_refinement import model as critical
from model.particle_primary_gravity_premises import model as modular
from model.particle_primary_local_fluctuation import model as local
from model.particle_local_wilson_density import density
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
 'model/particle_primary_work_spectrum/__init__.py',
 'model/particle_primary_work_spectrum/__main__.py',
 'model/particle_primary_work_spectrum/model.py',
 'model/particle_primary_work_spectrum/test_model.py',
 'model/particle_primary_work_spectrum/PROOF.md',
 'model/particle_primary_work_spectrum/README.md',
 'model/particle_primary_critical_refinement/PROOF.md',
 'model/particle_primary_current_tail/PROOF.md',
 'model/particle_primary_local_fluctuation/PROOF.md',
 'model/particle_primary_gravity_premises/PROOF.md',
 'model/particle_primary_charge_conjugation/PROOF.md',
 'model/particle_local_wilson_density/PROOF.md',
 'Lean/Thermodynamics/FluctuationTheorems.lean',
 'Lean/QFT/GroundRecordTail.lean',
)


def current_second_moment(member):
    """Full-CAR trace norm at each parent gauge coordinate, then compression."""
    contract.verify_member(member);source=member['source'];n=source['level']
    R=2**n;a=F(source['spacing']);V=source['side']**3;B=source['interior_hopping_bonds']
    sites=4 if n==0 else (2*R-1)**2*(2*R+5)
    geometry=critical.profile_geometry(n)
    d2=geometry['three_axis_gradient_square_sum']
    parent_upper=3000*sites*F(B,V)*d2
    if n==0:parent_upper=1800*sites*F(B,V)*4/(a*a)
    inherited=F(19620) if n==0 else F(324000*R**4)
    chosen=min(parent_upper,inherited)
    return {'active_source_sites_upper':sites,'one_particle_rank_upper':60*sites,
        'one_body_square_trace_formula':'120/a^2 * sum_b (1-cos(delta_b))',
        'origin_averaged_full_CAR_second_moment_upper':parent_upper,
        'inherited_second_moment_upper':inherited,
        'selected_second_moment_upper':chosen,
        'exact_work_relative_bound_formula':'M2 <= 1800*M*B*W/(a*V*h_n)',
        'uniform_second_moment_proved':False,
        'bound_growth_is_O_a_minus_four':True,
        'parent_anchor_substituted_for_selected_density':False}


def downward_tail(beta_lower,threshold):
    b=ground._q(beta_lower,'positive modular slope',True)
    E=ground._q(threshold,'positive source energy threshold',True)
    x=b*E;k=x.numerator//x.denominator
    return {'threshold':E,'exact_upper_formula':'1/(1+exp(beta_low*threshold))',
        'rational_upper_factored':{'numerator':1,'denominator_form':'1+2^exponent','exponent':k},
        'rounding_is_outward':True}


def polynomial_control(member):
    """A distinct two-level carrier uses literal source polynomial coefficients."""
    contract.verify_member(member);a=F(member['source']['spacing'])
    slope=modular.selected_modular_slope(member);beta=F(slope['beta']);low=F(slope['modular_slope_lower'])
    source=member['state'];upper=F(source['spectral_upper']);lower=F(source['spectral_lower'])
    gap=1/a;eps=F(source['target_ground_excess'])
    if gap>upper-lower or 1/low>eps or low<2*beta/3:
        raise ArithmeticError('polynomial control left its proved source-parameter domain')
    return {'same_carrier_or_source_as_primary':False,
        'copied_polynomial_parameters':{'beta':beta,'Taylor_order':source['Taylor_order'],
            'dyadic_scale':source['dyadic_scale'],'spectral_upper':upper},
        'two_energy_levels':[upper-gap,upper],
        'negative_log_excited_to_ground_probability_ratio_interval':[low/a,beta/a],
        'selected_two_level_ground_excess_upper':1/low,
        'primary_target_ground_excess':eps,
        'phase_generator':'sqrt(a)*sigma_x; C=sigma_z reverses its sign',
        'uniform_double_commutator_norm':2,
        'original_transfer_probabilities':'sin(sqrt(a))^2 * (p_ground at +1/a, p_excited at -1/a); remaining mass at0',
        'work_limit':1,'second_moment_asymptotic':'1/a',
        'modular_work_divided_by_beta_liminf_lower':F(2,3),
        'affinity_limit':1,'all_fixed_positive_energy_weight_limit':0,
        'fixed_finite_unitary_history_operator_limit':'identity',
        'countercontrol_is_primary_source_collapse_proof':False}


def selected_bound(member):
    contract.verify_member(member);source=member['source'];n=source['level']
    a=F(source['spacing']);eps=F(member['state']['target_ground_excess'])
    upstream=critical.bind_member(member);p=upstream['result']['profiles'][0]
    slope=modular.selected_modular_slope(member);low=F(slope['modular_slope_lower']);beta=F(slope['beta'])
    width=F(member['state']['spectral_upper'])-F(member['state']['spectral_lower'])
    norm=max(abs(F(member['state']['spectral_upper'])),abs(F(member['state']['spectral_lower'])))
    mean_lower=F(p['original_signed_work_lower']);mean_upper=F(p['original_signed_work_upper'])
    if n==0:mean_upper=min(mean_upper,F(180))
    clock_error=F(upstream['result']['complete_energy_clock_signed_mean_error_upper'])
    negative=min(eps,1/(3*low))
    if not negative<eps or mean_lower<=clock_error:raise ArithmeticError('actual improved tail or record contrast failed')
    second=current_second_moment(member)
    window=local.positive_window(mean_lower,second['selected_second_moment_upper'])
    clock=upstream['result']['complete_two_energy_clock'];eta=F(clock['energy_resolution'])
    clock_loss=2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper'])
    lo,hi=window['original_transfer_interval'];probability=window['original_weight_lower']
    if lo<=2*eta or clock_loss>=probability/2:
        raise ArithmeticError('common critical clock does not resolve the original positive window')
    reported={'reported_transfer_interval':[lo-2*eta,hi+2*eta],
        'original_probability_lower':probability-clock_loss,
        'two_complete_clock_probability_loss_upper':clock_loss,
        'same_critical_two_clock_plan_sha256':contract.digest(clock),
        'energy_resolution':eta,
        'all_fine_records_retained':2*source['side']**3*clock['clock_plan']['clock_size']**2,
        'conditioned_or_postselected_probability':False,
        'full_joint_clock_histogram_executed':False}
    return contract.serial({'schema':'oph.primary_work_spectrum_result.v1','member_id':member['member_id'],**member['ids'],
        'critical_work_bound_sha256':upstream['result_sha256'],
        'actual_modular_slope':slope,'original_signed_work_interval':[mean_lower,mean_upper],
        'same_source_state_clock_and_retained_reference':True,
        'actual_unitary_menu':'all unit-tent origins and both signs, each with probability1/(2V)',
        'modular_fluctuation_identity':'E exp(-[K(E_final)-K(E_initial)])=1',
        'source_affinity_identity':'A=average Tr(sqrt(rho) U sqrt(rho) U_dagger)=E exp(-DeltaK/2)',
        'mean_modular_work_interval':[low*mean_lower,beta*mean_upper],
        'negative_absolute_first_moment_upper':negative,
        'negative_tail_at_half':downward_tail(low,F(1,2)),
        'original_absolute_first_moment_upper':mean_upper+2*negative,
        'source_CQ_half_trace_distance_from_no_operation_lower':mean_lower/width,
        'source_CQ_affinity_deficit_lower':low*mean_lower/(4+low*width),
        'finite_complete_clock_record_total_variation_lower':(mean_lower-clock_error)/(56*norm),
        'bounds_deteriorate_with_capacity':True,
        'affinity_is_untouched_joint_purification_overlap':False,
        'classical_energy_kernel_symmetric_and_doubly_stochastic':True,
        'finite_fluctuation_theorem_reference':'uniform energy-index law',
        'finite_fluctuation_entropy_variable':'log p_initial-log q_final, q=push(p,T)',
        'modular_work_and_classical_entropy_production_are_equal':False,
        'source_entropy_balance':'mean DeltaK = H(q)-H(p)+D(q||p)',
        'second_moment':second,'ideal_original_positive_window':window,
        'finite_reported_positive_window':reported,
        'uniform_selected_energy_weighted_tail_proved':False,
        'nonelastic_fixed_energy_limit_proved':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_work_spectrum_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('work-spectrum source or bound changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*critical.SOURCE_INPUTS,*modular.SOURCE_INPUTS,*local.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    members=[contract.member(n) for n in (0,1,2,8)]
    result={'schema':'oph.primary_work_spectrum_certificate.v1','members':[bind_member(x) for x in members],
        'distinct_literal_polynomial_controls':[contract.serial(polynomial_control(x)) for x in members],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('work-spectrum inputs changed during replay')
    return result
