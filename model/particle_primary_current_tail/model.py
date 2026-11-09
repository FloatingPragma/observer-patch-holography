"""Actual source double commutators and common mixed-word energy bounds."""
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_local_fluctuation import model as fluctuation
from model.particle_primary_spectral_scale import scale
from model.particle_primary_gravity_premises import model as modular
from model.particle_spatial_ground_state import ground
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_current_tail/__init__.py',
    'model/particle_primary_current_tail/__main__.py',
    'model/particle_primary_current_tail/model.py',
    'model/particle_primary_current_tail/test_model.py',
    'model/particle_primary_current_tail/PROOF.md',
    'model/particle_primary_current_tail/README.md',
    'model/particle_primary_local_fluctuation/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'Lean/QFT/GroundRecordTail.lean',
    'Lean/QFT/ColimitGNSRepresentation.lean',
)


def source_form_bound(member):
    contract.verify_member(member)
    n=member['source']['level'];a=F(member['source']['spacing'])
    R=2**n;L=member['source']['side']
    if a*R!=1 or L!=2*R*R:raise ArithmeticError('the source scaling changed')
    D=fluctuation.tent_sums(R)['one_axis_three_dimensional_gradient_square_sum']
    # Full Fock hopping norm30/a, three axes, pi²<10. At R=1
    # each open-cube corner has3 incident bonds, not6 wrapping edges.
    double=F(900) if n==0 else 900*a*D
    uniform=F(2025,2)
    if double>uniform:raise ArithmeticError('source double commutator is not uniform')
    return {'double_commutator_norm_upper':double,
        'uniform_double_commutator_norm_upper':uniform,
        'energy_square_root_drift_upper':sqrt_upper(double/2),
        'uniform_energy_square_root_drift_upper':F(45,2),
        'same_source_energy_zero_used':True,
        'ground_subtraction_is_only_in_the_bound':True,
        'source_or_state_or_clock_changed':False}


def mixed_word_bound(member,weyl_count=2,filter_widths=('1',)):
    contract.verify_member(member)
    k=ground._i(weyl_count,'number of Weyl factors',0)
    widths=tuple(ground._q(x,'fixed triangular filter width',True) for x in filter_widths)
    if k+len(widths)==0:raise ValueError('at least one physical factor is required')
    form=source_form_bound(member);epsilon=F(member['state']['target_ground_excess'])
    root=sqrt_upper(epsilon)
    old_drifts=tuple(2*sqrt_upper(2/T) for T in widths)
    drift=k*form['energy_square_root_drift_upper']+sum(old_drifts,F())
    uniform_drift=k*F(45,2)+sum(old_drifts,F())
    final=(root+drift)**2;uniform_final=(F(1,4)+uniform_drift)**2
    absolute=final+epsilon
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_form_bound':form,'Weyl_factors':k,'triangular_filter_widths':widths,
        'source_time_translates_and_adjoints_allowed':True,
        'mixed_word_operator_norm_upper':1,
        'original_spectral_mass_upper':1,
        'original_spectral_mass_is_exactly_one':len(widths)==0,
        'selected_ground_excess_upper':epsilon,
        'fixed_word_energy_square_root_drift_upper':drift,
        'fixed_word_final_ground_relative_energy_mean_upper':final,
        'fixed_word_absolute_transfer_first_moment_upper':absolute,
        'uniform_final_ground_relative_energy_mean_upper':uniform_final,
        'uniform_absolute_transfer_first_moment_upper':uniform_final+F(1,16),
        'downward_probability_numerator_upper':epsilon,
        'cyclic_time_squared_distance_coefficient_upper':2*absolute,
        'uniform_cyclic_time_squared_distance_coefficient_upper':2*(uniform_final+F(1,16)),
        'time_bound_formula':'||[alpha_t(W)-W]||_rho^2 <= min(4, coefficient*abs(t))',
        'common_unfiltered_Weyl_and_filtered_probe_time_GNS_constructed':True,
        'canonical_GNS_energy_nonnegative':True,
        'point_norm_continuity_of_unfiltered_sections_proved':False,
        'spatial_continuity_or_momentum_for_unfiltered_Weyls_proved':False,
        'finite_same_source_preparation_and_clock_retained':True,
        'retained_purification_joint_stationarity_assumed':False,
        'full_growing_reference_trace_norm_limit_proved':False,
        'effective_common_subsequence_rate_proved':False,
        'uniform_fixed_positive_energy_weight_proved':False,
        'physical_mass_table_GeV':[]})


def random_sign_history(member,pulses=3):
    contract.verify_member(member);k=ground._i(pulses,'number of pulses',1)
    a=F(member['source']['spacing']);R=1/a;n=member['source']['level']
    axis=fluctuation.tent_sums(int(R))['one_axis_three_dimensional_gradient_square_sum']
    cap=F(180) if n==0 else 450*a*axis
    epsilon=F(member['state']['target_ground_excess'])
    return contract.serial({'pulses':k,'retained_sign_histories':2**k,
        'each_original_sign_history_probability':F(1,2**k),
        'signs_independent_of_source_input':True,
        'common_origin_retained':True,'source_time_evolutions_allowed_between_pulses':True,
        'averaged_final_ground_excess_upper':epsilon+k*cap,
        'uniform_averaged_final_ground_excess_upper':F(1,16)+k*F(2025,4),
        'all_reference_and_sign_outputs_retained':True,
        'history_histogram_executed':False})


def physical_anchor_residual(member):
    contract.verify_member(member);a=F(member['source']['spacing']);n=member['source']['level']
    L=member['source']['side'];V=L**3;B=member['source']['interior_hopping_bonds']
    rows=scale.vacuum_channels()
    if F(rows['bond_vacuum_variance'])!=15:raise ArithmeticError('the native CAR packet changed')
    D=fluctuation.tent_sums(2**n)['one_axis_three_dimensional_gradient_square_sum']
    k_low=2*a*a*D;k_high=5*a*a*D
    if n==0:k_low=k_high=F(4)
    low=30*B*k_low/(a*a*V);high=30*B*k_high/(a*a*V)
    if low<F(160)/a:raise ArithmeticError('the physical anchor ultraviolet lower bound failed')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'anchor':'filled negative-beta determinant with constant quantum links',
        'anchor_is_physical':True,'anchor_is_selected_state':False,
        'Weyl_unitary_fixes_anchor':True,
        'all_five_native_blocks_retained':True,'bond_variance':15,
        'exact_averaged_commutator_square_formula':'30*B*k_n/(a^2*V)',
        'averaged_anchor_commutator_square_interval':[low,high],
        'uniform_anchor_lower_formula':'160/a',
        'uniform_anchor_lower':F(160)/a,
        'maximum_address_operator_commutator_norm_squared_at_least':low,
        'uniform_absolute_operator_second_moment_bound_possible':False,
        'selected_state_weighted_second_moment_lower_inferred':False,
        'source_noncollapse_or_no_go_inferred':False})


def modular_phase_response(member):
    """Complete sign/origin relative entropy on the source marginal only."""
    contract.verify_member(member)
    slope=modular.selected_modular_slope(member)
    y=F(slope['scaled_argument_upper']);order=slope['Taylor_order']
    delta=F(slope['relative_modular_slope_defect_upper'])
    if not order>=1 or not 0<=y<=F(1,2):raise ArithmeticError('selected modular domain changed')
    if not delta<=y/(order+y)<=F(1,3):raise ArithmeticError('uniform source modular slope failed')
    response=fluctuation.selected_bound(member)
    mean_lower=F(response['original_signed_first_moment_lower'])
    mean_upper=F(response['original_signed_first_moment_upper'])
    uniform_lower=F(response['uniform_original_signed_first_moment_lower'])
    uniform_upper=F(response['uniform_original_signed_first_moment_upper'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'comparison':'complete classical sign/origin record plus source marginal',
        'unchanged_reference_state':'uniform sign/origin register tensor selected source rho',
        'retained_sign_origin_blocks':response['complete_origin_sign_outcomes'],
        'same_source_mean_response_digest':contract.digest(response),
        'beta':slope['beta'],'Taylor_order':order,'scaled_argument_upper':y,
        'relative_modular_slope_defect_upper':delta,
        'uniform_relative_modular_slope_defect_upper':F(1,3),
        'normalized_complete_CQ_relative_entropy_interval':[(1-delta)*mean_lower,mean_upper],
        'uniform_normalized_complete_CQ_relative_entropy_interval':[F(2,3)*uniform_lower,uniform_upper],
        'normalization':'D(rho_CQ_output || uniform_record tensor rho_source) / beta',
        'capacity_uniform_positive_source_marginal_modular_response_proved':True,
        'source_state_replaced_by_exact_Gibbs':False,
        'all_sign_origin_outcomes_included':True,
        'source_reference_joint_relative_entropy_identified':False,
        'regional_diamond_modular_kernel_identified':False,
        'area_response_or_Einstein_dynamics_proved':False,
        'response_mean_passes_to_common_GNS_limit_proved':False,
        'physical_mass_table_GeV':[]})


def bind_member(member,weyl_count=2,filter_widths=('1',)):
    result={'member_id':member['member_id'],**member['ids'],
        'mixed_word':mixed_word_bound(member,weyl_count,filter_widths),
        'pure_Weyl_word':mixed_word_bound(member,weyl_count+len(filter_widths),()),
        'complete_random_sign_history':random_sign_history(member),
        'physical_anchor':physical_anchor_residual(member),
        'source_marginal_modular_response':modular_phase_response(member)}
    return {'schema':'oph.primary_current_tail.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value,weyl_count=2,filter_widths=('1',)):
    if value!=bind_member(member,weyl_count,filter_widths):raise ValueError('source or mixed-word consequence changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*fluctuation.SOURCE_INPUTS,
                                      *scale.SOURCE_INPUTS,*modular.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_current_tail_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('energy-form inputs changed during replay')
    return result
