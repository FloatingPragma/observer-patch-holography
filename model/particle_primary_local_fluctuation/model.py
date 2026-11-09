"""Exact finite local twists of the selected physical source.

Uniform energy tightness includes the elastic atom. A positive first
moment alone does not prove a nonzero limiting excitation probability.
"""
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_energy_tightness import model as energy
from model.particle_primary_origin_momentum import model as origin
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
    'model/particle_primary_local_fluctuation/__init__.py',
    'model/particle_primary_local_fluctuation/__main__.py',
    'model/particle_primary_local_fluctuation/model.py',
    'model/particle_primary_local_fluctuation/test_model.py',
    'model/particle_primary_local_fluctuation/PROOF.md',
    'model/particle_primary_energy_tightness/PROOF.md',
    'model/particle_primary_origin_momentum/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def tent_sums(radius):
    R=ground._i(radius,'tent radius',1)
    return {'radius':R,'one_dimensional_sum':F(R),
        'one_dimensional_square_sum':F(2*R*R+1,3*R),
        'one_dimensional_gradient_square_sum':F(2,R),
        'one_axis_three_dimensional_gradient_square_sum':F(2*(2*R*R+1)**2,9*R**3),
        'three_axis_gradient_absolute_sum':6*R*R,
        'nonzero_profile_sites':(2*R-1)**3}


def positive_window(mean_lower,second_moment_upper):
    """Original window from a signed mean and a full second moment."""
    mean=ground._q(mean_lower,'positive mean',True)
    moment=ground._q(second_moment_upper,'second moment',True)
    if moment<mean*mean:raise ValueError('second moment is smaller than squared mean')
    return {'original_transfer_interval':[mean/4,4*moment/mean],
        'original_weight_lower':mean*mean/(4*moment)}


def phase_implementation_bound(radius,per_site_radian_error):
    R=ground._i(radius,'tent radius',1)
    delta=ground._q(per_site_radian_error,'phase error')
    # Every centered number has norm30. This is the original complete
    # random-origin/sign instrument, including an arbitrary reference.
    return min(F(1),30*(2*R-1)**3*delta)


def selected_bound(member):
    contract.verify_member(member)
    source=member['source'];n=source['level'];L=source['side'];V=L**3
    a=F(source['spacing']);R=2**n;eps=F(member['state']['target_ground_excess'])
    if a*R!=1 or L!=2*R*R:raise ArithmeticError('primary spatial scaling changed')
    inherited=energy.selected_density_sum_rule(member)
    h=F(inherited['selected_negative_dimensionless_hopping_density_lower'])
    hstar=F(inherited['uniform_negative_dimensionless_hopping_density_lower'])
    sums=tent_sums(R);axis=sums['one_axis_three_dimensional_gradient_square_sum']
    k_lower=2*a*a*axis;k_upper=5*a*a*axis
    if n==0:k_lower=k_upper=F(4)
    mean_lower=k_lower*h/a
    mean_upper=k_upper*30*source['interior_hopping_bonds']/(a*V)
    uniform_lower=F(16,9)*hstar;uniform_upper=F(14175,32)
    local_upper=F(180) if n==0 else 90*k_upper/a
    if not 0<uniform_lower<=mean_lower<=mean_upper<=uniform_upper:
        raise ArithmeticError('same-state Weyl work bounds failed')
    if local_upper>F(2025,4):raise ArithmeticError('uniform local work cap failed')
    # On the eight-cell source this is the actual full-Fock corner star.
    # 16*15²*(3+sqrt6)<19620. Else use the exact tent L1-gradient
    # and pi²<10, without an electric-cutoff norm.
    second=F(19620) if n==0 else F(324000*R**4)
    band=positive_window(mean_lower,second)
    weight=band['original_weight_lower'];lo,hi=band['original_transfer_interval']
    clock=spectroscopy.clock_plan(member,lo/4,weight/16,weight/16)
    clock_loss=weight/8+F(clock['two_complete_clock_error_upper'])
    if clock_loss>weight/4:raise ArithmeticError('complete Weyl clock budget failed')
    geometry=origin.origin_geometry(n);fraction=geometry['retained_fraction']
    discard=geometry['discarded_fraction']
    retained_lower=max(F(),(mean_lower-discard*local_upper)/fraction)
    retained_upper=min(local_upper,mean_upper/fraction)
    return contract.serial({'schema':'oph.primary_local_fluctuation_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'hopping_density_producer_digest':contract.digest(contract.serial(inherited)),
        'probe':'exp(sign*i*pi*a*sum_y f_R(y-X)*(N_y-30I))',
        'profile':'product of three periodic discrete tents; only apparatus coefficients wrap',
        'profile_root_of_unity_order':2*R**4,'tent_sums':sums,
        'source_H_or_boundary_conditions_changed':False,
        'source_input_state_and_idle_reference_unchanged':True,
        'after_intervention_state_unchanged':False,
        'complete_origin_sign_outcomes':2*V,
        'same_origin_and_sign_retained_through_history':True,
        'unitary_twist_has_no_discarded_success_or_failure_sector':True,
        'original_total_spectral_weight':1,
        'cosine_gradient_sum_enclosure':[k_lower,k_upper],
        'original_signed_first_moment_lower':mean_lower,
        'original_signed_first_moment_upper':mean_upper,
        'uniform_original_signed_first_moment_lower':uniform_lower,
        'uniform_original_signed_first_moment_upper':uniform_upper,
        'original_negative_absolute_first_moment_upper':eps,
        'original_absolute_first_moment_upper':mean_upper+2*eps,
        'uniform_original_absolute_first_moment_upper':uniform_upper+F(1,8),
        'energy_tail_formula':'mu(|omega|>E) <= min(1, absolute_first_moment_upper/E)',
        'uniform_original_probability_tightness_proved':True,
        'subsequential_continuous_positive_energy_probability_measure':True,
        'limit_may_be_purely_elastic':True,
        'uniform_fixed_positive_band_noncollapse_proved':False,
        'retained_origin':{'geometry':geometry,
            'profile_unwrapped_on_retained_origins':n>=1,
            'original_record_distance_after_common_label_embedding_upper':discard,
            'signed_first_moment_lower':retained_lower,
            'signed_first_moment_upper':retained_upper,
            'uniform_absolute_first_moment_upper':F(2025,4)+F(1,8),
            'selected_source_density_replaced_by_averaged_state':False,
            'regular_momentum_dictionary_transferred_to_this_probe':False},
        'second_moment_upper':second,'finite_original_window':band,
        'complete_two_energy_clock':clock,
        'complete_fine_record_count':2*V*clock['clock_plan']['clock_size']**2,
        'finite_reported_transfer_interval':[lo/2,hi+lo/2],
        'finite_reported_original_weight_lower':weight-clock_loss,
        'initial_energy_record_preserves_source_marginal':True,
        'initial_energy_record_preserves_joint_purification':False,
        'full_primary_histogram_executed':False,
        'spectral_line_or_rest_momentum_or_particle_identification':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    value=selected_bound(member)
    return {'schema':'oph.primary_local_fluctuation_bound.v1','result':value,
        'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('Weyl response or its source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*energy.SOURCE_INPUTS,*origin.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_local_fluctuation_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2)],
        'retained_origin_member':bind_member(contract.member(32)),
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('Weyl inputs changed during replay')
    return result
