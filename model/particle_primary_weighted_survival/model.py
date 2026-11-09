"""Actual-source obstruction to subcritical normalized-density histories."""
from fractions import Fraction as F
from math import isqrt

from model.particle_common_contract import contract, provenance
from model.particle_primary_current_tail import model as forms
from model.particle_primary_local_fluctuation import model as fluctuation
from model.particle_spatial_ground_state import ground
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_weighted_survival/__init__.py',
    'model/particle_primary_weighted_survival/__main__.py',
    'model/particle_primary_weighted_survival/model.py',
    'model/particle_primary_weighted_survival/test_model.py',
    'model/particle_primary_weighted_survival/PROOF.md',
    'model/particle_primary_weighted_survival/README.md',
    'model/particle_primary_current_tail/PROOF.md',
    'model/particle_primary_energy_tightness/PROOF.md',
    'model/particle_primary_local_fluctuation/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def normalization(member):
    contract.verify_member(member)
    n=member['source']['level'];a=F(member['source']['spacing']);R=2**n
    sums=fluctuation.tent_sums(R)
    D=sums['one_axis_three_dimensional_gradient_square_sum']
    C=F(1,10) if n==0 else D/(10*a*R**6)
    if C>F(9,80)*a**4:raise ArithmeticError('normalized source gradient changed')
    # A relative rational bound, so no fixed absolute rounding floor masks
    # the vanishing source cooling energy at arbitrary refinement level.
    root_eps=F(1,2**((n+1)//2+2))*(F(3,2) if n%2 else 1)
    eps=F(member['state']['target_ground_excess'])
    if not eps<=root_eps**2<=F(9,8)*eps:raise ArithmeticError('cooling root bound changed')
    return {'normalized_generator':'Q_X=sum_y f_R(y-X)*(N_y-30)/(30*R^3)',
        'operator_norm_upper':1,'radius':R,'profile_sum':R**3,
        'double_commutator_norm_upper':C,
        'uniform_double_commutator_coefficient':F(9,80),
        'uniform_double_commutator_formula':'(9/80)*a^4',
        'phase_gate':'exp(i*pi*r*Q_X tensor A), ||A||<=1',
        'unit_strength_energy_root_drift_upper':F(3,4)*a*a,
        'selected_mean_ground_excess_upper':eps,
        'selected_energy_root_upper':root_eps,
        'exact_gate_control_spectrum':'rational eigenvalues in [-1,1]',
        'Gauss_and_word_cutoff_preserved':True}


def fixed_history(member,strengths=('1','1','1'),threshold='1',waiting_changes=('1',)):
    data=normalization(member)
    turns=tuple(ground._signed(x,'phase strength') for x in strengths)
    waits=tuple(ground._signed(x,'source waiting-time change') for x in waiting_changes)
    if not turns:raise ValueError('at least one phase call is required')
    eta=ground._q(threshold,'positive source-energy threshold',True)
    total=sum(map(abs,turns),F())
    root=data['selected_energy_root_upper']+data['unit_strength_energy_root_drift_upper']*total
    energy=root*root
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'normalization':data,'phase_strengths':turns,'total_absolute_phase_strength':total,
        'positive_energy_threshold':eta,'waiting_time_changes':waits,
        'final_ground_relative_energy_mean_upper':energy,
        'original_positive_transfer_probability_upper':min(F(1),energy/eta),
        'final_positive_excitation_probability_upper':min(F(1),energy/eta),
        'complete_waiting_comparison_half_trace_upper':min(F(1),
            root*sum((sqrt_upper(2*abs(t)) for t in waits),F())),
        'waiting_comparison_domain':'uncontrolled source waiting times; every subsequent output retained',
        'controlled_clock_ground_phase_removed':False,
        'controlled_waiting_time_indistinguishability_claimed':False,
        'complete_reference_and_pointer_outputs_retained':True,
        'coherent_control_and_adaptive_pointer_operations_allowed':True,
        'controlled_source_evolution_calls_have_zero_energy_drift':True,
        'filter_success_probability_not_renormalized':True,
        'full_spatial_histogram_executed':False,
        'subcritical_strength_condition':'a^2*S_n -> 0',
        'subcritical_complete_positive_energy_weight_vanishes':True,
        'subcritical_fixed_waiting_time_distinguishability_vanishes':True,
        'all_pointer_records_become_deterministic':False,
        'critical_strength_Weyl_noncollapse_or_collapse_proved':False,
        'physical_mass_table_GeV':[]})


def minimum_strength(member,probability='1/100',threshold='1'):
    data=normalization(member)
    p=ground._q(probability,'required original probability',True)
    eta=ground._q(threshold,'positive threshold',True)
    if p>1:raise ValueError('original probability cannot exceed one')
    x=p*eta;scale=1<<80
    low=F(isqrt(x.numerator*scale*scale//x.denominator),scale)
    if not low*low<=x<(low+F(1,scale))**2:raise ArithmeticError('rational root enclosure failed')
    lower=max(F(),low-data['selected_energy_root_upper'])/data['unit_strength_energy_root_drift_upper']
    return contract.serial({'required_original_probability':p,'positive_threshold':eta,
        'total_phase_strength_necessary_lower':lower,
        'unit_bounded_phase_calls_necessary_lower':-(-lower.numerator//lower.denominator),
        'formula':'S >= (4/(3*a^2))*(sqrt(p*eta)-sqrt(epsilon))_+',
        'gate_count_domain':'each normalized Q phase call has abs(r)<=1',
        'native_hardware_gate_lower_bound_claimed':False})


def finite_work(member,strength='1'):
    data=normalization(member);r=ground._signed(strength,'phase strength')
    a=F(member['source']['spacing']);R=data['radius']
    if abs(r)>30*R**4:raise ValueError('the sine-chord work domain is exceeded')
    source=member['source'];density=fluctuation.energy.selected_density_sum_rule(member)
    D=fluctuation.tent_sums(R)['one_axis_three_dimensional_gradient_square_sum']
    k_low=2*r*r*D/(900*R**6);k_high=5*r*r*D/(900*R**6)
    low=k_low*F(density['selected_negative_dimensionless_hopping_density_lower'])/a
    high=k_high*30*source['interior_hopping_bonds']/(a*source['side']**3)
    cap=5*r*r*data['double_commutator_norm_upper']
    return contract.serial({'complete_instrument':'uniform origin and retained +/- phases exp(+/-i*pi*r*Q_X)',
        'phase_strength':r,'original_signed_energy_mean_interval':[low,min(high,cap)],
        'uniform_work_upper':F(9,16)*r*r*a**4,
        'work_upper_formula':'(9/16)*r^2*a^4',
        'finite_nonzero_work_for_nonzero_strength':r!=0,
        'uniform_positive_work_at_fixed_strength_proved':False})


def critical_weyl(member):
    data=normalization(member);a=F(member['source']['spacing']);r=30/(a*a)
    if data['unit_strength_energy_root_drift_upper']*r!=F(45,2):
        raise ArithmeticError('the critical Weyl normalization changed')
    forms.source_form_bound(member)
    return contract.serial({'original_Weyl_strength':r,'identity':'pi*r*Q_X=pi*a*sum_y f_R(y-X)*(N_y-30)',
        'same_original_source_Weyl_unitary':True,'energy_root_drift_upper':F(45,2),
        'subcritical_obstruction_applies':False,
        'selected_energy_weighted_tail_bound_supplied':False,
        'fixed_positive_energy_survival_decided':False})


def bind_member(member):
    result={'member_id':member['member_id'],**member['ids'],
        'fixed_history':fixed_history(member),'necessary_strength':minimum_strength(member),
        'balanced_finite_work':finite_work(member),'critical_Weyl':critical_weyl(member)}
    return {'schema':'oph.primary_weighted_survival.v1','result':result,
        'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('source or phase-strength consequence changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*forms.SOURCE_INPUTS,*fluctuation.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_weighted_survival_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('weighted survival inputs changed during replay')
    return result
