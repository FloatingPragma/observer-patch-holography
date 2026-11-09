"""Actual source odd filtering, selector sum rules and required dressing depth."""
from fractions import Fraction as F
from math import isqrt

from model.particle_common_contract import contract, provenance
from model.particle_primary_selected_odd_response import model as odd
from model.particle_primary_current_tail import model as form
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_spatial_ground_state import ground
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
 'model/particle_primary_odd_low_energy/__init__.py',
 'model/particle_primary_odd_low_energy/__main__.py',
 'model/particle_primary_odd_low_energy/model.py',
 'model/particle_primary_odd_low_energy/test_model.py',
 'model/particle_primary_odd_low_energy/PROOF.md',
 'model/particle_primary_odd_low_energy/README.md',
 'model/particle_primary_selected_odd_response/PROOF.md',
 'model/particle_primary_current_tail/PROOF.md',
 'model/particle_filtered_probe/PROOF.md',
 'model/particle_primary_joint_spectroscopy/PROOF.md',
 'model/particle_common_contract/PROOF.md',
)


def _ceil(x):
    x=F(x);return -(-x.numerator//x.denominator)


def _sqrt_lower(x):
    x=F(x)
    return x/sqrt_upper(x) if x else F()


def _clock_steps(time,norm,allowance):
    square=time**3*norm**3/(12*allowance)
    n=isqrt(_ceil(square))
    if n*n<square:n+=1
    n=max(n,1)
    value=apparatus.cayley_service_bound(time,time,norm,n)
    if value['unitary_operator_error']>allowance:raise ArithmeticError('Cayley plan exceeds its allowance')
    return {'subdivisions':n,'operator_error_upper':value['unitary_operator_error']}


def low_pass_bound(member,half_width='1'):
    """Original zero-centered triangle success has an actual source upper bound."""
    contract.verify_member(member);T=ground._q(half_width,'positive source-time half-width',True)
    n=member['source']['level'];a=F(member['source']['spacing']);R=F(max(4,2**(n//2)))
    if not a*R<9:raise ArithmeticError('the chosen low-energy comparison cut failed')
    numerator=F(9,2) if n==0 else F(9)
    low=numerator/(9-a*R)**2;tail=min(F(1),16/(T**4*R**4))
    if T==1 and low+tail>=F(1,4):raise ArithmeticError('fixed-width low-pass upper bound failed')
    return contract.serial({'half_width':T,'comparison_energy_radius':R,
        'raw_original_weight_inside_comparison_band_upper':min(F(1),low),
        'filter_multiplier_outside_comparison_band_upper':tail,
        'ideal_zero_centered_blank_success_upper':min(F(1),low+tail),
        'uniform_fixed_R_raw_weight_formula':'mu([-R,R]) <= 9/(9-a R)^2 for a R<9',
        'fixed_T_lim_sup_blank_success_upper':F(1,9),
        'half_width_one_uniform_blank_success_upper':F(1,4),
        'positive_low_energy_success_lower_proved':False,
        'low_energy_component_excluded':False,
        'ideal_filter_is_not_the_complete_selector_marginal':True})


def dressing_depth(level,source_energy_radius='1',tail_target='1/200'):
    """A necessary depth for the specified original post-creation tail target."""
    n=ground._i(level,'source level',0)
    R=ground._q(source_energy_radius,'positive source-energy radius',True)
    eta=ground._q(tail_target,'nonnegative original tail target')
    p=F(1,50);a=F(1,2**n);L=1/(2*a)
    if eta>=p:raise ValueError('the target tail must be smaller than the proved escaping weight')
    # Lower-round the necessary drift; never round a lower bound upward.
    drift=max(F(),_sqrt_lower(L)*max(F(),_sqrt_lower(p)-sqrt_upper(eta))-sqrt_upper(R))
    count=_ceil(2*drift/45)
    if R==1 and eta==F(1,200):
        exact_case=max(F(),(_sqrt_lower(1/a)-20)/450)
        count=max(count,_ceil(exact_case))
    return contract.serial({'source_level':n,'source_spacing':a,'source_energy_radius':R,
        'original_tail_target':eta,'proved_original_escaping_weight':p,
        'proved_excitation_tail_threshold':L,
        'critical_Weyl_word_length_necessary_lower':count,
        'bound_formula':'k >= (2/45)*[sqrt(L)*(sqrt(p)-sqrt(eta))-sqrt(R)]_+',
        'simple_radius_one_target_formula':'k >= [(a^(-1/2)-20)/450]_+',
        'word_length_o_a_minus_half_cannot_produce_energy_tightness':True,
        'arbitrary_source_waits_and_adjoint_Weyls_allowed':True,
        'bound_is_for_every_native_gate_or_every_dressing':False,
        'fixed_positive_energy_residue_proved':False})


def centered_filter(member,source_bound=None):
    contract.verify_member(member)
    source_bound=odd.bind_member(member) if source_bound is None else source_bound
    odd.verify_bound(member,source_bound)
    r=source_bound['result'];a=F(member['source']['spacing']);n=member['source']['level']
    T=a/6;nu=9/a;M=1024
    grid=apparatus.quadrature(T,M)
    # The demodulated derivative is the actual centered residual, not ||H||.
    derivative=3/a;quadrature_error=derivative*T/M
    limits=member['preparation_plan']['source'];norm=max(abs(F(limits['spectral_lower'])),abs(F(limits['spectral_upper'])))
    source_clock=_clock_steps(T,norm,F(1,65536))
    scalar_clock=_clock_steps(T,nu,F(1,32768))
    implementation=2*source_clock['operator_error_upper']+scalar_clock['operator_error_upper']
    amplitude=quadrature_error+implementation
    if implementation>F(1,16384) or amplitude>F(9,16384):
        raise ArithmeticError('filter source/scalar clock budget failed')
    lo,hi=map(F,r['original_transfer_interval'])
    radius=max(abs(lo-nu),abs(hi-nu));argument=radius*T/2
    sinc_lower=1-argument**2/6
    if not 0<sinc_lower<=1 or argument>1:raise ArithmeticError('filter multiplier chart failed')
    raw=F(r['original_transfer_weight_lower']);ideal=raw*sinc_lower**4
    energy=spectroscopy.clock_plan(member,lo/4,'1/16384','1/16384')
    energy_loss=2*F(energy['clock_plan']['complete_clock_tail_upper'])+F(energy['two_complete_clock_error_upper'])
    original=ideal-2*amplitude;reported=original-energy_loss
    if reported<F(1,80):raise ArithmeticError('complete selected odd filter lost its uniform original band')
    sharp=None
    if n==0:
        narrow=r['member_zero_sharper_band'];weight=F(narrow['original_transfer_weight_lower'])
        value=weight*F(23,24)**4-2*amplitude-energy_loss
        if value<F(3,5):raise ArithmeticError('member-zero selected filtered band failed')
        eta=F(energy['energy_resolution'])
        sharp={'ideal_transfer_band':[3,15],'reported_transfer_band':[3-2*eta,15+2*eta],
            'reported_original_blank_success_weight_lower':value,'simple_reported_weight_lower':F(3,5)}
    V=member['source']['side']**3
    return contract.serial({'half_width':T,'analyzer_carrier_frequency':nu,'selector_dimension':M,
        'triangle_weight_sum':grid['total_weight'],'quadrature_step':grid['step'],
        'quadrature_node_formula':'-T+(2j+1)T/M, j=0,...,M-1',
        'quadrature_weight_formula':'(2/M)*(1-|t_j|/T)',
        'demodulated_derivative_norm_upper':derivative,'quadrature_operator_error_upper':quadrature_error,
        'original_source_norm_upper':norm,'source_echo_Cayley_plan':source_clock,
        'scalar_analyzer_Cayley_plan':scalar_clock,
        'exact_node_instrument_operator_error_upper':implementation,
        'blank_operator_error_to_continuous_triangle_upper':amplitude,
        'blank_original_event_error_upper':2*amplitude,
        'ideal_sinc_argument_radius':argument,'ideal_sinc_lower':sinc_lower,
        'ideal_source_transfer_band':[lo,hi],'ideal_original_blank_band_weight_lower':ideal,
        'finite_blank_original_band_weight_lower':original,'complete_energy_clock':energy,
        'complete_energy_clock_event_error_upper':energy_loss,
        'reported_transfer_band':[lo/2,hi+lo/2],
        'reported_original_blank_band_weight_lower':reported,'uniform_reported_weight_lower':F(1,80),
        'member_zero_sharper_reported_band':sharp,
        'complete_filter_fine_outcomes':8*V*M,
        'complete_two_clock_filter_fine_outcomes':8*V*M*energy['clock_plan']['clock_size']**2,
        'source_creation_success_and_failure_memory_retained':True,
        'all_selector_outcomes_and_quantum_reference_retained':True,
        'analyzer_carrier_subtracted_from_source_energy':False,
        'source_clock_or_energy_zero_changed':False,
        'spectral_band_is_cutoff_scale':True,'fixed_energy_continuum_survival_proved':False,
        'full_source_filter_histogram_executed':False,'division_by_success_probability':False})


def bind_member(member):
    contract.verify_member(member)
    inherited=odd.bind_member(member);odd.verify_bound(member,inherited)
    phase=form.source_form_bound(member)
    value=contract.serial({'schema':'oph.primary_odd_low_energy_result.v1','member_id':member['member_id'],**member['ids'],
        'selected_odd_receipt_sha256':inherited['result_sha256'],
        'source_uniform_phase_root_drift_upper':phase['uniform_energy_square_root_drift_upper'],
        'populated_centered_filter':centered_filter(member,inherited),
        'zero_centered_filter_bound':low_pass_bound(member),
        'necessary_phase_dressing_depth':dressing_depth(member['source']['level']),
        'complete_selector_summed_creation_energy_measure_exactly_original':True,
        'selector_sum_identity_valid_for_exact_Cayley_source_clocks':True,
        'complete_selector_sum_rule_assumes_stationary_source_marginal':True,
        'stationarity_of_joint_canonical_purification_assumed':False,
        'static_number_phase_conjugation_rule':'Z A_(x,alpha,b) Z^dagger = exp(3 i theta_x) A_(x,alpha,b)',
        'arbitrarily_growing_static_phase_conjugations_change_classical_local_instrument':False,
        'coherent_spatial_analyzers_included_in_static_invariance_claim':False,
        'all_possible_useful_success_filtered_observables_excluded':False,
        'fixed_positive_energy_odd_residue_proved':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_odd_low_energy.v1','result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('odd-filter source, controls or original probabilities changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*odd.SOURCE_INPUTS,*form.SOURCE_INPUTS,
                                     *apparatus.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_odd_low_energy_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'depth_samples':[dressing_depth(n) for n in (0,8,20,40)],
        'source_sha256':pins,'fixed_energy_survival_proved':False,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during odd-filter replay')
    return value
