"""Actual dressed bond seeds, complete finite filters and canonical H/P kernels."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_candidate_charged import model as charged
from model.particle_primary_momentum_transfer import model as momentum
from model.particle_primary_origin_momentum import model as origin
from model.particle_primary_joint_spectroscopy import spectroscopy

SOURCE_INPUTS=(
 'model/particle_candidate_charged_propagation/__init__.py',
 'model/particle_candidate_charged_propagation/__main__.py',
 'model/particle_candidate_charged_propagation/model.py',
 'model/particle_candidate_charged_propagation/test_model.py',
 'model/particle_candidate_charged_propagation/PROOF.md',
 'model/particle_candidate_charged_propagation/README.md',
 'model/particle_candidate_charged/results.json',
 'model/particle_primary_field_algebra/PROOF.md',
 'model/particle_primary_spatial_algebra/PROOF.md',
 'model/particle_primary_origin_momentum/PROOF.md',
 'model/particle_primary_momentum_transfer/PROOF.md',
 'Lean/QFT/GroundRecordTail.lean','Lean/QFT/PeriodicMomentum.lean',
 *charged.SOURCE_INPUTS,*momentum.SOURCE_INPUTS,*origin.SOURCE_INPUTS,
 *spectroscopy.SOURCE_INPUTS,
)


@lru_cache(maxsize=1)
def seed_control():
    comp=charged.compensated_control();d=charged.dressing_control()
    C=tuple(charged.regional.creator(j,8) for j in range(8))
    numbers=tuple(x*x.H for x in C)
    NL=charged.clean(s.kronecker_product(s.eye(4),d['NL']))
    NE=charged.clean(s.kronecker_product(s.eye(4),d['NE'])+numbers[6]+numbers[7])
    parity=s.diag(*((-1)**i.bit_count() for i in range(256)))
    for B in comp['rows']:
        if charged.clean(NL*B-B*NL)!=s.zeros(256) or charged.clean(NE*B-B*NE)!=s.zeros(256):
            raise ArithmeticError('the neutral bond seed changes a global native occupation')
        if parity*B!=B*parity:raise ArithmeticError('the bond seed is not CAR even')
    return {'source_seed':'O_(x,j,p,s)=chi_6(U_(x,x+e_j)) P0 c_(x,p)^dagger D_(x+e_j,s)/3',
        'forward_axes':3,'spin_rows_per_axis':4,'include_adjoint_and_reversed_bond':True,
        'owner_support_cells':2,'extra_physical_support_radius':'a',
        'full_Fock_seed_norm_upper':1,'summed_four_seed_effect_norm':1,
        'inside_Q6_increment':-6,'retained_compensator_Q6_increment':6,
        'global_native_increment':[0,0,0,0,0],'global_CAR_parity':0,
        'actual_native_spin_factors':'G1 at each endpoint; complete pair coefficients A1 plus T1',
        'all_Gauss_constraints_preserved':True,'first_constant_link_action_exact':True,
        'global_neutral_pair_is_a_charged_odd_field':False,
        'zero_intertwiner_requires_preserved_native_charge_and_CAR_parity_dictionary':True,
        'unbroken_charge_parity_intertwiner_to_single_charged_fermion':'zero',
        'boundary_superselection_or_extended_field_algebra_excluded':False,
        'seed_zero_outside_actual_source_bonds':True,
        'bulk_aperture_collar_is_a_separate_admission_test':True,
        'native_C6_collective_seed_identified_with_dressed_seed':False}


def filter_compiler(member,width='1/128',allowance='1/65536'):
    contract.verify_member(member);T=F(width);e=F(allowance)
    if T<=0 or e<=0:raise ValueError('positive rational filter width and allowance required')
    lo=F(member['state']['spectral_lower']);hi=F(member['state']['spectral_upper'])
    D=hi-lo;K=max(abs(lo),abs(hi))
    M=2*max(1,(D*T/(2*e)).__ceil__())
    steps=max(1,momentum.ceil_sqrt(K**3*T**3/(12*e)))
    quadrature=D*T/M;cayley=K**3*T**3/(12*steps**2)
    return {'triangle_width':T,'inner_time_cells':M,'Cayley_subdivisions':steps,
        'source_spectral_width':D,'source_norm_upper':K,
        'continuous_success_column_quadrature_error_upper':quadrature,
        'one_source_evolution_operator_error_upper':cayley,
        'complete_finite_node_instrument_half_trace_error_upper':2*cayley,
        'distinguished_success_column_amplitude_error_upper':quadrature+2*cayley,
        'distinguished_original_event_probability_error_upper':2*quadrature+4*cayley,
        'energy_zero_changed':False,'source_clock_changed':False,
        'time_nodes':'even midpoint grid on [-T,T]; exact positive triangular cell weights',
        'analyzer':'exact weighted Householder with distinguished row sqrt(weights)',
        'local_complete_instrument_outcomes':97,'fine_outcomes':97*M,
        'all_temporal_analyzer_and_probe_failure_records_retained':True,
        'comparison_failure_channels_assumed_equal':False,
        'large_source_or_selector_executed':False}


def triangle_band_retention(width,upper):
    T=F(width);U=F(upper)
    if T<=0 or U<0 or T*U/2>1:raise ValueError('the rational lower uses |T omega/2|<=1')
    return (1-(T*U)**2/24)**4


def ultraviolet_filter_upper(spacing,width='1/128',lower='39/4'):
    a=F(spacing);T=F(width);l=F(lower)
    if min(a,T,l)<=0:raise ValueError('positive spacing, filter width and band lower required')
    return min(F(1),16*a**4/(T*l)**4)


def filtered_original_window(member,parent=None):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the positive dressed spectral window is member zero only')
    parent=charged.bind_member(member) if parent is None else parent
    charged.verify_bound(member,parent)
    original=parent['result']['original_compensated_dressed_energy_window']
    bare=charged.regional.charged.bind_member(member)
    charged.regional.charged.verify_bound(member,bare)
    p=F(bare['result']['original_transfer_weight_lower'])
    eps=F(original['complete_history_residual_squared_upper'])
    if p<=F(37,120) or eps>=F(129,400)**2 or F(523,1000)**2>=F(8,9)*F(37,120):
        raise ArithmeticError('the actual ideal spectral comparison lost its rational margin')
    ideal=(F(523,1000)-F(129,400))**2
    if ideal!=F(160801,4000000) or ideal<=F(1,25):raise ArithmeticError('the ideal dressed lower changed')
    band=tuple(F(x) for x in bare['result']['original_transfer_band'])
    if band!=(F(39,4),F(135,4)):raise ArithmeticError('the authenticated ideal charged band differs')
    T=F(1,128);retention=triangle_band_retention(T,band[1]);filtered=ideal*retention
    if filtered<=F(39,1000):raise ArithmeticError('the fixed-time filtered band lost its weight')
    compiler=filter_compiler(member,T)
    clocks=spectroscopy.clock_plan(member,'1/64','1/32768','1/32768')
    loss=2*F(clocks['clock_plan']['complete_clock_tail_upper'])+F(clocks['two_complete_clock_error_upper'])
    recorded=filtered-compiler['distinguished_original_event_probability_error_upper']-loss
    if recorded<=F(3,80):raise ArithmeticError('the complete filtered energy records lost their margin')
    return {'dressed_source_result_sha256':parent['result_sha256'],
        'bare_ideal_source_result_sha256':bare['result_sha256'],
        'ideal_transfer_interval':band,'bare_original_ideal_weight_lower':p,
        'original_residual_squared_upper':eps,
        'original_dressed_ideal_weight_strict_lower':ideal,
        'original_dressed_ideal_simple_strict_lower':F(1,25),
        'fixed_triangle_width':T,'spectral_amplitude_multiplier':'sinc^2(T omega/2)',
        'spectral_probability_multiplier':'sinc^4(T omega/2)',
        'triangle_band_retention_lower':retention,
        'original_time_filtered_ideal_weight_strict_lower':filtered,
        'original_time_filtered_ideal_simple_strict_lower':F(39,1000),
        'filter_compiler':compiler,'complete_source_clock_plan':clocks,
        'complete_two_clock_event_loss_upper':loss,
        'reported_transfer_interval':[band[0]-F(1,32),band[1]+F(1,32)],
        'original_filtered_recorded_weight_strict_lower':recorded,
        'original_filtered_recorded_simple_strict_lower':F(3,80),
        'full_clock_probe_records':compiler['fine_outcomes']*int(clocks['clock_plan']['clock_size'])**2,
        'all_other_temporal_analyzer_rows_and_reference_retained':True,
        'conditioned_on_temporal_success_or_probe_success':False,
        'time_filter_is_spatial_regularization':False,
        'positive_fixed_capacity_record_is_a_surviving_continuum_pole':False,
        'pair_energy_divided_into_endpoint_masses':False}


def joint_plan(member,spatial_width='1',time_width='1/128',spatial_aperture='1',
               time_aperture='1/1000',spatial_index=(1,0,0),time_index=0):
    contract.verify_member(member);n=member['source']['level'];a=F(member['source']['spacing'])
    R,T,Q,S=map(F,(spatial_width,time_width,spatial_aperture,time_aperture))
    if min(R,T,Q,S)<=0 or (R/a).denominator!=1 or (Q/a).denominator!=1:
        raise ValueError('positive grid-aligned spatial widths and positive times required')
    k=momentum._frequency(spatial_index)
    if type(time_index) is not int:raise ValueError('integer temporal Fourier index required')
    geo=origin.origin_geometry(n)
    # A forward bond consumes one extra cell; two owner cells enter the actual Duhamel cap.
    bulk=origin.bulk_bound(member,R+Q+a,T+S,2);active=bulk['all_origins_bulk_admitted']
    lo=F(member['state']['spectral_lower']);hi=F(member['state']['spectral_upper'])
    width=hi-lo;norm=max(abs(lo),abs(hi));target=F(1,2**(n+12))
    M=2*max(1,(width*T/(2*target)).__ceil__())
    steps=max(1,momentum.ceil_sqrt(norm**3*(T+S)**3/(12*target)))
    eU=norm**3*(T+S)**3/(12*steps**2);eq=width*T/M
    Nt=max(2,2**n);h=2*S/Nt
    grid=a*(6/R+F(22,7)*sum(abs(x) for x in k)/Q)+h*(2/T+F(22,7)*abs(time_index)/S)
    bits=bulk['operator_error']['negative_exponent']
    eb=F() if bulk['operator_error']['exact_zero'] else (F(2) if bits<0 else F(1,2**min(bits,n+12)))
    amp=min(F(2),eq+2*eU+grid+eb) if active else None
    inner=int(2*R/a)**3*M;outer=int(2*Q/a)**3*Nt
    return {'parameters':{'spatial_width':R,'time_width':T,'spatial_aperture':Q,
            'time_aperture':S,'spatial_index':k,'time_index':time_index},
        'seed':'actual Gauss-even scalar-dressed charged bond; each chosen spin/axis row has norm at most one',
        'seed_support_owners':2,'seed_extra_radius':a,'origin_geometry':geo,
        'active_source_aperture':active,'inactive_convention':'zero success; complete failure retained',
        'complete_analyzer':'positive inner triangle weights then Cartesian outer space/time Fourier rows',
        'inner_analyzer_dimension':inner,'outer_analyzer_dimension':outer,
        'selector_dimension':inner*outer,'fine_outcomes_per_origin':2*inner*outer,
        'complete_fine_outcomes_with_retained_origin':2*inner*outer*geo['origin_count'],
        'outer_Fourier_shape':[int(2*Q/a)]*3+[Nt],
        'distinguished_Fourier_bin':[x%int(2*Q/a) for x in k]+[time_index%Nt],
        'source_spectral_interval':[lo,hi],'source_norm_upper':norm,
        'inner_triangle_cells':M,'Cayley_subdivisions':steps,
        'inner_quadrature_operator_error_upper':eq,'one_source_evolution_operator_error_upper':eU,
        'complete_node_instrument_half_trace_error_upper':min(F(1),2*eU),
        'regular_outer_grid_operator_error_upper':min(F(2),grid),
        'same_source_bulk_bound':bulk,'bulk_error_usable_upper':eb,
        'success_amplitude_comparison_upper':amp,
        'original_probability_comparison_upper':min(F(1),2*amp) if active else None,
        'spatial_regular_Lipschitz_per_coordinate':2/R,'time_regular_Lipschitz':2/T,
        'canonical_joint_kernel':'sinc^2(S(E-nu)) product_i sinc^2(Q(P_i-q_i))',
        'frequency_centers':{'q_over_pi':[F(x)/Q for x in k],'nu_over_pi':F(time_index)/S},
        'one_common_origin_retained_across_history':True,'quantum_reference_retained':True,
        'same_source_clock_energy_zero_and_preparation':True,
        'common_countable_algebra_extension_and_subsequence':True,
        'old_GNS_cyclic_subspace_embeds_and_H_P_intertwine':True,
        'whole_extended_GNS_Hilbert_space_equals_old_space':False,
        'canonical_energy_nonnegative':True,'canonical_H_P_strongly_commute':True,
        'finite_open_source_translation_invariant':False,
        'finite_Fourier_bin_is_a_conserved_momentum_eigenvalue':False,
        'fixed_origin_state_replaced_by_translation_invariant_state':False,
        'limiting_joint_spectral_weight_lower_proved':False,
        'state_subsequence_convergence_modulus_supplied':False,
        'large_source_or_selector_executed':False}


def bind_member(member):
    contract.verify_member(member);seed=seed_control()
    parent=charged.bind_member(member);charged.verify_bound(member,parent)
    n=member['source']['level'];a=F(member['source']['spacing'])
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'charged_parent_result_sha256':parent['result_sha256'],'actual_seed_control':seed,
        'filtered_original_window':filtered_original_window(member,parent) if n==0 else None,
        'charged_joint_momentum_plan':joint_plan(member),
        'fixed_time_triangle_UV_component_upper':{
            'width':F(1,128),'scaled_positive_band_lower':F(39,4),
            'probability_multiplier_upper':ultraviolet_filter_upper(a),
            'formula':'min(1,16 a^4/(T l)^4) on omega>=l/a',
            'UV_band_population_at_all_members_assumed':False,
            'whole_charged_field_proved_zero':False},
        'physical_particle_map':[],'physical_mass_table_GeV':[],
        'next_needed':'positive original finite-energy and finite-momentum weight after spatial regularization; a charged-sector or asymptotic endpoint dictionary distinct from the neutral pair'})
    return {'schema':'oph.particle_candidate_charged_propagation_bound.v1',
        'result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,bound):
    if bound!=bind_member(member):raise ValueError('charged propagation source or result differs')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.particle_candidate_charged_propagation.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('charged propagation source changed during replay')
    return value
