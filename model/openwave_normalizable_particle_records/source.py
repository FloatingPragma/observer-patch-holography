"""Exact packet normalization, source-pole transfer and complete record bounds."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import sympy as sy
from model.openwave_baryon_poles import source as poles
from model.openwave_pole_weak_current import source as current
from model.openwave_uniform_pole_bands.source import member
from model.openwave_hadron_refinement_flow.analytic import constants
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_capacity_stress.source import enclosure
KAPPA=F(1,512);POSITION_RADIUS=8192


def moment(power,kappa=KAPPA):
    if type(power)is not int or power<0 or power%2:raise ValueError('nonnegative even moment required')
    return 15*kappa**power/F((power+1)*(power+3)*(power+5))


def cosine_mean(kappa=KAPPA,terms=8):
    partial=sum((-1)**j*moment(2*j,kappa)/factorial(2*j)for j in range(terms))
    following=(-1)**terms*moment(2*terms,kappa)/factorial(2*terms)
    return min(partial,partial+following),max(partial,partial+following)


@lru_cache(None)
def packet(proton=False):
    rest=poles.phase(proton,F(1),F(0));g=F(rest['pole']['complete_ground_gap_lower_source'])
    _,norm,_=poles.bands.repair_matrix(proton);_,_,_,p,_,unit=constants(*degree_coefficients())
    perturb=6*p['t']*KAPPA+norm*KAPPA*KAPPA
    assert 2*perturb<g
    distance=perturb/(g-perturb)
    trial_sine=sqrt_upper(1-F(rest['pole']['prepared_trial_pole_weight_lower']))+distance
    assert trial_sine<1
    weight=1-trial_sine**2
    A,G=poles.exact_symbol(F(1),F(0),proton);B,_=poles.exact_symbol(F(-1),F(0),proton)
    v=poles.vector(rest['trial_source_coordinates']);n=F(rest['trial_Gram_norm_squared'])
    coefficient=poles.real_fraction((v.conjugate().T*(B-A)*v)[0])/(2*n)
    assert coefficient>0
    coslo,coshi=cosine_mean();mean=F(rest['trial_source_energy'])
    energy=((mean+coefficient*(1-coshi))*unit,(mean+coefficient*(1-coslo))*unit)
    R=F(POSITION_RADIUS);tail=F(5)/(2*KAPPA*KAPPA*R*R);posttail=min(F(1),tail+sqrt_upper(1-weight))
    velocity=6*p['t']+2*norm
    return {'assignment':'proton'if proton else'neutron',
        'native_trial_coordinates':rest['trial_source_coordinates'],
        'full_source_band_projector_distance_to_rest_upper':str(distance),
        'uniform_band_projector_probability_lower':str(weight),
        'uniform_band_projector_probability_display':enclosure((weight,F(1))),
        'original_mixed_input_trial_then_band_probability_lower':str(weight/5),
        'uniform_trial_complement_probability_upper':str(1-weight),
        'complete_packet_trial_mean_energy_GeV_interval':list(map(str,energy)),
        'complete_packet_trial_mean_energy_GeV_display':enclosure(energy),
        'source_H_k_minus_H_zero_norm_upper_on_packet_support':str(perturb),
        'position_record':{'observable':'Leftmost occupied source cell, as a declared translation-covariant position readout.',
            'radius_cells':POSITION_RADIUS,'unfiltered_packet_outside_probability_upper':str(tail),
            'unfiltered_left_and_right_tail_probability_upper_each':str(tail/2),
            'band_conditioned_outside_probability_upper':str(posttail),
            'band_conditioned_inside_probability_lower':str(1-posttail),
            'complete_outcomes':['x < -R','-R <= x <= R','x > R'],
            'ground_band_filter_success_and_failure_both_retained':True,
            'joint_sharp_position_and_momentum_claimed':False},
        'propagation':{'source_velocity_norm_upper_cells_per_source_time':str(velocity),
            'complete_unfiltered_time_t_position_second_moment_bound':'(sqrt(5/(2*kappa^2))+abs(t)*velocity_bound)^2',
            'valid_all_real_times':True,'source_time_is_physical_seconds':False},
        'complete_position_readout_mean_work_absolute_bound_source':str(2*(6*p['t']+4*norm)),
        'complete_position_readout_mean_work_absolute_bound_GeV':str(2*(6*p['t']+4*norm)*unit)}


def run():
    k=sy.Symbol('k',real=True);kap=sy.Rational(KAPPA)
    density=sy.Rational(15,16)/kap*(1-k*k/kap**2)**2
    assert sy.integrate(density,(k,-kap,kap))==1
    assert sy.integrate(k*k*density,(k,-kap,kap))==sy.Rational(moment(2))
    A2=sy.Rational(15,16)/kap
    position_var=sy.integrate(A2*4*k*k/kap**4,(k,-kap,kap))
    assert position_var==5/(2*kap**2)
    members=[packet(p)for p in(False,True)];n,p=members
    dn,dp=F(n['full_source_band_projector_distance_to_rest_upper']),F(p['full_source_band_projector_distance_to_rest_upper'])
    rest=current.point(F(1),F(0))['currents'];rows={}
    for name,norm in(('vector',1),('axial',3)):
        lo,hi=map(F,rest[name]['actual_source_pole_amplitude_magnitude_interval'])
        error=F(norm)*sqrt_upper(2)*(dn+dp)
        low=max(F(0),lo-error);high=min(F(norm),hi+error);assert low>0
        rows[name]={'uniform_source_pole_amplitude_magnitude_interval':[str(low),str(high)],
            'uniform_source_pole_amplitude_display':enclosure((low,high)),
            'conditional_band_input_current_and_proton_band_probability_interval':[str(low*low/9),str(high*high/9)],
            'original_mixed_input_trial_neutron_band_current_proton_band_probability_lower':str(F(n['uniform_band_projector_probability_lower'])*low*low/45)}
    return {'schema':'oph.normalizable_native_particle_records.v1',
        'verdict':'PASS_NORMALIZABLE_INFINITE_CHAIN_PARTICLE_PACKET_AND_COMPLETE_RECORDS',
        'normalizable_packet':{'momentum_half_width_radians_exact':str(KAPPA),
            'momentum_amplitude':'sqrt(15/(16*kappa))*(1-(k/kappa)^2) for |k|<=kappa; zero otherwise, with Lebesgue dk normalization.',
            'exact_L2_norm_squared':'1','even_moment_formula':'15*kappa^(2m)/[(2m+1)(2m+3)(2m+5)]',
            'momentum_variance_exact':str(moment(2)),
            'position_mean_cells_exact':'0','position_variance_cells_squared_exact':str(F(position_var)),
            'position_wavefunction':'For x!=0: sqrt(15/(16*kappa))/sqrt(2*pi) * 4*(sin(kappa*x)-kappa*x*cos(kappa*x))/(kappa^2*x^3); at x=0 use its continuous limit.',
            'position_probability_at_origin_exact':'5*kappa/(6*pi)',
            'norm_argument':'Unitary Fourier transform L2([-pi,pi],dk) to l2(Z). The compactly supported continuous piecewise-polynomial amplitude lies in H1; Parseval gives finite position variance.',
            'carrier_is_infinite_chain':True,'finite_ring_or_delta_momentum_substitution':False},
        'members':members,'uniform_native_weak_packet_records':rows,
        'complete_protocol':{'input':'The normalized spatial packet tensor the native five-coordinate maximally mixed input P/5.',
            'trial_filter':'K_s=|v><v|, K_f=I-K_s in internal source coordinates, with actual inherited rational v.',
            'trial_success_exact':'1/5','trial_failure_exact':'4/5',
            'nonselective_trial_preparation_work_exact':'0 on this specified input',
            'band_filter':'The direct integral of the isolated full-source rank-one ground projections and its full orthogonal complement.',
            'band_filter_is_bounded_and_preserves_normalizability':True,
            'band_filter_nonselective_source_work_exact':'0; the complete band projector commutes with the source H.',
            'native_current_instrument':current.run()['complete_native_instrument'],
            'position_readout':'Complete position effects are executed as terminal alternatives, not a simultaneous sharp position/momentum record.',
            'position_readout_work':'Tr[H(sum_x P_x rho P_x-rho)] is retained; it is not asserted zero.',
            'arbitrary_idle_external_quantum_reference':True,
            'finite_local_compilation_of_exact_band_filter_proved':False},
        'finite_theorem_reuse':'Uniform all-momentum source isolation, exact finite preparation/residue bound, spectral-projector perturbation, native CAR current contraction and Fourier Parseval.',
        'new_empirical_parameter_or_mass_fit':False,
        'ancestry':'Electron anchor, strong/axial mass calibration, supplied speed/length; kappa and the observer position radius are declared packet/readout choices.',
        'physical_four_dimensional_particle_identification':False,'autonomous_beta_decay_or_lifetime_proved':False,
        'scope':'Actual normalizable states and complete bounded instruments on the infinite one-dimensional source. The exact band filter is a declared nonlocal spectral service; its finite local implementation and physical 3+1 transfer remain open.'}
