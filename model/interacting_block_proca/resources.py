"""Exact finite-window bounds, with the frozen native source retained."""
from fractions import Fraction as F
import math

from model.interacting_block_quadrature.model import rational_moment_bounds,ceil_fraction,ceil_sqrt_fraction
from model.interacting_block_quadrature.reference import apparatus_capacity
from .capacity import logical_capacity


def dyadic_sqrt_upper(value,bits=40):
    denominator=2**bits
    return F(ceil_sqrt_fraction(F(value)*denominator**2),denominator)


def filter_bounds(P):
    P=F(P)
    return {'momentum_radius_over_rest_frequency':str(P),
        'spatial_filter_norm_squared_upper':str(1+P*P),
        'time_filter_phi_pi_coefficient_norm_upper':str(P),
        'polarization_tensor_norm_upper':str(1+2*P*P),
        'time_filter_uses_actual_lattice_frequency_denominator':True,
        'mass_and_speed_units':'tau=Delta*t, xi=sqrt(6)*x; mass=1 in these source-derived units'}


def finite_source_resources(a,side,P=F(1,4),smearing_norm=F(1,8),kicks=2,accuracy=F(1,1000)):
    a,P,smearing_norm=map(F,(a,P,smearing_norm))
    assert 0<a<=1 and P>0 and side>=2 and kicks>=1
    D=1/a;t=(1-a*a)/(12*a);V=side**3;z=6;edges=3*V
    # Each four-vector test has this norm, so each of its spatial and
    # temporal parts is bounded by it. sqrt(1+P²)<=1+P²/2 is rational.
    f=(1+P*P/2)*smearing_norm;g=P*smearing_norm
    displacement_squared=(a*f*f+g*g/a)/2
    per_kick=dyadic_sqrt_upper(displacement_squared)
    R=kicks*per_kick/a
    amplification=1+8*R+6*R**2+4*R**3+R**4
    vacuum=rational_moment_bounds(a,side)
    base=rational_moment_bounds(a,side,True)
    M2=vacuum['second_moment'];M4=amplification**2*base['n_plus_two_fourth_moment']
    M4root=ceil_sqrt_fraction(M4)
    segments=kicks+1;tail_count=2*segments+4*kicks
    T=F(1);e_dynamics=F(accuracy);e_initial=e_dynamics/5;e_energy=e_initial**2/4
    C=D+12*t*edges
    K=max(2,ceil_fraction(8*T*C*M4root/e_dynamics),ceil_sqrt_fraction(2*M2))
    def cutoff(k):
        weight=M2/(k+1)**2
        return ((D+2*z*t)*M2/(k+1)+6*z*t*V*M2/(k+1)**2)/(1-weight)
    while cutoff(K)>e_energy/2:K*=2
    hp_target=min(e_energy/2,e_dynamics/(4*T))
    N=max(K+2,ceil_fraction(12*z*V*t*(K+1)*(K+2)/hp_target),
        ceil_sqrt_fraction(4*tail_count*M4root/e_dynamics))
    hp=12*z*V*t*(K+1)*(K+2)/N
    cut=cutoff(K)
    assert cut+hp<=e_energy and 2*(cut+hp)<=e_initial**2
    dynamical_bulk=T*hp+2*T*C*F(M4root,K+3)
    tails=tail_count*F(M4root,(N+3)**2)
    assert hp*T<=e_dynamics/4 and dynamical_bulk<=e_dynamics/2 and tails<=e_dynamics/4
    return {'spacing':str(a),'macro_side':side,'macro_sites':V,
        'source_physical_box_side':str(a*side),'unit_speed_box_side':'sqrt(6)*'+str(a*side),
        'momentum_radius_in_unit_speed_coordinates':str(P),
        'time_window_Delta_times_T':'1','Weyl_tests':kicks,'evolution_segments':segments,
        'requested_accuracy_for_records_not_measured_physics':str(e_dynamics),
        'each_four_vector_test_L2_norm_upper':str(smearing_norm),
        'filtered_spatial_coefficient_norm_upper':str(f),'filtered_momentum_coefficient_norm_upper':str(g),
        'canonical_displacement_squared_upper_per_kick':str(displacement_squared),
        'canonical_displacement_norm_upper_per_kick':str(per_kick),
        'sequence_displacement_R':str(R),'sequence_fourth_moment_upper_exact':str(M4),
        'comparison_total_number_cutoff_K':K,'actual_replicas_per_macrocell_N':N,
        'actual_native_source_cells':8*N*V,
        'actual_ground_energy_error_over_Delta_upper_exact':str(cut+hp),
        'actual_ground_vector_error_upper_exact':str(e_initial),
        'complete_instrument_dynamics_and_tails_upper_exact':str(dynamical_bulk+tails),
        'complete_source_vs_finite_free_record_TV_upper_exact':str(e_initial+dynamical_bulk+tails),
        'simple_complete_record_error_upper':str(19*e_dynamics/20),
        'full_native_complement_retained':True,
        'retained_apparatus':apparatus_capacity(kicks),
        'source_capacity_nats_expression':str(N*V)+'*(240 log(2)+8 log(5)+24 log(119))',
        'source_and_declared_readout_capacity_nats_expression':str(N*V)+'*(240 log(2)+8 log(5)+24 log(119))+log(3072)+'+str(2*kicks)+'*log(2)',
        'logical_capacity':logical_capacity(N,V,kicks),
        'readout_control_and_flag_qubits':2*kicks,
        'declared_readout_apparatus_dimension':3072*4**kicks,
        'finite_linear_filter_coefficient_count_upper':6*V*kicks,
        'admitted_matrix_dimension_expression':f'binomial({N}+3,3)^{V}',
        'exact_filter_gates_are_declared_finite_services':True,
        'gate_synthesis_and_ground_preparation_cost_included':False,
        'full_source_matrix_numerically_enumerated':False,
        'same_H_state_reference_clock_as_scalar_readouts':True}


def continuum_bounds(a,P=F(1,4),T=F(1),kicks=2,test_norm_sum=F(1,4)):
    a,P,T,S=map(F,(a,P,T,test_norm_sum))
    # q=p/sqrt6, tau=Delta*t. This is the inherited dispersion error
    # expressed in source-derived mass-one, speed-one coordinates.
    omega_error=a*a*(P*P/2+P**4/4)
    projector=1+2*P*P
    covariance=projector*(1+T)*omega_error/2
    root_histories=dyadic_sqrt_upper(F(2)**kicks)
    records=root_histories*covariance*S*S
    return {'spacing':str(a),'unit_speed_momentum_radius':str(P),
        'frequency_error_upper':str(omega_error),
        'polarization_tensor_norm_upper':str(projector),
        'covariance_error_coefficient_upper':str(covariance),
        'complete_Weyl_record_spacing_error_upper':str(records),
        'sum_of_four_vector_test_L2_norms_upper':str(S),
        'volume_error_formula':'(2P+2h)^3*(3h/2)*sup_i,p|partial_i F(p)|/(2pi)^3',
        'F':'hhat(p)^dag Pi(p) khat(p) exp(-i omega_a(p)t)/(2omega_a(p))',
        'uniform_derivative_bound':'L/2*(||dh|| ||k||+||h|| ||dk||)+(sqrt(2L)+L*(1+T)*P/2)||h|| ||k||, L=1+2P^2',
        'complete_record_covariance_transfer':'TV<=2^(r/2)*delta_cov*(sum_j||h_j||)^2 when |delta C_jk|<=delta_cov||h_j||||h_k||',
        'Fourier_tail_for_P_ge_1':'3*C_h*C_k*P^(5-2s)/[4pi^2(2s-5)], s>5/2',
        'microcausality_scope':'full limiting distribution; bounded momentum projections have controlled Fourier tails, not exact locality',
        'boosted_momentum_radius_upper':'exp(rapidity_bound)*sqrt(1+P^2)',
        'exact_finite_Lorentz_or_boost_symmetry_claim':False}
