"""Constructive small-hopping graph and full-H state/record error bounds."""
import sympy as s
from .lattice import gap_data,hopping_norm,local_norm


def parameters(side=2,kappa=s.Integer(1),coupling=None):
    atomic=gap_data(side,kappa);gap=atomic['gap'];tau=hopping_norm(side)
    maximum=s.simplify(gap/(5*tau))
    coupling=maximum if coupling is None else s.sympify(coupling)
    if coupling.has(s.Float) or coupling.is_real is not True or coupling.is_finite is not True:
        raise ValueError('an exact finite real hopping coefficient is required')
    v=s.simplify(abs(coupling)*tau)
    if v>gap/5:raise ValueError('hopping coefficient exceeds the certified finite-volume interval')
    separated=s.simplify(gap-v);radius=s.simplify(2*v/separated)
    contraction=s.simplify(4*v**2/separated**2)
    centered=s.simplify(atomic['volume']*(local_norm()+3*atomic['kappa']+12)+v)
    assert 0<=radius<=s.Rational(1,2) and 0<=contraction<=s.Rational(1,4)
    return {'side':side,'kappa':atomic['kappa'],'volume':atomic['volume'],
            'gap':gap,'tau':tau,'lambda_max':maximum,'lambda':coupling,'v':v,
            'separated':separated,'radius':radius,'contraction':contraction,'centered_norm':centered}


def iteration_bounds(iterations,side=2,kappa=s.Integer(1),coupling=None):
    if type(iterations) is not int or iterations<0:raise ValueError('nonnegative exact iteration count required')
    p=parameters(side,kappa,coupling);r=p['radius'];q=p['contraction'];h=p['centered_norm']
    graph=s.simplify(r*q**iterations)
    isometry=s.simplify((1+2*r)*graph)
    residual=s.simplify(4*h*isometry)
    return {'iterations':iterations,'graph_operator_error_upper_bound':str(graph),
            'isometry_operator_error_upper_bound':str(isometry),
            'full_H_residual_operator_upper_bound':str(residual),
            'state_or_complete_record_error_at_time_t':'min(1, |t| * full_H_residual_operator_upper_bound)',
            'stationary_band_density_trace_distance_upper_bound':str(min(s.Integer(1),isometry))}


def certificate(side=2,kappa=s.Integer(1)):
    p=parameters(side,kappa);assert p['v']==p['gap']/5
    assert p['radius']==s.Rational(1,2) and p['contraction']==s.Rational(1,4)
    return {'comparison_Hamiltonian':'H_lambda=H0+lambda*T on the same full PW carrier',
            'side':side,'volume':p['volume'],'link_penalty':str(p['kappa']),
            'complete_hopping_norm_upper_bound':str(p['tau']),
            'admitted_hopping_coefficient_interval':[str(-p['lambda_max']),str(p['lambda_max'])],
            'actual_unit_hopping_coefficient_is_in_interval':bool(1<=p['lambda_max']),
            'bare_atomic_full_H_residual_bound':'|lambda|*tau',
            'bare_atomic_state_and_complete_record_error':'min(1, |t|*|lambda|*tau)',
            'low_band_rank':2**p['volume'],
            'band_energy_window_relative_to_volume_E0':['-v','v'],
            'band_to_rest_gap_lower_bound_at_admitted_endpoint':str(p['gap']-2*p['v']),
            'all_Q_states_in_graph_iteration':True,
            'graph_domain':'complete physical lattice Gauss kernel; Q=I_physical-P0',
            'graph_block_definitions':'B=lambda QTP0; D=Q(H_lambda-volume*E0)Q',
            'exact_Q_block_lower_bound_at_endpoint':str(p['separated']),
            'graph_iteration':'X0=0; X_(n+1)=-D^(-1)B+D^(-1)X_n B^dagger X_n',
            'certified_graph_ball_radius':str(p['radius']),
            'contraction_factor':str(p['contraction']),
            'centered_full_H_norm_upper_bound_at_endpoint':str(p['centered_norm']),
            'normalized_isometry':'J_n=(I;X_n)(I+X_n^dagger X_n)^(-1/2)',
            'same_H_effective_generator':'K_n=J_n^dagger H_lambda J_n',
            'exact_stationary_dressed_density':'P_lambda/2^volume, complete continued low-energy band',
            'stationary_band_is_proved_to_be_full_H_ground_projector':False,
            'certified_iteration_bounds':[iteration_bounds(n,side,kappa) for n in (0,4,8,16)],
            'finite_volume_range_can_shrink_with_volume':True,
            'unit_hopping_ground_or_thermodynamic_phase_transferred':False}
