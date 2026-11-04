"""Exact identities and rational error bounds for actual-ground textures."""
from fractions import Fraction as F

import sympy as s

from ..strong_coupling_refinement.bounds import rational
from ..bipartite_spin_limit.bounds import encode


def static_error(exchange, derivative_bound, support_radius, scale):
    J, K, M, L = (rational(v, positive=True)
                   for v in (exchange, derivative_bound, support_radius, scale))
    if min(K, M, L) < 1:
        raise ValueError('K,M,L must be at least one')
    return 10*J*K*K*(2*M+4)**3/L


def exact_identities():
    a, b, q1, q2, q3 = s.symbols('a b q1 q2 q3', real=True)
    diagonal = (1-2*(q2*q2+q3*q3), 1-2*(q1*q1+q3*q3), 1-2*(q1*q1+q2*q2))
    excess = a*(diagonal[0]-1)+b*(diagonal[1]+diagonal[2]-2)
    lower = -(a+b)*(1-diagonal[0])
    if s.expand(excess-lower+4*b*q1*q1) != 0:
        raise ArithmeticError('rotation energy or positive residual changed')
    w1, w2, w3 = s.symbols('w1 w2 w3', real=True)
    A = s.Matrix([[0, -w3, w2], [w3, 0, -w1], [-w2, w1, 0]])
    C = s.diag(a, b, b)
    expected = -(a+b)*(w2*w2+w3*w3)-2*b*w1*w1
    if s.expand(s.trace(C*A*A)-expected) != 0 or s.trace(C*A) != 0:
        raise ArithmeticError('nonlinear static continuum density changed')
    f = s.symbols('f', real=True)
    degree = 2*s.integrate(s.sin(f)**2, (f, 0, s.pi))/s.pi
    if degree != 1:
        raise ArithmeticError('compact unit winding profile changed')
    x = s.symbols('x', real=True)
    profile = (1-x*x)**4
    square = s.integrate(profile**2, (x, -1, 1))
    derivative = s.integrate(s.diff(profile, x)**2, (x, -1, 1))
    return {
        'rotation_energy_residual': '-4*b*q1^2>=0',
        'one_site_pi_rotation_energy': '-24*J*b>=0',
        'static_density': 'J/2*(kappa*(w2^2+w3^2)-2*b*w1^2)',
        'compact_Hopf_lift_absolute_degree': int(degree),
        'coplanar_C3_test_profile': '(1-x1^2)^4_+*(1-x2^2)^4_+*(1-x3^2)^4_+',
        'coplanar_profile_gradient_integral': str(3*derivative*square**2),
    }


def certificate():
    return {
        'identities': exact_identities(),
        'rational_static_error_example': encode({
            'J': F(1), 'derivatives_up_to_order3_bounded_by_K': F(128),
            'support_radius': F(1), 'scale': F(2**48),
            'conditional_error_upper': static_error(1, 128, 1, 2**48)}),
        'actual_order_lower': 'm>1/1000',
        'actual_kappa_lower': 'kappa>=m^2/2>1/2000000',
        'rotation_orbit_direction_energy_lower': 'J*kappa/2 sum_edges |n_x-n_y|^2',
        'radial_hedgehog_energy_lower_per_radius': 'liminf DeltaE/L>=4*pi*J*kappa',
        'compact_Hopf_energy_in_fixed_J_clock': 'L*E_cont+O(1), E_cont>0',
        'compact_Hopf_energy_in_soft_spin_clock': 'order L^2 for K/(J/L)',
        'actual_quantum_ground_used': True,
        'product_Neel_ground_assumed': False,
        'local_unitary_orbit_only': True,
        'energy_minimized_over_all_states': False,
        'exact_microscopic_topological_superselection': False,
        'dynamical_sigma_model_derived': False,
        'finite_mass_topological_particle_proved': False,
        'divergent_mean_excludes_finite_energy_spectral_component': False,
        'all_possible_quantum_defects_excluded': False,
        'interacting_relativistic_QFT_complete': False,
    }
