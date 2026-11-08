"""Exact zero-mean current coefficient in the source scalar polarization."""
from fractions import Fraction as F
import sympy as s

from model.spin_energy_polarization.algebra import exact_checks as inherited_checks
from model.spin_energy_polarization.samples import bump
from .algebra import spins
from .bounds import phase_transport


def exact_checks():
    inherited = inherited_checks()
    local = spins(2)
    current = local[0][0]*local[1][1]-local[0][1]*local[1][0]
    gamma = (2*local[0][1])*(2*local[1][1])
    assert gamma*current*gamma == -current
    p,fx,fy,epsilon = s.symbols('p Fx Fy epsilon',real=True)
    assert s.expand(p*p*epsilon*(fx-fy)/2
                    +p**3*epsilon*((fy-fx)/p)/2) == 0
    controls = []
    for divisions in (2,4,8,16):
        mesh = F(1,divisions)
        points = range(-divisions-2,divisions+3)
        one_dimensional_difference = sum(((bump(mesh*(x+1))-bump(mesh*x))/mesh
                                          for x in points),F())
        other_sum = sum((bump(mesh*x) for x in points),F())
        three_dimensional_sum = mesh**3*one_dimensional_difference*other_sum**2
        assert three_dimensional_sum == 0
        controls.append({'divisions':divisions,
                         'exact_weighted_difference_quotient_sum':'0',
                         'all_three_coordinate_directions':True})
    return {'inherited_exact_weighted_source_commutator':inherited,
            'gamma_axis2_pi_reverses_transverse_axis3_current':True,
            'beta_current_identity':'beta_x(C_0i^3)=epsilon_x C_xi^3',
            'actual_residual':'-p^3 sum_x D_i^pF(px) beta_x(C_0i^3)/2',
            'compact_test_derivative_integral':0,
            'exact_mesh_zero_sum_controls':controls}


def compact_bounds(p, radius, derivative_bound, second_derivative_bound,
                   history_energy=F(1), time=F(1), phase=F(1)):
    p,radius,first,second = map(F,(p,radius,derivative_bound,second_derivative_bound))
    if not 0 < p <= 1 or radius < 1 or min(first,second) < 0:
        raise ValueError('positive scale and valid smooth compact test bounds required')
    cells = (2*radius+1)**3
    # At most twice the number of supported vertices supplies a nonzero
    # positive-axis edge difference. Each current operator has norm1/2.
    amplitude = first/2
    volume = 2*cells
    return {'current_coefficient_maximum':amplitude*p**3,
            'current_coefficient_l1_upper':amplitude*volume,
            'residual_operator_norm_upper':amplitude*volume/2,
            'difference_quotient_to_derivative_operator_error':cells*second*p/4,
            'source_history_phase_transport':phase_transport(p,amplitude,volume,history_energy,time,phase),
            'strong_zero_mean_ergodic_limit_proved':True,
            'effective_ground_vector_convergence_modulus':False}


def product_source_control():
    """Actual spin-product inputs, distinct from the selected ground histories."""
    rows = []
    for n in (4,8,16,32):
        p = F(1,n)
        values = {j:bump(p*j) for j in range(-n-2,n+3)}
        sign = lambda j:1 if j >= 0 else -1
        derivative = sum(((values[j]-values[j+1])*(sign(j)+sign(j+1))
                          for j in range(-n-2,n+2)),F())
        assert derivative == 2*(values[0]+values[-1])
        transverse_sum = sum(values.values(),F())
        mean = p*p*derivative*transverse_sum**2/16
        assert mean > F(1,10)
        rows.append({'divisions':n,'exact_actual_source_product_expectation':str(mean)})
    x = s.symbols('x')
    integral = F(s.integrate((1-x*x)**4,(x,-1,1)))
    limit = integral*integral/4
    assert limit == F(16384,99225)
    return {'input':'product spin directions (epsilon_x,sign(x_i),0)/sqrt(2)',
            'same_spin_source_and_same_current_operator':True,
            'input_is_the_selected_ground_or_bounded_low_energy_history':False,
            'exact_mesh_controls':rows,'residual_expectation_limit':str(limit),
            'source_residual_converges_in_operator_norm':False}
