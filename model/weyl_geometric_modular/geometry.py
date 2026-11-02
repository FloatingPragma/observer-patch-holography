"""Exact diamond geometry, normalization and matter small-ball remainder."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s


def positive(value, name):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError('positive exact rational '+name+' required')
    return F(value)


def matter_remainder_pi_squared(radius, hessian_bound, speed=1):
    """Bound divided by pi^2, for physical/model-time energy density."""
    R = positive(radius, 'radius'); v = positive(speed, 'speed')
    if type(hessian_bound) not in (int, F) or hessian_bound < 0:
        raise ValueError('nonnegative exact Hessian norm bound required')
    return F(4, 35)*F(hessian_bound)*R**6/v


@lru_cache(maxsize=1)
def diamond_certificate():
    t, x, y, z = s.symbols('t x y z', real=True)
    R, r, u, E = s.symbols('R r u E', positive=True)
    coords = (t, x, y, z); eta = (1, -1, -1, -1)
    xi = ((R*R-t*t-x*x-y*y-z*z)/(2*R), -t*x/R, -t*y/R, -t*z/R)
    div = sum(s.diff(xi[mu], coords[mu]) for mu in range(4))
    for mu in range(4):
        for nu in range(4):
            residual = eta[nu]*s.diff(xi[nu], coords[mu])+eta[mu]*s.diff(xi[mu], coords[nu])
            residual -= s.Rational(1, 2)*div*eta[mu]*int(mu == nu)
            if s.expand(residual) != 0:
                raise ArithmeticError('diamond vector is not conformal Killing')
    # E=e^{-2*pi*s}; differentiation of E produces -2*pi*E.
    flow = R*((R+u)-E*(R-u))/((R+u)+E*(R-u))
    if s.factor(-2*E*s.diff(flow, E)-(R*R-flow*flow)/R) != 0:
        raise ArithmeticError('null-coordinate modular flow normalization failed')
    weight = (R*R-r*r)/(2*R)
    integral = s.integrate(4*s.pi*r*r*weight, (r, 0, R))
    second = s.integrate(4*s.pi*r**4*weight, (r, 0, R))
    if integral != 4*s.pi*R**4/15 or second != 4*s.pi*R**6/35:
        raise ArithmeticError('small-ball kernel coefficients changed')
    return {'conformal_Killing_components': 16, 'divergence': '-4*x0/R',
            'inverse_modular_flow_null_coordinate': str(flow),
            'inverse_modular_flow_generator': '2*pi*xi',
            'weight': '(R^2-r^2)/(2R)',
            'integral_weight': '4*pi*R^4/15',
            'integral_weight_times_r_squared': '4*pi*R^6/35',
            'modular_energy_in_model_clock': '(2*pi/v)*integral_B w_R*energy_density',
            'leading_coefficient': '8*pi^2*R^4/(15*v)',
            'C2_matter_remainder_bound': '4*pi^2*M*R^6/(35*v)',
            'remainder_domain': 'operator norm of spatial Hessian of varied energy density <=M on B_R',
            'lattice_to_continuum_modular_convergence_rate_proved': False,
            'curved_metric_remainder_proved': False,
            'collar_or_edge_remainder_proved': False}


def flat_geometry_countercontrol(energy, coupling):
    energy = positive(energy, 'nonvacuum integrated energy')
    coupling = positive(coupling, 'coupling')
    # The vacuum fixes Lambda=0 on the same flat metric. Integrating the 00
    # equation in a positive-energy one-particle state then leaves this defect.
    return coupling*energy


def certificate():
    return {'diamond': diamond_certificate(),
            'example_remainder_over_pi_squared': str(matter_remainder_pi_squared(F(1, 10), 1, F(1, 8))),
            'same_flat_metric_vacuum_and_excitation_Einstein_defect': str(flat_geometry_countercontrol(1, 1)),
            'countercontrol_scope': 'rules out deriving state-responsive Einstein dynamics from the fixed free-background matter data alone',
            'Newton_constant_derived': False, 'area_density_derived': False}
